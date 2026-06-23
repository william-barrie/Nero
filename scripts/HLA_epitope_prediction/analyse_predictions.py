#!/usr/bin/env python3
"""
Analyse HLA epitope predictions to identify strong and weak binders.

This script aggregates predictions by allele and reports:
- Number of unique proteins with strong binders
- Percentage of strong and weak binders per HLA allele
"""

import sys
import argparse
from pathlib import Path
import pandas as pd
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def analyse_class(df_class, hla_class, strong_threshold, weak_threshold,
                  output_dir, organism_name):
    """
    Analyse predictions for a single HLA class.

    Args:
        df_class: DataFrame with predictions for this class
        hla_class: 'I' or 'II'
        strong_threshold: %Rank threshold for strong binders
        weak_threshold: %Rank threshold for weak binders
        output_dir: Output directory
        organism_name: Organism name for filenames

    Returns:
        Dictionary with summary statistics
    """
    if len(df_class) == 0:
        logger.warning(f"No Class {hla_class} predictions found")
        return None

    logger.info(f"\n=== Class {hla_class} Analysis ===")
    logger.info(f"Total Class {hla_class} predictions: {len(df_class)}")

    # Calculate per-allele statistics
    allele_stats = []

    for allele in sorted(df_class['hla_allele'].unique()):
        df_allele = df_class[df_class['hla_allele'] == allele]
        total_peptides = len(df_allele)

        # Count binders
        strong_mask = df_allele['rank'] < strong_threshold
        weak_mask = (df_allele['rank'] >= strong_threshold) & (df_allele['rank'] < weak_threshold)

        n_strong = strong_mask.sum()
        n_weak = weak_mask.sum()

        pct_strong = 100.0 * n_strong / total_peptides if total_peptides > 0 else 0
        pct_weak = 100.0 * n_weak / total_peptides if total_peptides > 0 else 0

        # Count unique proteins with strong binders
        unique_proteins_strong = df_allele[strong_mask]['protein_id'].nunique() if n_strong > 0 else 0

        allele_stats.append({
            'hla_allele': allele,
            'total_peptides': total_peptides,
            'strong_binders_n': n_strong,
            'strong_binders_pct': round(pct_strong, 2),
            'weak_binders_n': n_weak,
            'weak_binders_pct': round(pct_weak, 2),
            'unique_proteins_with_strong_binders': unique_proteins_strong
        })

    # Create DataFrame and save
    stats_df = pd.DataFrame(allele_stats)
    output_file = output_dir / f"{organism_name}_class_{hla_class.lower()}_allele_statistics.tsv"
    stats_df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Saved Class {hla_class} allele statistics to {output_file}")

    # Print top alleles by unique proteins with strong binders
    top_alleles = stats_df.nlargest(10, 'unique_proteins_with_strong_binders')
    if len(top_alleles) > 0 and top_alleles['unique_proteins_with_strong_binders'].max() > 0:
        logger.info(f"\nTop 10 Class {hla_class} alleles by unique proteins with strong binders:")
        for _, row in top_alleles.iterrows():
            logger.info(f"  {row['hla_allele']}: {row['unique_proteins_with_strong_binders']} proteins "
                       f"({row['strong_binders_pct']}% strong, {row['weak_binders_pct']}% weak)")

    # Overall summary
    total_strong = stats_df['strong_binders_n'].sum()
    total_weak = stats_df['weak_binders_n'].sum()
    total_all = stats_df['total_peptides'].sum()

    overall_pct_strong = 100.0 * total_strong / total_all if total_all > 0 else 0
    overall_pct_weak = 100.0 * total_weak / total_all if total_all > 0 else 0

    logger.info(f"\nOverall Class {hla_class} statistics:")
    logger.info(f"  Strong binders (rank < {strong_threshold}): {total_strong:,} ({overall_pct_strong:.2f}%)")
    logger.info(f"  Weak binders (rank {strong_threshold}-{weak_threshold}): {total_weak:,} ({overall_pct_weak:.2f}%)")
    logger.info(f"  Total predictions: {total_all:,}")

    return {
        'class': hla_class,
        'total_predictions': total_all,
        'strong_binders': total_strong,
        'strong_binders_pct': round(overall_pct_strong, 2),
        'weak_binders': total_weak,
        'weak_binders_pct': round(overall_pct_weak, 2),
        'unique_proteins': df_class['protein_id'].nunique(),
        'proteins_with_strong_binders': df_class[df_class['rank'] < strong_threshold]['protein_id'].nunique() if total_strong > 0 else 0,
        'n_alleles': df_class['hla_allele'].nunique(),
        'strong_threshold': strong_threshold,
        'weak_threshold': weak_threshold
    }


def analyse_predictions(
    input_file: Path,
    output_dir: Path,
    class_i_strong_threshold: float = 0.5,
    class_i_weak_threshold: float = 2.0,
    class_ii_strong_threshold: float = 1.0,
    class_ii_weak_threshold: float = 5.0,
    organism_name: str = None
):
    """
    Analyse HLA predictions and generate summary statistics.

    Args:
        input_file: Path to combined predictions CSV
        output_dir: Directory for output files
        class_i_strong_threshold: %Rank threshold for Class I strong binders
        class_i_weak_threshold: %Rank threshold for Class I weak binders
        class_ii_strong_threshold: %Rank threshold for Class II strong binders
        class_ii_weak_threshold: %Rank threshold for Class II weak binders
        organism_name: Optional organism name for output files
    """
    logger.info(f"Loading predictions from {input_file}")

    try:
        df = pd.read_csv(input_file)
    except FileNotFoundError:
        logger.error(f"File not found: {input_file}")
        return False
    except pd.errors.EmptyDataError:
        logger.error(f"File is empty: {input_file}")
        return False

    # Validate required columns
    required_cols = ['hla_class', 'hla_allele', 'protein_id', 'rank']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        logger.error(f"Missing required columns: {missing_cols}")
        logger.info(f"Available columns: {df.columns.tolist()}")
        return False

    logger.info(f"Loaded {len(df):,} predictions")
    logger.info(f"Classes found: {df['hla_class'].unique()}")

    output_dir.mkdir(parents=True, exist_ok=True)

    # Determine organism name from input file if not provided
    if organism_name is None:
        organism_name = input_file.stem.replace('_hla_predictions', '')

    # Analyse Class I
    df_i = df[df['hla_class'] == 'I'].copy()
    summary_i = analyse_class(df_i, 'I', class_i_strong_threshold, class_i_weak_threshold,
                               output_dir, organism_name)

    # Analyse Class II
    df_ii = df[df['hla_class'] == 'II'].copy()
    summary_ii = analyse_class(df_ii, 'II', class_ii_strong_threshold, class_ii_weak_threshold,
                                output_dir, organism_name)

    # Generate combined summary
    logger.info("\n=== Combined Summary ===")

    summary_data = []
    if summary_i:
        summary_data.append(summary_i)
    if summary_ii:
        summary_data.append(summary_ii)

    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_file = output_dir / f"{organism_name}_summary.tsv"
        summary_df.to_csv(summary_file, sep='\t', index=False)
        logger.info(f"\nSaved summary to {summary_file}")
        logger.info("\n" + summary_df.to_string(index=False))

    return True


def organism_name_from_predictions(predictions_file: Path) -> str:
    """Derive the organism/strain (or sequence identifier) name from a
    predictions file written by the pipeline.

    Files are named ``{name}_hla_predictions.csv`` (or, with the legacy
    convention, ``{name}_hla_binding_results.csv``), where ``{name}`` is the
    organism/strain or, for bare-identifier FASTA records, the sequence ID.
    """
    stem = predictions_file.stem
    for suffix in ('_hla_predictions', '_hla_binding_results'):
        if stem.endswith(suffix):
            return stem[: -len(suffix)]
    return stem


def discover_prediction_dirs(results_dir: Path,
                             select: list = None) -> list:
    """Find every directory under ``results_dir`` that holds a predictions file.

    This is the companion to running ``hla_epitope_predictor.py --fasta``: that
    run writes one ``<name>/<name>_hla_predictions.csv`` per organism/strain (or
    per sequence identifier when headers carry no organism), and this walks the
    base ``--output-dir`` to collect them all automatically.

    Args:
        results_dir: Base output directory from a prediction run.
        select: Optional list of names. When given, only directories whose
            organism/identifier name matches one of these (case-insensitive
            substring) are returned, so a subset of sequence identifiers can be
            analysed.

    Returns:
        Sorted list of directories, each containing a predictions CSV.
    """
    results_dir = Path(results_dir)
    if not results_dir.exists():
        raise FileNotFoundError(f"Results directory not found: {results_dir}")

    predictions_files = sorted(results_dir.rglob("*_hla_predictions.csv"))
    if not predictions_files:
        predictions_files = sorted(results_dir.rglob("*_hla_binding_results.csv"))

    if not predictions_files:
        raise ValueError(
            f"No '*_hla_predictions.csv' files found under {results_dir}. "
            f"Point --results-dir at the --output-dir used for the prediction run."
        )

    # De-duplicate directories while preserving order.
    dirs = []
    seen = set()
    for pf in predictions_files:
        if pf.parent not in seen:
            seen.add(pf.parent)
            dirs.append(pf.parent)

    if select:
        sel_lower = [s.lower() for s in select]
        filtered = []
        for d in dirs:
            pf = (sorted(d.glob("*_hla_predictions.csv")) or
                  sorted(d.glob("*_hla_binding_results.csv")))[0]
            name = organism_name_from_predictions(pf).lower()
            if any(s in name or s in d.name.lower() for s in sel_lower):
                filtered.append(d)

        if not filtered:
            available = ", ".join(
                organism_name_from_predictions(
                    (sorted(d.glob("*_hla_predictions.csv")) or
                     sorted(d.glob("*_hla_binding_results.csv")))[0]
                ) for d in dirs
            )
            raise ValueError(
                f"No predictions matched --select {select}. "
                f"Available: {available}"
            )
        dirs = filtered

    return dirs


def batch_analyse(analysis_dirs: list, output_dir: Path,
                  class_i_strong_threshold: float = 0.5,
                  class_i_weak_threshold: float = 2.0,
                  class_ii_strong_threshold: float = 1.0,
                  class_ii_weak_threshold: float = 5.0):
    """
    Analyse multiple organism predictions in batch.

    Args:
        analysis_dirs: List of analysis directories
        output_dir: Output directory for aggregated results
        class_i_strong_threshold: Class I strong threshold
        class_i_weak_threshold: Class I weak threshold
        class_ii_strong_threshold: Class II strong threshold
        class_ii_weak_threshold: Class II weak threshold
    """
    logger.info(f"Batch analysis of {len(analysis_dirs)} organisms")
    output_dir.mkdir(parents=True, exist_ok=True)

    for analysis_dir in analysis_dirs:
        analysis_path = Path(analysis_dir)

        if not analysis_path.exists():
            logger.warning(f"Directory not found: {analysis_path}")
            continue

        # Find predictions file
        predictions_files = list(analysis_path.glob("*_hla_predictions.csv"))

        if not predictions_files:
            # Try old naming convention
            predictions_files = list(analysis_path.glob("*_hla_binding_results.csv"))

        if not predictions_files:
            logger.warning(f"No predictions file found in {analysis_path}")
            continue

        predictions_file = predictions_files[0]
        # Prefer the name embedded in the predictions filename (the organism /
        # strain / sequence identifier), falling back to the directory name.
        organism_name = organism_name_from_predictions(predictions_file)
        if not organism_name:
            organism_name = analysis_path.name.replace('_analysis', '')

        logger.info(f"\n{'='*60}")
        logger.info(f"Analysing: {organism_name}")
        logger.info(f"{'='*60}")

        analyse_predictions(
            predictions_file,
            output_dir,
            class_i_strong_threshold,
            class_i_weak_threshold,
            class_ii_strong_threshold,
            class_ii_weak_threshold,
            organism_name
        )


def main():
    parser = argparse.ArgumentParser(
        description="Analyse HLA epitope predictions",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyse single organism
  python analyse_predictions.py --input smallpox_analysis/smallpox_hla_predictions.csv

  # Batch analyse multiple organisms
  python analyse_predictions.py --batch \\
      smallpox_analysis HIV_analysis plague_analysis

  # Analyse every organism/identifier produced by a --fasta prediction run
  # (auto-discovers each <name>/<name>_hla_predictions.csv under the dir)
  python analyse_predictions.py --results-dir ./results

  # Same, but only a subset of sequence identifiers / organisms
  python analyse_predictions.py --results-dir ./results --select seq_A seq_B

  # Custom thresholds
  python analyse_predictions.py --input results.csv \\
      --class-i-strong 0.05 --class-i-weak 0.5 \\
      --class-ii-strong 0.1 --class-ii-weak 1.0
        """
    )

    # Input mode
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        '--input', type=Path,
        help='Single predictions CSV file'
    )
    input_group.add_argument(
        '--batch', nargs='+',
        help='Multiple analysis directories to process'
    )
    input_group.add_argument(
        '--results-dir', type=Path,
        help='Base output directory from a prediction run (e.g. the --output-dir '
             'used with hla_epitope_predictor.py --fasta). Every '
             '<name>/<name>_hla_predictions.csv beneath it is analysed.'
    )

    parser.add_argument(
        '--select', nargs='+', metavar='NAME',
        help='With --results-dir, only analyse organisms/strains or sequence '
             'identifiers whose name matches one of these (case-insensitive '
             'substring).'
    )

    # Class I thresholds
    parser.add_argument(
        '--class-i-strong', type=float, default=0.5,
        help='%%Rank threshold for Class I strong binders (default: 0.5)'
    )
    parser.add_argument(
        '--class-i-weak', type=float, default=2.0,
        help='%%Rank threshold for Class I weak binders (default: 2.0)'
    )

    # Class II thresholds
    parser.add_argument(
        '--class-ii-strong', type=float, default=1.0,
        help='%%Rank threshold for Class II strong binders (default: 1.0)'
    )
    parser.add_argument(
        '--class-ii-weak', type=float, default=5.0,
        help='%%Rank threshold for Class II weak binders (default: 5.0)'
    )

    # Output
    parser.add_argument(
        '--output-dir', type=Path, default='./analysed_results',
        help='Output directory (default: ./analysed_results)'
    )

    args = parser.parse_args()

    if args.select and not args.results_dir:
        parser.error("--select can only be used together with --results-dir.")

    try:
        if args.input:
            # Single file analysis
            success = analyse_predictions(
                args.input,
                args.output_dir,
                args.class_i_strong,
                args.class_i_weak,
                args.class_ii_strong,
                args.class_ii_weak
            )
            sys.exit(0 if success else 1)

        elif args.batch:
            # Batch analysis
            batch_analyse(
                args.batch,
                args.output_dir,
                args.class_i_strong,
                args.class_i_weak,
                args.class_ii_strong,
                args.class_ii_weak
            )
            sys.exit(0)

        elif args.results_dir:
            # Auto-discover per-organism/identifier prediction dirs (e.g. from
            # a --fasta run) and batch-analyse them, optionally filtered.
            analysis_dirs = discover_prediction_dirs(args.results_dir, args.select)
            logger.info(
                f"Discovered {len(analysis_dirs)} prediction directory(ies) "
                f"under {args.results_dir}"
            )
            batch_analyse(
                analysis_dirs,
                args.output_dir,
                args.class_i_strong,
                args.class_i_weak,
                args.class_ii_strong,
                args.class_ii_weak
            )
            sys.exit(0)

    except Exception as e:
        logger.error(f"Analysis failed: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
