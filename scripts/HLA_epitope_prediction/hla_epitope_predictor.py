#!/usr/bin/env python3
"""
HLA Epitope Prediction Pipeline

Predicts HLA Class I and II binding for pathogen proteomes with support
for incremental processing and result caching.

Author: Refactored for better maintainability and incremental processing
"""

import os
import re
import sys
import subprocess
import json
import shutil
from pathlib import Path
from typing import List, Dict, Set, Tuple, Optional
from dataclasses import dataclass
import argparse
import logging

import pandas as pd
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class AnalysisConfig:
    """Configuration for HLA epitope analysis."""
    email: str
    organism_query: str
    output_dir: Path
    organism_name: str
    netmhcpan_path: Path
    netmhciipan_path: Path
    signalp_model_dir: Path
    n_parallel_jobs: int = 20
    class_i_peptide_lengths: List[int] = None
    class_i_step_size: int = 3
    class_ii_peptide_length: int = 15
    class_ii_step_size: int = 3
    batch_size_proteins: int = 100
    max_proteins: int = 15000
    temp_dir: Optional[Path] = None

    def __post_init__(self):
        if self.class_i_peptide_lengths is None:
            self.class_i_peptide_lengths = [9]  # Conservative default: 9-mers only
        self.output_dir = Path(self.output_dir)
        self.netmhcpan_path = Path(self.netmhcpan_path)
        self.netmhciipan_path = Path(self.netmhciipan_path)
        self.signalp_model_dir = Path(self.signalp_model_dir)

        # Set up temp directory
        if self.temp_dir is None:
            # Default: use a temp subdirectory within output_dir
            self.temp_dir = self.output_dir / "temp"
        else:
            self.temp_dir = Path(self.temp_dir)

        self.temp_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Temporary files will be stored in: {self.temp_dir}")


class ProteomeDownloader:
    """Handles downloading proteomes from NCBI."""

    def __init__(self, email: str, batch_size: int = 100):
        self.email = email
        self.batch_size = batch_size
        Entrez.email = email

    def download(self, organism_query: str, output_file: Path, max_proteins: int = 15000) -> Path:
        """
        Download proteome from NCBI.

        Args:
            organism_query: NCBI search query (e.g., "Variola virus[Organism]")
            output_file: Path to save FASTA file
            max_proteins: Maximum number of proteins to download

        Returns:
            Path to downloaded FASTA file
        """
        logger.info(f"Searching NCBI for: {organism_query}")

        # Search for protein sequences
        search_handle = Entrez.esearch(
            db="protein",
            term=organism_query,
            retmax=max_proteins
        )
        search_results = Entrez.read(search_handle)
        search_handle.close()

        protein_ids = search_results["IdList"]
        logger.info(f"Found {len(protein_ids)} proteins")

        if not protein_ids:
            raise ValueError(f"No proteins found for: {organism_query}")

        # Download in batches
        all_records = []
        n_batches = (len(protein_ids) - 1) // self.batch_size + 1

        for i in range(0, len(protein_ids), self.batch_size):
            batch_ids = protein_ids[i:i + self.batch_size]
            batch_num = i // self.batch_size + 1
            logger.info(f"Downloading batch {batch_num}/{n_batches}")

            fetch_handle = Entrez.efetch(
                db="protein",
                id=batch_ids,
                rettype="fasta",
                retmode="text"
            )

            batch_records = list(SeqIO.parse(fetch_handle, "fasta"))
            all_records.extend(batch_records)
            fetch_handle.close()

        # Write to file
        output_file.parent.mkdir(parents=True, exist_ok=True)
        with open(output_file, 'w') as f:
            SeqIO.write(all_records, f, "fasta")

        logger.info(f"Downloaded {len(all_records)} sequences to {output_file}")
        return output_file


def sanitize_name(name: str) -> str:
    """Convert an arbitrary organism/strain/sequence name into a safe filename."""
    safe = re.sub(r'[^\w\-.]+', '_', name.strip())
    safe = safe.strip('_')
    return safe or "unknown"


class FastaInputHandler:
    """Parse a local FASTA file and group sequences by organism/strain.

    Supports FASTA files that contain sequences from multiple organisms or
    strains. The organism/strain is inferred from each record's header line so
    that the pipeline can either analyse every organism in the file
    automatically (one analysis per organism) or just a user-selected subset.

    Recognised header conventions:
      - NCBI style: ``>ACC description [Organism name strain X]`` (organism in
        the trailing square brackets).
      - UniProt style: ``>db|ACC|NAME ... OS=Organism name OX=taxid ...``
        (organism after the ``OS=`` token).
    """

    UNKNOWN_ORGANISM = "unknown_organism"

    @staticmethod
    def parse_organism_from_header(description: str) -> str:
        """Infer the organism/strain name from a FASTA header description.

        Returns ``UNKNOWN_ORGANISM`` when no organism can be confidently
        extracted so those sequences are still grouped together rather than
        silently dropped.
        """
        if not description:
            return FastaInputHandler.UNKNOWN_ORGANISM

        # UniProt style: organism follows "OS=" up to the next KEY= token.
        os_match = re.search(r'\bOS=(.+?)(?:\s+[A-Z]{2}=|$)', description)
        if os_match:
            return os_match.group(1).strip()

        # NCBI style: organism (and often strain) in trailing square brackets.
        brackets = re.findall(r'\[([^\[\]]+)\]', description)
        if brackets:
            return brackets[-1].strip()

        return FastaInputHandler.UNKNOWN_ORGANISM

    @staticmethod
    def group_by_organism(fasta_file: Path) -> "dict":
        """Group the records of a FASTA file by inferred organism/strain.

        Returns an insertion-ordered dict mapping organism name to the list of
        SeqRecords belonging to it.
        """
        fasta_file = Path(fasta_file)
        if not fasta_file.exists():
            raise FileNotFoundError(f"FASTA file not found: {fasta_file}")

        groups: Dict[str, List[SeqRecord]] = {}
        n_records = 0
        with open(fasta_file) as f:
            for record in SeqIO.parse(f, "fasta"):
                n_records += 1
                organism = FastaInputHandler.parse_organism_from_header(record.description)
                groups.setdefault(organism, []).append(record)

        if n_records == 0:
            raise ValueError(f"No sequences found in FASTA file: {fasta_file}")

        logger.info(
            f"Parsed {n_records} sequences from {fasta_file} "
            f"spanning {len(groups)} organism/strain group(s)"
        )
        for organism, records in groups.items():
            logger.info(f"  - {organism}: {len(records)} sequence(s)")

        return groups

    @staticmethod
    def select_groups(groups: "dict", selections: Optional[List[str]]) -> "dict":
        """Filter grouped records down to a user-specified selection.

        Each selection string is matched (case-insensitively) first against
        organism/strain names (substring match) and then, if it matches no
        organism, against individual sequence IDs/descriptions. Selected
        sequences that are picked by ID are grouped under that organism so the
        normal per-organism workflow still applies.

        With no selections, all groups are returned unchanged.
        """
        if not selections:
            return groups

        selected: Dict[str, List[SeqRecord]] = {}
        for selection in selections:
            sel_lower = selection.lower()

            # 1) Match against organism/strain group names.
            org_matches = [org for org in groups if sel_lower in org.lower()]
            if org_matches:
                for org in org_matches:
                    selected.setdefault(org, [])
                    existing_ids = {r.id for r in selected[org]}
                    for rec in groups[org]:
                        if rec.id not in existing_ids:
                            selected[org].append(rec)
                continue

            # 2) Otherwise match against individual sequence IDs/descriptions.
            id_matched = False
            for org, records in groups.items():
                for rec in records:
                    if sel_lower == rec.id.lower() or sel_lower in rec.description.lower():
                        selected.setdefault(org, [])
                        if rec.id not in {r.id for r in selected[org]}:
                            selected[org].append(rec)
                        id_matched = True

            if not id_matched:
                available = ", ".join(groups.keys())
                raise ValueError(
                    f"Selection '{selection}' did not match any organism/strain "
                    f"or sequence ID in the FASTA file. "
                    f"Available organisms: {available}"
                )

        return selected

    @staticmethod
    def write_group(records: List[SeqRecord], output_file: Path) -> Path:
        """Write a group of records to a FASTA file."""
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        with open(output_file, 'w') as f:
            SeqIO.write(records, f, "fasta")
        logger.info(f"Wrote {len(records)} sequences to {output_file}")
        return output_file


class SignalPFilter:
    """Filters proteins using SignalP to identify secreted proteins."""

    def __init__(self, model_dir: Path):
        self.model_dir = model_dir

    def _create_short_id_fasta(self, input_fasta: Path, output_fasta: Path) -> Dict[str, str]:
        """Create FASTA with shortened IDs to avoid SignalP filename issues."""
        id_mapping = {}
        short_records = []

        with open(input_fasta) as f:
            for i, record in enumerate(SeqIO.parse(f, "fasta")):
                short_id = f"prot_{i+1:06d}"
                id_mapping[short_id] = record.id

                short_record = SeqRecord(
                    record.seq,
                    id=short_id,
                    description=""
                )
                short_records.append(short_record)

        with open(output_fasta, 'w') as f:
            SeqIO.write(short_records, f, "fasta")

        logger.info(f"Created short-ID FASTA with {len(short_records)} sequences")
        return id_mapping

    def run(self, input_fasta: Path, output_dir: Path, organism_type: str = "other") -> List[str]:
        """
        Run SignalP and return list of secreted protein IDs.

        Args:
            input_fasta: Input proteome FASTA
            output_dir: Directory for SignalP output
            organism_type: "eukarya", "gram+", "gram-", or "other"

        Returns:
            List of protein IDs with signal peptides
        """
        logger.info("Running SignalP analysis...")
        output_dir.mkdir(parents=True, exist_ok=True)

        # Create shortened FASTA
        short_fasta = output_dir / "proteins_short_ids.fasta"
        id_mapping = self._create_short_id_fasta(input_fasta, short_fasta)

        # Run SignalP
        cmd = [
            "signalp6",
            "--fastafile", str(short_fasta),
            "--output_dir", str(output_dir),
            "--format", "none",
            "--mode", "fast",
            "--organism", organism_type,
            "--model_dir", str(self.model_dir)
        ]

        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info("SignalP completed successfully")
        except subprocess.CalledProcessError as e:
            logger.error(f"SignalP failed: {e.stderr}")
            raise

        # Parse results
        secreted_ids = []
        prediction_file = output_dir / "prediction_results.txt"

        if prediction_file.exists():
            with open(prediction_file) as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 2 and parts[1] == "SP":
                        short_id = parts[0]
                        original_id = id_mapping.get(short_id, short_id)
                        secreted_ids.append(original_id)

        logger.info(f"Found {len(secreted_ids)} secreted proteins")
        return secreted_ids

    def filter_fasta(self, input_fasta: Path, output_fasta: Path, protein_ids: List[str]) -> Path:
        """Filter FASTA file to keep only specified protein IDs."""
        logger.info(f"Filtering FASTA to {len(protein_ids)} proteins...")

        protein_id_set = set(protein_ids)
        filtered_records = []

        with open(input_fasta) as f:
            for record in SeqIO.parse(f, "fasta"):
                # Handle different FASTA header formats
                record_id = record.id.split('|')[1] if '|' in record.id else record.id

                if record_id in protein_id_set:
                    filtered_records.append(record)

        output_fasta.parent.mkdir(parents=True, exist_ok=True)
        with open(output_fasta, 'w') as f:
            SeqIO.write(filtered_records, f, "fasta")

        logger.info(f"Filtered to {len(filtered_records)} sequences")
        return output_fasta


class PeptideGenerator:
    """Generates peptides from protein sequences."""

    @staticmethod
    def generate_peptides(fasta_file: Path, peptide_length: int, step_size: int = 1) -> List[Dict]:
        """
        Generate overlapping peptides of fixed length.

        Args:
            fasta_file: Input FASTA file
            peptide_length: Length of peptides to generate
            step_size: Step size for sliding window (default: 1 = fully overlapping)

        Returns:
            List of dicts with protein_id, peptide_sequence, position
        """
        logger.info(f"Generating {peptide_length}-mer peptides (step={step_size})...")

        peptides = []
        with open(fasta_file) as f:
            for record in SeqIO.parse(f, "fasta"):
                protein_id = record.id
                sequence = str(record.seq)

                for i in range(0, len(sequence) - peptide_length + 1, step_size):
                    peptide = sequence[i:i + peptide_length]
                    # Skip peptides with ambiguous amino acids
                    if 'X' not in peptide and '*' not in peptide:
                        peptides.append({
                            'protein_id': protein_id,
                            'peptide_sequence': peptide,
                            'position': i + 1,
                            'peptide_length': peptide_length
                        })

        logger.info(f"Generated {len(peptides)} peptides")
        return peptides

    @staticmethod
    def generate_variable_length_peptides(fasta_file: Path, peptide_lengths: List[int], step_size: int = 1) -> List[Dict]:
        """
        Generate peptides of multiple lengths (for Class I).

        Args:
            fasta_file: Input FASTA file
            peptide_lengths: List of peptide lengths to generate
            step_size: Step size for sliding window (default: 1 = fully overlapping)

        Returns:
            List of dicts with protein_id, peptide_sequence, position, peptide_length
        """
        logger.info(f"Generating peptides of lengths: {peptide_lengths} (step={step_size})")

        peptides = []
        with open(fasta_file) as f:
            for record in SeqIO.parse(f, "fasta"):
                protein_id = record.id
                sequence = str(record.seq)

                for length in peptide_lengths:
                    for i in range(0, len(sequence) - length + 1, step_size):
                        peptide = sequence[i:i + length]
                        if 'X' not in peptide and '*' not in peptide:
                            peptides.append({
                                'protein_id': protein_id,
                                'peptide_sequence': peptide,
                                'position': i + 1,
                                'peptide_length': length
                            })

        logger.info(f"Generated {len(peptides)} variable-length peptides")
        return peptides


class HLAPredictor:
    """Base class for HLA binding prediction."""

    def __init__(self, tool_path: Path, output_dir: Path, temp_dir: Optional[Path] = None):
        self.tool_path = tool_path
        self.output_dir = output_dir
        self.results_dir = output_dir / "allele_results"
        self.results_dir.mkdir(parents=True, exist_ok=True)
        self.temp_dir = temp_dir

        if not self.tool_path.exists():
            raise FileNotFoundError(f"Tool not found: {self.tool_path}")

    def get_completed_alleles(self) -> Set[str]:
        """Get set of alleles that have already been processed."""
        completed = set()
        status_file = self.output_dir / "completed_alleles.json"

        if status_file.exists():
            with open(status_file) as f:
                completed = set(json.load(f))

        return completed

    def mark_allele_complete(self, allele: str):
        """Mark an allele as completed."""
        completed = self.get_completed_alleles()
        completed.add(allele)

        status_file = self.output_dir / "completed_alleles.json"
        with open(status_file, 'w') as f:
            json.dump(sorted(list(completed)), f, indent=2)

    def predict_allele(self, allele: str, peptide_file: Path) -> Path:
        """Run prediction for a single allele. Must be implemented by subclass."""
        raise NotImplementedError

    def run_predictions(self, peptides: List[Dict], alleles: List[str],
                       n_jobs: int = 1, force: bool = False) -> Path:
        """
        Run predictions for all alleles with incremental processing.

        Args:
            peptides: List of peptide dictionaries
            alleles: List of HLA alleles
            n_jobs: Number of parallel jobs
            force: If True, rerun even completed alleles

        Returns:
            Path to combined results file
        """
        logger.info(f"Running predictions for {len(alleles)} alleles")

        # Write peptides to file
        peptide_file = self.output_dir / "peptides.txt"
        with open(peptide_file, 'w') as f:
            for peptide in peptides:
                f.write(f"{peptide['peptide_sequence']}\n")

        # Filter to uncompleted alleles
        completed = set() if force else self.get_completed_alleles()
        alleles_to_run = [a for a in alleles if a not in completed]

        if alleles_to_run:
            logger.info(f"Processing {len(alleles_to_run)} alleles ({len(completed)} already completed)")
        else:
            logger.info(f"All {len(alleles)} alleles already completed")
            return self._combine_results(alleles)

        # Create commands for parallel execution
        commands_file = self.output_dir / "prediction_commands.txt"
        with open(commands_file, 'w') as f:
            for allele in alleles_to_run:
                result_file = self.results_dir / f"{self._sanitize_allele_name(allele)}.txt"
                cmd = self._build_command(allele, peptide_file, result_file)
                f.write(f"{cmd}\n")

        # Run in parallel with TMPDIR set to redirect temporary files
        logger.info(f"Running {len(alleles_to_run)} predictions with {n_jobs} parallel jobs...")
        parallel_cmd = f"cat {commands_file} | parallel -j {n_jobs}"

        # Set environment with custom TMPDIR to control temporary file location
        env = os.environ.copy()
        if self.temp_dir:
            env['TMPDIR'] = str(self.temp_dir)
            logger.info(f"Setting TMPDIR={self.temp_dir} for netMHC tools")

        try:
            subprocess.run(parallel_cmd, shell=True, check=True, capture_output=True, text=True, env=env)
            logger.info("Predictions completed successfully")

            # Mark all as complete
            for allele in alleles_to_run:
                self.mark_allele_complete(allele)

        except subprocess.CalledProcessError as e:
            logger.warning(f"Some predictions failed. Check individual results.")
            # Mark successful ones as complete
            for allele in alleles_to_run:
                result_file = self.results_dir / f"{self._sanitize_allele_name(allele)}.txt"
                if result_file.exists() and result_file.stat().st_size > 0:
                    self.mark_allele_complete(allele)

        # Combine all results
        return self._combine_results(alleles)

    def _sanitize_allele_name(self, allele: str) -> str:
        """Convert allele name to safe filename."""
        return allele.replace(':', '_').replace('*', '_').replace('/', '_')

    def _build_command(self, allele: str, peptide_file: Path, result_file: Path) -> str:
        """Build command line for prediction tool."""
        raise NotImplementedError

    def _combine_results(self, alleles: List[str]) -> Path:
        """Combine individual allele results into single file."""
        combined_file = self.output_dir / "combined_results.txt"

        logger.info(f"Combining results from {len(alleles)} alleles...")
        with open(combined_file, 'w') as outf:
            for allele in alleles:
                result_file = self.results_dir / f"{self._sanitize_allele_name(allele)}.txt"

                if not result_file.exists():
                    logger.warning(f"Missing results for {allele}")
                    continue

                try:
                    with open(result_file) as inf:
                        outf.write(inf.read())
                        outf.write('\n')
                except Exception as e:
                    logger.warning(f"Failed to read {allele}: {e}")

        logger.info(f"Combined results saved to {combined_file}")
        return combined_file


class NetMHCpanPredictor(HLAPredictor):
    """HLA Class I predictor using netMHCpan."""

    def _build_command(self, allele: str, peptide_file: Path, result_file: Path) -> str:
        return (f"{self.tool_path} -f {peptide_file} -a {allele} "
                f"-inptype 1 -xls -xlsfile {result_file}")


class NetMHCIIpanPredictor(HLAPredictor):
    """HLA Class II predictor using netMHCIIpan."""

    def _build_command(self, allele: str, peptide_file: Path, result_file: Path) -> str:
        return (f"{self.tool_path} -f {peptide_file} -a {allele} "
                f"-inptype 1 -xls -xlsfile {result_file}")


class ResultParser:
    """Parses prediction results into DataFrame."""

    @staticmethod
    def parse_netmhc_results(results_file: Path, peptides: List[Dict],
                            hla_class: str) -> pd.DataFrame:
        """
        Parse netMHCpan or netMHCIIpan results.

        Args:
            results_file: Combined results file
            peptides: Original peptide list
            hla_class: 'I' or 'II'

        Returns:
            DataFrame with columns: protein_id, peptide_sequence, position,
                                   peptide_length, hla_class, hla_allele,
                                   score, rank
        """
        logger.info(f"Parsing HLA Class {hla_class} results...")

        # Create peptide lookup
        peptide_lookup = {p['peptide_sequence']: p for p in peptides}

        all_results = []
        current_allele = None
        lines_processed = 0
        alleles_found = []
        data_lines_attempted = 0

        with open(results_file) as f:
            for line in f:
                lines_processed += 1
                line = line.strip()

                if not line:
                    continue

                # Detect allele name (header line without tabs)
                if '\t' not in line and (line.startswith('HLA-') or
                                        line.startswith('DRB') or
                                        line.startswith('DQA') or
                                        line.startswith('DPA') or
                                        line.startswith('DPB') or
                                        line.startswith('DQB')):
                    if 'Pos' not in line:
                        current_allele = line
                        alleles_found.append(line)
                        logger.debug(f"Found allele: {current_allele}")
                    continue

                # Skip header lines
                if line.startswith('Pos') or 'Peptide' in line:
                    continue

                # Parse data lines
                parts = line.split('\t')
                if len(parts) >= 7 and parts[0].isdigit():
                    data_lines_attempted += 1
                    try:
                        peptide_seq = parts[1]

                        # Skip invalid peptides
                        if not peptide_seq or pd.isna(peptide_seq):
                            continue

                        # Get score and rank (columns may vary)
                        score = float(parts[5]) if parts[5] != 'NA' else None
                        rank = float(parts[6]) if parts[6] != 'NA' else None

                        # For Class II, might be in different columns
                        if hla_class == 'II' and len(parts) >= 8:
                            score = float(parts[6]) if parts[6] != 'NA' else None
                            rank = float(parts[7]) if parts[7] != 'NA' else None

                        # Get protein metadata
                        peptide_info = peptide_lookup.get(peptide_seq, {})

                        all_results.append({
                            'protein_id': peptide_info.get('protein_id', 'Unknown'),
                            'peptide_sequence': peptide_seq,
                            'position': peptide_info.get('position', 0),
                            'peptide_length': peptide_info.get('peptide_length', len(peptide_seq)),
                            'hla_class': hla_class,
                            'hla_allele': current_allele,
                            'score': score,
                            'rank': rank
                        })

                    except (ValueError, IndexError) as e:
                        logger.debug(f"Could not parse line: {line[:100]}")
                        continue

        # Diagnostic logging
        logger.info(f"Parse diagnostics: {lines_processed} lines processed, "
                   f"{len(alleles_found)} alleles found, "
                   f"{data_lines_attempted} data lines attempted, "
                   f"{len(all_results)} results parsed")

        if alleles_found:
            logger.info(f"Alleles detected: {', '.join(alleles_found[:5])}..." if len(alleles_found) > 5 else f"Alleles detected: {', '.join(alleles_found)}")

        if not all_results:
            logger.warning(f"No valid results found in {results_file}")
            logger.warning(f"Diagnostics: {lines_processed} total lines, {len(alleles_found)} alleles found, {data_lines_attempted} data lines attempted")
            if data_lines_attempted == 0:
                logger.warning("No data lines matched the expected format (7+ columns, first column digit)")
                # Sample first few non-empty lines for debugging
                with open(results_file) as f:
                    sample_lines = [line.strip() for line in f if line.strip()][:10]
                    logger.warning(f"First few lines of file for debugging:")
                    for i, line in enumerate(sample_lines, 1):
                        logger.warning(f"  Line {i}: {line[:200]}")
            return pd.DataFrame(columns=['protein_id', 'peptide_sequence', 'position',
                                        'peptide_length', 'hla_class', 'hla_allele',
                                        'score', 'rank'])

        df = pd.DataFrame(all_results)
        logger.info(f"Parsed {len(df)} Class {hla_class} predictions")
        return df


class HLAEpitopePipeline:
    """Main pipeline orchestrating the full analysis."""

    def __init__(self, config: AnalysisConfig):
        self.config = config
        self.output_dir = config.output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Initialize components
        self.downloader = ProteomeDownloader(config.email)
        self.signalp = SignalPFilter(config.signalp_model_dir)
        self.peptide_gen = PeptideGenerator()

    def run_class_i_analysis(self, proteome_file: Path,
                            alleles: List[str]) -> pd.DataFrame:
        """Run HLA Class I analysis on all proteins."""
        logger.info("=== HLA Class I Analysis (all proteins) ===")
        logger.info(f"Peptide lengths: {self.config.class_i_peptide_lengths}")
        logger.info(f"Step size: {self.config.class_i_step_size}")

        # Generate peptides
        peptides = self.peptide_gen.generate_variable_length_peptides(
            proteome_file,
            self.config.class_i_peptide_lengths,
            self.config.class_i_step_size
        )

        # Run predictions
        class_i_dir = self.output_dir / "class_i"
        predictor = NetMHCpanPredictor(self.config.netmhcpan_path, class_i_dir, self.config.temp_dir)
        results_file = predictor.run_predictions(
            peptides, alleles, self.config.n_parallel_jobs
        )

        # Parse results
        df = ResultParser.parse_netmhc_results(results_file, peptides, 'I')
        return df

    def run_class_ii_analysis(self, proteome_file: Path,
                             alleles: List[str]) -> pd.DataFrame:
        """Run HLA Class II analysis on secreted proteins only."""
        logger.info("=== HLA Class II Analysis (secreted proteins) ===")
        logger.info(f"Peptide length: {self.config.class_ii_peptide_length}")
        logger.info(f"Step size: {self.config.class_ii_step_size}")

        # Filter for secreted proteins
        signalp_dir = self.output_dir / "signalp"
        secreted_ids = self.signalp.run(proteome_file, signalp_dir, organism_type="other")

        # Check if any secreted proteins were found
        if not secreted_ids:
            logger.warning("No secreted proteins found - skipping Class II analysis")
            # Return empty DataFrame with correct structure
            return pd.DataFrame(columns=['protein_id', 'peptide_sequence', 'position',
                                        'peptide_length', 'hla_class', 'hla_allele',
                                        'score', 'rank'])

        secreted_file = self.output_dir / f"{self.config.organism_name}_secreted.fasta"
        self.signalp.filter_fasta(proteome_file, secreted_file, secreted_ids)

        # Generate peptides
        peptides = self.peptide_gen.generate_peptides(
            secreted_file,
            self.config.class_ii_peptide_length,
            self.config.class_ii_step_size
        )

        # Check if any peptides were generated
        if not peptides:
            logger.warning("No peptides generated from secreted proteins - skipping Class II analysis")
            return pd.DataFrame(columns=['protein_id', 'peptide_sequence', 'position',
                                        'peptide_length', 'hla_class', 'hla_allele',
                                        'score', 'rank'])

        # Run predictions
        class_ii_dir = self.output_dir / "class_ii"
        predictor = NetMHCIIpanPredictor(self.config.netmhciipan_path, class_ii_dir, self.config.temp_dir)
        results_file = predictor.run_predictions(
            peptides, alleles, self.config.n_parallel_jobs
        )

        # Parse results
        df = ResultParser.parse_netmhc_results(results_file, peptides, 'II')
        return df

    def run_full_analysis(self, class_i_alleles: List[str],
                         class_ii_alleles: List[str],
                         proteome_file: Optional[Path] = None) -> Path:
        """
        Run complete analysis pipeline.

        Args:
            class_i_alleles: List of HLA Class I alleles
            class_ii_alleles: List of HLA Class II alleles
            proteome_file: Optional pre-existing proteome FASTA to analyse. When
                provided (e.g. from a user-supplied FASTA file), the NCBI
                download step is skipped. When None, the proteome is downloaded
                from NCBI using the configured organism query.

        Returns:
            Path to combined results CSV
        """
        logger.info("=== Starting HLA Epitope Prediction Pipeline ===")
        logger.info(f"Organism: {self.config.organism_query}")
        logger.info(f"Output directory: {self.output_dir}")
        logger.info("")
        logger.info("Configuration:")
        logger.info(f"  Class I: {len(class_i_alleles)} alleles, lengths={self.config.class_i_peptide_lengths}, step={self.config.class_i_step_size}")
        logger.info(f"  Class II: {len(class_ii_alleles)} alleles, length={self.config.class_ii_peptide_length}, step={self.config.class_ii_step_size}")
        logger.info(f"  Parallel jobs: {self.config.n_parallel_jobs}")
        logger.info("")

        # Obtain proteome: use a supplied FASTA if given, otherwise download.
        if proteome_file is not None:
            proteome_file = Path(proteome_file)
            if not proteome_file.exists():
                raise FileNotFoundError(f"Proteome FASTA not found: {proteome_file}")
            logger.info(f"Using provided proteome FASTA: {proteome_file}")
        else:
            proteome_file = self.output_dir / f"{self.config.organism_name}_proteome.fasta"
            if not proteome_file.exists():
                self.downloader.download(
                    self.config.organism_query,
                    proteome_file,
                    self.config.max_proteins
                )
            else:
                logger.info(f"Using existing proteome: {proteome_file}")

        # Run Class I analysis
        class_i_df = self.run_class_i_analysis(proteome_file, class_i_alleles)

        # Run Class II analysis
        class_ii_df = self.run_class_ii_analysis(proteome_file, class_ii_alleles)

        # Combine results
        logger.info("=== Combining Results ===")
        combined_df = pd.concat([class_i_df, class_ii_df], ignore_index=True)

        # Save results
        results_file = self.output_dir / f"{self.config.organism_name}_hla_predictions.csv"
        combined_df.to_csv(results_file, index=False)

        logger.info("=== Analysis Complete ===")
        logger.info(f"Class I: {len(class_i_df)} predictions")
        logger.info(f"Class II: {len(class_ii_df)} predictions")
        logger.info(f"Total: {len(combined_df)} predictions")
        logger.info(f"Results: {results_file}")

        # Clean up temporary files
        self._cleanup_temp_files()

        return results_file

    def _cleanup_temp_files(self):
        """Clean up temporary files created during analysis."""
        if self.config.temp_dir and self.config.temp_dir.exists():
            try:
                # Count files before cleanup
                temp_files = list(self.config.temp_dir.glob("*"))
                if temp_files:
                    logger.info(f"Cleaning up {len(temp_files)} temporary files from {self.config.temp_dir}")
                    shutil.rmtree(self.config.temp_dir)
                    # Recreate empty temp dir for potential future use
                    self.config.temp_dir.mkdir(parents=True, exist_ok=True)
                    logger.info("Temporary files cleaned up successfully")
                else:
                    logger.info("No temporary files to clean up")
            except Exception as e:
                logger.warning(f"Failed to clean up temporary files: {e}")
                logger.warning("You may need to manually remove files from: {self.config.temp_dir}")


def load_alleles_from_file(file_path: Path) -> List[str]:
    """Load alleles from a text file (one per line)."""
    if not file_path.exists():
        raise FileNotFoundError(f"Allele file not found: {file_path}")

    with open(file_path) as f:
        alleles = [line.strip() for line in f if line.strip() and not line.startswith('#')]

    logger.info(f"Loaded {len(alleles)} alleles from {file_path}")
    return alleles


def get_default_class_i_alleles() -> List[str]:
    """Return default HLA Class I alleles."""
    return [
        'HLA-A01:01', 'HLA-A02:01', 'HLA-A03:01', 'HLA-A11:01',
        'HLA-A23:01', 'HLA-A24:02', 'HLA-A25:01', 'HLA-A29:02',
        'HLA-A30:01', 'HLA-A31:01', 'HLA-A32:01', 'HLA-A68:01',
        'HLA-B07:02', 'HLA-B08:01', 'HLA-B13:02', 'HLA-B14:01',
        'HLA-B14:02', 'HLA-B15:01', 'HLA-B18:01', 'HLA-B27:05',
        'HLA-B37:01', 'HLA-B38:01', 'HLA-B40:01', 'HLA-B44:02',
        'HLA-B44:03', 'HLA-B49:01', 'HLA-B50:01', 'HLA-B51:01',
        'HLA-B55:01', 'HLA-B57:01', 'HLA-C01:02', 'HLA-C02:02',
        'HLA-C03:03', 'HLA-C03:04', 'HLA-C04:01', 'HLA-C05:01',
        'HLA-C06:02', 'HLA-C07:01', 'HLA-C07:02', 'HLA-C07:04',
        'HLA-C08:02', 'HLA-C12:03', 'HLA-C14:02', 'HLA-C15:02',
        'HLA-C16:01'
    ]


def get_default_class_ii_alleles() -> List[str]:
    """Return default HLA Class II alleles."""
    dr_alleles = [
        'DRB5_9901', 'DRB1_1301', 'DRB1_0701', 'DRB1_0801',
        'DRB1_0301', 'DRB1_0101', 'DRB1_1601', 'DRB5_0202',
        'DRB1_1401', 'DRB1_1501', 'DRB5_0101', 'DRB1_0102',
        'DRB1_1201', 'DRB1_1302', 'DRB3_0301', 'DRB1_1303',
        'DRB3_0101', 'DRB3_9901', 'DRB4_0101', 'DRB4_0103',
        'DRB4_9901'
    ]

    # DP/DQ combinations
    dp_dq_combinations = [
        'HLA-DPA10103-DPB10101', 'HLA-DPA10201-DPB10101', 'HLA-DPA10202-DPB10101',
        'HLA-DQA10101-DPB10101', 'HLA-DQA10102-DPB10101', 'HLA-DQA10103-DPB10101',
        'HLA-DQA10201-DPB10101', 'HLA-DQA10301-DPB10101', 'HLA-DQA10401-DPB10101',
        'HLA-DQA10501-DPB10101', 'HLA-DPA10103-DPB10401', 'HLA-DPA10201-DPB10401',
        'HLA-DPA10202-DPB10401', 'HLA-DQA10101-DPB10401', 'HLA-DQA10102-DPB10401',
        'HLA-DQA10103-DPB10401', 'HLA-DQA10201-DPB10401', 'HLA-DQA10301-DPB10401',
        'HLA-DQA10401-DPB10401', 'HLA-DQA10501-DPB10401', 'HLA-DPA10103-DPB10402',
        'HLA-DPA10201-DPB10402', 'HLA-DPA10202-DPB10402', 'HLA-DQA10101-DPB10402',
        'HLA-DQA10102-DPB10402', 'HLA-DQA10103-DPB10402', 'HLA-DQA10201-DPB10402',
        'HLA-DQA10301-DPB10402', 'HLA-DQA10401-DPB10402', 'HLA-DQA10501-DPB10402',
        'HLA-DPA10103-DPB10501', 'HLA-DPA10201-DPB10501', 'HLA-DPA10202-DPB10501',
        'HLA-DQA10101-DPB10501', 'HLA-DQA10102-DPB10501', 'HLA-DQA10103-DPB10501',
        'HLA-DQA10201-DPB10501', 'HLA-DQA10301-DPB10501', 'HLA-DQA10401-DPB10501',
        'HLA-DQA10501-DPB10501', 'HLA-DPA10103-DPB11001', 'HLA-DPA10201-DPB11001',
        'HLA-DPA10202-DPB11001', 'HLA-DQA10101-DPB11001', 'HLA-DQA10102-DPB11001',
        'HLA-DQA10103-DPB11001', 'HLA-DQA10201-DPB11001', 'HLA-DQA10301-DPB11001',
        'HLA-DQA10401-DPB11001', 'HLA-DQA10501-DPB11001', 'HLA-DPA10103-DPB11101',
        'HLA-DPA10201-DPB11101', 'HLA-DPA10202-DPB11101', 'HLA-DQA10101-DPB11101',
        'HLA-DQA10102-DPB11101', 'HLA-DQA10103-DPB11101', 'HLA-DQA10201-DPB11101',
        'HLA-DQA10301-DPB11101', 'HLA-DQA10401-DPB11101', 'HLA-DQA10501-DPB11101',
        'HLA-DPA10103-DPB11301', 'HLA-DPA10201-DPB11301', 'HLA-DPA10202-DPB11301',
        'HLA-DQA10101-DPB11301', 'HLA-DQA10102-DPB11301', 'HLA-DQA10103-DPB11301',
        'HLA-DQA10201-DPB11301', 'HLA-DQA10301-DPB11301', 'HLA-DQA10401-DPB11301',
        'HLA-DQA10501-DPB11301', 'HLA-DPA10103-DPB11701', 'HLA-DPA10201-DPB11701',
        'HLA-DPA10202-DPB11701', 'HLA-DQA10101-DPB11701', 'HLA-DQA10102-DPB11701',
        'HLA-DQA10103-DPB11701', 'HLA-DQA10201-DPB11701', 'HLA-DQA10301-DPB11701',
        'HLA-DQA10401-DPB11701', 'HLA-DQA10501-DPB11701', 'HLA-DPA10103-DQB10501',
        'HLA-DPA10201-DQB10501', 'HLA-DPA10202-DQB10501', 'HLA-DQA10101-DQB10501',
        'HLA-DQA10102-DQB10501', 'HLA-DQA10103-DQB10501', 'HLA-DQA10201-DQB10501',
        'HLA-DQA10301-DQB10501', 'HLA-DQA10401-DQB10501', 'HLA-DQA10501-DQB10501',
        'HLA-DPA10103-DQB10603', 'HLA-DPA10201-DQB10603', 'HLA-DPA10202-DQB10603',
        'HLA-DQA10101-DQB10603', 'HLA-DQA10102-DQB10603', 'HLA-DQA10103-DQB10603',
        'HLA-DQA10201-DQB10603', 'HLA-DQA10301-DQB10603', 'HLA-DQA10401-DQB10603',
        'HLA-DQA10501-DQB10603', 'HLA-DPA10103-DQB10202', 'HLA-DPA10201-DQB10202',
        'HLA-DPA10202-DQB10202', 'HLA-DQA10101-DQB10202', 'HLA-DQA10102-DQB10202',
        'HLA-DQA10103-DQB10202', 'HLA-DQA10201-DQB10202', 'HLA-DQA10301-DQB10202',
        'HLA-DQA10401-DQB10202', 'HLA-DQA10501-DQB10202', 'HLA-DPA10103-DQB10402',
        'HLA-DPA10201-DQB10402', 'HLA-DPA10202-DQB10402', 'HLA-DQA10101-DQB10402',
        'HLA-DQA10102-DQB10402', 'HLA-DQA10103-DQB10402', 'HLA-DQA10201-DQB10402',
        'HLA-DQA10301-DQB10402', 'HLA-DQA10401-DQB10402', 'HLA-DQA10501-DQB10402',
        'HLA-DPA10103-DQB10201', 'HLA-DPA10201-DQB10201', 'HLA-DPA10202-DQB10201',
        'HLA-DQA10101-DQB10201', 'HLA-DQA10102-DQB10201', 'HLA-DQA10103-DQB10201',
        'HLA-DQA10201-DQB10201', 'HLA-DQA10301-DQB10201', 'HLA-DQA10401-DQB10201',
        'HLA-DQA10501-DQB10201', 'HLA-DPA10103-DQB10301', 'HLA-DPA10201-DQB10301',
        'HLA-DPA10202-DQB10301', 'HLA-DQA10101-DQB10301', 'HLA-DQA10102-DQB10301',
        'HLA-DQA10103-DQB10301', 'HLA-DQA10201-DQB10301', 'HLA-DQA10301-DQB10301',
        'HLA-DQA10401-DQB10301', 'HLA-DQA10501-DQB10301', 'HLA-DPA10103-DQB10302',
        'HLA-DPA10201-DQB10302', 'HLA-DPA10202-DQB10302', 'HLA-DQA10101-DQB10302',
        'HLA-DQA10102-DQB10302', 'HLA-DQA10103-DQB10302', 'HLA-DQA10201-DQB10302',
        'HLA-DQA10301-DQB10302', 'HLA-DQA10401-DQB10302', 'HLA-DQA10501-DQB10302',
        'HLA-DPA10103-DQB10303', 'HLA-DPA10201-DQB10303', 'HLA-DPA10202-DQB10303',
        'HLA-DQA10101-DQB10303', 'HLA-DQA10102-DQB10303', 'HLA-DQA10103-DQB10303',
        'HLA-DQA10201-DQB10303', 'HLA-DQA10301-DQB10303', 'HLA-DQA10401-DQB10303',
        'HLA-DQA10501-DQB10303', 'HLA-DPA10103-DQB10502', 'HLA-DPA10201-DQB10502',
        'HLA-DPA10202-DQB10502', 'HLA-DQA10101-DQB10502', 'HLA-DQA10102-DQB10502',
        'HLA-DQA10103-DQB10502', 'HLA-DQA10201-DQB10502', 'HLA-DQA10301-DQB10502',
        'HLA-DQA10401-DQB10502', 'HLA-DQA10501-DQB10502', 'HLA-DPA10103-DQB10503',
        'HLA-DPA10201-DQB10503', 'HLA-DPA10202-DQB10503', 'HLA-DQA10101-DQB10503',
        'HLA-DQA10102-DQB10503', 'HLA-DQA10103-DQB10503', 'HLA-DQA10201-DQB10503',
        'HLA-DQA10301-DQB10503', 'HLA-DQA10401-DQB10503', 'HLA-DQA10501-DQB10503',
        'HLA-DPA10103-DQB10602', 'HLA-DPA10201-DQB10602', 'HLA-DPA10202-DQB10602',
        'HLA-DQA10101-DQB10602', 'HLA-DQA10102-DQB10602', 'HLA-DQA10103-DQB10602',
        'HLA-DQA10201-DQB10602', 'HLA-DQA10301-DQB10602', 'HLA-DQA10401-DQB10602',
        'HLA-DQA10501-DQB10602', 'HLA-DPA10103-DQB10604', 'HLA-DPA10201-DQB10604',
        'HLA-DPA10202-DQB10604', 'HLA-DQA10101-DQB10604', 'HLA-DQA10102-DQB10604',
        'HLA-DQA10103-DQB10604', 'HLA-DQA10201-DQB10604', 'HLA-DQA10301-DQB10604',
        'HLA-DQA10401-DQB10604', 'HLA-DQA10501-DQB10604', 'HLA-DPA10103-DQB10609',
        'HLA-DPA10201-DQB10609', 'HLA-DPA10202-DQB10609', 'HLA-DQA10101-DQB10609',
        'HLA-DQA10102-DQB10609', 'HLA-DQA10103-DQB10609', 'HLA-DQA10201-DQB10609',
        'HLA-DQA10301-DQB10609', 'HLA-DQA10401-DQB10609', 'HLA-DQA10501-DQB10609'
    ]

    return dr_alleles + dp_dq_combinations


def main():
    parser = argparse.ArgumentParser(
        description="HLA Epitope Prediction Pipeline with incremental processing",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with conservative defaults (9-mers, step=3, ~10x faster)
  python hla_epitope_predictor.py --email user@example.com --organism "Variola virus[Organism]"

  # Run directly from a local FASTA file (no NCBI download)
  python hla_epitope_predictor.py --fasta proteome.fasta --output-dir ./my_analysis

  # FASTA with multiple organisms/strains: list what was detected, then run all
  python hla_epitope_predictor.py --fasta multi_organism.fasta --list-organisms
  python hla_epitope_predictor.py --fasta multi_organism.fasta --output-dir ./results

  # FASTA input but only analyse a chosen organism/strain or sequence
  python hla_epitope_predictor.py --fasta multi_organism.fasta --select "Variola virus"
  python hla_epitope_predictor.py --fasta multi_organism.fasta --select NP_042094.1

  # Run with custom allele files
  python hla_epitope_predictor.py --email user@example.com --organism "HIV[Organism]" \\
      --class-i-alleles class_i_alleles.txt --class-ii-alleles class_ii_alleles.txt

  # Add new alleles (will skip already-completed ones)
  python hla_epitope_predictor.py --email user@example.com --organism "Variola virus[Organism]" \\
      --class-i-alleles new_alleles.txt --output-dir ./existing_analysis

  # Comprehensive mode (all lengths, step=1) for thorough analysis
  python hla_epitope_predictor.py --email user@example.com --organism "HIV[Organism]" \\
      --class-i-lengths 8 9 10 11 --class-i-step 1 --class-ii-step 1

  # Ultra-fast screening mode (non-overlapping peptides)
  python hla_epitope_predictor.py --email user@example.com --organism "Large organism[Organism]" \\
      --class-i-step 9 --class-ii-step 15
        """
    )

    # Input source: either an NCBI query (--organism) or a local FASTA file
    # (--fasta). Exactly one of these must be provided.
    parser.add_argument("--email",
                       help="Email for NCBI Entrez API (required with --organism)")
    parser.add_argument("--organism",
                       help="NCBI organism query (e.g., 'Variola virus[Organism]'). "
                            "Mutually exclusive with --fasta.")
    parser.add_argument("--fasta", type=Path,
                       help="Local FASTA file to analyse instead of downloading from "
                            "NCBI. May contain multiple organisms/strains (detected "
                            "from header lines); by default every organism is analysed.")
    parser.add_argument("--select", nargs='+', metavar="NAME",
                       help="With --fasta, only analyse the given organism/strain "
                            "name(s) or sequence ID(s). May be given multiple values.")
    parser.add_argument("--combine", action="store_true",
                       help="With --fasta, treat the whole file as a single proteome "
                            "instead of splitting by organism/strain.")
    parser.add_argument("--list-organisms", action="store_true",
                       help="With --fasta, list the organisms/strains detected in the "
                            "file (and sequence counts) then exit without running.")

    # Optional arguments
    parser.add_argument("--output-dir", default="./hla_analysis",
                       help="Output directory (default: ./hla_analysis)")
    parser.add_argument("--name",
                       help="Name for output files (default: auto-generated from organism)")

    # Allele specification
    parser.add_argument("--class-i-alleles", type=Path,
                       help="File with Class I alleles (one per line). Uses defaults if not provided.")
    parser.add_argument("--class-ii-alleles", type=Path,
                       help="File with Class II alleles (one per line). Uses defaults if not provided.")

    # Tool paths
    parser.add_argument("--netmhcpan-path",
                       default="/home/jdn321/software/netMHCpan-4.2/netMHCpan",
                       help="Path to netMHCpan executable")
    parser.add_argument("--netmhciipan-path",
                       default="/home/jdn321/software/netMHCIIpan-4.3/netMHCIIpan",
                       help="Path to netMHCIIpan executable")
    parser.add_argument("--signalp-model-dir",
                       default="/home/jdn321/software/signalp6_fast/signalp-6-package/models/",
                       help="Path to SignalP6 models directory")

    # Temporary file management
    parser.add_argument("--temp-dir", type=Path,
                       help="Directory for temporary files (default: output_dir/temp). "
                            "netMHC tools create temporary files that can fill /tmp on servers.")

    # Processing parameters
    parser.add_argument("--n-jobs", type=int, default=20,
                       help="Number of parallel jobs (default: 20)")
    parser.add_argument("--max-proteins", type=int, default=15000,
                       help="Maximum proteins to download (default: 15000)")

    # Peptide generation parameters
    parser.add_argument("--class-i-lengths", type=int, nargs='+', default=[9],
                       help="Class I peptide lengths (default: 9). Example: --class-i-lengths 9 10 11")
    parser.add_argument("--class-i-step", type=int, default=3,
                       help="Step size for Class I peptides (default: 3 = ~10x speedup)")
    parser.add_argument("--class-ii-length", type=int, default=15,
                       help="Class II peptide length (default: 15)")
    parser.add_argument("--class-ii-step", type=int, default=3,
                       help="Step size for Class II peptides (default: 3 = ~3x speedup)")

    args = parser.parse_args()

    # Validate input source: exactly one of --organism / --fasta.
    if bool(args.organism) == bool(args.fasta):
        parser.error("Provide exactly one input source: --organism (NCBI) or --fasta (local file).")
    if args.organism and not args.email:
        parser.error("--email is required when using --organism (NCBI download).")
    if not args.fasta and (args.select or args.combine or args.list_organisms):
        parser.error("--select/--combine/--list-organisms can only be used with --fasta.")

    # Build the helper that constructs a config + allele lists for one analysis.
    def build_config(organism_query: str, output_dir: Path, organism_name: str) -> AnalysisConfig:
        return AnalysisConfig(
            email=args.email or "",
            organism_query=organism_query,
            output_dir=output_dir,
            organism_name=organism_name,
            netmhcpan_path=args.netmhcpan_path,
            netmhciipan_path=args.netmhciipan_path,
            signalp_model_dir=args.signalp_model_dir,
            n_parallel_jobs=args.n_jobs,
            class_i_peptide_lengths=args.class_i_lengths,
            class_i_step_size=args.class_i_step,
            class_ii_peptide_length=args.class_ii_length,
            class_ii_step_size=args.class_ii_step,
            max_proteins=args.max_proteins,
            temp_dir=args.temp_dir
        )

    # Load or use default alleles (shared across all analyses in this run).
    if args.class_i_alleles:
        class_i_alleles = load_alleles_from_file(args.class_i_alleles)
    else:
        class_i_alleles = get_default_class_i_alleles()
        logger.info(f"Using {len(class_i_alleles)} default Class I alleles")

    if args.class_ii_alleles:
        class_ii_alleles = load_alleles_from_file(args.class_ii_alleles)
    else:
        class_ii_alleles = get_default_class_ii_alleles()
        logger.info(f"Using {len(class_ii_alleles)} default Class II alleles")

    try:
        if args.fasta:
            # --- Local FASTA input (one analysis per organism/strain) ---
            groups = FastaInputHandler.group_by_organism(args.fasta)

            if args.list_organisms:
                print(f"\nDetected {len(groups)} organism/strain group(s) in {args.fasta}:")
                for organism, records in groups.items():
                    print(f"  - {organism}: {len(records)} sequence(s)")
                print("")
                return

            if args.combine:
                # Treat the entire file as a single proteome.
                combined_name = args.name or sanitize_name(Path(args.fasta).stem)
                groups = {combined_name: [r for recs in groups.values() for r in recs]}
            else:
                groups = FastaInputHandler.select_groups(groups, args.select)

            base_output_dir = Path(args.output_dir)
            multiple = len(groups) > 1
            logger.info(f"Running analysis for {len(groups)} organism/strain group(s)")

            result_files = []
            for organism, records in groups.items():
                organism_name = args.name if (args.name and not multiple) else sanitize_name(organism)
                # Nest per-organism outputs only when there is more than one group.
                output_dir = base_output_dir / organism_name if multiple else base_output_dir
                output_dir.mkdir(parents=True, exist_ok=True)

                # Write this group's sequences to its own proteome FASTA.
                proteome_file = output_dir / f"{organism_name}_proteome.fasta"
                FastaInputHandler.write_group(records, proteome_file)

                config = build_config(
                    organism_query=f"FASTA input: {organism}",
                    output_dir=output_dir,
                    organism_name=organism_name,
                )
                pipeline = HLAEpitopePipeline(config)
                result_file = pipeline.run_full_analysis(
                    class_i_alleles, class_ii_alleles, proteome_file=proteome_file
                )
                result_files.append(result_file)

            print(f"\n{'='*80}")
            print("Analysis completed successfully!")
            for result_file in result_files:
                print(f"Results saved to: {result_file}")
            print(f"{'='*80}\n")

        else:
            # --- NCBI download input (original behaviour) ---
            organism_name = args.name or sanitize_name(args.organism)
            config = build_config(
                organism_query=args.organism,
                output_dir=Path(args.output_dir),
                organism_name=organism_name,
            )
            pipeline = HLAEpitopePipeline(config)
            result_file = pipeline.run_full_analysis(class_i_alleles, class_ii_alleles)

            print(f"\n{'='*80}")
            print("Analysis completed successfully!")
            print(f"Results saved to: {result_file}")
            print(f"{'='*80}\n")

    except Exception as e:
        logger.error(f"Analysis failed: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
