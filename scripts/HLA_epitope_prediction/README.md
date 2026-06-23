# HLA Epitope Prediction Pipeline

Comprehensive pipeline for predicting HLA Class I and II epitopes from pathogen proteomes with support for incremental processing and result caching.

## Features

- **HLA Class I & II prediction** using netMHCpan and netMHCIIpan
- **RefSeq reference sequences**: Uses NCBI-curated reference genomes for consistency
- **Incremental processing**: Add new alleles without re-running completed ones
- **Progress tracking**: Resume from failures automatically
- **Configurable alleles**: Load from external files
- **Memory-efficient**: Adjustable parallelization
- **Modular design**: Clean, maintainable code

## NCBI Query Format

This pipeline uses **RefSeq[Filter]** in all organism queries to ensure:

- ✅ **High-quality sequences**: NCBI-curated reference genomes only
- ✅ **Non-redundant**: Eliminates duplicate sequences from the same organism
- ✅ **Consistent annotation**: Standardized gene/protein nomenclature
- ✅ **Reproducibility**: Same reference strain across analyses
- ✅ **No strain name issues**: Avoids "strain not found" errors

**Example query format:**
```bash
"Variola virus[Organism] AND RefSeq[Filter]"
"Mycobacterium tuberculosis[Organism] AND RefSeq[Filter]"
```

This automatically fetches the official RefSeq reference genome without needing to specify exact strain names.

## Input Sources

The pipeline accepts **either** of two mutually exclusive input sources:

1. **NCBI download** (`--organism` + `--email`): downloads the proteome from
   NCBI as described above.
2. **Local FASTA file** (`--fasta`): analyses sequences from a FASTA file you
   already have, with no network access required.

### Using a Local FASTA File

A single FASTA file may contain sequences for **multiple organisms or strains**.
The organism/strain is inferred from each header line, supporting both common
conventions:

- **NCBI style** — organism (and strain) in trailing square brackets:
  `>NP_042094.1 hypothetical protein [Variola virus]`
- **UniProt style** — organism after the `OS=` token:
  `>sp|P12345|NAME ... OS=Homo sapiens OX=9606 GN=ABC PE=1 SV=2`

Sequences whose organism cannot be determined are grouped under
`unknown_organism` rather than dropped.

```bash
# Inspect which organisms/strains are in the file (no analysis is run)
python3 hla_epitope_predictor.py --fasta proteome.fasta --list-organisms

# Analyse EVERY organism/strain automatically (one analysis each).
# When more than one organism is present, results are written to a
# per-organism subdirectory of --output-dir.
python3 hla_epitope_predictor.py --fasta proteome.fasta --output-dir ./results

# Analyse only a chosen organism/strain (name substring, case-insensitive)
python3 hla_epitope_predictor.py --fasta proteome.fasta --select "Variola virus"

# Analyse only a chosen sequence by its FASTA ID
python3 hla_epitope_predictor.py --fasta proteome.fasta --select NP_042094.1

# Select multiple organisms and/or sequence IDs at once
python3 hla_epitope_predictor.py --fasta proteome.fasta --select "Measles" NP_042094.1

# Treat the whole file as a single combined proteome (ignore organism splitting)
python3 hla_epitope_predictor.py --fasta proteome.fasta --combine --name my_dataset
```

Parallel execution (`--n-jobs`, GNU `parallel`) works exactly as with NCBI
input; each organism is analysed in turn with its alleles processed in parallel.

## Quick Start

### 1. Install Dependencies

```bash
# Python packages
pip install biopython pandas

# External tools (adjust paths in scripts):
# - netMHCpan-4.2
# - netMHCIIpan-4.3
# - SignalP 6.0
# - GNU parallel
```

### 2. Basic Usage

```bash
# Make scripts executable
chmod +x hla_epitope_predictor.py analyse_predictions.py

# Run analysis for one organism
python3 hla_epitope_predictor.py \
    --email your@email.com \
    --organism "Variola virus[Organism] AND RefSeq[Filter]" \
    --output-dir ./smallpox_analysis
```

### 3. Batch Analysis

```bash
# Edit run_pathogen_analysis.sh to set your email and paths
nano run_pathogen_analysis.sh

# Run for all organisms
chmod +x run_pathogen_analysis.sh
./run_pathogen_analysis.sh
```

### 4. Analyse Results

```bash
python3 analyse_predictions.py --batch \
    smallpox_analysis HIV_analysis plague_analysis \
    --output-dir ./analysed_results
```

## File Structure

```
scripts/HLA_epitope_prediction/
├── hla_epitope_predictor.py       # Main prediction pipeline
├── analyse_predictions.py          # Results analysis
├── run_pathogen_analysis.sh        # Batch runner for 29 pathogens
├── alleles_class_i_default.txt     # Default Class I alleles (45)
├── alleles_class_ii_default.txt    # Default Class II alleles (189)
├── README.md                        # This file
├── PERFORMANCE_GUIDE.md             # Performance tuning guide
└── requirements.txt                 # Python dependencies
```

## Usage Examples

### Example 1: Run with Default Alleles

```bash
python3 hla_epitope_predictor.py \
    --email user@example.com \
    --organism "HIV[Organism] AND RefSeq[Filter]" \
    --output-dir ./HIV_analysis \
    --n-jobs 40
```

### Example 2: Custom Allele Files

```bash
# Create custom allele list
cat > my_alleles.txt <<EOF
HLA-A02:01
HLA-A24:02
HLA-B07:02
EOF

python3 hla_epitope_predictor.py \
    --email user@example.com \
    --organism "Variola virus[Organism]" \
    --class-i-alleles my_alleles.txt \
    --output-dir ./test_analysis
```

### Example 3: Adding New Alleles (Incremental)

```bash
# Initial run
python3 hla_epitope_predictor.py \
    --email user@example.com \
    --organism "Variola virus[Organism]" \
    --output-dir ./smallpox_analysis

# Add new alleles to file
echo "HLA-A24:03" >> alleles_class_i_default.txt
echo "HLA-B35:01" >> alleles_class_i_default.txt

# Re-run: only new alleles are processed
python3 hla_epitope_predictor.py \
    --email user@example.com \
    --organism "Variola virus[Organism]" \
    --output-dir ./smallpox_analysis  # Same directory!
```

### Example 4: Analyse Results

The analysis script now reports detailed per-allele statistics including:
- Number of unique proteins with strong binders
- Percentage of strong binders (rank < threshold) for each HLA allele
- Percentage of weak binders (rank between strong and weak thresholds) for each HLA allele

```bash
# Single organism
python3 analyse_predictions.py \
    --input smallpox_analysis/smallpox_hla_predictions.csv \
    --class-i-strong 0.5 \
    --class-i-weak 2.0 \
    --class-ii-strong 1.0 \
    --class-ii-weak 5.0

# Multiple organisms
python3 analyse_predictions.py --batch \
    smallpox_analysis \
    HIV_analysis \
    plague_analysis \
    --output-dir ./analysed_results
```

Output files:
- `{organism}_class_i_allele_statistics.tsv` - Per-allele statistics for Class I
- `{organism}_class_ii_allele_statistics.tsv` - Per-allele statistics for Class II
- `{organism}_summary.tsv` - Overall summary statistics

## Output Files

After running the pipeline:

```
output_dir/
├── organism_proteome.fasta              # Downloaded proteome
├── organism_secreted.fasta              # Filtered for Class II
├── signalp/                             # SignalP results
├── class_i/                             # Class I analysis
│   ├── peptides.txt
│   ├── completed_alleles.json           # Progress tracking!
│   ├── allele_results/                  # Per-allele results
│   │   ├── HLA-A01_01.txt
│   │   ├── HLA-A02_01.txt
│   │   └── ...
│   └── combined_results.txt
├── class_ii/                            # Class II analysis
│   ├── peptides.txt
│   ├── completed_alleles.json
│   ├── allele_results/
│   └── combined_results.txt
└── organism_hla_predictions.csv         # Final results!
```

## Result Format

The `*_hla_predictions.csv` contains:

| Column | Description |
|--------|-------------|
| `protein_id` | Source protein identifier from NCBI |
| `peptide_sequence` | Peptide sequence (8-11 AA for Class I, 15 AA for Class II) |
| `position` | Position in source protein |
| `peptide_length` | Length of peptide |
| `hla_class` | 'I' or 'II' |
| `hla_allele` | HLA allele name |
| `score` | Binding score |
| `rank` | Percentile rank (lower = stronger binding) |

## Interpreting Results

### Binding Thresholds

Common thresholds for `rank` (%Rank):

| Category | Class I | Class II |
|----------|---------|----------|
| Very strong | < 0.05 | < 0.1 |
| Strong | < 0.5 | < 1.0 |
| Weak | < 2.0 | < 5.0 |

### Analysis Workflow

```bash
# 1. Filter strong binders
python3 -c "
import pandas as pd
df = pd.read_csv('organism_hla_predictions.csv')
strong = df[(df['hla_class']=='I') & (df['rank']<0.5)]
print(f'Strong Class I binders: {len(strong)}')
"

# 2. Find proteins with broad HLA coverage
python3 -c "
import pandas as pd
df = pd.read_csv('organism_hla_predictions.csv')
strong_i = df[(df['hla_class']=='I') & (df['rank']<0.5)]
coverage = strong_i.groupby('protein_id')['hla_allele'].nunique()
top_proteins = coverage.sort_values(ascending=False).head(10)
print('Top proteins by HLA coverage:')
print(top_proteins)
"
```

## Memory Management

Adjust `--n-jobs` based on available memory and proteome size:

| Proteome Size | Proteins | Recommended --n-jobs | Memory |
|---------------|----------|---------------------|--------|
| Small         | < 500    | 40-60               | ~40GB  |
| Medium        | 500-2000 | 30-40               | ~30GB  |
| Large         | > 2000   | 20-30               | ~20GB  |

**Rule of thumb**: Each parallel job needs ~1GB RAM.

## Troubleshooting

### Issue: Script fails partway through

**Solution**: Just re-run the same command. Completed alleles are automatically skipped.

```bash
# Check progress
cat output_dir/class_i/completed_alleles.json
cat output_dir/class_ii/completed_alleles.json
```

### Issue: Out of memory errors

**Solution**: Reduce `--n-jobs`:

```bash
python3 hla_epitope_predictor.py \
    --organism "Large organism[Organism]" \
    --n-jobs 15  # Reduced from default 20
```

### Issue: Too many open files

**Solution**: Increase system limits:

```bash
ulimit -n 4096
```

### Issue: Can't find netMHC tools

**Solution**: Specify paths explicitly:

```bash
python3 hla_epitope_predictor.py \
    --netmhcpan-path /path/to/netMHCpan \
    --netmhciipan-path /path/to/netMHCIIpan \
    --signalp-model-dir /path/to/signalp/models/
```

## Comparison with Previous Version

| Feature | Old Version | New Version |
|---------|-------------|-------------|
| Organism queries | Exact strain names | RefSeq filter (automatic) |
| Alleles | Hardcoded in Python | External text files |
| Adding alleles | Re-run all (hours) | Only new ones (minutes) |
| Progress tracking | None | `completed_alleles.json` |
| Resume capability | No | Yes (automatic) |
| Memory usage | ~80GB (80 jobs) | ~20-40GB (20-40 jobs) |
| Error handling | Fail entire run | Continue on errors |
| Code organization | Monolithic (1119 lines) | Modular (950 lines) |
| Documentation | Limited | Comprehensive |

## Advanced Usage

### Check Which Alleles Are Supported

```bash
# List netMHCpan supported alleles
/path/to/netMHCpan -listMHC

# List netMHCIIpan supported alleles
/path/to/netMHCIIpan -list
```

### Force Re-run Specific Allele

```bash
# Remove from completed list
python3 -c "
import json
with open('output_dir/class_i/completed_alleles.json') as f:
    completed = set(json.load(f))
completed.discard('HLA-A02:01')
with open('output_dir/class_i/completed_alleles.json', 'w') as f:
    json.dump(sorted(list(completed)), f)
"

# Re-run
python3 hla_epitope_predictor.py --organism "..." --output-dir output_dir
```

### Export for Downstream Analysis

```bash
# Export strong binders for vaccine design
python3 -c "
import pandas as pd
df = pd.read_csv('organism_hla_predictions.csv')
strong = df[(df['rank'] < 0.5)]
strong.to_csv('strong_binders.csv', index=False)
print(f'Exported {len(strong)} strong binders')
"
```

## Citation

If you use this pipeline, please cite the underlying tools:

- **netMHCpan**: Reynisson et al. (2020) Nucleic Acids Res.
- **netMHCIIpan**: Reynisson et al. (2020) Nucleic Acids Res.
- **SignalP 6.0**: Teufel et al. (2022) Nat Biotechnol.

## Support

For issues or questions:
1. Check `PERFORMANCE_GUIDE.md` for performance tuning
2. Review log output for error messages
3. Verify tool paths are correct
4. Check `completed_alleles.json` for progress

## License

This pipeline is provided as-is for academic research.
