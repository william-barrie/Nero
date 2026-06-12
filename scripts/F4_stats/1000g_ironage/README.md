# HLA African-admixture f4 / f4-ratio analysis (1000G + Iron Age)

Scans ancient African source populations for a localised excess of "African"
affinity in the HLA region of chromosome 6 of modern European 1000G
populations, relative to a chromosome-1 genome-wide baseline.

## Pipeline

| Step | Script | Output |
|------|--------|--------|
| 1. Build genotypes | `prep_data.sh` | `{1,6}.1000g.sea.nocommas.{geno,snp,ind,ind_}` |
| ↳ VCF→EIGENSTRAT  | `vcf_to_eigenstrat.sh` (called by step 1) | `.geno/.snp/.ind` |
| ↳ population labels | `make_ind_labels.py` (called by step 1) | `.ind_` |
| 2. f4 analysis | `f4_HLA_analysis_v5.R` | CSVs, PDFs, manifest, logs |

`prep_data.sh` is idempotent (skips steps whose output already exists; set
`FORCE=1` to rebuild) and writes atomically (`.tmp` + rename).

## Running the analysis

```bash
bash prep_data.sh                 # produces the eigenstrat files
Rscript f4_sanity_checks.R        # positive/negative controls (run this first)
Rscript f4_HLA_analysis_v5.R      # the WHG-reference scan
```

`f4_common.R` holds the shared streaming reader + f4/jackknife machinery, so
the sanity checks and the scan compute f4 identically. The final scan is
restricted to the **WHG reference** (`REFERENCES = "WHG"`; override with
`NERO_REFERENCES="WHG,Yamnaya,..."`).

### Sanity checks (`f4_sanity_checks.R`)

A panel of well-known admixture events with textbook expected directions,
computed on chr1 with the same block jackknife, written to `sanity_checks.csv`:
Steppe ancestry N>S Europe, CHG in Yamnaya, hunter-gatherer resurgence
(Late>Early Neolithic), Anatolian-farmer ancestry in S Europe, African gene
flow Iberia>Finland, plus a CEU-vs-GBR null. Each is reported PASS / CHECK
against a `|z| > 3` threshold. A failing control means the f4 set-up is not
reproducing a known result — investigate before trusting the HLA scan.

The R step is configurable via environment variables (defaults reproduce the
previous in-place behaviour):

| Variable | Default | Meaning |
|----------|---------|---------|
| `NERO_DATADIR`   | `.`     | directory holding the `.geno/.snp/.ind_` files and the metadata TSV |
| `NERO_OUTDIR`    | `.`     | directory for all CSV/PDF/log outputs |
| `NERO_CHUNK_SNPS`| `50000` | SNPs per streaming chunk — the main memory/speed knob |

Example: `NERO_OUTDIR=results NERO_CHUNK_SNPS=20000 Rscript f4_HLA_analysis_v5.R`

## Why v5 (the memory fix)

The genotype files are huge — e.g. chr1 is ~2.7M SNPs × ~15,400 samples
across ~1,760 populations. v4 read the entire `.geno` into a dense
`n_SNPs × n_samples` integer matrix (~169 GB for chr1, plus a far larger
transient from `strsplit()` on single characters), which is what produced

```
Error: cannot allocate vector of size 400.5 Gb
```

The analysis only ever uses a few dozen populations (the outgroup, the YRI
sister, 5 European targets, 8 ancestry references, and the ancient African
candidate sources). v5 therefore:

1. resolves that small population set from the `.ind_` labels **before**
   touching genotypes;
2. streams the `.geno` file in chunks of `NERO_CHUNK_SNPS` rows, parsing only
   the sample columns for those populations (via `utf8ToInt`, avoiding
   per-character string overhead);
3. collapses straight to a per-SNP, per-population allele-frequency matrix
   (`n_SNPs × n_pops`, on the order of ~1 GB), which is all the f4 jackknife
   machinery needs.

Allele frequencies are computed identically to v4, so **results are
unchanged**; only the memory footprint differs. Equivalence is exercised by a
small unit check during development (frequencies match to < 1e-12).

## Statistic and sign convention

`compute_f4(A, B; C, D) = mean_SNP (pA − pB)(pC − pD)`, where `p` is the
REF-allele frequency. Standard errors use a block (leave-one-out) jackknife
over physical-position blocks; the f4-ratio jackknife re-estimates the ratio
per deleted block.

The headline contrast is `f4(OUTGROUP, source; target, reference)` evaluated in
the HLA window vs the chr1 baseline:

```
z_diff_chr1 = (f4_HLA − f4_chr1) / sqrt(se_HLA² + se_chr1²)
```

A **negative** `z_diff_chr1` means the HLA region makes the European `target`
look *more* like the ancient African `source` (relative to `reference`) than
the genome-wide chr1 baseline does — i.e. localised excess source affinity at
HLA. Sources are ranked most-negative-first.

> Note: as with any f4 statistic, the direction of "more source-like" is only
> interpretable under the assumed population topology. Sanity-check the sign
> against a known positive/negative control before drawing conclusions.

## Outputs

- `populations_used.csv` — every population loaded, its role, and per-chr sample counts (run manifest).
- `f4_contrast_results.csv` — per (source, target, reference): f4 in HLA / chr6-rest / chr1, the contrasts, z, p, and BH-FDR.
- `f4_ratio_results.csv` — f4-ratio α in HLA vs chr1, Δα, z, p, and BH-FDR.
- `source_summary.csv` — per-source aggregates and counts of significant pairs (raw / Bonferroni / FDR).
- `chr6_f4_sliding_window.{csv,pdf}`, `chr6_alpha_sliding_window.{csv,pdf}` — best-triple scans (data + figure).
- `f4_hla_excess_heatmap_bestsource.pdf`, `source_x_target_heatmap.pdf`, `source_x_reference_heatmap.pdf`.
- `run_params.txt`, `sessionInfo.txt` — parameters and R session for reproducibility.
