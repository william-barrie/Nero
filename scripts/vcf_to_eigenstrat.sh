#!/usr/bin/env bash
# =============================================================================
# vcf_to_eigenstrat.sh
# Direct VCF -> EIGENSTRAT conversion, single streaming pass.
#
# Replaces the vcftools/convertf-based convertVCFtoEigenstrat.sh for systems
# where the open-fd hard limit (`ulimit -Hn`) is below the number of samples.
# vcftools' PED writer opens one temp file per sample; bcftools+awk do not.
#
# Outputs (in CWD), drop-in compatible with the v3/v4 R analysis script:
#   <prefix>.geno   one row per SNP, one char per sample (count of REF alleles:
#                   0=hom alt, 1=het, 2=hom ref, 9=missing)
#   <prefix>.snp    snp_id, chrom, gd (Morgans, 2 cM/Mb), pos, ref, alt
#   <prefix>.ind    sample, sex='?', pop='?'   (pop assigned later)
#
# Filters: biallelic SNPs only (bcftools -m2 -M2 -v snps), drop sites
# monomorphic among non-missing genotypes (MAC>=1, applied in awk).
#
# Atomic outputs: writes to <prefix>.{ind,snp,geno}.tmp.$$ and renames on
# success. If killed mid-run, no partial final files are left behind.
# =============================================================================
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <vcf[.gz]> [rec_cM_per_Mb=2]" >&2
    exit 1
fi

VCF="$1"
REC="${2:-2}"

[[ -s "$VCF" ]] || { echo "VCF not found or empty: $VCF" >&2; exit 1; }

PREFIX=$(basename "${VCF%.gz}")
PREFIX="${PREFIX%.vcf}"

IND_TMP="${PREFIX}.ind.tmp.$$"
SNP_TMP="${PREFIX}.snp.tmp.$$"
GENO_TMP="${PREFIX}.geno.tmp.$$"

# Always remove .tmp files on exit; final files are only created on success
trap 'rm -f "$IND_TMP" "$SNP_TMP" "$GENO_TMP"' EXIT

# 1. .ind (sample, sex='?', pop='?')
bcftools query -l "$VCF" \
    | awk 'BEGIN{OFS="\t"} {print $1, "?", "?"}' \
    > "$IND_TMP"
N_IND=$(wc -l < "$IND_TMP")
echo "[$(date +%H:%M:%S)] ${N_IND} samples" >&2

# 2. Stream filter + encode in a single pass to .snp.tmp + .geno.tmp
bcftools view -m2 -M2 -v snps -Ou "$VCF" \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
  | awk -v rec="$REC" -v snp="$SNP_TMP" -v geno="$GENO_TMP" '
    BEGIN {
        OFS = "\t"
        last_chrom = ""; last_pos = 0; cum_gd = 0
        n_kept = 0; n_mono = 0
    }
    {
        chrom = $1; pos = $2 + 0; ref = $3; alt = $4

        # Encode genotypes (count of REF alleles) and track polymorphism
        line = ""; c0 = 0; c1 = 0; c2 = 0
        for (i = 5; i <= NF; i++) {
            gt = $i
            if      (gt == "0/0" || gt == "0|0") { c = "2"; c2++ }
            else if (gt == "1/1" || gt == "1|1") { c = "0"; c0++ }
            else if (gt == "0/1" || gt == "0|1" ||
                     gt == "1/0" || gt == "1|0") { c = "1"; c1++ }
            else                                   c = "9"
            line = line c
        }

        # MAC >= 1: keep if any het OR both homozygous classes present
        if (c1 == 0 && (c0 == 0 || c2 == 0)) { n_mono++; next }

        if (chrom != last_chrom) {
            last_chrom = chrom; last_pos = pos; cum_gd = 0
        }
        cum_gd += 1e-8 * rec * (pos - last_pos)
        last_pos = pos

        printf "%s:%s\t%s\t%.8f\t%d\t%s\t%s\n", \
               chrom, pos, chrom, cum_gd, pos, ref, alt > snp
        print line > geno
        n_kept++
    }
    END {
        printf "  kept=%d  dropped_monomorphic=%d\n", n_kept, n_mono \
               > "/dev/stderr"
    }
'

# 3. Atomic promotion (only reached on success)
mv "$IND_TMP"  "${PREFIX}.ind"
mv "$SNP_TMP"  "${PREFIX}.snp"
mv "$GENO_TMP" "${PREFIX}.geno"

N_SNPS=$(wc -l < "${PREFIX}.geno")
echo "[$(date +%H:%M:%S)] ${N_SNPS} SNPs -> ${PREFIX}.{geno,snp,ind}" >&2
