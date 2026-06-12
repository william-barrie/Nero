#!/usr/bin/env bash
# =============================================================================
# prep_data.sh
# Prepare eigenstrat genotype data (chr1 + chr6) for HLA f4 analysis.
#
# Pipeline (per chromosome, idempotent):
#   2b. Merge raw ancient VCF with the matching 1000G freeze9 VCF.
#   2c. Strip stray commas inside quoted fields.
#   2d. Convert directly to EIGENSTRAT (no vcftools / convertf).
#   2e. Add population labels to the .ind file -> .ind_.
#
# Idempotency: each step is skipped if its output already exists non-empty.
# Set FORCE=1 to re-run everything regardless of cached outputs.
#
# Atomic writes: every step writes to a .tmp first and only renames on
# success, so a killed/interrupted run does not leave a partial file masked
# as completed.
#
# Assumes vcf_to_eigenstrat.sh and make_ind_labels.py are alongside.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---- Inputs -----------------------------------------------------------------
ANCIENT_DIR="/projects/lundbeck/scratch/vrb229/project_sea/data/250601_impute_sea/merge/vcf"
THOUSANDG_DIR="/projects/lundbeck/data/1000genomes_2015_nature"
META="${SCRIPT_DIR}/sampleInfo.240504_impute_neosea.v07.240627_impute_iadk_cluster.map.ascii.250601_impute_sea_sampleId.tsv"
NAME2ID="${SCRIPT_DIR}/name2id"
# Ancient reference panel (grouped labels: WHG/EHG/CHG/Farmer*/Yamnaya).
REF_POP="/datasets/ukb-AUDIT/ref_pop_ids_mapped_grouped"
# Modern 1000G population codes (CEU/GBR/TSI/IBS/FIN/YRI/...). Download the
# sample table from the IGSR data portal (https://www.internationalgenome.org)
# and place it here, or point IGSR at its location.
IGSR="${IGSR:-${SCRIPT_DIR}/igsr_samples.tsv}"

CHRS=(1 6)
OUT_DIR="${SCRIPT_DIR}"
CONVERT_SH="${SCRIPT_DIR}/vcf_to_eigenstrat.sh"
LABEL_PY="${SCRIPT_DIR}/make_ind_labels.py"

FORCE="${FORCE:-0}"

# Returns 0 (run) if FORCE=1 or any named output is missing/empty,
# returns 1 (skip) if every named output is present and non-empty.
should_run() {
    [[ "$FORCE" == "1" ]] && return 0
    for f in "$@"; do
        [[ -s "$f" ]] || return 0
    done
    return 1
}

for CHR in "${CHRS[@]}"; do
    echo
    echo "===================  CHR ${CHR}  ==================="

    ANCIENT_VCF="${ANCIENT_DIR}/${CHR}.250601_impute_sea.glimpse.vcf.gz"
    THOUSANDG_VCF="${THOUSANDG_DIR}/${CHR}.1000g.freeze9.umich.GRCh37.snps.biallelic.pass.vcf.gz"
    MERGED_VCF="${OUT_DIR}/${CHR}.1000g.sea.vcf.gz"
    NOCOMMAS_VCF="${OUT_DIR}/${CHR}.1000g.sea.nocommas.vcf.gz"
    PREFIX="${CHR}.1000g.sea.nocommas"

    # Sanity: ancient VCF must be indexed (bcftools merge requires it)
    if [[ ! -f "${ANCIENT_VCF}.tbi" && ! -f "${ANCIENT_VCF}.csi" ]]; then
        echo "ERROR: ancient VCF not indexed: ${ANCIENT_VCF}" >&2
        echo "       If you can't write to that directory, symlink it locally and" >&2
        echo "       run 'bcftools index -t' on the symlink before retrying." >&2
        exit 1
    fi

    # 2b. Merge ancient + 1000G
    if should_run "${MERGED_VCF}"; then
        echo "[$(date +%H:%M:%S)] Merging ancient + 1000G..."
        [[ -f "${THOUSANDG_VCF}.tbi" ]] || bcftools index -f -t "${THOUSANDG_VCF}"
        bcftools merge "${THOUSANDG_VCF}" "${ANCIENT_VCF}" \
            -Oz -o "${MERGED_VCF}.tmp"
        mv "${MERGED_VCF}.tmp" "${MERGED_VCF}"
    else
        echo "  [skip] $(basename "${MERGED_VCF}") already present"
    fi

    # 2c. Strip commas inside quoted fields
    if should_run "${NOCOMMAS_VCF}"; then
        echo "[$(date +%H:%M:%S)] Stripping commas inside quoted fields..."
        zcat "${MERGED_VCF}" \
            | perl -pe 's/(".+?[^\\]")/($ret = $1) =~ (s#,##g); $ret/ge' \
            | bgzip -c > "${NOCOMMAS_VCF}.tmp"
        mv "${NOCOMMAS_VCF}.tmp" "${NOCOMMAS_VCF}"
    else
        echo "  [skip] $(basename "${NOCOMMAS_VCF}") already present"
    fi

    # 2d. Convert directly to EIGENSTRAT
    if should_run "${OUT_DIR}/${PREFIX}.geno" \
                  "${OUT_DIR}/${PREFIX}.snp"  \
                  "${OUT_DIR}/${PREFIX}.ind"; then
        echo "[$(date +%H:%M:%S)] Converting VCF -> EIGENSTRAT..."
        ( cd "${OUT_DIR}" && bash "${CONVERT_SH}" "${NOCOMMAS_VCF}" )
    else
        echo "  [skip] eigenstrat files present for ${PREFIX}"
    fi

    # 2e. Build labelled .ind_
    if should_run "${OUT_DIR}/${PREFIX}.ind_"; then
        echo "[$(date +%H:%M:%S)] Building labelled .ind_..."
        python3 "${LABEL_PY}" \
            --ind     "${OUT_DIR}/${PREFIX}.ind" \
            --meta    "${META}" \
            --name2id "${NAME2ID}" \
            --ref-pop "${REF_POP}" \
            --igsr    "${IGSR}" \
            --out     "${OUT_DIR}/${PREFIX}.ind_.tmp"
        mv "${OUT_DIR}/${PREFIX}.ind_.tmp" "${OUT_DIR}/${PREFIX}.ind_"
    else
        echo "  [skip] ${PREFIX}.ind_ already present"
    fi

    echo "[$(date +%H:%M:%S)] Done: chr ${CHR}"
done

echo
echo "All chromosomes processed. To re-run any step, delete its output file"
echo "or invoke with FORCE=1 (e.g. 'FORCE=1 bash prep_data.sh')."
echo "Run f4_HLA_analysis_v5.R next (memory-bounded; v4 retained for provenance)."
