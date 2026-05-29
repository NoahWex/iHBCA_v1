#!/bin/bash
#SBATCH --job-name=risk_l1_voom_de
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH --output=/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519/logs/02b_voom_de_l1_%j.out
#SBATCH --error=/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519/logs/02b_voom_de_l1_%j.err

# L1-resolution voom DE: 2 cohorts x 1 formula x 22 L1 types = 44 model fits.
# Reuses scripts/02_voom_de_sweep.R (after the 2026-05-20 brca_genotype filter
# extension) with --counts pointing to L1 pseudobulk and --inquiry-yaml at the
# new inquiry_l1.yaml config.
#
# Pattern source: 02_voom_de_sweep.sh (this stage). Same engine, L1 substrate.
# Tier 3 sizing (32G, 4 CPUs, 1:30) — fewer groups (22 vs 42) at L1 grain so
# wall time is shorter than the L2 sweep.

set -euo pipefail

PROJECT_ROOT="/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519"
GENE_MAPPING="/share/crsp/lab/dalawson/nwechter/iHBCAv1_upload/iHBCA_upload/publication/mappings/gene_symbol_to_ensembl.tsv"

COUNTS="${PROJECT_ROOT}/outputs/pseudobulk/risk_l1_voom_pseudobulk_counts.csv"
PB_META="${PROJECT_ROOT}/outputs/pseudobulk/risk_l1_voom_pseudobulk_meta.csv"
DONOR_META="${PROJECT_ROOT}/outputs/pseudobulk/enriched_donor_metadata.csv"
INQUIRY_YAML="${PROJECT_ROOT}/config/inquiry_l1.yaml"

module load singularity
CONTAINER="/dfs7/singularity_containers/rcic/JHUB3/Rocky8_jupyter_base_R4.3.3_Spatial.sif"

for f in "${COUNTS}" "${PB_META}" "${DONOR_META}" "${INQUIRY_YAML}" \
         "${GENE_MAPPING}" "${PROJECT_ROOT}/scripts/02_voom_de_sweep.R"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: missing input: ${f}" >&2
        exit 1
    fi
done

singularity exec \
    --bind /share/crsp/lab/dalawson:/share/crsp/lab/dalawson:rw \
    --bind /dfs7:/dfs7:ro \
    --bind /pub/nwechter:/pub/nwechter:ro \
    --env "R_LIBS_USER=/pub/nwechter/biojhub4_dir/Rocky8_jupyter_base_R4.3.3_Spatial.sif/R/library" \
    "${CONTAINER}" \
    Rscript "${PROJECT_ROOT}/scripts/02_voom_de_sweep.R" \
        --project-root  "${PROJECT_ROOT}" \
        --counts        "${COUNTS}" \
        --pb-meta       "${PB_META}" \
        --donor-meta    "${DONOR_META}" \
        --gene-mapping  "${GENE_MAPPING}" \
        --inquiry-yaml  "${INQUIRY_YAML}"

# 02_voom_de_sweep.R writes to outputs/de_results/{cohort_id}/{formula_id}/
# and overwrites outputs/design_audit.csv. To avoid clobbering the L2 audit,
# the inquiry name in inquiry_l1.yaml is "risk_l1_voom_20260520" but the script
# writes to a fixed path. Move the L1 audit aside post-run.
mv "${PROJECT_ROOT}/outputs/design_audit.csv" \
   "${PROJECT_ROOT}/outputs/design_audit_l1.csv" || true

echo
echo "L1 DE complete. Outputs:"
ls -la "${PROJECT_ROOT}/outputs/de_results/A_full/A2_full/" 2>/dev/null | head
echo "---"
ls -la "${PROJECT_ROOT}/outputs/de_results/A_tested/A2_tested/" 2>/dev/null | head
