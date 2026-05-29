#!/bin/bash
#SBATCH --job-name=risk_l2_voom
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519/logs/02_voom_%j.out
#SBATCH --error=/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519/logs/02_voom_%j.err

# Step 2: limma voom DE sweep across 8 (cohort, formula) cells x ~25 L2 types.
# Tier 3 sizing per ~/.claude/rules/hpc-resource-rules.md.
#
# Pattern source: risk_main_effects_20260507/run/run_stageF3.sh (the existing
# inquiry's within-parent marker voom job — same container, same R libs).

set -euo pipefail

PROJECT_ROOT="/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519"
PB_DIR="${PROJECT_ROOT}/outputs/pseudobulk"
COUNTS="${PB_DIR}/risk_l2_voom_pseudobulk_counts.csv"
PB_META="${PB_DIR}/risk_l2_voom_pseudobulk_meta.csv"
# Use cohort_BR1_vs_AR_design.csv as donor meta — it has string-coded covariates
# (parity_binary = "parous"/"nulliparous", not 0/1) matching the inquiry.yaml
# filter values. Pre-filtered to 159 BR1+AR donors of the 287-donor atlas.
DONOR_META="/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/da_pipeline/inquiries/risk_main_effects_20260507/outputs/stageA_cohorts/cohort_BR1_vs_AR_design.csv"
GENE_MAPPING="/share/crsp/lab/dalawson/nwechter/iHBCAv1_upload/iHBCA_upload/publication/mappings/gene_symbol_to_ensembl.tsv"
INQUIRY_YAML="${PROJECT_ROOT}/config/inquiry.yaml"

module load singularity

CONTAINER="/dfs7/singularity_containers/rcic/JHUB3/Rocky8_jupyter_base_R4.3.3_Spatial.sif"
R_LIBS_USER_HOST="/pub/nwechter/biojhub4_dir/Rocky8_jupyter_base_R4.3.3_Spatial.sif/R/library"

echo "=== risk_l2_voom_20260519 / Step 2 voom sweep ==="
echo "Host: $(hostname)"
echo "Start: $(date)"

for f in "${COUNTS}" "${PB_META}" "${DONOR_META}" "${GENE_MAPPING}" "${INQUIRY_YAML}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: missing input: ${f}" >&2
        exit 1
    fi
done

singularity exec \
    --bind /share/crsp/lab/dalawson:/share/crsp/lab/dalawson:rw \
    --bind /dfs7:/dfs7:ro \
    --bind /dfs8:/dfs8:ro \
    --bind "${R_LIBS_USER_HOST}:/home/jovyan/R/library:ro" \
    --env "R_LIBS_USER=/home/jovyan/R/library" \
    "${CONTAINER}" \
    Rscript "${PROJECT_ROOT}/scripts/02_voom_de_sweep.R" \
        --project-root "${PROJECT_ROOT}" \
        --counts       "${COUNTS}" \
        --pb-meta      "${PB_META}" \
        --donor-meta   "${DONOR_META}" \
        --gene-mapping "${GENE_MAPPING}" \
        --inquiry-yaml "${INQUIRY_YAML}"

echo "Done: $(date)"
echo
echo "Outputs:"
ls -la "${PROJECT_ROOT}/outputs/design_audit.csv" || true
echo "DE CSVs:"
find "${PROJECT_ROOT}/outputs/de_results" -name '*.csv' | head -10
echo "Total DE CSVs: $(find "${PROJECT_ROOT}/outputs/de_results" -name '*.csv' | wc -l)"
