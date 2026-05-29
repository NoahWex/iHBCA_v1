#!/bin/bash
#SBATCH --job-name=risk_l2_prep
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519/logs/00_prep_%j.out
#SBATCH --error=/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519/logs/00_prep_%j.err

# Step 0: build cell_to_l2.csv + augment donor metadata with sample_type_coarse.
# Tier 2 sizing per ~/.claude/rules/hpc-resource-rules.md (3-CSV concat + groupby).
#
# Pattern source: load_paths convention from
#   risk_main_effects_20260507/run/*.sh (paths read from config/paths.yaml
#   via small shell parser).

set -euo pipefail

PROJECT_ROOT="/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519"

COMPONENTS="/share/crsp/lab/dalawson/nwechter/iHBCAv1_upload/.dev/assembly_rebuild_20260405/outputs/components"
LABELS_FULL="/share/crsp/lab/dalawson/nwechter/iHBCA_publication/publication/analysis/annotation/labels_full.csv"
EXISTING_DONOR_META="/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/da_pipeline/inquiries/risk_main_effects_20260507/config/derived_donor_metadata.csv"

OUT_DIR="${PROJECT_ROOT}/outputs/pseudobulk"
mkdir -p "${OUT_DIR}"

module load singularity

CONTAINER="/dfs7/singularity_containers/rcic/JHUB3/Rocky8_jupyter_base_R4.3.3_Spatial.sif"

echo "=== risk_l2_voom_20260519 / Step 0 prep ==="
echo "Host: $(hostname)"
echo "Start: $(date)"

singularity exec \
    --bind /share/crsp/lab/dalawson:/share/crsp/lab/dalawson:rw \
    --bind /dfs7:/dfs7:ro \
    --bind /dfs8:/dfs8:ro \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
    "${CONTAINER}" \
    python3 "${PROJECT_ROOT}/scripts/00_build_cell_to_l2.py" \
        --epi-meta   "${COMPONENTS}/epi_metadata_enriched.csv" \
        --imm-meta   "${COMPONENTS}/imm_metadata_enriched.csv" \
        --str-meta   "${COMPONENTS}/str_metadata_enriched.csv" \
        --labels-full "${LABELS_FULL}" \
        --output-csv         "${OUT_DIR}/cell_to_l2.csv" \
        --donor-summary-csv  "${OUT_DIR}/donor_x_l2_cellcounts.csv" \
        --donor-meta-out     "${OUT_DIR}/enriched_donor_metadata.csv" \
        --existing-donor-meta "${EXISTING_DONOR_META}" \
        --stc-heterogeneity-csv "${OUT_DIR}/sample_type_coarse_heterogeneity.csv"

echo "Done: $(date)"
