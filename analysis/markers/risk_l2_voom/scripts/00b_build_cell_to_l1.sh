#!/bin/bash
#SBATCH --job-name=risk_l1_build
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519/logs/00b_cell_to_l1_%j.out
#SBATCH --error=/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519/logs/00b_cell_to_l1_%j.err

# Pattern source: 00_build_cell_to_l2.sh — identical wrapper, swaps the Python
# script and output paths to write cell_to_l1.csv instead.
# Tier 2 sizing (8G, 2 CPUs, 20 min) for the cell→L1 alignment.

set -euo pipefail

PROJECT_ROOT="/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519"
COMPONENTS="/share/crsp/lab/dalawson/nwechter/iHBCAv1_upload/.dev/assembly_rebuild_20260405/outputs/components"
LABELS_FULL="/share/crsp/lab/dalawson/nwechter/iHBCA_publication/publication/analysis/annotation/labels_full.csv"

OUT_DIR="${PROJECT_ROOT}/outputs/pseudobulk"
mkdir -p "${OUT_DIR}"

module load singularity
CONTAINER="/dfs7/singularity_containers/rcic/JHUB3/Rocky8_jupyter_base_R4.3.3_Spatial.sif"

singularity exec \
    --bind /share/crsp/lab/dalawson:/share/crsp/lab/dalawson:rw \
    --bind /dfs7:/dfs7:ro \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
    "${CONTAINER}" \
    python3 "${PROJECT_ROOT}/scripts/00b_build_cell_to_l1.py" \
        --epi-meta          "${COMPONENTS}/epi_metadata_enriched.csv" \
        --imm-meta          "${COMPONENTS}/imm_metadata_enriched.csv" \
        --str-meta          "${COMPONENTS}/str_metadata_enriched.csv" \
        --labels-full       "${LABELS_FULL}" \
        --output-csv        "${OUT_DIR}/cell_to_l1.csv" \
        --donor-summary-csv "${OUT_DIR}/donor_x_l1_cell_counts.csv"

echo "Outputs:"
ls -la "${OUT_DIR}/cell_to_l1.csv" "${OUT_DIR}/donor_x_l1_cell_counts.csv"
