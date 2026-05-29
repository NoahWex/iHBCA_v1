#!/bin/bash
#SBATCH --job-name=ihbca_anno_imm_structural
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=${HPC_LOGS_DIR}/%x_%j.out
#SBATCH --error=${HPC_LOGS_DIR}/%x_%j.err

# =============================================================================
# Pre-annotation structural evidence package for the immune compartment.
#
# Produces 13-resolution Leiden UMAP grid, native-label F1 concordance heatmaps
# per resolution, study composition, and per-resolution top-marker heatmaps.
# Inputs are canonical scANVI integration outputs from iHBCAv1_upload; output
# is a self-contained evidence package consumed by annotation review.
#
# Tier 3 (32G / 4 CPU / 2h) per hpc-resource-rules.md — accounts for ~162K-cell
# compartment counts matrix + 13 UMAP renderings + F1 grid.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../config/load_paths.sh"

COMPARTMENT="imm"
OUT_DIR="${PROJECT_PUBLICATION}/analysis/annotation/${COMPARTMENT}"
STAGE_OUT="${OUT_DIR}/structural/"
SCANVI_DIR="${SOURCE_IHBCAV1_ASSEMBLY_SCANVI}${COMPARTMENT}/n_latent_50"

LEIDEN="${SCANVI_DIR}/${COMPARTMENT}_leiden_barcoded.csv"
UMAP="${SCANVI_DIR}/${COMPARTMENT}_scanvi_umap_barcoded.csv"
META="${SOURCE_IHBCAV1_ASSEMBLY_COMPONENTS}${COMPARTMENT}_metadata_enriched.csv"
NATIVE="${SOURCE_IHBCAV1_ASSEMBLY_COMPONENTS}native_labels.csv"
CONFIG="${SCRIPT_DIR}/../config/compartment_config.yaml"
SCRIPT="${SCRIPT_DIR}/../scripts/00_build_structural_package.py"

mkdir -p "${STAGE_OUT}" "${PROJECT_PUBLICATION}/logs"

for f in "$LEIDEN" "$UMAP" "$META" "$NATIVE" "$CONFIG" "$SCRIPT"; do
    [ -f "$f" ] || { echo "ERROR: missing input: $f"; exit 1; }
done

module load singularity

singularity exec --no-home \
    ${BIND_MOUNTS} \
    --env "PYTHONUSERBASE=${PYTHON_LIBS_SCVI_GPU_2025Q3_USERBASE}" \
    --env "PYTHONUNBUFFERED=1" \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
    "${CONTAINER_PYTHON_SCVI_GPU_2025Q3}" \
    python "$SCRIPT" \
        --compartment "$COMPARTMENT" \
        --leiden-csv "$LEIDEN" \
        --umap-csv "$UMAP" \
        --metadata-csv "$META" \
        --native-labels-csv "$NATIVE" \
        --config "$CONFIG" \
        --out-dir "$OUT_DIR"

echo "=== Complete ==="
