#!/bin/bash
#SBATCH --job-name=ihbca_anno_str_markers
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=${HPC_LOGS_DIR}/%x_%j.out
#SBATCH --error=${HPC_LOGS_DIR}/%x_%j.err

# =============================================================================
# Per-resolution top-marker heatmaps for the stromal compartment.
#
# Tier 4 (64G / 4 CPU / 2h) per hpc-resource-rules.md — 974K-cell counts
# matrix (~3.8GB) dominates the memory footprint.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../config/load_paths.sh"

COMPARTMENT="str"
OUT_DIR="${PROJECT_PUBLICATION}/analysis/annotation/${COMPARTMENT}"
STAGE_OUT="${OUT_DIR}/structural/"
SCANVI_DIR="${SOURCE_IHBCAV1_ASSEMBLY_SCANVI}${COMPARTMENT}/n_latent_50"

LEIDEN="${SCANVI_DIR}/${COMPARTMENT}_leiden_barcoded.csv"
COUNTS="${SOURCE_IHBCAV1_ASSEMBLY_COMPONENTS}${COMPARTMENT}_counts.npz"
GENE_DATA="${SOURCE_IHBCAV1_ASSEMBLY_COMPONENTS}gene_data.csv"
LIMMA_DIR="${SOURCE_IHBCA_V1_ANNOTATION}markers/limma/stromal"
CONFIG="${SCRIPT_DIR}/../config/compartment_config.yaml"
SCRIPT="${SCRIPT_DIR}/../scripts/01_compute_markers.py"

mkdir -p "${STAGE_OUT}" "${PROJECT_PUBLICATION}/logs"

for f in "$LEIDEN" "$COUNTS" "$GENE_DATA" "$CONFIG" "$SCRIPT"; do
    [ -f "$f" ] || { echo "ERROR: missing input: $f"; exit 1; }
done
[ -d "$LIMMA_DIR" ] || { echo "ERROR: missing limma dir: $LIMMA_DIR"; exit 1; }

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
        --counts-npz "$COUNTS" \
        --gene-data "$GENE_DATA" \
        --limma-dir "$LIMMA_DIR" \
        --config "$CONFIG" \
        --out-dir "$OUT_DIR"

echo "=== Complete ==="
