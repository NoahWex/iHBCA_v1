#!/bin/bash
#SBATCH --job-name=build_sketch
#SBATCH --account=dalawson_lab
#SBATCH --partition=free
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=04:00:00

# =============================================================================
# Build geometric sketch of integrated iHBCA object
# =============================================================================
# Produces all-breast-cells-sketch.h5ad (~280K cells, 1000/donor via geosketch)
# from the fully enriched all-breast-cells.h5ad (2.1M cells).
#
# 256GB: full h5ad load (~95GB) + subsetting + write overhead.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
REPO_ROOT="${PROJECT_IHBCAV1_UPLOAD}"
SCRIPTS="${REPO_ROOT}/publication/scripts"
OUTDIR="${REPO_ROOT}/publication/outputs/integrated_objects"
CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"

INPUT="${OUTDIR}/all-breast-cells.h5ad"
OUTPUT="${OUTDIR}/all-breast-cells-sketch.h5ad"

echo "============================================="
echo "Build sketch object"
echo "  Date: $(date)"
echo "  Node: $(hostname)"
echo "  Job ID: ${SLURM_JOB_ID}"
echo "  Memory: ${SLURM_MEM_PER_NODE}MB"
echo "  Input: ${INPUT}"
echo "  Output: ${OUTPUT}"
echo "============================================="

if [ ! -f "${INPUT}" ]; then
    echo "ERROR: Input h5ad not found: ${INPUT}"
    exit 1
fi

echo "Input size: $(ls -lh "${INPUT}" | awk '{print $5}')"

module load singularity

singularity exec \
    --no-mount bind-paths \
    ${BIND_MOUNTS} \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
    "$CONTAINER" \
    python3 "${SCRIPTS}/build_sketch.py" \
        --input "$INPUT" \
        --output "$OUTPUT" \
        --cells-per-donor 1000 \
        --seed 42

echo ""
echo "============================================="
echo "Sketch build COMPLETE"
echo "Output size: $(ls -lh "${OUTPUT}" | awk '{print $5}')"
echo "Date: $(date)"
echo "============================================="
