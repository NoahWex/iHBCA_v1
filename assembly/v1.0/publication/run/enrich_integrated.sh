#!/bin/bash
#SBATCH --job-name=enrich_integrated
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=06:00:00

# =============================================================================
# Enrich integrated h5ad: all-breast-cells
# =============================================================================
# Enriches all-breast-cells.h5ad with var, obs, uns, and log_normalized layer.
#
# 256GB for all-breast-cells (2.1M cells + log_normalized layer).
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
REPO_ROOT="${PROJECT_IHBCAV1_UPLOAD}"
source "${REPO_ROOT}/publication/run/load_config.sh"
SCRIPTS="${REPO_ROOT}/publication/scripts"
# Env-var override for alternate builds (e.g., INTEGRATED_DIR=.../cxg_build/integrated_objects)
OUTDIR="${INTEGRATED_DIR:-${OUTPUT_INTEGRATED_OBJECTS}}"
CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"

echo "============================================="
echo "Enrich integrated h5ads"
echo "  Date: $(date)"
echo "  Node: $(hostname)"
echo "  Job ID: ${SLURM_JOB_ID}"
echo "  Memory: ${SLURM_MEM_PER_NODE}MB"
echo "============================================="

run_enrich() {
    local H5AD="$1"
    local SKIP_LAYERS="${2:-}"
    local FILENAME
    FILENAME=$(basename "$H5AD")

    echo ""
    echo "--- Processing: ${FILENAME} ---"
    if [ ! -f "${H5AD}" ]; then
        echo "ERROR: h5ad not found: ${H5AD}"
        return 1
    fi
    echo "Input size: $(ls -lh "${H5AD}" | awk '{print $5}')"

    singularity exec \
        --no-mount bind-paths \
        ${BIND_MOUNTS} \
        --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
        --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
        "$CONTAINER" \
        python3 "${SCRIPTS}/enrich_h5ads.py" \
            --h5ad "$H5AD" \
            --mode integrated \
            --repo-root "$REPO_ROOT" \
            ${SKIP_LAYERS}

    echo "Output size: $(ls -lh "${H5AD}" | awk '{print $5}')"
    echo "DONE: ${FILENAME}"
}

module load singularity

run_enrich "${OUTDIR}/all-breast-cells.h5ad"

echo ""
echo "============================================="
echo "Integrated enrichment COMPLETE"
echo "Date: $(date)"
echo "============================================="
