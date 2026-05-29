#!/bin/bash
#SBATCH --job-name=cl_crosswalk
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=02:00:00

# =============================================================================
# Build CL term crosswalk CSVs for non-CxG source datasets
# =============================================================================
# Reads integrated h5ad (backed mode) + source h5ads + cell_id_mapping.csv
# Outputs: publication/mappings/cl_term_crosswalk_{study}.csv
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
REPO_ROOT="${PROJECT_IHBCAV1_UPLOAD}"
PY_CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"

echo "=== CL Term Crosswalk Builder ==="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo ""

# Verify key inputs exist
INTEGRATED="${REPO_ROOT}/publication/outputs/integrated_objects/all-breast-cells.h5ad"
if [ ! -f "$INTEGRATED" ]; then
    echo "ERROR: Integrated h5ad not found: $INTEGRATED"
    exit 1
fi
echo "Input: $INTEGRATED ($(ls -lh "$INTEGRATED" | awk '{print $5}'))"

# Create output directories
mkdir -p "${REPO_ROOT}/publication/mappings"
mkdir -p "${REPO_ROOT}/publication/logs"

# Load singularity
module load singularity

echo ""
echo "Starting crosswalk build..."
singularity exec \
    --no-mount bind-paths \
    ${BIND_MOUNTS} \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
    "$PY_CONTAINER" \
    python3 "${REPO_ROOT}/publication/scripts/build_cl_crosswalk.py" \
        --repo-root "$REPO_ROOT"

echo ""
echo "=== Crosswalk Build Complete ==="
echo "Exit code: $?"
