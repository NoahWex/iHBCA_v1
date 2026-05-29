#!/bin/bash
#SBATCH --job-name=l15_mapping
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=01:00:00

# =============================================================================
# Build level1.5_annotation -> CL term mapping table
# =============================================================================
# Reads integrated h5ad (backed mode), extracts label -> CL distributions
# Output: publication/mappings/level15_to_cl_mapping_raw.csv
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
REPO_ROOT="${PROJECT_IHBCAV1_UPLOAD}"
PY_CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"

echo "=== Level1.5 -> CL Mapping Builder ==="
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
echo "Starting mapping build..."
singularity exec \
    --no-mount bind-paths \
    ${BIND_MOUNTS} \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
    "$PY_CONTAINER" \
    python3 "${REPO_ROOT}/publication/scripts/build_level15_mapping.py" \
        --repo-root "$REPO_ROOT"

echo ""
echo "=== Mapping Build Complete ==="
echo "Exit code: $?"
