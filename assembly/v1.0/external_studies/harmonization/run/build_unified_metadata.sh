#!/bin/bash
#SBATCH --job-name=build_unified
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00

# =============================================================================
# Build Unified Cell Metadata - SLURM Wrapper
# =============================================================================
# Estimated memory: ~8-10GB (64GB requested for safety)
# Estimated time: ~30-60 minutes for full run
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../../config/load_paths.sh"

# Setup
BASE_DIR="${SOURCE_IHBCAV1_HARMONIZATION%/}"
SCRIPT="${BASE_DIR}/scripts/build_unified_metadata.py"
CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"

echo "=============================================="
echo "Job: Build Unified Cell Metadata"
echo "=============================================="
echo "Start time: $(date)"
echo "Host: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Working directory: ${BASE_DIR}"
echo ""

# Load singularity
module load singularity

# Run the script
echo "Running build_unified_metadata.py..."
singularity exec \
    ${BIND_MOUNTS} \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
    "${CONTAINER}" \
    python "${SCRIPT}" \
        --base-dir "${BASE_DIR}" \
        --ihbca-source "${SOURCE_AUTHOR_SHARE}/ihbca_level1.5_annotations.csv"

# Completion
echo ""
echo "=============================================="
echo "Job completed: $(date)"
echo "=============================================="
