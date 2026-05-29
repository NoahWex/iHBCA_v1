#!/bin/bash
#SBATCH --job-name=prep_joint_refs
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2:00:00

# =============================================================================
# prepare_joint_references.sh
# =============================================================================
# Convert author_share CSV files to parquet format.
# Requires 64GB for loading 2.7GB CSV into memory.
#
# Usage:
#   sbatch prepare_joint_references.sh
# =============================================================================

set -euo pipefail

# Paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../../config/load_paths.sh"
BASE_PATH="${SOURCE_IHBCAV1_EXTERNAL_STUDIES%/}"
SCRIPT="${BASE_PATH}/harmonization/scripts/prepare_joint_references.py"
LOG_DIR="${BASE_PATH}/harmonization/logs"

# Container
CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"

# Create log directory
mkdir -p "$LOG_DIR"

echo "=============================================="
echo "Prepare Joint References"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start time: $(date)"
echo "=============================================="

# Load singularity
module load singularity

# Run Python script in container
singularity exec \
    ${BIND_MOUNTS} \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
    "$CONTAINER" \
    python3 "$SCRIPT"

echo ""
echo "=============================================="
echo "Completed"
echo "End time: $(date)"
echo "=============================================="
