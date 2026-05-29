#!/bin/bash
#SBATCH --job-name=ihbca_ref
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00

# =============================================================================
# Extract iHBCA Reference from author_share
# =============================================================================
# Extracts the iHBCA cell inventory from ihbca_level1.5_annotations.csv.
# This becomes the reference for mapping component study cell IDs.
#
# Outputs:
#   - outputs/ihbca_reference/ihbca_reference_summary.yaml
#   - outputs/ihbca_reference/ihbca_cell_inventory.csv
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../../config/load_paths.sh"

echo "=== Extract iHBCA Reference ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo

# Paths
PROJECT_DIR="${PROJECT_IHBCAV1_UPLOAD}"
PHASE_A_DIR="${PROJECT_DIR}/external_studies/harmonization"
SCRIPT="${PHASE_A_DIR}/scripts/extract_ihbca_reference.py"
INPUT="${SOURCE_AUTHOR_SHARE}/ihbca_level1.5_annotations.csv"
OUTPUT_DIR="${PHASE_A_DIR}/outputs/ihbca_reference"

# Container
CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"

# Verify input exists
if [ ! -f "${INPUT}" ]; then
    echo "ERROR: Input file not found: ${INPUT}"
    exit 1
fi

# Load singularity (REQUIRED - not available on compute nodes by default)
module load singularity

# Create output directory
mkdir -p "${OUTPUT_DIR}"

echo "Input: ${INPUT}"
echo "Output: ${OUTPUT_DIR}"
echo "Container: ${CONTAINER}"
echo

# Run extraction
singularity exec \
    ${BIND_MOUNTS} \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
    "${CONTAINER}" \
    python3 "${SCRIPT}" \
        --input "${INPUT}" \
        --output-dir "${OUTPUT_DIR}"

echo
echo "Completed: $(date)"
