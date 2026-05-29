#!/bin/bash
#SBATCH --job-name=map_gray
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:30:00

# =============================================================================
# Map Gray Component Study Cells to iHBCA Reference
# =============================================================================
# Maps Gray cells (52,681) to iHBCA cell IDs by barcode matching.
# Expected: 100% match rate (no collisions for Gray).
#
# Verifies: gray-rm-donor-mismatch (RM-* donors should appear with ~19K WT cells)
#
# Outputs:
#   - studies/gray/outputs/gray_cells.csv (all cells with gray_* prefixed columns)
#   - studies/gray/outputs/mapping_findings.yaml (batch→patient mapping discovered)
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../../config/load_paths.sh"

echo "=== Map Gray Cells ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo

# Paths
PROJECT_DIR="${PROJECT_IHBCAV1_UPLOAD}"
PHASE_A_DIR="${PROJECT_DIR}/external_studies/harmonization"
SCRIPT="${PHASE_A_DIR}/scripts/map_gray_cells.py"

GRAY_H5AD="${SOURCE_COMPONENT_STUDIES}gray.h5ad"
IHBCA_INVENTORY="${PHASE_A_DIR}/outputs/ihbca_reference/ihbca_cell_inventory.csv"
OUTPUT_DIR="${PHASE_A_DIR}/studies/gray/outputs"

# Container
CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"

# Verify inputs exist
for INPUT in "${GRAY_H5AD}" "${IHBCA_INVENTORY}"; do
    if [ ! -f "${INPUT}" ]; then
        echo "ERROR: Input file not found: ${INPUT}"
        exit 1
    fi
done

# Load singularity
module load singularity

# Create output directory
mkdir -p "${OUTPUT_DIR}"

echo "Gray h5ad: ${GRAY_H5AD}"
echo "iHBCA inventory: ${IHBCA_INVENTORY}"
echo "Output: ${OUTPUT_DIR}"
echo

# Run mapping
singularity exec \
    ${BIND_MOUNTS} \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
    "${CONTAINER}" \
    python3 "${SCRIPT}" \
        --gray-h5ad "${GRAY_H5AD}" \
        --ihbca-inventory "${IHBCA_INVENTORY}" \
        --output-dir "${OUTPUT_DIR}"

echo
echo "Completed: $(date)"
