#!/bin/bash
#SBATCH --job-name=track2_murrow
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../../../../config/load_paths.sh"

# =============================================================================
# MURROW TRACK 2 EXTRACTION
# =============================================================================
# Extract raw data schema from murrow.rds
# =============================================================================

echo "=== MURROW TRACK 2 EXTRACTION ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start time: $(date)"
echo ""

# Paths
BASE_DIR="${SOURCE_IHBCAV1_HARMONIZATION%/}/studies/murrow"
SCRIPT="${BASE_DIR}/scripts/extract_rds.R"
CONTAINER="${CONTAINER_R_SPATIAL_4_3_3}"
R_LIBS_USER="${R_LIBS_R_SPATIAL_4_3_3}"

# Create logs directory if needed
mkdir -p "${BASE_DIR}/logs"

# Verify source file exists
SOURCE_FILE="${SOURCE_COMPONENT_STUDIES}/murrow_l1_original.rds"
if [ ! -f "$SOURCE_FILE" ]; then
    echo "ERROR: Source file not found: $SOURCE_FILE"
    exit 1
fi
echo "Source file: $SOURCE_FILE"
echo "Source size: $(ls -lh "$SOURCE_FILE" | awk '{print $5}')"
echo ""

# Load singularity
module load singularity

# Run extraction
echo "Starting R extraction script..."
singularity exec \
    --bind ${SHARED_LAB}:${SHARED_LAB}:ro \
    --bind ${USER_ROOT}:${USER_ROOT}:rw \
    ${BIND_MOUNTS} \
    --bind "$R_LIBS_USER:/home/jovyan/R/library:ro" \
    --env "R_LIBS_USER=/home/jovyan/R/library" \
    "$CONTAINER" \
    Rscript "$SCRIPT"

EXIT_CODE=$?

echo ""
echo "Exit code: $EXIT_CODE"
echo "End time: $(date)"

if [ $EXIT_CODE -eq 0 ]; then
    echo "=== EXTRACTION COMPLETED SUCCESSFULLY ==="
else
    echo "=== EXTRACTION FAILED ==="
fi

exit $EXIT_CODE
