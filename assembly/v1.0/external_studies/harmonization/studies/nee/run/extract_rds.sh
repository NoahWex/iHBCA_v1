#!/bin/bash
#SBATCH --job-name=nee_extract
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../../../../config/load_paths.sh"

# =============================================================================
# NEE RDS DATA EXTRACTION
# =============================================================================
# UNSUPERVISED extraction from Nee RDS file
# Focus: Verify cell ID patterns, metadata columns, and donor distribution
# =============================================================================

echo "============================================================================="
echo "NEE RDS EXTRACTION - SLURM JOB"
echo "============================================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "============================================================================="

# Configuration
CONTAINER="${CONTAINER_R_SPATIAL_4_3_3}"
SCRIPT="${SOURCE_IHBCAV1_HARMONIZATION%/}/studies/nee/scripts/extract_rds_data.R"
R_LIBS_USER="${R_LIBS_R_SPATIAL_4_3_3}"

# Verify inputs exist
if [ ! -f "$CONTAINER" ]; then
    echo "ERROR: Container not found: $CONTAINER"
    exit 1
fi

if [ ! -f "$SCRIPT" ]; then
    echo "ERROR: R script not found: $SCRIPT"
    exit 1
fi

# Load singularity
module load singularity

echo ""
echo "Running R script in container..."
echo ""

# Execute R script
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
echo "============================================================================="
echo "JOB COMPLETED"
echo "============================================================================="
echo "Exit code: $EXIT_CODE"
echo "Finished: $(date)"
echo "============================================================================="

exit $EXIT_CODE
