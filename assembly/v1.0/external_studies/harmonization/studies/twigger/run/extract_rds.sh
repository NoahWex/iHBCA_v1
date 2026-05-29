#!/bin/bash
#SBATCH --job-name=twigger_rds_extract
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00

# ==============================================================================
# TWIGGER RDS EXTRACTION JOB
# ==============================================================================
# Purpose: Extract metadata from Twigger RDS to resolve critical mapping issues
# Priority: HIGH
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../../../../config/load_paths.sh"


echo "============================================================"
echo "TWIGGER RDS EXTRACTION"
echo "============================================================"
echo "Job ID:     ${SLURM_JOB_ID}"
echo "Node:       $(hostname)"
echo "Start time: $(date)"
echo "============================================================"

# Paths
PROJECT_BASE="${PROJECT_IHBCAV1_UPLOAD}"
STUDY_DIR="${PROJECT_BASE}/external_studies/harmonization/studies/twigger"
SCRIPT="${STUDY_DIR}/scripts/extract_rds_metadata.R"
CONTAINER="${CONTAINER_R_SPATIAL_4_3_3}"

# R library path for user packages
R_LIBS_USER="${R_LIBS_R_SPATIAL_4_3_3}"

# Ensure logs directory exists
mkdir -p "${STUDY_DIR}/logs"

# Ensure extracted directory exists
mkdir -p "${STUDY_DIR}/extracted"

echo ""
echo "Script:     ${SCRIPT}"
echo "Container:  ${CONTAINER}"
echo "R_LIBS:     ${R_LIBS_USER}"
echo ""

# Load singularity module
module load singularity

# Verify container exists
if [[ ! -f "${CONTAINER}" ]]; then
    echo "ERROR: Container not found: ${CONTAINER}"
    exit 1
fi

# Verify script exists
if [[ ! -f "${SCRIPT}" ]]; then
    echo "ERROR: Script not found: ${SCRIPT}"
    exit 1
fi

echo "Running extraction script..."
echo "============================================================"

singularity exec \
    --bind ${SHARED_LAB}:${SHARED_LAB}:ro \
    --bind ${USER_ROOT}:${USER_ROOT}:rw \
    ${BIND_MOUNTS} \
    --bind "${R_LIBS_USER}:/home/jovyan/R/library:ro" \
    --env "R_LIBS_USER=/home/jovyan/R/library" \
    "${CONTAINER}" \
    Rscript "${SCRIPT}"

EXIT_CODE=$?

echo ""
echo "============================================================"
echo "Job completed with exit code: ${EXIT_CODE}"
echo "End time: $(date)"
echo "============================================================"

if [[ ${EXIT_CODE} -eq 0 ]]; then
    echo "SUCCESS: Extraction complete"
    echo "Output: ${STUDY_DIR}/extracted/raw_data_extraction.yaml"
else
    echo "FAILED: Check logs for errors"
fi

exit ${EXIT_CODE}
