#!/bin/bash
#SBATCH --job-name=reed_extract
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00

# =============================================================================
# REED RDS EXTRACTION
# =============================================================================
# Extract and analyze Reed RDS file to define canonical iHBCA format
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../../../../config/load_paths.sh"


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

PROJECT_ROOT="${PROJECT_IHBCAV1_UPLOAD}"
STUDY_DIR="${PROJECT_ROOT}/external_studies/harmonization/studies/reed"
SCRIPT="${STUDY_DIR}/scripts/extract_rds.R"
LOG_DIR="${STUDY_DIR}/run/logs"

CONTAINER="${CONTAINER_R_SPATIAL_4_3_3}"

# R user library path
R_LIBS_USER="${R_LIBS_R_SPATIAL_4_3_3}"

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

echo "=============================================="
echo "REED RDS EXTRACTION"
echo "=============================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start time: $(date)"
echo "Host: $(hostname)"
echo "Working directory: ${STUDY_DIR}"
echo "=============================================="

# Ensure log directory exists
mkdir -p "${LOG_DIR}"

# Ensure output directory exists
mkdir -p "${STUDY_DIR}/extracted"

# -----------------------------------------------------------------------------
# Load modules
# -----------------------------------------------------------------------------

module load singularity

# -----------------------------------------------------------------------------
# Run extraction
# -----------------------------------------------------------------------------

echo ""
echo "Running extraction script..."
echo "Container: ${CONTAINER}"
echo "Script: ${SCRIPT}"
echo ""

singularity exec \
    --bind ${SHARED_LAB}:${SHARED_LAB}:ro \
    --bind ${USER_ROOT}:${USER_ROOT}:rw \
    ${BIND_MOUNTS} \
    --bind "${R_LIBS_USER}:/home/jovyan/R/library:ro" \
    --env "R_LIBS_USER=/home/jovyan/R/library" \
    "${CONTAINER}" \
    Rscript "${SCRIPT}"

# -----------------------------------------------------------------------------
# Completion
# -----------------------------------------------------------------------------

echo ""
echo "=============================================="
echo "Extraction complete!"
echo "End time: $(date)"
echo "Output: ${STUDY_DIR}/extracted/raw_data_extraction.yaml"
echo "=============================================="
