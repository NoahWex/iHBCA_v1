#!/bin/bash
# =============================================================================
# PAL NORM B1: UNSUPERVISED RDS EXTRACTION
# =============================================================================
# Extracts ALL cell IDs, metadata, and structure from the NormB1Total Seurat
# RDS file using unsupervised pattern analysis.
#
# Answers open questions from edge_cases.yaml:
#   1. cell_id_format_b1: Confirm BRCA1 cell ID format
#   2. cell_count_verification: Verify expected cell counts
#   3. donor_0123_exclusion_verification: Count N_0123_* cells
#   4. metadata_columns_in_rds: What metadata columns exist?
#   5. cell_distribution_by_donor: Cell distribution across donors
#   6. ihbca_donor_mapping_verification: Verify donor mappings
#
# Usage:
#   sbatch ${SOURCE_IHBCAV1_HARMONIZATION%/}/studies/pal_norm_b1/run/extract_rds.sh
# =============================================================================

#SBATCH --job-name=pal_normb1_extract
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

set -e

echo "=== Job Info ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start time: $(date)"
echo ""

# =============================================================================
# Configuration
# =============================================================================

PROJECT_ROOT="${PROJECT_IHBCAV1_UPLOAD}"
STUDY_DIR="${PROJECT_ROOT}/external_studies/harmonization/studies/pal_norm_b1"
SCRIPT="${STUDY_DIR}/scripts/extract_rds.R"
CONTAINER="${CONTAINER_R_SPATIAL_4_3_3}"
R_LIBS_USER="${R_LIBS_R_SPATIAL_4_3_3}"

# Ensure output directory exists
mkdir -p "${STUDY_DIR}/extracted"
mkdir -p "${STUDY_DIR}/logs"

# Load singularity
module load singularity

# =============================================================================
# Execution
# =============================================================================

echo "Running extraction script..."
echo "Script: ${SCRIPT}"
echo "Container: ${CONTAINER}"
echo ""

singularity exec \
    --bind ${SHARED_LAB}:${SHARED_LAB}:ro \
    --bind ${USER_ROOT}:${USER_ROOT}:rw \
    ${BIND_MOUNTS} \
    --bind "${R_LIBS_USER}:/home/jovyan/R/library:ro" \
    --env "R_LIBS_USER=/home/jovyan/R/library" \
    "${CONTAINER}" \
    Rscript "${SCRIPT}"

echo ""
echo "=== Job Complete ==="
echo "End time: $(date)"
echo "Output: ${STUDY_DIR}/extracted/raw_data_extraction.yaml"
