#!/bin/bash
#SBATCH --job-name=kumar_extract
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00

# =============================================================================
# KUMAR TRACK 2: UNSUPERVISED RDS EXTRACTION
# =============================================================================
# Extracts cell IDs, metadata, and structure from the ORIGINAL Kumar RDS file.
# Answers open questions from Track 1 (edge_cases.yaml):
#   1. missing_samples_explanation - c01-c65 gap
#   2. cell_id_unsupervised_analysis - actual format discovery
#   3. metadata_column_mapping - supplemental to RDS mapping
#   4. ethnicity_as_condition - is ethnicity present?
#   5. bmi_as_condition - is BMI present?
#   6. cell_distribution - 714K cells across samples
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../../../../config/load_paths.sh"


echo "=========================================="
echo "Kumar Track 2: Unsupervised RDS Extraction"
echo "Start time: $(date)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "=========================================="

# Paths
BASE_DIR="${SOURCE_IHBCAV1_HARMONIZATION%/}"
KUMAR_DIR="${BASE_DIR}/studies/kumar"
SCRIPT="${KUMAR_DIR}/scripts/extract_rds.R"
R_LIBS_USER="${R_LIBS_R_SPATIAL_4_3_3}"
CONTAINER="${CONTAINER_R_SPATIAL_4_3_3}"

# Create log directory if needed
LOG_DIR="${KUMAR_DIR}/logs"
mkdir -p "${LOG_DIR}"

# Load singularity
module load singularity

echo ""
echo "Source file: ${SOURCE_COMPONENT_STUDIES}/kumar.rds"
echo "Output: ${KUMAR_DIR}/extracted/raw_data_extraction.yaml"
echo ""

# Run extraction
singularity exec \
    --bind ${SHARED_LAB}:${SHARED_LAB}:ro \
    --bind ${USER_ROOT}:${USER_ROOT}:rw \
    ${BIND_MOUNTS} \
    --bind "${R_LIBS_USER}:/home/jovyan/R/library:ro" \
    --env "R_LIBS_USER=/home/jovyan/R/library" \
    "${CONTAINER}" \
    Rscript "${SCRIPT}"

echo ""
echo "=========================================="
echo "Extraction complete"
echo "End time: $(date)"
echo "=========================================="

# Show output file
OUTPUT_FILE="${KUMAR_DIR}/extracted/raw_data_extraction.yaml"
if [ -f "${OUTPUT_FILE}" ]; then
    echo ""
    echo "Output file created: ${OUTPUT_FILE}"
    echo "File size: $(du -h "${OUTPUT_FILE}" | cut -f1)"
else
    echo ""
    echo "WARNING: Output file not found!"
    exit 1
fi
