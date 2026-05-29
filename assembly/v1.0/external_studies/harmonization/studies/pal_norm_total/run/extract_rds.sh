#!/bin/bash
#SBATCH --job-name=extract_pal_normtotal
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00

# =============================================================================
# PAL NORM_TOTAL: RAW DATA SCHEMA EXTRACTION
# =============================================================================
# Track 2 Wave 1 - Unsupervised extraction from SeuratObject_NormTotal.rds
#
# Open Questions to Answer:
# 1. cell_id_type_field: Confirm cell IDs use 'mix' or 'total'?
# 2. cell_id_format_verification: What is the complete cell ID format?
# 3. cell_count_verification: Confirm total cell count is 54,332
# 4. metadata_columns_in_rds: What metadata columns exist?
# 5. cell_distribution: How do cells distribute across 13 samples?
# 6. menopausal_assignment_verification: Does RDS match R script assignments?
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../../../../config/load_paths.sh"


echo "=========================================="
echo "Pal NormTotal: Raw Data Schema Extraction"
echo "Start time: $(date)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "=========================================="

# Paths
BASE_DIR="${SOURCE_IHBCAV1_HARMONIZATION%/}"
SCRIPT="${BASE_DIR}/scripts/extract_raw_schema.R"
R_LIBS_USER="${R_LIBS_R_SPATIAL_4_3_3}"
CONTAINER="${CONTAINER_R_SPATIAL_4_3_3}"

# Create study-specific logs directory if needed
mkdir -p "${BASE_DIR}/studies/pal_norm_total/logs"

# Load singularity
module load singularity

echo ""
echo "Source file: ${SOURCE_PAL_ORIGINAL%/}/SeuratObject_NormTotal.rds"
echo "Output: ${BASE_DIR}/studies/pal_norm_total/extracted/raw_data_extraction.yaml"
echo ""

singularity exec \
    --bind ${SHARED_LAB}:${SHARED_LAB}:ro \
    --bind ${USER_ROOT}:${USER_ROOT}:rw \
    ${BIND_MOUNTS} \
    --bind "${R_LIBS_USER}:/home/jovyan/R/library:ro" \
    --env "R_LIBS_USER=/home/jovyan/R/library" \
    "${CONTAINER}" \
    Rscript "${SCRIPT}" --study pal_norm_total

echo ""
echo "=========================================="
echo "Extraction complete"
echo "End time: $(date)"
echo "=========================================="

# Show output file
OUTPUT_FILE="${BASE_DIR}/studies/pal_norm_total/extracted/raw_data_extraction.yaml"
if [ -f "${OUTPUT_FILE}" ]; then
    echo ""
    echo "Output file created:"
    ls -la "${OUTPUT_FILE}"
    echo ""
    echo "First 100 lines of extraction:"
    head -100 "${OUTPUT_FILE}"
fi
