#!/bin/bash
#SBATCH --job-name=pal_epi_extract
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00

# =============================================================================
# Pal NormEpi: Extract Raw Data from RDS
# =============================================================================
# Performs UNSUPERVISED extraction of cell IDs, metadata, and structure from
# the Pal NormEpi Seurat object to answer open questions in edge_cases.yaml.
#
# Output: studies/pal_norm_epi/extracted/raw_data_extraction.yaml
#
# Usage:
#   sbatch extract_rds.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../../../../config/load_paths.sh"


# Configuration
CONTAINER="${CONTAINER_R_SPATIAL_4_3_3}"
R_LIBS_USER="${R_LIBS_R_SPATIAL_4_3_3}"
SCRIPT="${SOURCE_IHBCAV1_HARMONIZATION%/}/studies/pal_norm_epi/scripts/extract_rds.R"
LOG_DIR="${SOURCE_IHBCAV1_HARMONIZATION%/}/studies/pal_norm_epi/logs"

echo "============================================"
echo "Pal NormEpi: RDS Extraction"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Time: $(date)"
echo "============================================"

# Create logs directory
mkdir -p "${LOG_DIR}"

# Load singularity module
module load singularity

# Run the R script in container
singularity exec \
    --bind ${SHARED_LAB}:${SHARED_LAB}:ro \
    --bind ${USER_ROOT}:${USER_ROOT}:rw \
    ${BIND_MOUNTS} \
    --bind "${R_LIBS_USER}:/home/jovyan/R/library:ro" \
    --env "R_LIBS_USER=/home/jovyan/R/library" \
    "${CONTAINER}" \
    Rscript "${SCRIPT}"

echo "============================================"
echo "Extraction complete"
echo "Time: $(date)"
echo "============================================"
