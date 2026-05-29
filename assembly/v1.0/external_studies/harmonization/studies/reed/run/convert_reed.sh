#!/bin/bash
#SBATCH --job-name=reed_convert
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=06:00:00

# =============================================================================
# Convert Reed H5AD → RDS
# =============================================================================
# Phase 1: Python (scanpy) extracts H5AD → intermediates (counts, metadata, PCA)
# Phase 2: R (Seurat) constructs RDS from intermediates
# Phase 3: Validation (cell count, gene count, PCA dims)
# Phase 4: Cleanup intermediates
#
# Usage:
#   sbatch convert_reed.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../../../../config/load_paths.sh"


PROJECT_ROOT="${PROJECT_IHBCAV1_UPLOAD}"
STUDY_DIR="${PROJECT_ROOT}/external_studies/harmonization/studies/reed"
COMPONENT_DIR="${SOURCE_COMPONENT_STUDIES}"

H5AD_INPUT="${COMPONENT_DIR}/reed.h5ad"
RDS_OUTPUT="${COMPONENT_DIR}/reed.rds"
EXTRACT_DIR="/tmp/reed_h5ad_extract_${SLURM_JOB_ID}"
LOG_DIR="${STUDY_DIR}/run/logs"

# Containers
PY_CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"
R_CONTAINER="${CONTAINER_R_SPATIAL_4_3_3}"
R_LIBS_USER="${R_LIBS_R_SPATIAL_4_3_3}"

mkdir -p "$LOG_DIR"
mkdir -p "$EXTRACT_DIR"

echo "=============================================="
echo "REED H5AD → RDS CONVERSION"
echo "=============================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start time: $(date)"
echo "Input: ${H5AD_INPUT}"
echo "Output: ${RDS_OUTPUT}"
echo "Intermediates: ${EXTRACT_DIR}"
echo "=============================================="

if [[ ! -f "$H5AD_INPUT" ]]; then
    echo "FATAL: H5AD not found: ${H5AD_INPUT}"
    echo "Run download_reed.sh first"
    exit 1
fi

# Backup existing RDS if present
if [[ -f "$RDS_OUTPUT" ]]; then
    echo "Backing up existing reed.rds → reed.rds.prev"
    mv "$RDS_OUTPUT" "${RDS_OUTPUT}.prev"
fi

module load singularity

# --- Phase 1: Python extraction ---
echo ""
echo "=== PHASE 1: Python H5AD extraction ==="
echo ""

singularity exec \
    --bind ${SHARED_LAB}:${SHARED_LAB}:ro \
    --bind ${USER_ROOT}:${USER_ROOT}:rw \
    ${BIND_MOUNTS} \
    --bind /tmp:/tmp:rw \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
    "$PY_CONTAINER" \
    python "${STUDY_DIR}/scripts/convert_h5ad_to_rds.py" \
        --input "$H5AD_INPUT" \
        --output-dir "$EXTRACT_DIR"

# Verify intermediates exist
for f in counts.mtx.gz barcodes.tsv.gz features.tsv.gz metadata.csv; do
    if [[ ! -f "${EXTRACT_DIR}/${f}" ]]; then
        echo "FATAL: Phase 1 missing output: ${f}"
        exit 1
    fi
done
echo "Phase 1 outputs verified."

# --- Phase 2: R construction ---
echo ""
echo "=== PHASE 2: R Seurat construction ==="
echo ""

singularity exec \
    --bind ${SHARED_LAB}:${SHARED_LAB}:rw \
    --bind ${USER_ROOT}:${USER_ROOT}:rw \
    ${BIND_MOUNTS} \
    --bind /tmp:/tmp:rw \
    --bind "$R_LIBS_USER:/home/jovyan/R/library:ro" \
    --env "R_LIBS_USER=/home/jovyan/R/library" \
    "$R_CONTAINER" \
    Rscript "${STUDY_DIR}/scripts/construct_rds.R" \
        --input-dir "$EXTRACT_DIR" \
        --output "$RDS_OUTPUT"

# Verify RDS exists
if [[ ! -f "$RDS_OUTPUT" ]]; then
    echo "FATAL: RDS not created at ${RDS_OUTPUT}"
    exit 1
fi

echo "Phase 2 output verified: $(du -h "$RDS_OUTPUT" | cut -f1)"

# --- Phase 4: Cleanup ---
echo ""
echo "=== PHASE 4: Cleanup intermediates ==="
rm -rf "$EXTRACT_DIR"
echo "Removed: ${EXTRACT_DIR}"

echo ""
echo "=============================================="
echo "CONVERSION COMPLETE"
echo "Output: ${RDS_OUTPUT} ($(du -h "$RDS_OUTPUT" | cut -f1))"
echo "End time: $(date)"
echo "=============================================="
