#!/bin/bash
#SBATCH --job-name=reassemble_src
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --array=0-6
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=02:00:00

# =============================================================================
# Phase 3: Re-assemble all 7 source h5ads after validation fixes
# =============================================================================
# Fixes applied: A1-A4, A5-partial, A6-A10 (script changes),
#                A8 (CxG raw counts for gray/twigger),
#                A12 (v24+v32 gene mapping for murrow/nee/pal)
#
# Array index → study mapping:
#   0=gray  1=kumar  2=murrow  3=nee  4=twigger  5=reed  6=pal
#
# 128GB for all tasks (reed needs ~76GB, others less but uniform is simpler)
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
REPO_ROOT="${PROJECT_IHBCAV1_UPLOAD}"
source "${REPO_ROOT}/publication/run/load_config.sh"
SCRIPTS="${REPO_ROOT}/publication/scripts"
CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"

# Study and filename arrays from pipeline.yaml (via load_config.sh)
INDEX="${SLURM_ARRAY_TASK_ID:?Must run as SLURM array job}"
STUDY="${STUDIES_ARRAY[$INDEX]}"
FILENAME="${FILENAMES_ARRAY[$INDEX]}"

# Env-var overrides for alternate builds (e.g., TARGET=cxg OUTPUT_DIR=.../cxg_build/)
OUTPUT_DIR="${OUTPUT_DIR:-${OUTPUT_SOURCE_DATASETS}}"

echo "============================================="
echo "Re-assemble source h5ad: ${STUDY}"
echo "  Array index: ${INDEX}"
echo "  Target: ${TARGET:-hca}"
echo "  Output dir: ${OUTPUT_DIR}"
echo "  Date: $(date)"
echo "  Node: $(hostname)"
echo "  Job ID: ${SLURM_JOB_ID}"
echo "  Memory: ${SLURM_MEM_PER_NODE}MB"
echo "============================================="

mkdir -p "$OUTPUT_DIR"

module load singularity

singularity exec \
    --no-mount bind-paths \
    ${BIND_MOUNTS} \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
    "$CONTAINER" \
    python3 "${SCRIPTS}/assemble_h5ad.py" \
        --study "$STUDY" \
        --repo-root "$REPO_ROOT" \
        --target "${TARGET:-hca}" \
        --output "${OUTPUT_DIR}/${FILENAME}"

echo ""
echo "Re-assembly COMPLETE: ${STUDY}"
echo "Output:"
ls -lh "${OUTPUT_DIR}/"*.h5ad 2>/dev/null | grep -i "$STUDY" || echo "  (check output directory)"
