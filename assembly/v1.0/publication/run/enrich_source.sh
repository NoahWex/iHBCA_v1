#!/bin/bash
#SBATCH --job-name=enrich_source
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --array=0-6
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=02:00:00

# =============================================================================
# Enrich source h5ads with var annotations, uns documentation, obs polish
# =============================================================================
# Array index -> study mapping:
#   0=gray  1=kumar  2=murrow  3=nee  4=twigger  5=reed  6=pal
#
# 128GB for all tasks (reed is 803K cells, needs headroom)
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

# Env-var override for alternate builds (e.g., SOURCE_DIR=.../cxg_build/source_datasets)
SOURCE_DIR="${SOURCE_DIR:-${OUTPUT_SOURCE_DATASETS}}"
H5AD="${SOURCE_DIR}/${FILENAME}"

echo "============================================="
echo "Enrich source h5ad: ${STUDY}"
echo "  File: ${FILENAME}"
echo "  Array index: ${INDEX}"
echo "  Date: $(date)"
echo "  Node: $(hostname)"
echo "  Job ID: ${SLURM_JOB_ID}"
echo "  Memory: ${SLURM_MEM_PER_NODE}MB"
echo "============================================="

if [ ! -f "${H5AD}" ]; then
    echo "ERROR: h5ad not found: ${H5AD}"
    exit 1
fi

echo "Input size: $(ls -lh "${H5AD}" | awk '{print $5}')"

module load singularity

singularity exec \
    --no-mount bind-paths \
    ${BIND_MOUNTS} \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
    "$CONTAINER" \
    python3 "${SCRIPTS}/enrich_h5ads.py" \
        --h5ad "$H5AD" \
        --mode source \
        --study "$STUDY" \
        --repo-root "$REPO_ROOT" \
        --skip-layers

echo ""
echo "Enrichment COMPLETE: ${STUDY}"
echo "Output size: $(ls -lh "${H5AD}" | awk '{print $5}')"
