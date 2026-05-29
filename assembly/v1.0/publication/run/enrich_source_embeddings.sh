#!/bin/bash
#SBATCH --job-name=enrich_emb
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --array=0-6
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=00:45:00

# =============================================================================
# Enrich source h5ads with joint scVI embeddings from integrated object
# =============================================================================
# Ports X_scvi_100 from all-breast-cells.h5ad to each source study.
# Computes per-study UMAP, adds obs['in_ihbca_integrated'] flag.
# Removes X_scVI_joint / X_scVI_native (superseded).
#
# Array index -> study mapping:
#   0=gray  1=kumar  2=murrow  3=nee  4=twigger  5=reed  6=pal
#
# 64GB needed — actual peak 52-64GB across studies (backed mode + obsm operations)
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
REPO_ROOT="${PROJECT_IHBCAV1_UPLOAD}"
SCRIPTS="${REPO_ROOT}/publication/scripts"
CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"

STUDIES=(gray kumar murrow nee twigger reed pal)
FILENAMES=(gray2022.h5ad kumar2023.h5ad murrow2022.h5ad nee2023.h5ad twigger2022.h5ad reed2024.h5ad pal2021.h5ad)
INDEX="${SLURM_ARRAY_TASK_ID:?Must run as SLURM array job}"
STUDY="${STUDIES[$INDEX]}"
FILENAME="${FILENAMES[$INDEX]}"

# Env-var override for alternate builds (e.g., SOURCE_DIR=.../cxg_build/source_datasets)
SOURCE_DIR="${SOURCE_DIR:-${REPO_ROOT}/publication/outputs/source_datasets}"

echo "============================================="
echo "Enrich source embeddings: ${STUDY}"
echo "  Array index: ${INDEX}"
echo "  Source dir: ${SOURCE_DIR}"
echo "  Date: $(date)"
echo "  Node: $(hostname)"
echo "  Job ID: ${SLURM_JOB_ID}"
echo "  Memory: ${SLURM_MEM_PER_NODE}MB"
echo "============================================="

module load singularity

singularity exec \
    --no-mount bind-paths \
    ${BIND_MOUNTS} \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
    "$CONTAINER" \
    python3 "${SCRIPTS}/enrich_source_embeddings.py" \
        --study "$STUDY" \
        --repo-root "$REPO_ROOT" \
        --source "${SOURCE_DIR}/${FILENAME}"

echo ""
echo "Enrichment COMPLETE: ${STUDY}"
echo "Output:"
ls -lh "${SOURCE_DIR}/"*.h5ad 2>/dev/null | grep -i "$STUDY" || echo "  (check output directory)"
