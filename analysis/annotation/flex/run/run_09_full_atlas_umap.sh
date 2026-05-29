#!/bin/bash
#SBATCH --job-name=trackc_09_full_atlas_umap
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/09_full_atlas_umap_%j.out
#SBATCH --error=.slurm_stubs/09_full_atlas_umap_%j.err

# Step 9: project the merged FLEX L2S labels onto the integrated atlas UMAP.
# Sanity panel confirming compartment splits are tight on the joint scvi_n100
# UMAP and lineage neighborhoods emerge spatially as expected.

set -euo pipefail

CONTAINER_TYPE=python_spatial_2025Q2
if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
source "${SCRIPT_DIR}/_common.sh"

mkdir -p "${HPC_LOGS_DIR}"
exec >> "${HPC_LOGS_DIR}/09_full_atlas_umap_${SLURM_JOB_ID:-local}.out" \
    2>> "${HPC_LOGS_DIR}/09_full_atlas_umap_${SLURM_JOB_ID:-local}.err"

SCRIPT="${CFG_SCRIPTS_DIR}/09_render_full_atlas_umap.py"
FULL_BUNDLE="${CFG_INTEGRATION_INTERMEDIATE}/full/scvi_n100"
LABELS_CSV="${CFG_OUTPUTS_ROOT}/flex_l2s_labels.csv"
OUT_DIR="${CFG_OUTPUTS_ROOT}"

if [ ! -f "${FULL_BUNDLE}/umap.csv" ]; then
    echo "ERROR: ${FULL_BUNDLE}/umap.csv missing"
    exit 1
fi
if [ ! -f "${LABELS_CSV}" ]; then
    echo "ERROR: ${LABELS_CSV} missing - run Step 7 (resolve_l2s_labels) first"
    exit 1
fi

_load_singularity

singularity exec --no-home \
    ${BIND_MOUNTS:-} \
    --bind "${JOB_TEMP_DIR}:${JOB_TEMP_DIR}:rw" \
    --env "PYTHONUNBUFFERED=1" \
    --env "NUMBA_CACHE_DIR=${JOB_TEMP_DIR}/numba_cache" \
    --env "MPLCONFIGDIR=${JOB_TEMP_DIR}/mpl_config" \
    "${CONTAINER_PATH}" \
    python "${SCRIPT}" \
        --full-bundle "${FULL_BUNDLE}" \
        --labels-csv "${LABELS_CSV}" \
        --out-dir "${OUT_DIR}"

echo "=== Step 9 complete ==="
