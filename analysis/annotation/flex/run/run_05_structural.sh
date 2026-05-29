#!/bin/bash
#SBATCH --job-name=trackc_05_structural
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --array=0-2%3
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/05_structural_%A_%a.out
#SBATCH --error=.slurm_stubs/05_structural_%A_%a.err

# Step 5: pre-review structural diagnostic package per compartment. Renders
# UMAP grids, F1 matrices, and per-cluster heatmaps used during human review
# to draft the per-compartment annotation_v2s yaml.

set -euo pipefail

CONTAINER_TYPE=python_spatial_2025Q2
if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
source "${SCRIPT_DIR}/_common.sh"

mkdir -p "${HPC_LOGS_DIR}"
exec >> "${HPC_LOGS_DIR}/05_structural_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}.out" \
    2>> "${HPC_LOGS_DIR}/05_structural_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}.err"

COMPARTMENTS=(Epithelial Immune Stromal)
COMPARTMENT="${COMPARTMENTS[$SLURM_ARRAY_TASK_ID]}"

RESOLUTIONS=(0.3 0.5 0.8 1.0 1.5 3.0 5.0)

SCRIPT="${CFG_SCRIPTS_DIR}/05_render_structural.py"
BUNDLE_DIR="${CFG_INTEGRATION_INTERMEDIATE}/${COMPARTMENT}/scvi_n100"
CLUSTERS="${CFG_OUTPUTS_ROOT}/clusters/${COMPARTMENT}_clusters.csv"
UMAP="${BUNDLE_DIR}/umap.csv"
UCELL="${CFG_OUTPUTS_ROOT}/ucell/${COMPARTMENT}_ucell_per_cell.csv"
LIMMA_DIR="${CFG_OUTPUTS_ROOT}/limma"
OUT_DIR="${CFG_OUTPUTS_ROOT}/${COMPARTMENT}"

mkdir -p "${OUT_DIR}"

_load_singularity

singularity exec --no-home \
    ${BIND_MOUNTS:-} \
    --bind "${JOB_TEMP_DIR}:${JOB_TEMP_DIR}:rw" \
    --env "PYTHONUNBUFFERED=1" \
    --env "NUMBA_CACHE_DIR=${JOB_TEMP_DIR}/numba_cache" \
    --env "MPLCONFIGDIR=${JOB_TEMP_DIR}/mpl_config" \
    "${CONTAINER_PATH}" \
    python "${SCRIPT}" \
        --compartment "${COMPARTMENT}" \
        --clusters-csv "${CLUSTERS}" \
        --umap-csv "${UMAP}" \
        --ucell-per-cell-csv "${UCELL}" \
        --limma-dir "${LIMMA_DIR}" \
        --resolutions "${RESOLUTIONS[@]}" \
        --out-dir "${OUT_DIR}"

echo "=== Step 5 complete (${COMPARTMENT}) ==="
