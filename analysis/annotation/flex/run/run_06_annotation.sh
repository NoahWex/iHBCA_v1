#!/bin/bash
#SBATCH --job-name=trackc_06_annotation
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --array=0-2%3
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/06_annotation_%A_%a.out
#SBATCH --error=.slurm_stubs/06_annotation_%A_%a.err

# Step 6: post-review annotation evidence package per compartment. Consumes
# the per-compartment annotation_yamls/{compartment}.yaml decision, the
# multi-resolution clusters table, and the integration bundle to render
# label-level heatmaps, label_summary.csv, and label_umap.pdf.

set -euo pipefail

CONTAINER_TYPE=python_spatial_2025Q2
if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
source "${SCRIPT_DIR}/_common.sh"

mkdir -p "${HPC_LOGS_DIR}"
exec >> "${HPC_LOGS_DIR}/06_annotation_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}.out" \
    2>> "${HPC_LOGS_DIR}/06_annotation_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}.err"

COMPARTMENTS=(Epithelial Immune Stromal)
COMPARTMENT="${COMPARTMENTS[$SLURM_ARRAY_TASK_ID]}"

case "${COMPARTMENT}" in
    Epithelial) YAML_SHORT=epi ; YAML_LOWER=epithelial ;;
    Immune)     YAML_SHORT=imm ; YAML_LOWER=immune ;;
    Stromal)    YAML_SHORT=str ; YAML_LOWER=stromal ;;
esac

SCRIPT="${CFG_SCRIPTS_DIR}/06_render_annotation.py"
OUT_DIR="${CFG_OUTPUTS_ROOT}/${COMPARTMENT}"
YAML_DRAFT="${CFG_ANNOTATION_YAMLS}/${YAML_LOWER}.yaml"
V1_YAML="${CFG_V1_ANNOTATION_YAMLS}/annotation_v2_${YAML_SHORT}.yaml"
CLUSTERS="${CFG_OUTPUTS_ROOT}/clusters/${COMPARTMENT}_clusters.csv"
BUNDLE_DIR="${CFG_INTEGRATION_INTERMEDIATE}/${COMPARTMENT}/scvi_n100"
LIMMA_DIR="${CFG_OUTPUTS_ROOT}/limma"

if [ ! -f "${YAML_DRAFT}" ]; then
    echo "ERROR: per-compartment yaml not found: ${YAML_DRAFT}"
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
        --bundle-dir "${BUNDLE_DIR}" \
        --clusters-csv "${CLUSTERS}" \
        --yaml-draft "${YAML_DRAFT}" \
        --v1-yaml "${V1_YAML}" \
        --limma-dir "${LIMMA_DIR}" \
        --compartment "${COMPARTMENT}" \
        --out-dir "${OUT_DIR}"

echo "=== Step 6 complete (${COMPARTMENT}) ==="
