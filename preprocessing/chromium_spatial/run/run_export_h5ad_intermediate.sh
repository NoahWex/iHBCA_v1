#!/bin/bash
#SBATCH --job-name=h5ad_export
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --array=0-2
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/h5ad_export_%A_%a.out
#SBATCH --error=.slurm_stubs/h5ad_export_%A_%a.err

# Stage 1 of compartment integration report rendering.
# Reads compartment integrated.h5ad and writes a flat intermediate directory
# (counts.mtx.gz + cells.tsv + genes.tsv + obs.csv + latent.csv + umap.csv +
# manifest.json) that the R Rmd in stage 2 consumes without reticulate.
#
# Container: scgpt_gpu (Python + anndata + scipy)
# Tier:      4cpu/24G/30min — 100K x 18K sparse h5ad load + write fits in <8G
#
# Tasks:
#   0  Epithelial  102,952 cells
#   1  Stromal     119,141 cells
#   2  Immune       47,942 cells

if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CONTAINER_TYPE=scgpt_gpu
source "${PIPELINE_ROOT}/run/_common.sh"


# Redirect stdout/stderr to canonical publication logs location.
# SBATCH directives above land tiny stubs in .slurm_stubs/; real output goes
# to the publication-tier ${HPC_LOGS_DIR}/ which is resolved by _common.sh.
mkdir -p "${HPC_LOGS_DIR}"
exec >> "${HPC_LOGS_DIR}/h5ad_export_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.out" 2>> "${HPC_LOGS_DIR}/h5ad_export_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.err"
COMPARTMENTS=(Epithelial Stromal Immune)
COMPARTMENT=${COMPARTMENTS[$SLURM_ARRAY_TASK_ID]}
CONFIG_LABEL=scvi_n100

INPUT_H5AD="${CFG_UPSTREAM_ROOT}/compartments_step15/${COMPARTMENT}/integration/${CONFIG_LABEL}/integrated.h5ad"
OUTPUT_DIR="${CFG_PROJECT_ROOT}/outputs/integration_intermediate/${COMPARTMENT}/${CONFIG_LABEL}"

mkdir -p "${OUTPUT_DIR}"

echo "Compartment:  ${COMPARTMENT}"
echo "Config:       ${CONFIG_LABEL}"
echo "Input h5ad:   ${INPUT_H5AD}"
echo "Output dir:   ${OUTPUT_DIR}"

run_python_singularity \
    "${CFG_SCRIPTS_DIR}/step_18_integration_preview/08a_export_h5ad_intermediate.py" \
    --h5ad "${INPUT_H5AD}" \
    --output-dir "${OUTPUT_DIR}" \
    --compartment "${COMPARTMENT}" \
    --config-label "${CONFIG_LABEL}" \
    ${DRY_RUN_FLAG}

echo "Exit: $?, End: $(date)"
