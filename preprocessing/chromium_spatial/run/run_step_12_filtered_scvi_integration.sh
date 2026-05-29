#!/bin/bash
#SBATCH --job-name=step_12_filtered_scvi
#SBATCH --partition=free-gpu
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=128G
#SBATCH --gres=gpu:1
#SBATCH --time=08:00:00
#SBATCH --output=.slurm_stubs/step_12_%j.out
#SBATCH --error=.slurm_stubs/step_12_%j.err

# Step 12 Filtered scVI Cross-Patient Integration (Tier 5 / GPU)
#
# Second-round scVI integration on the FOUR-WAY filtered cell set. Uses
# Step 11 post-artifact VFs and the frozen step 10b retention list.
#
# Canonical model reuse:
#   DRY_RUN=1 bash run_step_12_filtered_scvi_integration.sh
#   LOAD_MODEL=1 sbatch run_step_12_filtered_scvi_integration.sh
# Default retrains; set LOAD_MODEL=1 to reuse $CFG_PREPROCESSING_STEP_12_MODEL.

CONTAINER_TYPE=scgpt_gpu
if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${PIPELINE_ROOT}/run/_common.sh"


# Redirect stdout/stderr to canonical publication logs location.
# SBATCH directives above land tiny stubs in .slurm_stubs/; real output goes
# to the publication-tier ${HPC_LOGS_DIR}/ which is resolved by _common.sh.
mkdir -p "${HPC_LOGS_DIR}"
exec >> "${HPC_LOGS_DIR}/step_12_${SLURM_JOB_ID}.out" 2>> "${HPC_LOGS_DIR}/step_12_${SLURM_JOB_ID}.err"
WRAPPER="${CFG_SCRIPTS_DIR}/step_12_filtered_scvi_integration/wrapper.py"
if [ ! -f "$WRAPPER" ]; then
    echo "FATAL: wrapper missing: $WRAPPER" >&2
    exit 2
fi

EXTRA_FLAGS=""
if [ "${LOAD_MODEL:-0}" = "1" ]; then
    EXTRA_FLAGS="${EXTRA_FLAGS} --load-model"
fi

echo "Wrapper:          $WRAPPER"
echo "Load model:       ${LOAD_MODEL:-0}"
echo "Output root:      ${CFG_PREPROCESSING_STEP_12}"
echo "Canonical model:  ${CFG_PREPROCESSING_STEP_12_MODEL}"
echo "Retention list:   ${CFG_PREPROCESSING_CELL_RETENTION_LIST}"

run_python_singularity "$WRAPPER" \
    --project-root "${CFG_PROJECT_ROOT}" \
    --central-manifest "${CFG_PREPROCESSING_CENTRAL_CELL_STATUS}" \
    --cell-retention-list "${CFG_PREPROCESSING_CELL_RETENTION_LIST}" \
    --step-11-dir "${CFG_PREPROCESSING_STEP_11}" \
    --raw-data-manifest "${CFG_RAW_DATA_SAMPLE_MANIFEST}" \
    --module-configs "${CFG_PROJECT_ROOT}/config/module_configs.yaml" \
    --output-root "${CFG_PREPROCESSING_STEP_12}" \
    --model-path "${CFG_PREPROCESSING_STEP_12_MODEL}" \
    ${DRY_RUN_FLAG} \
    ${EXTRA_FLAGS}

echo "Exit: $?, End: $(date)"
