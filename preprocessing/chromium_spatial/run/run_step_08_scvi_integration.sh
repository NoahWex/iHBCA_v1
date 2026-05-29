#!/bin/bash
#SBATCH --job-name=step_08_scvi
#SBATCH --partition=free-gpu
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=128G
#SBATCH --gres=gpu:1
#SBATCH --time=08:00:00
#SBATCH --output=.slurm_stubs/step_08_%j.out
#SBATCH --error=.slurm_stubs/step_08_%j.err

# Step 08 scVI Cross-Patient Integration (Tier 5 / GPU)
#
# Trains the canonical cross-patient scVI model on the THREE-WAY QC gate
# (step 01 ∩ step 02 ∩ step 03) using step 05 post-QC variable features.
#
# Canonical model reuse:
#   DRY_RUN=1 bash run_step_08_scvi_integration.sh     # dry run only
#   LOAD_MODEL=1 sbatch run_step_08_scvi_integration.sh  # reuse frozen model
# Leave LOAD_MODEL unset for a full retrain.

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
exec >> "${HPC_LOGS_DIR}/step_08_${SLURM_JOB_ID}.out" 2>> "${HPC_LOGS_DIR}/step_08_${SLURM_JOB_ID}.err"
WRAPPER="${CFG_SCRIPTS_DIR}/step_08_scvi_integration/wrapper.py"
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
echo "Output root:      ${CFG_PREPROCESSING_STEP_08}"
echo "Canonical model:  ${CFG_PREPROCESSING_STEP_08_MODEL}"

run_python_singularity "$WRAPPER" \
    --project-root "${CFG_PROJECT_ROOT}" \
    --central-manifest "${CFG_PREPROCESSING_CENTRAL_CELL_STATUS}" \
    --step-05-dir "${CFG_PREPROCESSING_STEP_05}" \
    --raw-data-manifest "${CFG_RAW_DATA_SAMPLE_MANIFEST}" \
    --module-configs "${CFG_PROJECT_ROOT}/config/module_configs.yaml" \
    --output-root "${CFG_PREPROCESSING_STEP_08}" \
    --model-path "${CFG_PREPROCESSING_STEP_08_MODEL}" \
    ${DRY_RUN_FLAG} \
    ${EXTRA_FLAGS}

echo "Exit: $?, End: $(date)"
