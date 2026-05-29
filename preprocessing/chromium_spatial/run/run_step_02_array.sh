#!/bin/bash
#SBATCH --job-name=step_02_cell_filtering
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --array=0-3
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/step_02_%A_%a.out
#SBATCH --error=.slurm_stubs/step_02_%A_%a.err

# Step 02: Per-patient MAD cell filtering (scanpy).
#
# Container:   scgpt_gpu (Jupyter R 4.4.2 Giotto Spatial Python image)
# Tier:        Tier 3 (observed peak 30-50G; 32G budget fits the patient
#              object plus scanpy overhead).
# Pattern src: run_recluster_integration.sh (PATIENT_ID env + python dispatch)

set -euo pipefail

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
exec >> "${HPC_LOGS_DIR}/step_02_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.out" 2>> "${HPC_LOGS_DIR}/step_02_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.err"
PATIENT_ARRAY=(Pat1 Pat2 UCI604 UCI220228)
if [ "${SLURM_ARRAY_TASK_ID}" -ge "${#PATIENT_ARRAY[@]}" ]; then
    echo "FATAL: task id ${SLURM_ARRAY_TASK_ID} out of range (max ${#PATIENT_ARRAY[@]})" >&2
    exit 2
fi
PATIENT_ID="${PATIENT_ARRAY[${SLURM_ARRAY_TASK_ID}]}"
export PATIENT_ID

WRAPPER="${CFG_SCRIPTS_DIR}/step_02_cell_filtering/wrapper.py"

echo "Step:          step_02_cell_filtering"
echo "Patient:       ${PATIENT_ID}"
echo "Wrapper:       ${WRAPPER}"
echo "Output root:   ${CFG_PREPROCESSING_STEP_02}"

run_python_singularity "${WRAPPER}" \
    --patient-id "${PATIENT_ID}" \
    --central-cell-status "${CFG_PREPROCESSING_CENTRAL_CELL_STATUS}" \
    --raw-data-manifest "${CFG_RAW_DATA_SAMPLE_MANIFEST}" \
    --step-01-root "${CFG_PREPROCESSING_STEP_01}" \
    --output-root "${CFG_PREPROCESSING_STEP_02}" \
    ${DRY_RUN_FLAG}

echo "Exit: $?, End: $(date)"
