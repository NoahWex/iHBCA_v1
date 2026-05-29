#!/bin/bash
#SBATCH --job-name=step_04_classify_singler
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --array=0-61
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/step_04_%A_%a.out
#SBATCH --error=.slurm_stubs/step_04_%A_%a.err

# Step 04 per-sample classify: run 03b_classify_singler.R for one sample.
# Submit AFTER run_step_04_train.sh completes (models must exist first).
#
# Required env vars:
#   SC_REFERENCE_DIR  — Kumar SC reference dir (MTX format); used to resolve
#                       h5_path via the raw data manifest
# Resources: Tier 2 (8 CPU / 64G / 4h, array of 62 samples)
#
# Sequence:
#   1. sbatch run_step_04_train.sh          (one-shot, Tier 5, ~4h)
#   2. sbatch run_step_04_array.sh          (array 0-61, Tier 2, ~30min each)
#   3. sbatch run_step_04_harmonize.sh      (one-shot, Tier 2, LEIDEN_CSV=...)

set -euo pipefail

if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CONTAINER_TYPE=r_spatial
source "${PIPELINE_ROOT}/run/_common.sh"


# Redirect stdout/stderr to canonical publication logs location.
# SBATCH directives above land tiny stubs in .slurm_stubs/; real output goes
# to the publication-tier ${HPC_LOGS_DIR}/ which is resolved by _common.sh.
mkdir -p "${HPC_LOGS_DIR}"
exec >> "${HPC_LOGS_DIR}/step_04_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.out" 2>> "${HPC_LOGS_DIR}/step_04_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.err"
SCRIPT="${CFG_SCRIPTS_DIR}/step_04_label_transfer/03b_classify_singler.R"
MODELS_DIR="${CFG_PREPROCESSING_STEP_04}/models/singler_models"

# Resolve sample_id and h5_path from the raw data manifest by task index
MANIFEST="${CFG_RAW_DATA_SAMPLE_MANIFEST}"
SAMPLE_LINE=$(awk -v idx="${SLURM_ARRAY_TASK_ID}" 'NR == idx+2 {print}' "${MANIFEST}")
SAMPLE_ID=$(echo "${SAMPLE_LINE}" | cut -f1)
H5_PATH=$(echo "${SAMPLE_LINE}" | cut -f2)

if [ -z "$SAMPLE_ID" ] || [ -z "$H5_PATH" ]; then
  echo "FATAL: Could not resolve sample for task ${SLURM_ARRAY_TASK_ID} from ${MANIFEST}" >&2
  exit 1
fi

echo "Step:         step_04_classify_singler"
echo "Sample index: ${SLURM_ARRAY_TASK_ID}"
echo "Sample ID:    ${SAMPLE_ID}"
echo "H5 path:      ${H5_PATH}"
echo "Models dir:   ${MODELS_DIR}"

run_r_singularity "${SCRIPT}" \
    --sample_id   "${SAMPLE_ID}" \
    --h5_path     "${H5_PATH}" \
    --models_dir  "${MODELS_DIR}" \
    --output_dir  "${CFG_PREPROCESSING_STEP_04}/${SAMPLE_ID}" \
    ${DRY_RUN_FLAG}

echo "Exit: $?, End: $(date)"
