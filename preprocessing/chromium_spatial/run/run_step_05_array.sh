#!/bin/bash
#SBATCH --job-name=step_05_post_qc_vfs
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:45:00
#SBATCH --array=0-61
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/step_05_%A_%a.out
#SBATCH --error=.slurm_stubs/step_05_%A_%a.err

# Step 05: Post-QC BigSur VFs (per-sample array).
#
# Container:   r_spatial
# Tier:        Tier 2 (smaller than step 01 because the input cell set
#              is already filtered; observed peak ~10G).
# Pattern src: run_step_01_array.sh

set -euo pipefail

CONTAINER_TYPE=r_spatial
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
exec >> "${HPC_LOGS_DIR}/step_05_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.out" 2>> "${HPC_LOGS_DIR}/step_05_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.err"
WRAPPER="${CFG_SCRIPTS_DIR}/step_05_post_qc_vfs/wrapper.R"

echo "Step:          step_05_post_qc_vfs"
echo "Sample index:  ${SLURM_ARRAY_TASK_ID}"
echo "Wrapper:       ${WRAPPER}"
echo "Output root:   ${CFG_PREPROCESSING_STEP_05}"

run_r_singularity "${WRAPPER}" \
    --sample-index "${SLURM_ARRAY_TASK_ID}" \
    --raw-data-manifest "${CFG_RAW_DATA_SAMPLE_MANIFEST}" \
    --central-cell-status "${CFG_PREPROCESSING_CENTRAL_CELL_STATUS}" \
    --step-01-root "${CFG_PREPROCESSING_STEP_01}" \
    --output-root "${CFG_PREPROCESSING_STEP_05}" \
    ${DRY_RUN_FLAG}

echo "Exit: $?, End: $(date)"
