#!/bin/bash
#SBATCH --job-name=step_07_final_report
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/step_07_%j.out
#SBATCH --error=.slurm_stubs/step_07_%j.err

# Step 07: Final HTML report (single job).
#
# Container:   scgpt_gpu (pandas + jinja2 + matplotlib)
# Tier:        Tier 1 (aggregation + rendering; observed peak ~4G).
# Pattern src: run_step_02_array.sh (scaled down to single job)

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
exec >> "${HPC_LOGS_DIR}/step_07_${SLURM_JOB_ID}.out" 2>> "${HPC_LOGS_DIR}/step_07_${SLURM_JOB_ID}.err"
WRAPPER="${CFG_SCRIPTS_DIR}/step_07_final_report/wrapper.py"

echo "Step:          step_07_final_report"
echo "Wrapper:       ${WRAPPER}"
echo "Output root:   ${CFG_PREPROCESSING_STEP_07}"

run_python_singularity "${WRAPPER}" \
    --preprocessing-root "${CFG_PREPROCESSING_ROOT}" \
    --central-cell-status "${CFG_PREPROCESSING_CENTRAL_CELL_STATUS}" \
    --raw-data-manifest "${CFG_RAW_DATA_SAMPLE_MANIFEST}" \
    --output-root "${CFG_PREPROCESSING_STEP_07}" \
    ${DRY_RUN_FLAG}

echo "Exit: $?, End: $(date)"
