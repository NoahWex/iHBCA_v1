#!/bin/bash
#SBATCH --job-name=step_01_baseline_vfs
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --array=0-62
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/step_01_%A_%a.out
#SBATCH --error=.slurm_stubs/step_01_%A_%a.err

# Step 01: Baseline VF / BigSur variance selection (per-sample array).
#
# Container:   r_spatial (Rocky8 R4.3.3 Spatial image with Seurat + BigSur)
# Tier:        Tier 2 (observed peak ~12G; budget 24G gives 2x headroom).
# Pattern src: run_singler_l1.sh (argparse + _common.sh container dispatch)
#
# Array: SLURM_ARRAY_TASK_ID is the 0-based row index into
#        $CFG_RAW_DATA_SAMPLE_MANIFEST. Positions with no chromium data
#        (e.g. Pat2_P1) are handled by the wrapper's graceful-skip path.

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
exec >> "${HPC_LOGS_DIR}/step_01_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.out" 2>> "${HPC_LOGS_DIR}/step_01_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.err"
WRAPPER="${CFG_SCRIPTS_DIR}/step_01_baseline_vfs/wrapper.R"

if [ ! -f "${WRAPPER}" ]; then
    echo "FATAL: wrapper missing: ${WRAPPER}" >&2
    exit 2
fi

echo "Step:          step_01_baseline_vfs"
echo "Wrapper:       ${WRAPPER}"
echo "Sample index:  ${SLURM_ARRAY_TASK_ID}"
echo "Output root:   ${CFG_PREPROCESSING_STEP_01}"

run_r_singularity "${WRAPPER}" \
    --sample-index "${SLURM_ARRAY_TASK_ID}" \
    --raw-data-manifest "${CFG_RAW_DATA_SAMPLE_MANIFEST}" \
    --output-root "${CFG_PREPROCESSING_STEP_01}" \
    --qc-status-dir "${CFG_PREPROCESSING_QC_STATUS}" \
    ${DRY_RUN_FLAG}

echo "Exit: $?, End: $(date)"
