#!/bin/bash
#
# Step 16 — BigSur VFs (FIVE-WAY filtered) — SLURM array, one task per sample.
#
# Tier 2 per task (empirical peak ~13 GB for BigSur on ~5K cells per
# sample). Array bounds 0-61 match the 62 samples carried through the
# FIVE-WAY filter. Step 16 also performs an in-place atomic merge of the
# step_15_pass column into central_cell_status.csv the first time any
# array task runs — subsequent tasks see the column and skip the merge.

#SBATCH --job-name=step_16_bigsur_vfs
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:45:00
#SBATCH --array=0-61
#SBATCH --output=.slurm_stubs/step_16_%A_%a.out
#SBATCH --error=.slurm_stubs/step_16_%A_%a.err

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
exec >> "${HPC_LOGS_DIR}/step_16_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.out" 2>> "${HPC_LOGS_DIR}/step_16_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.err"
mkdir -p "${CFG_LOGS_DIR}"

SCRIPT_DIR="${CFG_SCRIPTS_DIR}/step_16_bigsur_vfs"
WRAPPER="${SCRIPT_DIR}/wrapper.R"

if [ ! -f "${WRAPPER}" ]; then
    echo "FATAL: wrapper not found: ${WRAPPER}" >&2
    exit 2
fi

export CFG_PROJECT_ROOT
export CFG_PREPROCESSING_STEP_16
export CFG_PREPROCESSING_STEP16_VF_DIR
export CFG_PREPROCESSING_STEP_01
export CFG_PREPROCESSING_CENTRAL_CELL_STATUS
export CFG_PREPROCESSING_CONTAMINATION_CELLS
export CFG_RAW_DATA_ROOT
export CFG_RAW_DATA_SAMPLE_MANIFEST
export PROJECT_ROOT="${CFG_PROJECT_ROOT}"

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"

echo "--- Step 16 BigSur VFs (task ${TASK_ID}) ---"
echo "Wrapper:           ${WRAPPER}"
echo "Sample index:      ${TASK_ID}"
echo "Step 16 VF dir:    ${CFG_PREPROCESSING_STEP16_VF_DIR}"
echo "Central manifest:  ${CFG_PREPROCESSING_CENTRAL_CELL_STATUS}"
echo "--------------------------------------------"

if [ "${DRY_RUN}" = "1" ]; then
    # CP5.7.9 pattern: invoke wrapper inside the container so Rscript resolves
    # on SLURM compute nodes (where Rscript is not on the default PATH).
    run_r_singularity "${WRAPPER}" --dry-run \
        --project-root "${CFG_PROJECT_ROOT}" \
        --step-16-dir "${CFG_PREPROCESSING_STEP_16}" \
        --step-01-dir "${CFG_PREPROCESSING_STEP_01}" \
        --central-manifest "${CFG_PREPROCESSING_CENTRAL_CELL_STATUS}" \
        --contamination-cells "${CFG_PREPROCESSING_CONTAMINATION_CELLS}" \
        --raw-data-root "${CFG_RAW_DATA_ROOT}"
    exit $?
fi

run_r_singularity "${WRAPPER}" \
    --sample-index "${TASK_ID}" \
    --project-root "${CFG_PROJECT_ROOT}" \
    --step-16-dir "${CFG_PREPROCESSING_STEP_16}" \
    --step-01-dir "${CFG_PREPROCESSING_STEP_01}" \
    --central-manifest "${CFG_PREPROCESSING_CENTRAL_CELL_STATUS}" \
    --contamination-cells "${CFG_PREPROCESSING_CONTAMINATION_CELLS}" \
    --raw-data-root "${CFG_RAW_DATA_ROOT}"
