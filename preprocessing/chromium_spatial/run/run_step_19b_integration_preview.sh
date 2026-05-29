#!/bin/bash
#
# Step 19b — Integration Preview (per-compartment array).
#
# Tier 4: Seurat load + winner UMAP injection + marker analysis + Rmd render
# per compartment. Runs after the step 17 sweep winner sidecars exist.
#
# Tasks:
#   0  Epithelial
#   1  Stromal
#   2  Immune

#SBATCH --job-name=step_19b_integration_preview
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --array=0-2
#SBATCH --output=.slurm_stubs/step_19b_%A_%a.out
#SBATCH --error=.slurm_stubs/step_19b_%A_%a.err

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
exec >> "${HPC_LOGS_DIR}/step_19b_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.out" 2>> "${HPC_LOGS_DIR}/step_19b_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.err"
COMPARTMENTS=(Epithelial Stromal Immune)
COMPARTMENT=${COMPARTMENTS[$SLURM_ARRAY_TASK_ID]}

SCRIPT_DIR="${CFG_SCRIPTS_DIR}/step_19b_integration_preview"
WRAPPER="${SCRIPT_DIR}/wrapper.R"
RMD="${SCRIPT_DIR}/integration_preview.Rmd"

OUTDIR="${CFG_PREPROCESSING_STEP_19B}/${COMPARTMENT}"
OUTPUT_FILE="integration_preview_${COMPARTMENT}.html"

mkdir -p "${OUTDIR}/reports"

export CFG_PROJECT_ROOT
export CFG_PREPROCESSING_STEP_19B
export CFG_CANONICAL_SWEEP_WINNER
export PROJECT_ROOT="${CFG_PROJECT_ROOT}"

echo "--- Step 19b Integration Preview ---"
echo "Compartment:     ${COMPARTMENT}"
echo "Rmd:             ${RMD}"
echo "Sweep winner:    ${CFG_CANONICAL_SWEEP_WINNER}"
echo "Output dir:      ${OUTDIR}"
echo "------------------------------------"

if [ "${DRY_RUN}" = "1" ]; then
    run_r_singularity "${WRAPPER}" --dry-run \
        --project-root "${CFG_PROJECT_ROOT}" \
        --step-19b-dir "${CFG_PREPROCESSING_STEP_19B}" \
        --sweep-winner "${CFG_CANONICAL_SWEEP_WINNER}" \
        --step-14-dir "${CFG_PREPROCESSING_STEP_14}" \
        --step-15-contamination "${CFG_PREPROCESSING_CONTAMINATION_CELLS}"
    exit $?
fi

run_r_singularity -e "rmarkdown::render('${RMD}', \
    output_dir = '${OUTDIR}/reports', \
    output_file = '${OUTPUT_FILE}', \
    intermediates_dir = '${JOB_TEMP_DIR}', \
    knit_root_dir = '${CFG_PROJECT_ROOT}', \
    params = list(PROJECT_ROOT = '${CFG_PROJECT_ROOT}', COMPARTMENT = '${COMPARTMENT}'))"

echo "Exit: $?, End: $(date)"
