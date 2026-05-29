#!/bin/bash
#
# Step 15 — Milo Contamination (canonical port, frozen-handoff validator)
#
# Step 15 outputs are a frozen handoff in the canonical pipeline. This
# wrapper does not rebuild the Milo graph or re-run DA testing — it simply
# verifies that `contamination_specific_cells.csv` exists and has the
# schema downstream consumers require.
#
# Tier 1 resources: no compute happens here, the job runs a quick path
# and schema check. If an actual re-execution is ever needed, see
# scripts/step_15_milo_contamination/FROZEN_NOTEBOOKS.md and submit the
# original Rmd stages from the dev tree.

#SBATCH --job-name=step_15_milo_contamination
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=.slurm_stubs/step_15_%j.out
#SBATCH --error=.slurm_stubs/step_15_%j.err

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
exec >> "${HPC_LOGS_DIR}/step_15_${SLURM_JOB_ID}.out" 2>> "${HPC_LOGS_DIR}/step_15_${SLURM_JOB_ID}.err"
mkdir -p "${CFG_LOGS_DIR}"

SCRIPT_DIR="${CFG_SCRIPTS_DIR}/step_15_milo_contamination"
WRAPPER="${SCRIPT_DIR}/wrapper.R"

if [ ! -f "${WRAPPER}" ]; then
    echo "FATAL: wrapper not found: ${WRAPPER}" >&2
    exit 2
fi

export CFG_PROJECT_ROOT
export CFG_PREPROCESSING_STEP_15
export CFG_PREPROCESSING_CONTAMINATION_CELLS
export CFG_PREPROCESSING_EPIDERMAL_SAMPLES
export CFG_PREPROCESSING_STEP_14
export PROJECT_ROOT="${CFG_PROJECT_ROOT}"

echo "--- Step 15 frozen-handoff validator ---"
echo "Wrapper:               ${WRAPPER}"
echo "Contamination CSV:     ${CFG_PREPROCESSING_CONTAMINATION_CELLS}"
echo "Epidermal samples YAML:${CFG_PREPROCESSING_EPIDERMAL_SAMPLES}"
echo "-----------------------------------------"

# Dry-run: invoke wrapper inside the container so Rscript resolves on SLURM
# compute nodes (where Rscript is not on the default PATH). CP5.7.9 pattern.
if [ "${DRY_RUN}" = "1" ]; then
    run_r_singularity "${WRAPPER}" --dry-run \
        --project-root "${CFG_PROJECT_ROOT}" \
        --step-15-dir "${CFG_PREPROCESSING_STEP_15}"
    exit $?
fi

run_r_singularity "${WRAPPER}" \
    --project-root "${CFG_PROJECT_ROOT}" \
    --step-15-dir "${CFG_PREPROCESSING_STEP_15}"
