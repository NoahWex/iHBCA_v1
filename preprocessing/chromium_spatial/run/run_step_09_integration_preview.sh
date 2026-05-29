#!/bin/bash
#SBATCH --job-name=step_09_preview
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=.slurm_stubs/step_09_%j.out
#SBATCH --error=.slurm_stubs/step_09_%j.err

# Step 09 Integration Preview (Tier 4 / r_spatial)
# Renders the step 08 QC + marker preview report. wrapper.R invokes
# rmarkdown::render itself, so the wrapper is called directly rather
# than via a helper.

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
exec >> "${HPC_LOGS_DIR}/step_09_${SLURM_JOB_ID}.out" 2>> "${HPC_LOGS_DIR}/step_09_${SLURM_JOB_ID}.err"
WRAPPER="${CFG_SCRIPTS_DIR}/step_09_integration_preview/wrapper.R"
if [ ! -f "$WRAPPER" ]; then
    echo "FATAL: wrapper missing: $WRAPPER" >&2
    exit 2
fi

echo "Wrapper:     $WRAPPER"
echo "Step 08 dir: ${CFG_PREPROCESSING_STEP_08}"
echo "Output root: ${CFG_PREPROCESSING_STEP_09}"

run_r_singularity "$WRAPPER" \
    --project-root "${CFG_PROJECT_ROOT}" \
    --step-08-dir "${CFG_PREPROCESSING_STEP_08}" \
    --central-manifest "${CFG_PREPROCESSING_CENTRAL_CELL_STATUS}" \
    --module-configs "${CFG_PROJECT_ROOT}/config/module_configs.yaml" \
    --raw-data-manifest "${CFG_RAW_DATA_SAMPLE_MANIFEST}" \
    --output-root "${CFG_PREPROCESSING_STEP_09}" \
    ${DRY_RUN_FLAG}

echo "Exit: $?, End: $(date)"
