#!/bin/bash
#SBATCH --job-name=step_13_filtered_preview
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=.slurm_stubs/step_13_%j.out
#SBATCH --error=.slurm_stubs/step_13_%j.err

# Step 13 Filtered Integration Preview (Tier 4 / r_spatial)
# Structural twin of run_step_09; inputs come from step 12 outputs.

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
exec >> "${HPC_LOGS_DIR}/step_13_${SLURM_JOB_ID}.out" 2>> "${HPC_LOGS_DIR}/step_13_${SLURM_JOB_ID}.err"
WRAPPER="${CFG_SCRIPTS_DIR}/step_13_filtered_integration_preview/wrapper.R"
if [ ! -f "$WRAPPER" ]; then
    echo "FATAL: wrapper missing: $WRAPPER" >&2
    exit 2
fi

echo "Wrapper:     $WRAPPER"
echo "Step 12 dir: ${CFG_PREPROCESSING_STEP_12}"
echo "Output root: ${CFG_PREPROCESSING_STEP_13}"

run_r_singularity "$WRAPPER" \
    --project-root "${CFG_PROJECT_ROOT}" \
    --step-12-dir "${CFG_PREPROCESSING_STEP_12}" \
    --module-configs "${CFG_PROJECT_ROOT}/config/module_configs.yaml" \
    --raw-data-manifest "${CFG_RAW_DATA_SAMPLE_MANIFEST}" \
    --output-root "${CFG_PREPROCESSING_STEP_13}" \
    ${DRY_RUN_FLAG}

echo "Exit: $?, End: $(date)"
