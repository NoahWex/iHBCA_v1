#!/bin/bash
#SBATCH --job-name=step_10b_validator
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=.slurm_stubs/step_10b_%j.out
#SBATCH --error=.slurm_stubs/step_10b_%j.err

# Step 10b Frozen Validator (Tier 1 / r_spatial)
#
# Per the revised migration model, step 10b is FROZEN: the canonical
# step10b_manifest_update.csv at $CFG_PREPROCESSING_CELL_RETENTION_LIST is
# treated as an input by downstream wrappers. This SLURM wrapper invokes
# the validator wrapper.R, which only confirms the retention list exists
# and has the expected schema. It does NOT render the interactive Rmd.
#
# To manually re-run the artifact-discovery notebook, open RStudio on HPC3
# against source/artifact_discovery.Rmd — not via SLURM.

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
exec >> "${HPC_LOGS_DIR}/step_10b_${SLURM_JOB_ID}.out" 2>> "${HPC_LOGS_DIR}/step_10b_${SLURM_JOB_ID}.err"
WRAPPER="${CFG_SCRIPTS_DIR}/step_10b_interactive_exploration/wrapper.R"
if [ ! -f "$WRAPPER" ]; then
    echo "FATAL: wrapper missing: $WRAPPER" >&2
    exit 2
fi

echo "Wrapper:            $WRAPPER"
echo "Cell retention:     ${CFG_PREPROCESSING_CELL_RETENTION_LIST}"
echo "Output root:        ${CFG_PREPROCESSING_STEP_10B}"
echo "Mode:               frozen validator (no batch compute)"

run_r_singularity "$WRAPPER" \
    --project-root "${CFG_PROJECT_ROOT}" \
    --cell-retention-list "${CFG_PREPROCESSING_CELL_RETENTION_LIST}" \
    --output-root "${CFG_PREPROCESSING_STEP_10B}" \
    ${DRY_RUN_FLAG}

echo "Exit: $?, End: $(date)"
