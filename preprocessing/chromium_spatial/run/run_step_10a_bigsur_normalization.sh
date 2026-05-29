#!/bin/bash
#SBATCH --job-name=step_10a_bigsur
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=.slurm_stubs/step_10a_%j.out
#SBATCH --error=.slurm_stubs/step_10a_%j.err

# Step 10a Add BigSur Pearson-Residual Normalization (Tier 3-4 / r_spatial)
# Attaches step 05 per-sample BigSur residuals as a new assay on the step 09
# integrated Seurat object and computes 5x4 discordance metrics used by
# step 10b manifold artifact review.

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
exec >> "${HPC_LOGS_DIR}/step_10a_${SLURM_JOB_ID}.out" 2>> "${HPC_LOGS_DIR}/step_10a_${SLURM_JOB_ID}.err"
WRAPPER="${CFG_SCRIPTS_DIR}/step_10a_add_residual_normalized_data/wrapper.R"
if [ ! -f "$WRAPPER" ]; then
    echo "FATAL: wrapper missing: $WRAPPER" >&2
    exit 2
fi

echo "Wrapper:     $WRAPPER"
echo "Step 05 dir: ${CFG_PREPROCESSING_STEP_05}"
echo "Step 09 dir: ${CFG_PREPROCESSING_STEP_09}"
echo "Output root: ${CFG_PREPROCESSING_STEP_10A}"

run_r_singularity "$WRAPPER" \
    --project-root "${CFG_PROJECT_ROOT}" \
    --step-09-dir "${CFG_PREPROCESSING_STEP_09}" \
    --step-05-dir "${CFG_PREPROCESSING_STEP_05}" \
    --raw-data-manifest "${CFG_RAW_DATA_SAMPLE_MANIFEST}" \
    --module-configs "${CFG_PROJECT_ROOT}/config/module_configs.yaml" \
    --output-root "${CFG_PREPROCESSING_STEP_10A}" \
    ${DRY_RUN_FLAG}

echo "Exit: $?, End: $(date)"
