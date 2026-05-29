#!/bin/bash
#SBATCH --job-name=step_11_post_vfs
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:45:00
#SBATCH --array=0-61
#SBATCH --output=.slurm_stubs/step_11_%A_%a.out
#SBATCH --error=.slurm_stubs/step_11_%A_%a.err

# Step 11 Post-Artifact Filtered VFs (Tier 2 / r_spatial, SLURM array 0-61)
#
# Per-sample recomputation of BigSur VFs on the FOUR-WAY filtered cell
# population. Pattern source: 01_Preprocessing run_step_11_array.sh for
# tier sizing; run_step_01 wrapper (P1) for the CFG_*-driven dispatch.
#
# NOTE: I3 audit confirmed step 11 IS populated in the current canonical
# output dir (285 files, 16.7 GB), contrary to an earlier sidecar_inventory
# entry that claimed step 11 was a no-op. This wrapper is included in the
# V1 dry-run matrix.

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
exec >> "${HPC_LOGS_DIR}/step_11_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.out" 2>> "${HPC_LOGS_DIR}/step_11_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.err"
WRAPPER="${CFG_SCRIPTS_DIR}/step_11_filtered_vfs/wrapper.R"
if [ ! -f "$WRAPPER" ]; then
    echo "FATAL: wrapper missing: $WRAPPER" >&2
    exit 2
fi

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"

echo "Wrapper:            $WRAPPER"
echo "Sample index:       $TASK_ID"
echo "Raw manifest:       ${CFG_RAW_DATA_SAMPLE_MANIFEST}"
echo "Central manifest:   ${CFG_PREPROCESSING_CENTRAL_CELL_STATUS}"
echo "Retention list:     ${CFG_PREPROCESSING_CELL_RETENTION_LIST}"
echo "Step 01 dir:        ${CFG_PREPROCESSING_STEP_01}"
echo "Output root:        ${CFG_PREPROCESSING_STEP_11}"

run_r_singularity "$WRAPPER" \
    --project-root "${CFG_PROJECT_ROOT}" \
    --sample-index "$TASK_ID" \
    --raw-data-manifest "${CFG_RAW_DATA_SAMPLE_MANIFEST}" \
    --central-manifest "${CFG_PREPROCESSING_CENTRAL_CELL_STATUS}" \
    --cell-retention-list "${CFG_PREPROCESSING_CELL_RETENTION_LIST}" \
    --step-01-dir "${CFG_PREPROCESSING_STEP_01}" \
    --module-configs "${CFG_PROJECT_ROOT}/config/module_configs.yaml" \
    --output-root "${CFG_PREPROCESSING_STEP_11}" \
    ${DRY_RUN_FLAG}

echo "Exit: $?, End: $(date)"
