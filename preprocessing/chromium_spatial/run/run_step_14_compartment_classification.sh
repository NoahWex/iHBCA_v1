#!/bin/bash
#
# Step 14 — Compartment Classification (canonical port)
#
# Tier 4: single-node Seurat reassembly + presto markers + Rmd render.
# Compute profile matches the dev-tree wrapper (6 cpus / 64G / 2h) and
# leaves headroom for the marker-discovery chunks on the full atlas.

#SBATCH --job-name=step_14_compartment_classification
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=.slurm_stubs/step_14_%j.out
#SBATCH --error=.slurm_stubs/step_14_%j.err

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
exec >> "${HPC_LOGS_DIR}/step_14_${SLURM_JOB_ID}.out" 2>> "${HPC_LOGS_DIR}/step_14_${SLURM_JOB_ID}.err"
# Canonical SLURM log destination (overrides the SBATCH directives above
# when CFG_LOGS_DIR differs from the submission cwd). _common.sh expands
# CFG_LOGS_DIR from paths.yaml.
mkdir -p "${CFG_LOGS_DIR}"

# Resolve absolute paths for wrapper + Rmd
SCRIPT_DIR="${CFG_SCRIPTS_DIR}/step_14_compartment_classification"
WRAPPER="${SCRIPT_DIR}/wrapper.R"

if [ ! -f "${WRAPPER}" ]; then
    echo "FATAL: wrapper not found: ${WRAPPER}" >&2
    exit 2
fi

# Export canonical CFG_* so the R wrapper reads them directly
export CFG_PROJECT_ROOT
export CFG_PREPROCESSING_STEP_14
export CFG_PREPROCESSING_STEP_12
export CFG_PREPROCESSING_STEP_08
export CFG_PREPROCESSING_STEP_04
# PROJECT_ROOT kept as an alias for the Rmd's `params$project_root` default
export PROJECT_ROOT="${CFG_PROJECT_ROOT}"

echo "--- Step 14 launch ---"
echo "Wrapper:          ${WRAPPER}"
echo "Output root:      ${CFG_PREPROCESSING_STEP_14}"
echo "L2 primary input: ${CFG_PREPROCESSING_STEP_12}"
echo "L2 fallback:      ${CFG_PREPROCESSING_STEP_08}"
echo "Step 04 labels:   ${CFG_PREPROCESSING_STEP_04}"
echo "----------------------"

# Dry-run: invoke wrapper inside the container so Rscript resolves on SLURM
# compute nodes (where Rscript is not on the default PATH). CP5.7.9 pattern.
if [ "${DRY_RUN}" = "1" ]; then
    run_r_singularity "${WRAPPER}" --dry-run \
        --project-root "${CFG_PROJECT_ROOT}" \
        --step-14-dir "${CFG_PREPROCESSING_STEP_14}"
    exit $?
fi

run_r_singularity "${WRAPPER}" \
    --project-root "${CFG_PROJECT_ROOT}" \
    --step-14-dir "${CFG_PREPROCESSING_STEP_14}"
