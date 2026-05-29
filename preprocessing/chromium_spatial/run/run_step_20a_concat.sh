#!/bin/bash
#
# Step 20a — Assemble full-object integration intermediate.
#
# Pattern src: run_step_20_finalize.sh (single-task, lightweight Python)
#
# Concatenates per-compartment counts matrices + obs metadata and pulls
# full-object latent/UMAP from step 17 sweep outputs.
# Run before step 20 (manifest builder).
#
# Tier 3: 270K cells x 18K genes sparse concat fits comfortably in 32G.

#SBATCH --job-name=step_20a_concat
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=.slurm_stubs/step_20a_concat_%j.out
#SBATCH --error=.slurm_stubs/step_20a_concat_%j.err

if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CONTAINER_TYPE=scgpt_gpu
source "${PIPELINE_ROOT}/run/_common.sh"


# Redirect stdout/stderr to canonical publication logs location.
# SBATCH directives above land tiny stubs in .slurm_stubs/; real output goes
# to the publication-tier ${HPC_LOGS_DIR}/ which is resolved by _common.sh.
mkdir -p "${HPC_LOGS_DIR}"
exec >> "${HPC_LOGS_DIR}/step_20a_concat_${SLURM_JOB_ID}.out" 2>> "${HPC_LOGS_DIR}/step_20a_concat_${SLURM_JOB_ID}.err"
SCRIPT="${CFG_SCRIPTS_DIR}/step_20_finalize/20a_concat_full_object.py"
STEP17_ROOT="${CFG_PREPROCESSING_STEP_17}"
OUT_DIR="${CFG_INTEGRATION_INTERMEDIATE}/full/scvi_n100"

mkdir -p "${OUT_DIR}"

echo "--- Step 20a Full Object Concat ---"
echo "Integration root: ${CFG_INTEGRATION_INTERMEDIATE}"
echo "Step 17 root:     ${STEP17_ROOT}"
echo "Output dir:       ${OUT_DIR}"
echo "-----------------------------------"

run_python_singularity python "${SCRIPT}" \
    --integration-root "${CFG_INTEGRATION_INTERMEDIATE}" \
    --step17-root      "${STEP17_ROOT}" \
    --out-dir          "${OUT_DIR}"

echo "Exit: $?, End: $(date)"
