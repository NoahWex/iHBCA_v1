#!/bin/bash
#
# Step 18 — DESeq2 pseudobulk per cluster resolution (compartment x resolution array).
#
# Pattern src: run_step_19b_integration_preview.sh (compartment array)
#
# Design: ~ patient_id + position + cluster (LRT vs ~ patient_id + position)
# Per-cluster LFC via lfcShrink(normal). Filters: min 10 cells per
# (cluster x sample) group; clusters present in < 2 patients excluded.
#
# 15 tasks: 3 compartments x 5 resolutions, fully parallel.
#
# Task mapping (COMPARTMENT_IDX * 5 + RESOLUTION_IDX):
#   0-4:   Epithelial @ 0.3, 0.5, 0.8, 1.0, 5.0
#   5-9:   Stromal    @ 0.3, 0.5, 0.8, 1.0, 5.0
#   10-14: Immune     @ 0.3, 0.5, 0.8, 1.0, 5.0
#
# Tier 3: one resolution per task, DESeq2 LRT + normal shrinkage.

#SBATCH --job-name=step_18_deseq
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --array=0-14
#SBATCH --output=.slurm_stubs/step_18_deseq_%A_%a.out
#SBATCH --error=.slurm_stubs/step_18_deseq_%A_%a.err

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
exec >> "${HPC_LOGS_DIR}/step_18_deseq_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.out" 2>> "${HPC_LOGS_DIR}/step_18_deseq_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.err"
COMPARTMENTS=(Epithelial Stromal Immune)
RESOLUTIONS=(0.3 0.5 0.8 1.0 5.0)

COMPARTMENT=${COMPARTMENTS[$((SLURM_ARRAY_TASK_ID / 5))]}
RESOLUTION=${RESOLUTIONS[$((SLURM_ARRAY_TASK_ID % 5))]}

SCRIPT="${CFG_SCRIPTS_DIR}/step_18_deseq/18_deseq_per_resolution.R"
INTEGRATION_DIR="${CFG_INTEGRATION_INTERMEDIATE}/${COMPARTMENT}/scvi_n100"
OUTDIR="${CFG_PREPROCESSING_STEP_18_DESEQ}/${COMPARTMENT}"

mkdir -p "${OUTDIR}"

echo "--- Step 18 DESeq2 Pseudobulk ---"
echo "Task:             ${SLURM_ARRAY_TASK_ID}"
echo "Compartment:      ${COMPARTMENT}"
echo "Resolution:       ${RESOLUTION}"
echo "Integration dir:  ${INTEGRATION_DIR}"
echo "Output dir:       ${OUTDIR}"
echo "---------------------------------"

run_r_singularity "${SCRIPT}" \
    --compartment     "${COMPARTMENT}" \
    --integration-dir "${INTEGRATION_DIR}" \
    --out-dir         "${OUTDIR}" \
    --resolution      "${RESOLUTION}" \
    --min-cells       10

echo "Exit: $?, End: $(date)"
