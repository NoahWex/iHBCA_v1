#!/bin/bash
#SBATCH --job-name=s17_build_merged
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --array=0-3
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/s17_build_merged_%A_%a.out
#SBATCH --error=.slurm_stubs/s17_build_merged_%A_%a.err

# ============================================================================
# step_17 pre-sweep: build per-target merged raw count matrices.
#
# Loads all 62 CellRanger H5 files ONCE per target and saves a merged
# H5AD so each sweep task reads a single file rather than re-loading 62.
#
# Array layout (task → target):
#   0 → full (all 258K step_15-passing cells)
#   1 → epi  (Epithelial compartment)
#   2 → str  (Stromal compartment)
#   3 → imm  (Immune compartment)
#
# MUST complete before submitting run_step_17_sweep.sh.
#
# Resources (Tier 5): 8 CPU / 128G / 2h
#   Outer-join concat of 62 sparse matrices (~258K × 33K genes) peaks ~80G
#   for the full target. Compartment tasks use proportionally less but share
#   the same profile for scheduling simplicity.
#
# To build a single target instead of all four:
#   sbatch --array=3 run_step_17_build_merged.sh   # immune only
#
# Pattern source: run/run_step_17_sweep.sh (SLURM header + _common.sh pattern)
# ============================================================================

set -euo pipefail

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
exec >> "${HPC_LOGS_DIR}/s17_build_merged_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.out" 2>> "${HPC_LOGS_DIR}/s17_build_merged_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.err"
TARGETS=(full epi str imm)
if [ "${SLURM_ARRAY_TASK_ID}" -ge "${#TARGETS[@]}" ]; then
    echo "FATAL: task ${SLURM_ARRAY_TASK_ID} out of range (0-$((${#TARGETS[@]}-1)))" >&2
    exit 2
fi
TARGET="${TARGETS[${SLURM_ARRAY_TASK_ID}]}"

SCRIPT="${CFG_SCRIPTS_DIR}/step_17_integration_sweep/source/build_merged_counts.py"

echo "Task:   ${SLURM_ARRAY_TASK_ID}"
echo "Target: ${TARGET}"

run_python_singularity "${SCRIPT}" \
    --target         "${TARGET}" \
    --project-root   "${CFG_PROJECT_ROOT}" \
    --temp-dir       "${JOB_TEMP_DIR}" \
    ${DRY_RUN_FLAG}

echo "Exit: $?, End: $(date)"
