#!/bin/bash
#
# 4b-joint step 2.1: pool per-position counts into cohort matrices.
# Concatenates per-position counts_{nuclear,cytoplasmic}.mtx.gz across the
# 65-sample clean cohort, projects onto the canonical 280-gene panel, and
# applies the Phase 1 noise-floor filter.
# Tier 2: 4 CPUs / 16G / 1h per task; array tasks 1=nuclear, 2=cytoplasmic.

#SBATCH --job-name=xen_pool_clean
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --array=1-2
#SBATCH --output=.slurm_stubs/pool_clean_%A_%a.out
#SBATCH --error=.slurm_stubs/pool_clean_%A_%a.err

if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
CONTAINER_TYPE=python_scvi_gpu_2025Q3
source "${PIPELINE_ROOT}/run/_common.sh"

mkdir -p "${CFG_HPC_LOGS_DIR}"
exec >> "${CFG_HPC_LOGS_DIR}/pool_clean_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.out" \
    2>> "${CFG_HPC_LOGS_DIR}/pool_clean_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.err"

SRCS=(nuclear cytoplasmic)
SRC=${SRCS[$((SLURM_ARRAY_TASK_ID - 1))]}

if [ "$SRC" = "nuclear" ]; then
    OUT="${CFG_XENIUM_POOLED_NUCLEAR}"
else
    OUT="${CFG_XENIUM_POOLED_CYTOPLASMIC}"
fi
mkdir -p "${OUT}"

run_python_singularity \
    "${PIPELINE_ROOT}/4b-joint/scripts/02_joint_embedding/01_pool_xenium_source.py" \
        --source       "${SRC}" \
        --manifest     "${CFG_XENIUM_MANIFEST}" \
        --legacy-root  "${CFG_XENIUM_COUNTS_PER_POSITION}" \
        --v2-root      "${CFG_XENIUM_COUNTS_V2_UCI604}" \
        --v2-patient   UCI604 \
        --phase1-qc    "${CFG_XENIUM_QC}" \
        --panel-genes  "${CFG_XENIUM_PANEL_GENES}" \
        --out-dir      "${OUT}" \
        --apply-filter

echo "Done task=${SLURM_ARRAY_TASK_ID} src=${SRC}: $(date)"
