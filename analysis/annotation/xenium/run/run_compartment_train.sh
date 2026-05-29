#!/bin/bash
#
# 4c-annotation step 9.1: per-compartment Concord retrain (array: epi/stro/imm).
#
# Subsets the joint cohort to one compartment (Epithelial / Stromal / Immune),
# re-trains a compartment-specific Concord embedding for finer L1.5 resolution.
#
# Tier 5: A100 GPU / 8 CPUs / 64G / 4h per task.

#SBATCH --job-name=xen_comp_train
#SBATCH --partition=free-gpu
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --gres=gpu:A100:1
#SBATCH --time=04:00:00
#SBATCH --array=0-2
#SBATCH --output=.slurm_stubs/comp_train_%A_%a.out
#SBATCH --error=.slurm_stubs/comp_train_%A_%a.err

if [ -n "${SLURM_JOB_ID:-}" ]; then SCRIPT_DIR="$(pwd)"; else SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"; fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CONTAINER_TYPE=python_scvi_gpu_2025Q3
source "${PIPELINE_ROOT}/run/_common.sh"

mkdir -p "${CFG_HPC_LOGS_DIR}"
exec >> "${CFG_HPC_LOGS_DIR}/comp_train_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.out" \
    2>> "${CFG_HPC_LOGS_DIR}/comp_train_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.err"

COMPARTMENTS=(Epithelial Stromal Immune)
COMPARTMENT=${COMPARTMENTS[$SLURM_ARRAY_TASK_ID]}
COMP_LC=$(echo "$COMPARTMENT" | tr '[:upper:]' '[:lower:]')
OUT="${CFG_XENIUM_COMPARTMENT_DIR}/${COMP_LC}"
mkdir -p "${OUT}/model"

run_python_singularity \
    "${PIPELINE_ROOT}/4c-annotation/scripts/01_compartment_train.py" \
        --compartment        "${COMPARTMENT}" \
        --flex-dir           "${CFG_FLEX_CONCORD_DIR}" \
        --flex-qc            "${CFG_FLEX_QC}" \
        --xenium-bundle      "${CFG_XENIUM_POOLED_NUCLEAR}" \
        --xenium-filter      "${CFG_XENIUM_THREE_AXIS_FILTER}" \
        --l0p5-xenium        "${CFG_XENIUM_L0P5}" \
        --panel-intersection "${CFG_XENIUM_LEGACY_PANEL_INTERSECTION}/panel_flex_intersection_${COMPARTMENT}.txt" \
        --out-dir            "${OUT}" \
        --n-latent           100

echo "Done comp_train ${COMPARTMENT}: $(date)"
