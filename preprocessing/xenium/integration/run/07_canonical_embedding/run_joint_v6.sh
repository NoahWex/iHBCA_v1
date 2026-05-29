#!/bin/bash
#
# 4b-joint phase 7: canonical Concord retrain on the filtered cohort (joint_v6).
#
# Re-runs joint FLEX+Xenium Concord training restricted to cells passing the
# three-axis filter (phase 6). Produces the canonical joint_v6 latent space
# used for all downstream annotation, DA, and figure work.
#
# Tier 5: A100 GPU / 8 CPUs / 64G / 4h.

#SBATCH --job-name=xen_joint_v6
#SBATCH --partition=free-gpu
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --gres=gpu:A100:1
#SBATCH --time=04:00:00
#SBATCH --output=.slurm_stubs/joint_v6_%j.out
#SBATCH --error=.slurm_stubs/joint_v6_%j.err

if [ -n "${SLURM_JOB_ID:-}" ]; then SCRIPT_DIR="$(pwd)"; else SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"; fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
CONTAINER_TYPE=python_scvi_gpu_2025Q3
source "${PIPELINE_ROOT}/run/_common.sh"

mkdir -p "${CFG_HPC_LOGS_DIR}"
exec >> "${CFG_HPC_LOGS_DIR}/joint_v6_${SLURM_JOB_ID}.out" 2>> "${CFG_HPC_LOGS_DIR}/joint_v6_${SLURM_JOB_ID}.err"

OUT="${CFG_XENIUM_INTEGRATION_ROOT}/embedding/joint_v6"
mkdir -p "${OUT}/model"

run_python_singularity \
    "${PIPELINE_ROOT}/4b-joint/scripts/07_canonical_embedding/01_joint_train_v2.py" \
        --flex-dir            "${CFG_FLEX_CONCORD_DIR}" \
        --flex-qc             "${CFG_FLEX_QC}" \
        --xenium-bundle       "${CFG_XENIUM_POOLED_NUCLEAR}" \
        --xenium-filter       "${CFG_XENIUM_THREE_AXIS_FILTER}" \
        --panel-intersection  "${CFG_XENIUM_LEGACY_PANEL_INTERSECTION}/panel_flex_intersection_Immune.txt" \
        --out-dir             "${OUT}" \
        --n-latent            100

echo "Done joint_v6: $(date)"
