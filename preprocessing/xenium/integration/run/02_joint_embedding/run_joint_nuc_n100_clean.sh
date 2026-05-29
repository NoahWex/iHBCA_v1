#!/bin/bash
#
# 4b-joint step 2.2: train joint FLEX+Xenium Concord embedding (canonical config).
# Nuclear source, n_latent=100. Selected by FLEX kNN purity benchmark.
# Production training run. Per-compartment Xenium compartment lookup is
# enabled via --compartment-dir (consumes phase 09 outputs once available;
# falls back to "unknown" on initial runs before phase 09 has populated the
# per-compartment cell lists).
# Tier 5: A100 GPU / 8 CPUs / 64G / 4h.

#SBATCH --job-name=xen_joint_nuc_clean
#SBATCH --partition=free-gpu
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --gres=gpu:A100:1
#SBATCH --time=04:00:00
#SBATCH --output=.slurm_stubs/joint_nuc_clean_%j.out
#SBATCH --error=.slurm_stubs/joint_nuc_clean_%j.err

if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
CONTAINER_TYPE=python_scvi_gpu_2025Q3
source "${PIPELINE_ROOT}/run/_common.sh"

mkdir -p "${CFG_HPC_LOGS_DIR}"
exec >> "${CFG_HPC_LOGS_DIR}/joint_nuc_clean_${SLURM_JOB_ID}.out" \
    2>> "${CFG_HPC_LOGS_DIR}/joint_nuc_clean_${SLURM_JOB_ID}.err"

OUT="${CFG_XENIUM_INTEGRATION_ROOT}/embedding/joint_nuc_n100"
mkdir -p "${OUT}/model"

run_python_singularity \
    "${PIPELINE_ROOT}/4b-joint/scripts/02_joint_embedding/02_joint_train.py" \
        --flex-dir            "${CFG_FLEX_CONCORD_DIR}" \
        --flex-qc             "${CFG_FLEX_QC}" \
        --xenium-bundle       "${CFG_XENIUM_POOLED_NUCLEAR}" \
        --source              nuclear \
        --panel-intersection  "${CFG_XENIUM_LEGACY_PANEL_INTERSECTION}/panel_flex_intersection_Immune.txt" \
        --compartment-dir     "${CFG_XENIUM_COMPARTMENT_DIR}" \
        --out-dir             "${OUT}" \
        --n-latent            100

echo "Done joint_nuc_n100_clean: $(date)"
