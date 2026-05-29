#!/bin/bash
#
# 4b-joint phase 4 (step a): kNN vote for L0.5 transfer to Xenium.
#
# For each Xenium cell, finds k=15 nearest FLEX neighbors in the joint latent
# space and votes L0.5 label by majority. Produces per-Xenium-cell L0.5
# prediction with confidence (vote fraction).
#
# Tier 3: 4 CPUs / 24G / 1h.

#SBATCH --job-name=xen_p4_knn
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=.slurm_stubs/p4_knn_%j.out
#SBATCH --error=.slurm_stubs/p4_knn_%j.err

if [ -n "${SLURM_JOB_ID:-}" ]; then SCRIPT_DIR="$(pwd)"; else SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"; fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
CONTAINER_TYPE=python_scvi_gpu_2025Q3
source "${PIPELINE_ROOT}/run/_common.sh"

mkdir -p "${CFG_HPC_LOGS_DIR}"
exec >> "${CFG_HPC_LOGS_DIR}/p4_knn_${SLURM_JOB_ID}.out" 2>> "${CFG_HPC_LOGS_DIR}/p4_knn_${SLURM_JOB_ID}.err"

L0P5="${CFG_XENIUM_LEGACY_PANEL_INTERSECTION}/../l0p5_per_cell.csv"
OUT_DIR="$(dirname "${CFG_XENIUM_L0P5}")"
mkdir -p "${OUT_DIR}"

run_python_singularity \
    "${PIPELINE_ROOT}/4b-joint/scripts/04_label_transfer/01_knn_vote_joint.py" \
        --latent   "${CFG_XENIUM_JOINT_LATENT}" \
        --obs      "${CFG_XENIUM_JOINT_OBS}" \
        --l0p5-csv "${L0P5}" \
        --out-csv  "${OUT_DIR}/xenium_knn_vote_clean.csv" \
        -k 30

echo "Done phase4_knn: $(date)"
