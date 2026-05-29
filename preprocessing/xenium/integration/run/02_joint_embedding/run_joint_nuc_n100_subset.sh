#!/bin/bash
#
# 4b-joint phase 2.2 — wrapper-resolution subset gate.
#
# CP3 subset test: validates the post-CP5.7.10 wrapper stack and the patched
# `--compartment-dir` flag flow on `02_joint_train.py` (intake gap #7,
# unexercised at promotion time).
#
# Tier 1: 4 CPUs / 16G / 30min, NO GPU. Validates: (a) walk-up `_common.sh`
# resolution under SLURM, (b) `--compartment-dir` flag is parsed by the
# patched script, (c) per-compartment leiden_assignments.csv lookup attempt
# fires in the script (logged as "Compartment lookup: ..." or "missing
# {path}" lines in stdout), (d) script reaches the Concord training step
# (full training would time out at 30min — that's acceptable for the
# wrapper-resolution gate; the compartment-dir code path fires before
# Concord training begins).
#
# Per coordinator's CP3 disposition (2026-04-28): the patched compartment
# lookup SEMANTIC test (full subset training with --n-cells-subsample +
# --skip-concord) is deferred to Track B v1.1 — needs script flags added
# to 02_joint_train.py first.

#SBATCH --job-name=xen_joint_subset
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=.slurm_stubs/joint_subset_%j.out
#SBATCH --error=.slurm_stubs/joint_subset_%j.err

if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
CONTAINER_TYPE=python_scvi_gpu_2025Q3
source "${PIPELINE_ROOT}/run/_common.sh"

mkdir -p "${CFG_HPC_LOGS_DIR}"
exec >> "${CFG_HPC_LOGS_DIR}/joint_subset_${SLURM_JOB_ID}.out" \
    2>> "${CFG_HPC_LOGS_DIR}/joint_subset_${SLURM_JOB_ID}.err"

# Pre-write-gate: subset outputs land in staging tree (NOT publication/).
OUT="${PIPELINE_ROOT}/run/.logs/subset_${SLURM_JOB_ID}"
mkdir -p "${OUT}/model"

run_python_singularity \
    "${PIPELINE_ROOT}/4b-joint/scripts/02_joint_embedding/02_joint_train.py" \
        --flex-dir            "${CFG_FLEX_CONCORD_UPSTREAM}" \
        --flex-qc             "${CFG_XENIUM_FLEX_QC_UPSTREAM}" \
        --xenium-bundle       "${CFG_XENIUM_POOLED_NUCLEAR_UPSTREAM}" \
        --source              nuclear \
        --panel-intersection  "${CFG_XENIUM_LEGACY_PANEL_INTERSECTION}/panel_flex_intersection_Immune.txt" \
        --compartment-dir     "${CFG_XENIUM_COMPARTMENT_DIR}" \
        --out-dir             "${OUT}" \
        --n-latent            8 \
        --seed                42

echo "Done joint subset wrapper-resolution gate: $(date)"
