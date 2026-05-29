#!/bin/bash
#
# 4b-joint phase 3: cluster alignment — FLEX Leiden on joint latent → L0.5 vocabulary.
#
# Runs FLEX-only Leiden clustering on the joint Concord latent space at multiple
# resolutions, then computes contingency between FLEX Leiden clusters and the
# L0.5 vocabulary. Defines which joint clusters correspond to canonical L0.5
# cell types.
#
# Tier 3: 4 CPUs / 32G / 2h.

#SBATCH --job-name=xen_phase3_clean
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=.slurm_stubs/phase3_%j.out
#SBATCH --error=.slurm_stubs/phase3_%j.err

if [ -n "${SLURM_JOB_ID:-}" ]; then SCRIPT_DIR="$(pwd)"; else SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"; fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
CONTAINER_TYPE=python_scvi_gpu_2025Q3
source "${PIPELINE_ROOT}/run/_common.sh"

mkdir -p "${CFG_HPC_LOGS_DIR}"
exec >> "${CFG_HPC_LOGS_DIR}/phase3_${SLURM_JOB_ID}.out" 2>> "${CFG_HPC_LOGS_DIR}/phase3_${SLURM_JOB_ID}.err"

OUT="$(dirname "${CFG_XENIUM_JOINT_LATENT}")"
mkdir -p "${OUT}"

# L0.5 reference (legacy FLEX L0.5-per-cell from the FLEX panel work)
L0P5="${CFG_XENIUM_LEGACY_PANEL_INTERSECTION}/../l0p5_per_cell.csv"

# Step 3.1: FLEX Leiden on joint latent
run_python_singularity \
    "${PIPELINE_ROOT}/4b-joint/scripts/03_cluster_alignment/01_flex_leiden_joint.py" \
        --latent  "${CFG_XENIUM_JOINT_LATENT}" \
        --obs     "${CFG_XENIUM_JOINT_OBS}" \
        --out-dir "${OUT}"

# Step 3.2: contingency table FLEX Leiden × L0.5
run_python_singularity \
    "${PIPELINE_ROOT}/4b-joint/scripts/03_cluster_alignment/02_contingency.py" \
        --leiden-csv "${OUT}/flex_leiden_joint.csv" \
        --l0p5-csv   "${L0P5}" \
        --out-dir    "${OUT}"

echo "Done phase3_clean: $(date)"
