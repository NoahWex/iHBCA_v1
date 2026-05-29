#!/bin/bash
#
# 4b-joint phase 6: three-axis Xenium filter.
#
# Applies three filters: (1) FLEX proximity in joint latent (cell must have a
# FLEX neighbor within q95 distance), (2) NMP rule (nuclear AND cyto NMP below
# threshold), (3) QC pass (pass_qc_whole). Output flags each cell pass/fail per
# axis and provides combined pass column.
#
# Tier 2: 4 CPUs / 16G / 1h.

#SBATCH --job-name=xen_p6_filter
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=.slurm_stubs/p6_filter_%j.out
#SBATCH --error=.slurm_stubs/p6_filter_%j.err

if [ -n "${SLURM_JOB_ID:-}" ]; then SCRIPT_DIR="$(pwd)"; else SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"; fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
CONTAINER_TYPE=python_spatial_2025Q4E2
source "${PIPELINE_ROOT}/run/_common.sh"

mkdir -p "${CFG_HPC_LOGS_DIR}"
exec >> "${CFG_HPC_LOGS_DIR}/p6_filter_${SLURM_JOB_ID}.out" 2>> "${CFG_HPC_LOGS_DIR}/p6_filter_${SLURM_JOB_ID}.err"

OUT="$(dirname "${CFG_XENIUM_THREE_AXIS_FILTER}")"
mkdir -p "${OUT}"

# Step 6.1: FLEX proximity (kNN xenium→flex + flex→flex self)
run_python_singularity \
    "${PIPELINE_ROOT}/4b-joint/scripts/06_filter/01_flex_proximity.py" \
        --latent  "${CFG_XENIUM_JOINT_LATENT}" \
        --obs     "${CFG_XENIUM_JOINT_OBS}" \
        --k       30 \
        --out-dir "${OUT}"

# Step 6.2: three-axis combine (proximity Q95 + NMP rule + QC whole)
run_python_singularity \
    "${PIPELINE_ROOT}/4b-joint/scripts/06_filter/02_three_axis_filter.py" \
        --proximity   "${OUT}/xenium_to_flex_knn.csv" \
        --flex-self   "${OUT}/flex_self_knn.csv" \
        --proximity-q 95.0 \
        --nmp-csv     "${CFG_XENIUM_NMP_SCORED}" \
        --qc-csv      "${CFG_XENIUM_QC}" \
        --out-csv     "${CFG_XENIUM_THREE_AXIS_FILTER}"

echo "Done phase6_filter: $(date)"
