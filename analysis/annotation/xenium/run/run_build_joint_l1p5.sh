#!/bin/bash
#
# 4c-annotation step 9.6: build cohort-wide L1.5 table from per-compartment outputs.
#
# Joins joint_obs.csv (cohort cell list + platform/compartment) with per-
# compartment cell_annotations.csv (L1.5 labels) and leiden_assignments.csv
# (multi-resolution leiden), plus library_covariates.csv (per-library + per-
# patient metadata). Produces joint_l1p5.csv: the master per-cell L1.5 table
# consumed by Phase D Xenium DA, LBridge bridge, and Figs 3-5.
#
# Tier 1: 2 CPUs / 8G / 15min.

#SBATCH --job-name=xen_build_joint_l1p5
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=.slurm_stubs/build_joint_l1p5_%j.out
#SBATCH --error=.slurm_stubs/build_joint_l1p5_%j.err

if [ -n "${SLURM_JOB_ID:-}" ]; then SCRIPT_DIR="$(pwd)"; else SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"; fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CONTAINER_TYPE=python_spatial_2025Q4E2
source "${PIPELINE_ROOT}/run/_common.sh"

mkdir -p "${CFG_HPC_LOGS_DIR}"
exec >> "${CFG_HPC_LOGS_DIR}/build_joint_l1p5_${SLURM_JOB_ID}.out" 2>> "${CFG_HPC_LOGS_DIR}/build_joint_l1p5_${SLURM_JOB_ID}.err"

run_python_singularity \
    "${PIPELINE_ROOT}/4c-annotation/scripts/06_build_joint_l1p5.py" \
        --joint-obs           "${CFG_XENIUM_JOINT_OBS}" \
        --compartment-dir     "${CFG_XENIUM_COMPARTMENT_DIR}" \
        --library-covariates  "${CFG_XENIUM_LIBRARY_COVARIATES}" \
        --out                 "${CFG_XENIUM_JOINT_L1P5}"

echo "Done build_joint_l1p5: $(date)"
