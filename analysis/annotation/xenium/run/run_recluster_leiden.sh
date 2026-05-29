#!/bin/bash
#
# 4c-annotation step 9.2b: re-run Leiden at additional resolutions per compartment.
#
# Adds finer-grained Leiden resolutions (used during cascade override decisions
# to split clusters that show heterogeneity at 0.5 but cleaner separation at
# 0.7 or 1.0).
#
# Tier 3: 4 CPUs / 32G / 1h per task.

#SBATCH --job-name=xen_recluster
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --array=0-2
#SBATCH --output=.slurm_stubs/recluster_%A_%a.out
#SBATCH --error=.slurm_stubs/recluster_%A_%a.err

if [ -n "${SLURM_JOB_ID:-}" ]; then SCRIPT_DIR="$(pwd)"; else SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"; fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CONTAINER_TYPE=python_spatial_2025Q4E2
source "${PIPELINE_ROOT}/run/_common.sh"

mkdir -p "${CFG_HPC_LOGS_DIR}"
exec >> "${CFG_HPC_LOGS_DIR}/recluster_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.out" \
    2>> "${CFG_HPC_LOGS_DIR}/recluster_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.err"

COMPARTMENTS=(Epithelial Stromal Immune)
COMPARTMENT=${COMPARTMENTS[$SLURM_ARRAY_TASK_ID]}
COMP_LC=$(echo "$COMPARTMENT" | tr '[:upper:]' '[:lower:]')
COMP_DIR="${CFG_XENIUM_COMPARTMENT_DIR}/${COMP_LC}"

run_python_singularity \
    "${PIPELINE_ROOT}/4c-annotation/scripts/02b_recluster_leiden.py" \
        --compartment "${COMPARTMENT}" \
        --latent      "${COMP_DIR}/joint_latent.csv" \
        --obs         "${COMP_DIR}/joint_obs.csv" \
        --out-dir     "${COMP_DIR}" \
        -k 30 \
        --resolutions 0.1 0.3 0.5 0.7 1.0 1.5

echo "Done recluster ${COMPARTMENT}: $(date)"
