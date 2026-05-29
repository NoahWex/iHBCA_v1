#!/bin/bash
#SBATCH --job-name=trackc_03a_pseudobulk
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --array=0-20%7
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/03a_pseudobulk_%A_%a.out
#SBATCH --error=.slurm_stubs/03a_pseudobulk_%A_%a.err

# Step 3a: aggregate counts to (cluster x library_id) pseudobulk per
# (compartment, resolution). Largest single load is the Stromal counts.mtx.gz
# at ~119K cells x 18K genes; sum-aggregation in pandas after sparse load.

set -euo pipefail

CONTAINER_TYPE=python_spatial_2025Q2
if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
source "${SCRIPT_DIR}/_common.sh"

mkdir -p "${HPC_LOGS_DIR}"
exec >> "${HPC_LOGS_DIR}/03a_pseudobulk_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}.out" \
    2>> "${HPC_LOGS_DIR}/03a_pseudobulk_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}.err"

COMPARTMENTS=(Epithelial Immune Stromal)
RESOLUTIONS=(0.3 0.5 0.8 1.0 1.5 3.0 5.0)
N_RES=${#RESOLUTIONS[@]}

C_IDX=$((SLURM_ARRAY_TASK_ID / N_RES))
R_IDX=$((SLURM_ARRAY_TASK_ID % N_RES))
COMPARTMENT="${COMPARTMENTS[$C_IDX]}"
RES="${RESOLUTIONS[$R_IDX]}"

SCRIPT="${CFG_SCRIPTS_DIR}/03a_pseudobulk.py"
BUNDLE_DIR="${CFG_INTEGRATION_INTERMEDIATE}/${COMPARTMENT}/scvi_n100"
CLUSTERS="${CFG_OUTPUTS_ROOT}/clusters/${COMPARTMENT}_clusters.csv"
OUT_DIR="${CFG_OUTPUTS_ROOT}/pseudobulk"

mkdir -p "${OUT_DIR}"

_load_singularity

echo "=== Step 3a pseudobulk: ${COMPARTMENT} leiden_${RES} ==="
echo "Start: $(date)"

singularity exec --no-home \
    ${BIND_MOUNTS:-} \
    --bind "${JOB_TEMP_DIR}:${JOB_TEMP_DIR}:rw" \
    --env "PYTHONUNBUFFERED=1" \
    --env "NUMBA_CACHE_DIR=${JOB_TEMP_DIR}/numba_cache" \
    --env "MPLCONFIGDIR=${JOB_TEMP_DIR}/mpl_config" \
    "${CONTAINER_PATH}" \
    python "${SCRIPT}" \
        --compartment "${COMPARTMENT}" \
        --bundle-dir "${BUNDLE_DIR}" \
        --clusters-csv "${CLUSTERS}" \
        --resolution "${RES}" \
        --out-dir "${OUT_DIR}"

echo "Done: $(date)"
