#!/bin/bash
#SBATCH --job-name=trackc_01_cluster
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --array=0-2%3
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:30:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/01_cluster_%A_%a.out
#SBATCH --error=.slurm_stubs/01_cluster_%A_%a.err

# Step 1: extract precomputed Leiden + recompute missing resolutions per compartment.
# Pulls {0.3, 0.5, 0.8, 1.0, 5.0} from the integration bundle's obs.csv and
# computes {1.5, 3.0} on the 50-D scVI latent via scanpy.pp.neighbors +
# leidenalg + igraph (scanpy's built-in Leiden is bypassed in favor of the
# native libraries to avoid a known scanpy regression on igraph >= 0.10).

set -euo pipefail

CONTAINER_TYPE=python_spatial_2025Q2
if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
source "${SCRIPT_DIR}/_common.sh"

mkdir -p "${HPC_LOGS_DIR}"
exec >> "${HPC_LOGS_DIR}/01_cluster_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}.out" \
    2>> "${HPC_LOGS_DIR}/01_cluster_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}.err"

COMPARTMENTS=(Epithelial Immune Stromal)
COMPARTMENT="${COMPARTMENTS[$SLURM_ARRAY_TASK_ID]}"

RESOLUTIONS=(0.3 0.5 0.8 1.0 1.5 3.0 5.0)
SCRIPT="${CFG_SCRIPTS_DIR}/01_cluster_compartment.py"
BUNDLE_DIR="${CFG_INTEGRATION_INTERMEDIATE}/${COMPARTMENT}/scvi_n100"
OUT_DIR="${CFG_OUTPUTS_ROOT}/clusters"

mkdir -p "${OUT_DIR}"

_load_singularity

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
        --resolutions "${RESOLUTIONS[@]}" \
        --out-dir "${OUT_DIR}"

echo "=== Step 1 complete (${COMPARTMENT}) ==="
