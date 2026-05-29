#!/bin/bash
#SBATCH --job-name=motif_05_type
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=slurm-%x_%j.out
#SBATCH --error=slurm-%x_%j.err

# Cell-type null: within-instance L1.5 label permutation (xy fixed); project
# shuffled features through canonical H; complement to spatial null.
# Tier 2 per hpc-resource-rules.md.

set -euo pipefail

if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$PWD"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
source "${PIPELINE_ROOT}/config/load_paths.sh"

CONTAINER="${CONTAINER_PYTHON_SCVI_GPU_2025Q3}"
PYTHONUSERBASE="${PYTHON_LIBS_SCGPT_GPU_2025Q3}"

JOB_TMP="/tmp/motif_05_${SLURM_JOB_ID:-local}"
mkdir -p "${JOB_TMP}" "${LOGS_ROOT}"
trap 'rm -rf "${JOB_TMP}"' EXIT

module load singularity
singularity exec \
    --cleanenv --containall --no-home \
    ${BIND_MOUNTS} \
    --bind /pub/nwechter:/pub/nwechter:rw \
    --bind /tmp:/tmp:rw \
    --env "PYTHONUSERBASE=${PYTHONUSERBASE}" \
    --env "NUMBA_CACHE_DIR=${JOB_TMP}/numba" \
    --env "MPLCONFIGDIR=${JOB_TMP}/mpl" \
    "${CONTAINER}" \
    python3 "${PIPELINE_ROOT}/scripts/05_validate_type_shuffle.py" \
        --project-root "${PIPELINE_ROOT}" \
        --cells-parquet "${SOURCE_MOTIF_JOINT_CELLS_PARQUET}" \
        --basis "${SOURCE_MOTIF_BASIS_CSV}" \
        --out-dir "${MOTIF_VALIDATION_DIR}/s05_type_shuffle"
echo "=== Done @ $(date) ==="
