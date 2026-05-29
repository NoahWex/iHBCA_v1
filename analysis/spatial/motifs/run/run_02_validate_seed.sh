#!/bin/bash
#SBATCH --job-name=motif_02_seed
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=slurm-%x_%j.out
#SBATCH --error=slurm-%x_%j.err

# Seed stability validation: 20 NMF refits at k=9 with independent random seeds.
# Reports per-program best-match cosine vs canonical basis + pairwise cohort matrix.
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

JOB_TMP="/tmp/motif_02_${SLURM_JOB_ID:-local}"
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
    python3 "${PIPELINE_ROOT}/scripts/02_validate_seed.py" \
        --project-root "${PIPELINE_ROOT}" \
        --cells-parquet "${SOURCE_MOTIF_JOINT_CELLS_PARQUET}" \
        --reference-basis "${SOURCE_MOTIF_BASIS_CSV}" \
        --out-dir "${MOTIF_VALIDATION_DIR}/s36_seed_stability"
echo "=== Done @ $(date) ==="
