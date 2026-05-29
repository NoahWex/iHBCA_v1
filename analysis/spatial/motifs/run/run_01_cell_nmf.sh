#!/bin/bash
#SBATCH --job-name=motif_01_cell_nmf
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=slurm-%x_%j.out
#SBATCH --error=slurm-%x_%j.err

# Cell-level density-weighted NMF on non-epithelial L1.5 neighborhoods (k=9).
# Tier 3 per hpc-resource-rules.md (32G / 4 CPU / 2h).

set -euo pipefail

if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$PWD"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Canonical resolver -- populates SOURCE_*, CONTAINER_*, PYTHON_LIBS_*, BIND_MOUNTS
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"

# Module-local resolver -- MOTIF_OUTPUTS_ROOT, MOTIF_VALIDATION_DIR, ...
source "${PIPELINE_ROOT}/config/load_paths.sh"

CONTAINER="${CONTAINER_PYTHON_SCVI_GPU_2025Q3}"
PYTHONUSERBASE="${PYTHON_LIBS_SCGPT_GPU_2025Q3}"

JOB_TMP="/tmp/motif_01_${SLURM_JOB_ID:-local}"
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
    python3 "${PIPELINE_ROOT}/scripts/01_cell_nmf.py" \
        --project-root "${PIPELINE_ROOT}" \
        --joint-l1p5 "${SOURCE_XENIUM_JOINT_L1P5}" \
        --annotated-dir "${SOURCE_XENIUM_SPACEFLOW_DIR}" \
        --manifest "${SOURCE_XENIUM_MANIFEST}" \
        --out-dir "${MOTIF_OUTPUTS_ROOT}"
echo "=== Done @ $(date) ==="
