#!/bin/bash
#SBATCH --job-name=motif_32b_hull_comp
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=slurm-%x_%j.out
#SBATCH --error=slurm-%x_%j.err

# Per-unit point-in-polygon over alpha hulls: counts ALL L1.5 cell types
# within each motif unit's geometry (vs motif_units.parquet's argmax-only
# composition). Captures the full local cellular context, including epi.
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

XENIUM_OBS_CSV="${SOURCE_XENIUM_NUCLEAR_BUNDLE}/xenium_obs.csv"

JOB_TMP="/tmp/motif_32b_${SLURM_JOB_ID:-local}"
mkdir -p "${JOB_TMP}" "${LOGS_ROOT}"
trap 'rm -rf "${JOB_TMP}"' EXIT

module load singularity
singularity exec \
    --cleanenv --containall --no-home \
    ${BIND_MOUNTS} \
    --bind /pub/nwechter:/pub/nwechter:rw \
    --bind /tmp:/tmp:rw \
    --env "PYTHONUSERBASE=${PYTHONUSERBASE}" \
    --env "MPLCONFIGDIR=${JOB_TMP}/mpl" \
    "${CONTAINER}" \
    python3 "${PIPELINE_ROOT}/scripts/32b_unit_hull_composition.py" \
        --motif-units "${SOURCE_MOTIF_UNITS_PARQUET}" \
        --xenium-obs "${XENIUM_OBS_CSV}" \
        --joint-l1p5 "${SOURCE_XENIUM_JOINT_L1P5}" \
        --out-dir "${MOTIF_UNITS_DIR}"
echo "=== Done @ $(date) ==="
