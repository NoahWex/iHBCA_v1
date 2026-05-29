#!/bin/bash
#SBATCH --job-name=trackc_05d_canonical
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --array=0-2%3
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/05d_canonical_%A_%a.out
#SBATCH --error=.slurm_stubs/05d_canonical_%A_%a.err

# Step 5d: V1 canonical + identity marker by cluster heatmap, one PDF per
# resolution. Cluster-mean expression normalized within marker; rows ordered
# by V1 label diagonal alignment.

set -euo pipefail

CONTAINER_TYPE=python_spatial_2025Q2
if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
source "${SCRIPT_DIR}/_common.sh"

mkdir -p "${HPC_LOGS_DIR}"
exec >> "${HPC_LOGS_DIR}/05d_canonical_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}.out" \
    2>> "${HPC_LOGS_DIR}/05d_canonical_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}.err"

COMPARTMENTS=(Epithelial Immune Stromal)
COMPARTMENT="${COMPARTMENTS[$SLURM_ARRAY_TASK_ID]}"

case "${COMPARTMENT}" in
    Epithelial) YAML_SHORT=epi ;;
    Immune)     YAML_SHORT=imm ;;
    Stromal)    YAML_SHORT=str ;;
esac

SCRIPT="${CFG_SCRIPTS_DIR}/05d_canonical_marker_heatmap.py"
BUNDLE_DIR="${CFG_INTEGRATION_INTERMEDIATE}/${COMPARTMENT}/scvi_n100"
YAML="${CFG_V1_ANNOTATION_YAMLS}/annotation_v2_${YAML_SHORT}.yaml"
CLUSTERS="${CFG_OUTPUTS_ROOT}/clusters/${COMPARTMENT}_clusters.csv"
ASSIGNMENTS_DIR="${CFG_OUTPUTS_ROOT}/assignments"
OUT_DIR="${CFG_OUTPUTS_ROOT}/${COMPARTMENT}"

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
        --bundle-dir "${BUNDLE_DIR}" \
        --clusters-csv "${CLUSTERS}" \
        --yaml "${YAML}" \
        --resolutions 0.3 0.5 0.8 1.0 1.5 3.0 \
        --assignments-dir "${ASSIGNMENTS_DIR}" \
        --compartment "${COMPARTMENT}" \
        --out-dir "${OUT_DIR}"

echo "=== Step 5d complete (${COMPARTMENT}) ==="
