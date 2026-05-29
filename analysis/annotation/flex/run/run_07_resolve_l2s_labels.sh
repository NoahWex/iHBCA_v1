#!/bin/bash
#SBATCH --job-name=trackc_07_resolve
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/07_resolve_%j.out
#SBATCH --error=.slurm_stubs/07_resolve_%j.err

# Step 7: resolve per-cell L2S labels for all 3 compartments from the
# per-compartment annotation yamls (v2_flex schema). Reads each
# annotation_yamls/{compartment}.yaml + clusters CSV, applies anchor labels,
# overlays fine_cluster_overrides, concatenates across compartments, and
# writes the canonical flex_l2s_labels.csv.

set -euo pipefail

CONTAINER_TYPE=python_spatial_2025Q2
if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
source "${SCRIPT_DIR}/_common.sh"

mkdir -p "${HPC_LOGS_DIR}"
exec >> "${HPC_LOGS_DIR}/07_resolve_${SLURM_JOB_ID:-local}.out" \
    2>> "${HPC_LOGS_DIR}/07_resolve_${SLURM_JOB_ID:-local}.err"

SCRIPT="${CFG_SCRIPTS_DIR}/07_resolve_l2s_labels.py"
YAMLS_DIR="${CFG_ANNOTATION_YAMLS}"
CLUSTERS_DIR="${CFG_OUTPUTS_ROOT}/clusters"
OUT_CSV="${CFG_OUTPUTS_ROOT}/flex_l2s_labels.csv"

for required in "${SCRIPT}" "${YAMLS_DIR}/epithelial.yaml" "${YAMLS_DIR}/immune.yaml" "${YAMLS_DIR}/stromal.yaml"; do
    [ -f "${required}" ] || { echo "ERROR: missing required input ${required}"; exit 1; }
done
[ -d "${CLUSTERS_DIR}" ] || { echo "ERROR: missing clusters dir ${CLUSTERS_DIR}"; exit 1; }

mkdir -p "$(dirname "${OUT_CSV}")"

_load_singularity

singularity exec --no-home \
    ${BIND_MOUNTS:-} \
    --bind "${JOB_TEMP_DIR}:${JOB_TEMP_DIR}:rw" \
    --env "PYTHONUNBUFFERED=1" \
    "${CONTAINER_PATH}" \
    python "${SCRIPT}" \
        --yamls-dir "${YAMLS_DIR}" \
        --clusters-dir "${CLUSTERS_DIR}" \
        --out-csv "${OUT_CSV}"

echo "=== Step 7 complete: ${OUT_CSV} ==="
ls -la "${OUT_CSV}"
