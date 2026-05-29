#!/bin/bash
#SBATCH --job-name=trackc_04b_concordance
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --array=0-2%3
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:45:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/04b_concordance_%A_%a.out
#SBATCH --error=.slurm_stubs/04b_concordance_%A_%a.err

# Step 4b: cross-resolution F1 concordance per V1 label. Thresholds per-cell
# UCell scores at the resolution-specific best cluster, computes F1 per label
# at each resolution, and identifies the per-label resolution within
# epsilon = 0.05 of the best F1 (anchor-resolution recommendation).

set -euo pipefail

CONTAINER_TYPE=python_spatial_2025Q2
if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
source "${SCRIPT_DIR}/_common.sh"

mkdir -p "${HPC_LOGS_DIR}"
exec >> "${HPC_LOGS_DIR}/04b_concordance_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}.out" \
    2>> "${HPC_LOGS_DIR}/04b_concordance_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}.err"

COMPARTMENTS=(Epithelial Immune Stromal)
COMPARTMENT="${COMPARTMENTS[$SLURM_ARRAY_TASK_ID]}"
case "${COMPARTMENT}" in
    Epithelial) YAML_SHORT=epi ;;
    Immune)     YAML_SHORT=imm ;;
    Stromal)    YAML_SHORT=str ;;
esac

RESOLUTIONS=(0.3 0.5 0.8 1.0 1.5 3.0 5.0)

SCRIPT="${CFG_SCRIPTS_DIR}/04b_resolution_concordance.py"
UCELL_PC="${CFG_OUTPUTS_ROOT}/ucell/${COMPARTMENT}_ucell_per_cell.csv"
CLUSTERS="${CFG_OUTPUTS_ROOT}/clusters/${COMPARTMENT}_clusters.csv"
YAML="${CFG_V1_ANNOTATION_YAMLS}/annotation_v2_${YAML_SHORT}.yaml"
OUT_DIR="${CFG_OUTPUTS_ROOT}/concordance"

mkdir -p "${OUT_DIR}"

_load_singularity

singularity exec --no-home \
    ${BIND_MOUNTS:-} \
    --bind "${JOB_TEMP_DIR}:${JOB_TEMP_DIR}:rw" \
    --env "PYTHONUNBUFFERED=1" \
    "${CONTAINER_PATH}" \
    python "${SCRIPT}" \
        --ucell-per-cell-csv "${UCELL_PC}" \
        --clusters-csv "${CLUSTERS}" \
        --resolutions "${RESOLUTIONS[@]}" \
        --yaml "${YAML}" \
        --compartment "${COMPARTMENT}" \
        --out-dir "${OUT_DIR}"

echo "=== Step 4b complete (${COMPARTMENT}) ==="
