#!/bin/bash
#SBATCH --job-name=trackc_04_assignment
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --array=0-20%7
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/04_assignment_%A_%a.out
#SBATCH --error=.slurm_stubs/04_assignment_%A_%a.err

# Step 4: best-match L2S assignment per cluster + novelty heuristics, per
# (compartment, resolution). Combines UCell per-cluster scores (canonical x 1 +
# identity x 2) with limma top markers. Novelty flags: low_confidence (combined
# score < 0.30), tight_margin (top1 - top2 < 0.05), marker_mismatch (<= 2
# canonical markers in cluster's top 50 limma genes).

set -euo pipefail

CONTAINER_TYPE=python_spatial_2025Q2
if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
source "${SCRIPT_DIR}/_common.sh"

mkdir -p "${HPC_LOGS_DIR}"
exec >> "${HPC_LOGS_DIR}/04_assignment_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}.out" \
    2>> "${HPC_LOGS_DIR}/04_assignment_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}.err"

COMPARTMENTS=(Epithelial Immune Stromal)
RESOLUTIONS=(0.3 0.5 0.8 1.0 1.5 3.0 5.0)
N_RES=${#RESOLUTIONS[@]}

C_IDX=$((SLURM_ARRAY_TASK_ID / N_RES))
R_IDX=$((SLURM_ARRAY_TASK_ID % N_RES))
COMPARTMENT="${COMPARTMENTS[$C_IDX]}"
RES="${RESOLUTIONS[$R_IDX]}"

case "${COMPARTMENT}" in
    Epithelial) YAML_SHORT=epi ;;
    Immune)     YAML_SHORT=imm ;;
    Stromal)    YAML_SHORT=str ;;
esac

SCRIPT="${CFG_SCRIPTS_DIR}/04_assignment_novelty.py"
UCELL="${CFG_OUTPUTS_ROOT}/ucell/${COMPARTMENT}_ucell_per_cluster.csv"
LIMMA_TOP="${CFG_OUTPUTS_ROOT}/limma/limma_top_markers_${COMPARTMENT}_leiden_${RES}.csv"
YAML="${CFG_V1_ANNOTATION_YAMLS}/annotation_v2_${YAML_SHORT}.yaml"
CLUSTERS="${CFG_OUTPUTS_ROOT}/clusters/${COMPARTMENT}_clusters.csv"
OUT_DIR="${CFG_OUTPUTS_ROOT}/assignments"

mkdir -p "${OUT_DIR}"

_load_singularity

singularity exec --no-home \
    ${BIND_MOUNTS:-} \
    --bind "${JOB_TEMP_DIR}:${JOB_TEMP_DIR}:rw" \
    --env "PYTHONUNBUFFERED=1" \
    "${CONTAINER_PATH}" \
    python "${SCRIPT}" \
        --ucell-csv "${UCELL}" \
        --limma-top-csv "${LIMMA_TOP}" \
        --yaml "${YAML}" \
        --clusters-csv "${CLUSTERS}" \
        --compartment "${COMPARTMENT}" \
        --resolution "${RES}" \
        --out-dir "${OUT_DIR}"

echo "=== Step 4 complete (${COMPARTMENT} leiden_${RES}) ==="
