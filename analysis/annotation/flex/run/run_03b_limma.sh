#!/bin/bash
#SBATCH --job-name=trackc_03b_limma
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --array=0-20%7
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/03b_limma_%A_%a.out
#SBATCH --error=.slurm_stubs/03b_limma_%A_%a.err

# Step 3b: limma-voom one-vs-rest de novo markers per cluster, per
# (compartment, resolution). Design ~ library_id + condition with dupCor on
# patient_id (FLEX has n_patients = 4). Memory bounded by voom on the full
# pseudobulk count matrix.

set -euo pipefail

CONTAINER_TYPE=r_spatial_4.3.3
if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
source "${SCRIPT_DIR}/_common.sh"

mkdir -p "${HPC_LOGS_DIR}"
exec >> "${HPC_LOGS_DIR}/03b_limma_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}.out" \
    2>> "${HPC_LOGS_DIR}/03b_limma_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}.err"

COMPARTMENTS=(Epithelial Immune Stromal)
RESOLUTIONS=(0.3 0.5 0.8 1.0 1.5 3.0 5.0)
N_RES=${#RESOLUTIONS[@]}

C_IDX=$((SLURM_ARRAY_TASK_ID / N_RES))
R_IDX=$((SLURM_ARRAY_TASK_ID % N_RES))
COMPARTMENT="${COMPARTMENTS[$C_IDX]}"
RES="${RESOLUTIONS[$R_IDX]}"

SCRIPT="${CFG_SCRIPTS_DIR}/03b_limma_markers.R"
PB_DIR="${CFG_OUTPUTS_ROOT}/pseudobulk"
OUT_DIR="${CFG_OUTPUTS_ROOT}/limma"
mkdir -p "${OUT_DIR}"

COUNTS="${PB_DIR}/${COMPARTMENT}_leiden_${RES}_pseudobulk_counts.csv"
META="${PB_DIR}/${COMPARTMENT}_leiden_${RES}_pseudobulk_meta.csv"

for f in "${COUNTS}" "${META}" "${SCRIPT}"; do
    [ -f "${f}" ] || { echo "ERROR: missing input: ${f}"; exit 1; }
done

_load_singularity

echo "=== Step 3b limma: ${COMPARTMENT} leiden_${RES} ==="
echo "Start: $(date)"

singularity exec --no-home \
    ${BIND_MOUNTS:-} \
    --bind "${R_LIBS_R_SPATIAL_4_3_3}:/home/jovyan/R/userlib:ro" \
    --bind "${JOB_TEMP_DIR}:${JOB_TEMP_DIR}:rw" \
    --env "R_LIBS_USER=/home/jovyan/R/userlib" \
    --env "TMPDIR=${JOB_TEMP_DIR}" \
    "${CONTAINER_PATH}" \
    Rscript "${SCRIPT}" \
        --counts "${COUNTS}" \
        --meta "${META}" \
        --output-dir "${OUT_DIR}" \
        --level-tag "${COMPARTMENT}_leiden_${RES}" \
        --n-top 50

echo "Done: $(date)"
