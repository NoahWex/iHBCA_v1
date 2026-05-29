#!/bin/bash
#SBATCH --job-name=trackc_02_ucell
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --array=0-2%3
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/02_ucell_%A_%a.out
#SBATCH --error=.slurm_stubs/02_ucell_%A_%a.err

# Step 2: UCell rank-based scoring per cell + per-cluster aggregation across all
# resolutions. Two score types per V1 label (canonical markers, identity markers)
# yield ~60 signatures over 50K-120K cells per compartment.

set -euo pipefail

CONTAINER_TYPE=r_spatial_4.3.3
if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
source "${SCRIPT_DIR}/_common.sh"

mkdir -p "${HPC_LOGS_DIR}"
exec >> "${HPC_LOGS_DIR}/02_ucell_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}.out" \
    2>> "${HPC_LOGS_DIR}/02_ucell_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}.err"

COMPARTMENTS=(Epithelial Immune Stromal)
COMPARTMENT="${COMPARTMENTS[$SLURM_ARRAY_TASK_ID]}"
case "${COMPARTMENT}" in
    Epithelial) YAML_SHORT=epi ;;
    Immune)     YAML_SHORT=imm ;;
    Stromal)    YAML_SHORT=str ;;
    *) echo "ERROR: unknown compartment ${COMPARTMENT}"; exit 1 ;;
esac

RESOLUTIONS=(0.3 0.5 0.8 1.0 1.5 3.0 5.0)
SCRIPT="${CFG_SCRIPTS_DIR}/02_ucell_scoring.R"
BUNDLE_DIR="${CFG_INTEGRATION_INTERMEDIATE}/${COMPARTMENT}/scvi_n100"
CLUSTERS="${CFG_OUTPUTS_ROOT}/clusters/${COMPARTMENT}_clusters.csv"
YAML="${CFG_V1_ANNOTATION_YAMLS}/annotation_v2_${YAML_SHORT}.yaml"
OUT_DIR="${CFG_OUTPUTS_ROOT}/ucell"

mkdir -p "${OUT_DIR}"

_load_singularity

singularity exec --no-home \
    ${BIND_MOUNTS:-} \
    --bind "${R_LIBS_R_SPATIAL_4_3_3}:/home/jovyan/R/userlib:ro" \
    --bind "${JOB_TEMP_DIR}:${JOB_TEMP_DIR}:rw" \
    --env "R_LIBS_USER=/home/jovyan/R/userlib" \
    --env "TMPDIR=${JOB_TEMP_DIR}" \
    "${CONTAINER_PATH}" \
    Rscript "${SCRIPT}" \
        --bundle-dir "${BUNDLE_DIR}" \
        --clusters-csv "${CLUSTERS}" \
        --resolutions "${RESOLUTIONS[@]}" \
        --yaml "${YAML}" \
        --out-dir "${OUT_DIR}" \
        --compartment "${COMPARTMENT}" \
        --ncores 4

echo "=== Step 2 complete (${COMPARTMENT}) ==="
