#!/bin/bash
#SBATCH --job-name=joint_milo
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --array=0-8

set -euo pipefail

STUDIES=(
    "gray"           # 0
    "kumar"          # 1
    "murrow"         # 2
    "nee"            # 3
    "twigger"        # 4
    "pal_norm_epi"   # 5
    "pal_norm_total" # 6
    "pal_norm_b1"    # 7
    "reed"           # 8
)

STUDY="${STUDIES[$SLURM_ARRAY_TASK_ID]}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../../config/load_paths.sh"
BASE_PATH="${SOURCE_IHBCAV1_EXTERNAL_STUDIES%/}"
SCRIPT="${BASE_PATH}/harmonization/scripts/build_joint_milo.R"
LOG_DIR="${BASE_PATH}/harmonization/logs"
STUDY_DIR="${BASE_PATH}/harmonization/outputs/study_objects/${STUDY}"

CONTAINER="${CONTAINER_R_SPATIAL_4_3_3}"
R_LIBS_USER="${R_LIBS_R_SPATIAL_4_3_3}"

mkdir -p "$LOG_DIR"

echo "=============================================="
echo "Build Joint Milo: $STUDY"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task: $SLURM_ARRAY_TASK_ID"
echo "Node: $(hostname)"
echo "Start time: $(date)"
echo "=============================================="

if [[ ! -f "${STUDY_DIR}/seurat.rds" ]]; then
    echo "ERROR: seurat.rds not found in ${STUDY_DIR}"
    exit 1
fi

if [[ ! -f "${STUDY_DIR}/embedding_joint.csv" ]]; then
    echo "ERROR: embedding_joint.csv not found in ${STUDY_DIR}"
    exit 1
fi

module load singularity

singularity exec \
    ${BIND_MOUNTS} \
    --bind "$R_LIBS_USER:/home/jovyan/R/library:ro" \
    --env "R_LIBS_USER=/home/jovyan/R/library" \
    "$CONTAINER" \
    Rscript "$SCRIPT" --study "$STUDY"

echo ""
echo "=============================================="
echo "Completed: $STUDY"
echo "End time: $(date)"
echo "=============================================="
