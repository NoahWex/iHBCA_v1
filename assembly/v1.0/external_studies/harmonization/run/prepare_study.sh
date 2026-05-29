#!/bin/bash
#SBATCH --job-name=prep_study
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --array=0-8

# =============================================================================
# prepare_study.sh — Phase 1 of parallel Milo build pipeline
# =============================================================================
# Prepares study objects (steps 1-6, no Milo build).
# Lighter than build_study_objects.sh — completes faster.
#
# Usage:
#   sbatch prepare_study.sh                                    # All studies
#   sbatch --array=0 prepare_study.sh                          # gray only
#   sbatch --mem=280G --time=4:00:00 --array=1,8 prepare_study.sh  # XL studies
# =============================================================================

set -euo pipefail

STUDIES=(
    "gray"           # 0 - 52K cells
    "kumar"          # 1 - 714K cells (XL)
    "murrow"         # 2 - 86K cells
    "nee"            # 3 - 230K cells
    "twigger"        # 4 - 111K cells
    "pal_norm_epi"   # 5 - 46K cells
    "pal_norm_total" # 6 - 48K cells
    "pal_norm_b1"    # 7 - 53K cells
    "reed"           # 8 - 803K cells (XL)
)

STUDY="${STUDIES[$SLURM_ARRAY_TASK_ID]}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../../config/load_paths.sh"
BASE_PATH="${SOURCE_IHBCAV1_EXTERNAL_STUDIES%/}"
SCRIPT="${BASE_PATH}/harmonization/scripts/prepare_study.R"
LOG_DIR="${BASE_PATH}/harmonization/logs"

CONTAINER="${CONTAINER_R_SPATIAL_4_3_3}"
R_LIBS_USER="${R_LIBS_R_SPATIAL_4_3_3}"

mkdir -p "$LOG_DIR"

echo "=============================================="
echo "Prepare Study (Phase 1): $STUDY"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task: $SLURM_ARRAY_TASK_ID"
echo "Node: $(hostname)"
echo "Start time: $(date)"
echo "Memory: $SLURM_MEM_PER_NODE"
echo "=============================================="

module load singularity

# Build extra args
EXTRA_ARGS=""
[ "${SKIP_JOINT:-0}" = "1" ] && EXTRA_ARGS="--skip-joint"

singularity exec \
    ${BIND_MOUNTS} \
    --bind "$R_LIBS_USER:/home/jovyan/R/library:ro" \
    --env "R_LIBS_USER=/home/jovyan/R/library" \
    "$CONTAINER" \
    Rscript "$SCRIPT" --study "$STUDY" $EXTRA_ARGS

echo ""
echo "=============================================="
echo "Completed: $STUDY (Phase 1)"
echo "End time: $(date)"
echo "=============================================="
