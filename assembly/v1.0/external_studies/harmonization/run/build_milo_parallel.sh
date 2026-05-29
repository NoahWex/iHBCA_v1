#!/bin/bash
#SBATCH --job-name=build_milo
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --array=0-17

# =============================================================================
# build_milo_parallel.sh — Phase 2 of parallel Milo build pipeline
# =============================================================================
# Builds Milo objects with BiocParallel multi-core k-NN.
# Array 0-17: two tasks per study (even=native, odd=joint).
#
# Usage:
#   sbatch build_milo_parallel.sh                              # All (18 tasks)
#   sbatch --array=0-1 build_milo_parallel.sh                  # gray only
#   sbatch --mem=280G --time=3-00:00:00 --cpus-per-task=16 \
#          --array=2-3,16-17 build_milo_parallel.sh            # XL: kumar + reed
#
# Index mapping:
#   0=gray/native  1=gray/joint   2=kumar/native  3=kumar/joint
#   4=murrow/native 5=murrow/joint 6=nee/native    7=nee/joint
#   8=twigger/native 9=twigger/joint
#   10=pal_epi/native 11=pal_epi/joint
#   12=pal_total/native 13=pal_total/joint
#   14=pal_b1/native 15=pal_b1/joint
#   16=reed/native  17=reed/joint
# =============================================================================

set -euo pipefail

STUDIES=(
    "gray"           # 0-1
    "kumar"          # 2-3
    "murrow"         # 4-5
    "nee"            # 6-7
    "twigger"        # 8-9
    "pal_norm_epi"   # 10-11
    "pal_norm_total" # 12-13
    "pal_norm_b1"    # 14-15
    "reed"           # 16-17
)

STUDY_IDX=$((SLURM_ARRAY_TASK_ID / 2))
MODE_IDX=$((SLURM_ARRAY_TASK_ID % 2))
STUDY="${STUDIES[$STUDY_IDX]}"
MODE=$( [ $MODE_IDX -eq 0 ] && echo "native" || echo "joint" )

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../../config/load_paths.sh"
BASE_PATH="${SOURCE_IHBCAV1_EXTERNAL_STUDIES%/}"
SCRIPT="${BASE_PATH}/harmonization/scripts/build_milo_worker.R"
LOG_DIR="${BASE_PATH}/harmonization/logs"

CONTAINER="${CONTAINER_R_SPATIAL_4_3_3}"
R_LIBS_USER="${R_LIBS_R_SPATIAL_4_3_3}"

mkdir -p "$LOG_DIR"

echo "=============================================="
echo "Build Milo (Phase 2): $STUDY [$MODE]"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task: $SLURM_ARRAY_TASK_ID (study=$STUDY_IDX, mode=$MODE_IDX)"
echo "Node: $(hostname)"
echo "Start time: $(date)"
echo "Memory: $SLURM_MEM_PER_NODE"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "=============================================="

# Dry-run support
DRY_RUN_FLAG=""
[ "${DRY_RUN:-0}" = "1" ] && DRY_RUN_FLAG="--dry-run"

module load singularity

singularity exec \
    ${BIND_MOUNTS} \
    --bind "$R_LIBS_USER:/home/jovyan/R/library:ro" \
    --env "R_LIBS_USER=/home/jovyan/R/library" \
    "$CONTAINER" \
    Rscript "$SCRIPT" --study "$STUDY" --mode "$MODE" \
        --cpus "$SLURM_CPUS_PER_TASK" $DRY_RUN_FLAG

echo ""
echo "=============================================="
echo "Completed: $STUDY [$MODE] (Phase 2)"
echo "End time: $(date)"
echo "=============================================="
