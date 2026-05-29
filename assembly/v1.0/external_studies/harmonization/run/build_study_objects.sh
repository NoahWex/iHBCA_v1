#!/bin/bash
#SBATCH --job-name=build_study
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --array=0-8

# =============================================================================
# build_study_objects.sh
# =============================================================================
# Build study objects for Phase B analysis.
# Runs as array job: one task per study.
#
# Usage:
#   sbatch build_study_objects.sh           # Run all studies
#   sbatch --array=0 build_study_objects.sh # Run gray only (index 0)
#   sbatch --array=1 build_study_objects.sh # Run kumar only (index 1)
#
# Memory requirements:
#   - gray, murrow, pal_*: 32GB sufficient
#   - nee, twigger: 64GB recommended
#   - kumar, reed: 128GB recommended (override with --mem=128G)
#
# For large studies:
#   sbatch --mem=128G --array=1 build_study_objects.sh  # kumar
#   sbatch --mem=128G --array=8 build_study_objects.sh  # reed
# =============================================================================

set -euo pipefail

# Study array (order matches SLURM_ARRAY_TASK_ID)
STUDIES=(
    "gray"           # 0 - 52K cells
    "kumar"          # 1 - 714K cells (large)
    "murrow"         # 2 - 86K cells
    "nee"            # 3 - 230K cells
    "twigger"        # 4 - 111K cells
    "pal_norm_epi"   # 5 - 46K cells
    "pal_norm_total" # 6 - 48K cells
    "pal_norm_b1"    # 7 - 53K cells
    "reed"           # 8 - 803K cells (large)
)

# Get study for this task
STUDY="${STUDIES[$SLURM_ARRAY_TASK_ID]}"

# Checkpoint: resume from step N (default: full pipeline)
START_STEP="${START_STEP:-1}"

# Paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../../config/load_paths.sh"
BASE_PATH="${SOURCE_IHBCAV1_EXTERNAL_STUDIES%/}"
SCRIPT="${BASE_PATH}/harmonization/scripts/prepare_study.R"
LOG_DIR="${BASE_PATH}/harmonization/logs"

# Container and R library paths
CONTAINER="${CONTAINER_R_SPATIAL_4_3_3}"
R_LIBS_USER="${R_LIBS_R_SPATIAL_4_3_3}"

# Create log directory
mkdir -p "$LOG_DIR"

echo "=============================================="
echo "Build Study Objects: $STUDY"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task: $SLURM_ARRAY_TASK_ID"
echo "Node: $(hostname)"
echo "Start time: $(date)"
echo "Memory: $SLURM_MEM_PER_NODE"
echo "Start step: $START_STEP"
echo "=============================================="

# Load singularity
module load singularity

# Run R script in container
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
