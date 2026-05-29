#!/bin/bash
#SBATCH --job-name=native_milo
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --array=0-8

# =============================================================================
# build_milo_native.sh
# =============================================================================
# Build Milo objects on native (study-specific) embeddings.
# Phase 2a of the parallel build pipeline.
#
# Usage:
#   sbatch build_milo_native.sh                              # All studies
#   sbatch --array=8 build_milo_native.sh                    # Reed only
#   sbatch --mem=128G --array=1 build_milo_native.sh         # Kumar (large)
#   sbatch --mem=128G --array=8 build_milo_native.sh         # Reed (large)
#
# Memory requirements:
#   - gray, murrow, pal_*: 64GB sufficient
#   - nee, twigger: 64GB recommended
#   - kumar, reed: 128GB recommended (override with --mem=128G)
# =============================================================================

set -euo pipefail

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

STUDY="${STUDIES[$SLURM_ARRAY_TASK_ID]}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../../config/load_paths.sh"
BASE_PATH="${SOURCE_IHBCAV1_EXTERNAL_STUDIES%/}"
SCRIPT="${BASE_PATH}/harmonization/scripts/build_milo_worker.R"
LOG_DIR="${BASE_PATH}/harmonization/logs"
STUDY_DIR="${BASE_PATH}/harmonization/outputs/study_objects/${STUDY}"

CONTAINER="${CONTAINER_R_SPATIAL_4_3_3}"
R_LIBS_USER="${R_LIBS_R_SPATIAL_4_3_3}"

mkdir -p "$LOG_DIR"

echo "=============================================="
echo "Build Native Milo: $STUDY"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task: $SLURM_ARRAY_TASK_ID"
echo "Node: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Start time: $(date)"
echo "=============================================="

if [[ ! -f "${STUDY_DIR}/seurat.rds" ]]; then
    echo "ERROR: seurat.rds not found in ${STUDY_DIR}"
    echo "Run prepare_study.R --study ${STUDY} first (Phase 1)"
    exit 1
fi

module load singularity

singularity exec \
    ${BIND_MOUNTS} \
    --bind "$R_LIBS_USER:/home/jovyan/R/library:ro" \
    --env "R_LIBS_USER=/home/jovyan/R/library" \
    "$CONTAINER" \
    Rscript "$SCRIPT" --study "$STUDY" --mode native --cpus "$SLURM_CPUS_PER_TASK"

echo ""
echo "=============================================="
echo "Completed: $STUDY"
echo "End time: $(date)"
echo "=============================================="
