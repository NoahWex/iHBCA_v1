#!/bin/bash
#SBATCH --job-name=patch_uns_title
#SBATCH --account=dalawson_lab
#SBATCH --partition=free
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00

# Set uns["title"] = "Integrated Human Breast Cell Atlas V1"
# in both integrated objects in upload-staging.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
REPO_ROOT="${PROJECT_IHBCAV1_UPLOAD}"
SCRIPTS="${REPO_ROOT}/publication/scripts"
STAGING="${REPO_ROOT}/upload-staging/integrated-objects"
CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"

TITLE="iHBCA v1"

echo "============================================="
echo "Patch uns[title] — integrated objects"
echo "  Date:    $(date)"
echo "  Node:    $(hostname)"
echo "  Job ID:  ${SLURM_JOB_ID}"
echo "  Title:   ${TITLE}"
echo "============================================="

module load singularity

singularity exec \
    --no-mount bind-paths \
    ${BIND_MOUNTS} \
    "$CONTAINER" \
    python3 "${SCRIPTS}/patch_uns_title.py" \
        --h5ad \
            "${STAGING}/all-breast-cells.h5ad" \
            "${STAGING}/all-breast-cells-sketch.h5ad" \
        --title "${TITLE}"

echo ""
echo "============================================="
echo "Patch COMPLETE: $(date)"
echo "Next: re-upload integrated-objects to S3"
echo "============================================="
