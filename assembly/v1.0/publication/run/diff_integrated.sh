#!/bin/bash
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --time=02:00:00
#SBATCH --job-name=diff_integrated

# Diff integrated object vs CxG published h5ad (OOM'd at 180GB, retrying at 256GB)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
REPO_ROOT="${PROJECT_IHBCAV1_UPLOAD}"
SCRIPTS="$REPO_ROOT/publication/scripts"
REPORT_DIR="$REPO_ROOT/publication/outputs/validation_reports"
CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"

NEW="$REPO_ROOT/publication/outputs/integrated_objects/all-breast-cells.h5ad"
REF="${SOURCE_INTEGRATED_ATLAS_REED}integration_iHBCA.h5ad"

echo "DIFF: integrated (all-breast-cells vs CxG published)"
echo "  New: $(ls -lh "$NEW" | awk '{print $5}')"
echo "  Ref: $(ls -lh "$REF" | awk '{print $5}')"
echo "  Memory: ${SLURM_MEM_PER_NODE}MB"
echo "  Date: $(date)"

module load singularity

singularity exec \
    --no-mount bind-paths \
    ${BIND_MOUNTS} \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
    "$CONTAINER" \
    python3 "$SCRIPTS/diff_h5ads.py" \
        --new "$NEW" \
        --reference "$REF" \
        --output "$REPORT_DIR/integrated_regen_diff.yaml"

echo ""
echo "Done: $(date)"
