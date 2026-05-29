#!/bin/bash
#SBATCH --job-name=patch_fif_raw
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00

# Remove feature_is_filtered from raw/var in restructured h5ads.
# HCA schema prohibits this column in raw.var.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
PROJECT="${PROJECT_IHBCAV1_UPLOAD}"
SCRIPT="$PROJECT/publication/scripts/patch_remove_feature_is_filtered_raw.py"
INTEGRATED="$PROJECT/publication/outputs/integrated_objects/all-breast-cells.h5ad"
SKETCH="$PROJECT/publication/outputs/integrated_objects/all-breast-cells-sketch.h5ad"

CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"
BINDS="${BIND_MOUNTS}"

module load singularity

echo "=== Patch: remove feature_is_filtered from raw/var ==="
echo "Date: $(date)"
echo "Job: $SLURM_JOB_ID"
echo ""

singularity exec $BINDS \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    "$CONTAINER" \
    python3 "$SCRIPT" "$INTEGRATED" "$SKETCH"

echo ""
echo "=== Done ==="
echo "Date: $(date)"
