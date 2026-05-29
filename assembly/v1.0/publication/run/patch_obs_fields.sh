#!/bin/bash
#SBATCH --job-name=patch_obs_fields
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --time=00:30:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1

# Patch obs columns in integrated + sketch h5ads:
#   1. facs_status: use FACS_status for murrow/pal/reed (sample-level source is correct)
#   2. n_genes, n_counts, percent_mito: categorical-of-strings → float32

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
PROJECT="${PROJECT_IHBCAV1_UPLOAD}"
SCRIPT=$PROJECT/publication/scripts/patch_obs_fields.py
INTEGRATED=$PROJECT/publication/outputs/integrated_objects/all-breast-cells.h5ad
SKETCH=$PROJECT/publication/outputs/integrated_objects/all-breast-cells-sketch.h5ad

CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"
BINDS="${BIND_MOUNTS}"

module load singularity

echo "=== Patch obs fields: facs_status + numeric dtypes ==="
echo "Date: $(date)"
echo "Job: $SLURM_JOB_ID"

singularity exec $BINDS \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    "$CONTAINER" \
    python3 "$SCRIPT" "$INTEGRATED" "$SKETCH"

echo ""
echo "=== Done ==="
echo "Date: $(date)"
