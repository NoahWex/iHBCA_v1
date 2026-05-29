#!/bin/bash
#SBATCH --job-name=patch_obs
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00

# Patch obs in integrated h5ads: add 21 L1 columns via h5py (no full load)
# Memory: ~2-5 GB per file (obs only, not X/layers)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
REPO_ROOT="${PROJECT_IHBCAV1_UPLOAD}"
SCRIPTS="${REPO_ROOT}/publication/scripts"
OUTDIR="${REPO_ROOT}/publication/outputs/integrated_objects"
CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"

echo "============================================="
echo "Patch obs: all 21 L1 columns via h5py"
echo "  Date: $(date)"
echo "  Node: $(hostname)"
echo "  Job ID: ${SLURM_JOB_ID}"
echo "============================================="

module load singularity

FILES=(all-breast-cells.h5ad breast-epithelial-lineage.h5ad breast-stromal-lineage.h5ad breast-immune-lineage.h5ad)

for F in "${FILES[@]}"; do
    echo ""
    echo "=== ${F} ==="
    H5AD="${OUTDIR}/${F}"
    if [ ! -f "${H5AD}" ]; then
        echo "ERROR: not found: ${H5AD}"
        exit 1
    fi
    echo "Size: $(ls -lh "${H5AD}" | awk '{print $5}')"

    singularity exec \
        --no-mount bind-paths \
        ${BIND_MOUNTS} \
        "$CONTAINER" \
        python3 "${SCRIPTS}/patch_obs_h5py.py" \
            --h5ad "$H5AD" \
            --repo-root "$REPO_ROOT" \
            --reorder

    echo "DONE: ${F}"
done

echo ""
echo "============================================="
echo "All patches COMPLETE"
echo "Date: $(date)"
echo "============================================="
