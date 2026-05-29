#!/bin/bash
#SBATCH --job-name=patch_var_uns
#SBATCH --account=dalawson_lab
#SBATCH --partition=free
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00

# Fix CxG 5.3.2 validation blockers on 8 h5ads (7 source + all-breast-cells)
# - Delete uns["layer_descriptions"] (deprecated)
# - Rename var["feature_biotype"] -> var["feature_biotype_gencode"] (reserved)
# Lineage splits excluded — will be regenerated from parent.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
REPO_ROOT="${PROJECT_IHBCAV1_UPLOAD}"
SCRIPTS="${REPO_ROOT}/publication/scripts"
SOURCE_DIR="${REPO_ROOT}/publication/outputs/source_datasets"
INTEGRATED_DIR="${REPO_ROOT}/publication/outputs/integrated_objects"
CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"

echo "============================================="
echo "Patch var/uns: CxG 5.3.2 fixes"
echo "  Date: $(date)"
echo "  Node: $(hostname)"
echo "  Job ID: ${SLURM_JOB_ID}"
echo "============================================="

module load singularity

FILES=(
    "${SOURCE_DIR}/gray2022.h5ad"
    "${SOURCE_DIR}/kumar2023.h5ad"
    "${SOURCE_DIR}/murrow2022.h5ad"
    "${SOURCE_DIR}/nee2023.h5ad"
    "${SOURCE_DIR}/twigger2022.h5ad"
    "${SOURCE_DIR}/reed2024.h5ad"
    "${SOURCE_DIR}/pal2021.h5ad"
    "${INTEGRATED_DIR}/all-breast-cells.h5ad"
)

singularity exec \
    --no-mount bind-paths \
    ${BIND_MOUNTS} \
    "$CONTAINER" \
    python3 "${SCRIPTS}/patch_var_uns_h5py.py" \
        --h5ad "${FILES[@]}"

echo ""
echo "============================================="
echo "All patches COMPLETE"
echo "Date: $(date)"
echo "============================================="
