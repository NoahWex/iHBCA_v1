#!/bin/bash
#SBATCH --job-name=plot_umaps
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --array=0-5
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=00:30:00

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
REPO_ROOT="${PROJECT_IHBCAV1_UPLOAD}"
SCRIPTS="${REPO_ROOT}/publication/scripts"
CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"

STUDIES=(gray kumar murrow nee twigger reed)
INDEX="${SLURM_ARRAY_TASK_ID:?Must run as SLURM array job}"
STUDY="${STUDIES[$INDEX]}"

echo "Plotting enriched UMAP: ${STUDY}"

module load singularity

singularity exec \
    --no-mount bind-paths \
    ${BIND_MOUNTS} \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
    "$CONTAINER" \
    python3 "${SCRIPTS}/plot_enriched_umaps.py" \
        --study "$STUDY" \
        --repo-root "$REPO_ROOT"

echo "Plot COMPLETE: ${STUDY}"
