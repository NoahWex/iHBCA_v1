#!/bin/bash
#SBATCH --job-name=diff_baseline
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=180G
#SBATCH --time=02:00:00

# =============================================================================
# Diff current h5ads against baseline manifest (structural comparison)
# =============================================================================
# Compares cell counts, gene counts, obs columns/dtypes, obsm keys, uns keys
# against the Phase 0 baseline manifest.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
REPO_ROOT="${PROJECT_IHBCAV1_UPLOAD}"
PUB_ROOT="${REPO_ROOT}/publication"
PY_CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"
MANIFEST="${PUB_ROOT}/outputs/baseline_manifest.yaml"

echo "=== Diff Against Baseline ==="
echo "Date: $(date)"
echo "Job ID: ${SLURM_JOB_ID}"
echo ""

module load singularity

# Diff integrated objects only (source datasets were not rebuilt)
singularity exec \
    --no-mount bind-paths \
    ${BIND_MOUNTS} \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    "$PY_CONTAINER" \
    python3 "${PUB_ROOT}/scripts/diff_against_baseline.py" \
        --manifest "$MANIFEST" \
        --repo-root "$REPO_ROOT" \
        --datasets all-breast-cells breast-epithelial-lineage breast-stromal-lineage breast-immune-lineage

echo ""
echo "Exit code: $?"
