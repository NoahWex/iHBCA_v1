#!/bin/bash
#SBATCH --job-name=snapshot_baseline
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=01:00:00

# =============================================================================
# Snapshot baseline manifest for pipeline_rebuild regression verification
# =============================================================================
# Reads all 11 h5ads (7 source + 4 integrated) in backed mode, records
# SHA256, cell/gene counts, obs columns/dtypes, obsm keys.
# Output: publication/outputs/baseline_manifest.yaml
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
REPO_ROOT="${PROJECT_IHBCAV1_UPLOAD}"
PUB_ROOT="${REPO_ROOT}/publication"
PY_CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"
OUTPUT="${PUB_ROOT}/outputs/baseline_manifest.yaml"

echo "=== Snapshot Baseline Manifest ==="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo ""

# Verify h5ads exist
echo "Checking source datasets..."
for f in gray2022 kumar2023 murrow2022 nee2023 twigger2022 reed2024 pal2021; do
    FILE="${PUB_ROOT}/outputs/source_datasets/${f}.h5ad"
    if [ -f "$FILE" ]; then
        echo "  OK: ${f}.h5ad ($(ls -lh "$FILE" | awk '{print $5}'))"
    else
        echo "  MISSING: ${f}.h5ad"
    fi
done

echo "Checking integrated objects..."
for f in all-breast-cells breast-epithelial-lineage breast-stromal-lineage breast-immune-lineage; do
    FILE="${PUB_ROOT}/outputs/integrated_objects/${f}.h5ad"
    if [ -f "$FILE" ]; then
        echo "  OK: ${f}.h5ad ($(ls -lh "$FILE" | awk '{print $5}'))"
    else
        echo "  MISSING: ${f}.h5ad"
    fi
done
echo ""

# Load singularity
module load singularity

# Run snapshot
echo "Running snapshot_baseline.py..."
singularity exec \
    --no-mount bind-paths \
    ${BIND_MOUNTS} \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    "$PY_CONTAINER" \
    python3 "${PUB_ROOT}/scripts/snapshot_baseline.py" \
        --repo-root "$REPO_ROOT" \
        --output "$OUTPUT"

echo ""
echo "=== Snapshot Complete ==="
echo "Output: $OUTPUT"
echo "Exit code: $?"
