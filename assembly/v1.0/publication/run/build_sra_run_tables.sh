#!/usr/bin/env bash
# =============================================================================
# Download SRA metadata and build standardized run tables
# =============================================================================
# Two-step pipeline:
#   Step 1: Download raw ENA/SDRF metadata (needs internet → login node)
#   Step 2: Parse into per-study CSVs (needs pandas → container)
#
# Usage:
#   # Run both steps on login node (or interactive session with internet):
#   bash publication/run/build_sra_run_tables.sh
#
#   # Or as SLURM job (Step 2 only, after manual download):
#   sbatch publication/run/build_sra_run_tables.sh
#
# =============================================================================
#SBATCH --job-name=build_sra_tables
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
REPO_ROOT="${PROJECT_IHBCAV1_UPLOAD}"
RAW_DIR="$REPO_ROOT/publication/mappings/sra_raw"
CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"

# Step 1: Download (skip if raw files already exist)
if [ ! -d "$RAW_DIR" ] || [ -z "$(ls -A "$RAW_DIR" 2>/dev/null)" ]; then
    echo "=== Step 1: Downloading raw SRA/SDRF metadata ==="
    bash "$REPO_ROOT/publication/run/download_sra_metadata.sh" "$REPO_ROOT"
else
    echo "=== Step 1: Raw files already exist in $RAW_DIR, skipping download ==="
    ls -lh "$RAW_DIR/"
fi

echo ""
echo "=== Step 2: Parsing into standardized run tables ==="

module load singularity 2>/dev/null || true

singularity exec \
    ${BIND_MOUNTS} \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
    "$CONTAINER" \
    python "$REPO_ROOT/publication/scripts/build_sra_run_tables.py" \
        --repo-root "$REPO_ROOT"

echo ""
echo "=== Complete ==="
