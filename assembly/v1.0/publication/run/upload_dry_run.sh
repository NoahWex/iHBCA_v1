#!/bin/bash
#SBATCH --job-name=upload_dry_run
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00

# =============================================================================
# Upload Dry Run: metadata verification + staging directory setup
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
PROJECT_ROOT="${PROJECT_IHBCAV1_UPLOAD}"
PYTHON_CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"

echo "=== Upload Dry Run ==="
echo "Date: $(date)"
echo "Host: $(hostname)"
echo ""

# ---------------------------------------------------------------------------
# Step 2: Verify tracker metadata alignment
# ---------------------------------------------------------------------------
echo ">>> Step 2: Verify tracker metadata alignment"

module load singularity

# Run metadata check — capture exit code but don't abort on mismatch
singularity exec \
    ${BIND_MOUNTS} \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
    "$PYTHON_CONTAINER" \
    python "$PROJECT_ROOT/publication/scripts/verify_tracker_metadata.py" \
        --project-root "$PROJECT_ROOT" || echo "(metadata check reported mismatches — see above)"

echo ""

# ---------------------------------------------------------------------------
# Step 6: Create upload staging directory with absolute symlinks
# ---------------------------------------------------------------------------
echo ">>> Step 6: Create upload staging directory"

STAGING="$PROJECT_ROOT/upload-staging"
mkdir -p "$STAGING/source-datasets"
mkdir -p "$STAGING/integrated-objects"

# Source datasets — absolute symlinks
for study in gray2022 kumar2023 murrow2022 nee2023 twigger2022 reed2024 pal2021; do
    src="$PROJECT_ROOT/publication/outputs/source_datasets/${study}.h5ad"
    dst="$STAGING/source-datasets/${study}.h5ad"
    if [ -L "$dst" ]; then
        echo "  Symlink exists: $dst → $(readlink "$dst")"
    elif [ -e "$src" ]; then
        ln -s "$src" "$dst"
        echo "  Created: $dst → $src"
    else
        echo "  ERROR: Source missing: $src"
    fi
done

# Integrated objects — absolute symlinks
for obj in all-breast-cells.h5ad all-breast-cells-sketch.h5ad; do
    src="$PROJECT_ROOT/publication/outputs/integrated_objects/$obj"
    dst="$STAGING/integrated-objects/$obj"
    if [ -L "$dst" ]; then
        echo "  Symlink exists: $dst → $(readlink "$dst")"
    elif [ -e "$src" ]; then
        ln -s "$src" "$dst"
        echo "  Created: $dst → $src"
    else
        echo "  ERROR: Source missing: $src"
    fi
done

echo ""

# ---------------------------------------------------------------------------
# Step 7: Dry run verification — check all files resolve
# ---------------------------------------------------------------------------
echo ">>> Step 7: Dry run verification"

echo ""
echo "Source datasets:"
ls -lh "$STAGING/source-datasets/"

echo ""
echo "Integrated objects:"
ls -lh "$STAGING/integrated-objects/"

echo ""
echo "File count verification:"
SRC_COUNT=$(ls -1 "$STAGING/source-datasets/"*.h5ad 2>/dev/null | wc -l)
INT_COUNT=$(ls -1 "$STAGING/integrated-objects/"*.h5ad 2>/dev/null | wc -l)
TOTAL=$((SRC_COUNT + INT_COUNT))
echo "  Source datasets: $SRC_COUNT (expected 7)"
echo "  Integrated objects: $INT_COUNT (expected 2)"
echo "  Total: $TOTAL (expected 9)"

echo ""
echo "Symlink resolution check:"
BROKEN=0
for f in "$STAGING"/source-datasets/*.h5ad "$STAGING"/integrated-objects/*.h5ad; do
    if [ -L "$f" ] && [ ! -e "$f" ]; then
        echo "  BROKEN: $f → $(readlink "$f")"
        BROKEN=$((BROKEN + 1))
    fi
done
if [ "$BROKEN" -eq 0 ]; then
    echo "  All symlinks resolve OK"
else
    echo "  ERROR: $BROKEN broken symlinks"
fi

echo ""
echo "File sizes (dereferenced):"
du -Lh "$STAGING/source-datasets/"*.h5ad "$STAGING/integrated-objects/"*.h5ad 2>/dev/null | sort -k1 -h
echo ""
TOTAL_SIZE=$(du -Lsh "$STAGING" 2>/dev/null | cut -f1)
echo "Total staging size: $TOTAL_SIZE"

echo ""
echo "=== Dry run complete ==="
echo "Date: $(date)"
