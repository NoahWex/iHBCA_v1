#!/bin/bash
#SBATCH --job-name=hca_upload
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00

# =============================================================================
# HCA Atlas Upload — 9 objects to s3://hca-atlas-tracker-data/breast/breast-v1/
#
# Objects:
#   source-datasets (7):  gray2022 kumar2023 murrow2022 nee2023 twigger2022 reed2024 pal2021
#   integrated-objects (2): all-breast-cells.h5ad, all-breast-cells-sketch.h5ad
#
# Total: ~165.5 GB via s5cmd
# =============================================================================

set -euo pipefail

# Ensure user-installed tools (aws, hca-smart-sync, s5cmd) are in PATH
export PATH="$HOME/.local/bin:$PATH"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
STAGING="${SOURCE_IHBCAV1_UPLOAD_STAGING}"

echo "=== HCA Atlas Upload ==="
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "Staging: $STAGING"
echo ""

module load python/3.10.2

# Verify tooling
echo "--- Tooling check ---"
which hca-smart-sync
which s5cmd
aws sts get-caller-identity --profile hca-tracker-upload | head -5
hca-smart-sync config show
echo ""

# --- Source datasets (66 GB) ---
echo "=========================================="
echo "Uploading source-datasets (7 files, ~66 GB)"
echo "Started: $(date)"
echo "=========================================="

cd "$STAGING"
echo "y" | hca-smart-sync sync source-datasets --local-path source-datasets --verbose

echo ""
echo "Source datasets upload complete: $(date)"
echo ""

# --- Integrated objects (99 GB) ---
echo "=========================================="
echo "Uploading integrated-objects (2 files, ~99 GB)"
echo "Started: $(date)"
echo "=========================================="

echo "y" | hca-smart-sync sync integrated-objects --local-path integrated-objects --verbose

echo ""
echo "Integrated objects upload complete: $(date)"
echo ""

echo "=== Upload complete ==="
echo "Finished: $(date)"
echo ""
echo "Next: verify all 9 objects on HCA Atlas Tracker"
