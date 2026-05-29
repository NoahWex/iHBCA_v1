#!/bin/bash
#SBATCH --job-name=reed_download
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00

# =============================================================================
# Download Reed H5AD from CELLxGENE
# =============================================================================
# Dataset: 44c8e245-4a1c-48d2-9ba1-59e3d94359bc
# Expected: ~803K cells, ~3-5GB H5AD
#
# Usage:
#   sbatch download_reed.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../../../../config/load_paths.sh"


DEST_DIR="${SOURCE_COMPONENT_STUDIES}"
DEST_FILE="${DEST_DIR}/reed.h5ad"
URL="https://datasets.cellxgene.cziscience.com/44c8e245-4a1c-48d2-9ba1-59e3d94359bc.h5ad"

LOG_DIR="${SOURCE_IHBCAV1_HARMONIZATION%/}/studies/reed/run/logs"
mkdir -p "$LOG_DIR"

echo "=============================================="
echo "REED H5AD DOWNLOAD"
echo "=============================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start time: $(date)"
echo "URL: ${URL}"
echo "Destination: ${DEST_FILE}"
echo "=============================================="

if [[ -f "$DEST_FILE" ]]; then
    echo "WARNING: ${DEST_FILE} already exists ($(du -h "$DEST_FILE" | cut -f1))"
    echo "Renaming to reed.h5ad.prev before download"
    mv "$DEST_FILE" "${DEST_FILE}.prev"
fi

echo ""
echo "Downloading..."
curl -L -o "$DEST_FILE" "$URL"

echo ""
echo "Download complete."
echo "File size: $(du -h "$DEST_FILE" | cut -f1)"
echo "MD5: $(md5sum "$DEST_FILE" | cut -d' ' -f1)"
echo ""
echo "=============================================="
echo "End time: $(date)"
echo "=============================================="
