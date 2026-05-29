#!/bin/bash
#SBATCH --job-name=c1_assemble
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=180G
#SBATCH --time=02:00:00

# =============================================================================
# C1: Assemble iHBCA integrated h5ad from author share primitives
# =============================================================================
# Constructs AnnData from author share files (NPZ counts, gene_data.csv,
# annotations CSV, embedding CSVs). Enriches ethnicity, populates Tier 1
# fields, writes all-breast-cells.h5ad.
# Output: ~6 h5ad files in integrated_objects/
#
# Build mode: author-share (default) or legacy (set MODE=legacy below)
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
REPO_ROOT="${PROJECT_IHBCAV1_UPLOAD}"
source "${REPO_ROOT}/publication/run/load_config.sh"
PUB_ROOT="${REPO_ROOT}/publication"
IHBCA_ROOT="${SOURCE_IHBCA_REED_ROOT%/}"
PY_CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"

# Build mode: "author-share" (default) or "legacy"
MODE="${1:-author-share}"

# --- Author share inputs ---
ANNOTATIONS="${IHBCA_ROOT}/author_share/ihbca_level1.5_annotations.csv"
EMBEDDINGS="${IHBCA_ROOT}/author_share/X_scVI100.csv"
UMAP="${IHBCA_ROOT}/author_share/X_scVI100_UMAP.csv"
GENE_DATA="${IHBCA_ROOT}/author_share/scVI_100Dims/gene_data.csv"
COUNTS_NPZ="${IHBCA_ROOT}/author_share/preintegration_HBCA_inner_ENSEMBL_simple_counts_matrix.npz"

# --- Legacy input ---
H5AD="${IHBCA_ROOT}/integrated_atlas/integration_iHBCA.h5ad"

# Output (env-var override for alternate builds, e.g., OUTPUT_DIR=.../cxg_build/)
OUTPUT_DIR="${OUTPUT_DIR:-${OUTPUT_INTEGRATED_OBJECTS}}"

echo "=== C1 Integrated Object Assembly ==="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Memory: ${SLURM_MEM_PER_NODE}MB"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Mode: ${MODE}"
echo ""

# Verify inputs based on mode
if [ "$MODE" = "author-share" ]; then
    for INPUT_FILE in "$ANNOTATIONS" "$EMBEDDINGS" "$GENE_DATA"; do
        if [ ! -f "$INPUT_FILE" ]; then
            echo "ERROR: Input not found: $INPUT_FILE"
            exit 1
        fi
        echo "Input: $INPUT_FILE ($(ls -lh "$INPUT_FILE" | awk '{print $5}'))"
    done
    # Optional inputs
    for INPUT_FILE in "$UMAP" "$COUNTS_NPZ"; do
        if [ -f "$INPUT_FILE" ]; then
            echo "Input: $INPUT_FILE ($(ls -lh "$INPUT_FILE" | awk '{print $5}'))"
        else
            echo "Optional: $INPUT_FILE (not found, skipping)"
        fi
    done
else
    for INPUT_FILE in "$H5AD" "$ANNOTATIONS" "$EMBEDDINGS"; do
        if [ ! -f "$INPUT_FILE" ]; then
            echo "ERROR: Input not found: $INPUT_FILE"
            exit 1
        fi
        echo "Input: $INPUT_FILE ($(ls -lh "$INPUT_FILE" | awk '{print $5}'))"
    done
fi
echo ""

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Load singularity
module load singularity

# Build command based on mode
if [ "$MODE" = "author-share" ]; then
    EXTRA_ARGS="--gene-data $GENE_DATA"
    [ -f "$COUNTS_NPZ" ] && EXTRA_ARGS="$EXTRA_ARGS --counts-npz $COUNTS_NPZ"
    [ -f "$UMAP" ] && EXTRA_ARGS="$EXTRA_ARGS --umap $UMAP"
else
    EXTRA_ARGS="--h5ad $H5AD"
fi

# Run assembly
echo "Starting assembly..."
singularity exec \
    --no-mount bind-paths \
    ${BIND_MOUNTS} \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
    "$PY_CONTAINER" \
    python3 "${PUB_ROOT}/scripts/assemble_integrated.py" \
        $EXTRA_ARGS \
        --annotations "$ANNOTATIONS" \
        --embeddings "$EMBEDDINGS" \
        --repo-root "$REPO_ROOT" \
        --output-dir "$OUTPUT_DIR" \
        --schema-version "5.3.2" \
        --target "${TARGET:-hca}" \
        --skip-splits

echo ""
echo "=== Assembly Complete ==="
echo "Output dir: $OUTPUT_DIR"
ls -lh "$OUTPUT_DIR/"
echo "Exit code: $?"
