#!/bin/bash
#SBATCH --job-name=ihbca_anno_str_annotation
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --output=${HPC_LOGS_DIR}/%x_%j.out
#SBATCH --error=${HPC_LOGS_DIR}/%x_%j.err

# =============================================================================
# Render the post-assignment annotation package for the stromal compartment.
#
# Tier 5 (128G / 8 CPU / 4h) per hpc-resource-rules.md — 974K cells x ~37K
# genes during one-vs-rest DE; runtime scales with label count.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../config/load_paths.sh"

COMPARTMENT="str"
OUT_DIR="${PROJECT_PUBLICATION}/analysis/annotation/${COMPARTMENT}"
STAGE_OUT="${OUT_DIR}/evidence/"
SCANVI_DIR="${SOURCE_IHBCAV1_ASSEMBLY_SCANVI}${COMPARTMENT}/n_latent_50"

V2_YAML="${SCRIPT_DIR}/../yamls/annotation_v2_${COMPARTMENT}.yaml"
LABELS="${OUT_DIR}/labels.csv"
LEIDEN="${SCANVI_DIR}/${COMPARTMENT}_leiden_barcoded.csv"
UMAP="${SCANVI_DIR}/${COMPARTMENT}_scanvi_umap_barcoded.csv"
COUNTS="${SOURCE_IHBCAV1_ASSEMBLY_COMPONENTS}${COMPARTMENT}_counts.npz"
GENE_DATA="${SOURCE_IHBCAV1_ASSEMBLY_COMPONENTS}gene_data.csv"
META="${SOURCE_IHBCAV1_ASSEMBLY_COMPONENTS}${COMPARTMENT}_metadata_enriched.csv"
CONFIG="${SCRIPT_DIR}/../config/compartment_config.yaml"
SCRIPT="${SCRIPT_DIR}/../scripts/03_render_annotation_package.py"

mkdir -p "${STAGE_OUT}" "${PROJECT_PUBLICATION}/logs"

for f in "$V2_YAML" "$LABELS" "$LEIDEN" "$UMAP" "$COUNTS" "$GENE_DATA" "$META" "$CONFIG" "$SCRIPT"; do
    [ -f "$f" ] || { echo "ERROR: missing input: $f"; exit 1; }
done

module load singularity

singularity exec --no-home \
    ${BIND_MOUNTS} \
    --env "PYTHONUSERBASE=${PYTHON_LIBS_SCVI_GPU_2025Q3_USERBASE}" \
    --env "PYTHONUNBUFFERED=1" \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
    "${CONTAINER_PYTHON_SCVI_GPU_2025Q3}" \
    python "$SCRIPT" \
        --v2-yaml "$V2_YAML" \
        --labels-csv "$LABELS" \
        --leiden-csv "$LEIDEN" \
        --umap-csv "$UMAP" \
        --counts-npz "$COUNTS" \
        --gene-data "$GENE_DATA" \
        --metadata-csv "$META" \
        --config "$CONFIG" \
        --out-dir "$OUT_DIR"

echo "=== Complete ==="
