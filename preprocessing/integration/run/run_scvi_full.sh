#!/bin/bash
#SBATCH --job-name=scvi_full
#SBATCH --account=dalawson_lab
#SBATCH --partition=free-gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --gres=gpu:1
#SBATCH --time=36:00:00
#SBATCH --array=0-2
#SBATCH --output=${HPC_LOGS_DIR}/%x_%A_%a.out
#SBATCH --error=${HPC_LOGS_DIR}/%x_%A_%a.err

# =============================================================================
# Per-compartment scVI training on the iHBCA v1 atlas (3 GPUs in parallel).
# =============================================================================
# Trains a scVI model with negative binomial likelihood, n_latent=N_LATENT,
# n_layers=3, max_epochs=300, on the 4000 most highly variable genes per
# compartment (Immune / Epithelial / Stromal). The trained model + 50D latent
# representation feed scANVI fine-tuning (see run_scanvi_full.sh).
#
# Resources: Tier 5+GPU (128G / 8 CPU / 36h / 1 A100) — the 36h budget covers
# the largest compartment (Stromal, ~974K cells); Immune (~162K) and
# Epithelial (~991K) finish well within budget.
#
# Submit pattern:
#   source publication/config/load_paths.sh
#   sbatch publication/preprocessing/integration/run/run_scvi_full.sh
#
# (load_paths.sh must be sourced BEFORE sbatch so #SBATCH log paths resolve.)
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../config/load_paths.sh"

INTEGRATION_ROOT="${PROJECT_PUBLICATION}/preprocessing/integration"
SCVI_SCRIPT="${INTEGRATION_ROOT}/scripts/01_compartment_scvi_full.py"

COMPONENTS="${SOURCE_IHBCAV1_ASSEMBLY_COMPONENTS}"
GENE_DATA="${COMPONENTS}gene_data.csv"

COMPARTMENTS=("Immune" "Epithelial" "Stromal")
SHORTS=("imm" "epi" "str")
COMP="${COMPARTMENTS[$SLURM_ARRAY_TASK_ID]}"
SHORT="${SHORTS[$SLURM_ARRAY_TASK_ID]}"

COUNTS_NPZ="${COMPONENTS}${SHORT}_counts.npz"
# Prefer enriched metadata (with native_* columns) when available; the scVI
# script handles the plain fallback by skipping native label columns.
ENRICHED_META="${COMPONENTS}${SHORT}_metadata_enriched.csv"
PLAIN_META="${COMPONENTS}${SHORT}_metadata.csv"
if [ -f "$ENRICHED_META" ]; then
    METADATA="$ENRICHED_META"
    echo "Using enriched metadata (with native_* columns)"
else
    METADATA="$PLAIN_META"
    echo "Using plain metadata (no native_* columns)"
fi

N_LATENT="${N_LATENT:-50}"
OUTPUT_DIR="${SOURCE_IHBCAV1_ASSEMBLY_ROOT}scvi/${SHORT}/n_latent_${N_LATENT}"
mkdir -p "$OUTPUT_DIR"

echo "=== scVI: $COMP n_latent=$N_LATENT (task $SLURM_ARRAY_TASK_ID) ==="
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "GPU: ${CUDA_VISIBLE_DEVICES:-unset}"
echo "Counts: $COUNTS_NPZ"
echo "Gene data: $GENE_DATA"
echo "Metadata: $METADATA"
echo "Output: $OUTPUT_DIR"
echo ""

if [ ! -f "$COUNTS_NPZ" ]; then
    echo "ERROR: $COUNTS_NPZ not found"
    echo "Available NPZs in ${COMPONENTS}:"
    ls -lh "${COMPONENTS}"*_counts.npz 2>&1 || echo "(none found)"
    exit 1
fi

for f in "$GENE_DATA" "$METADATA"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: $f not found"
        exit 1
    fi
done

# scVI reads gene_df["symbol"] for human-readable gene names in the output h5ad.
if ! head -1 "$GENE_DATA" | grep -q "symbol"; then
    echo "ERROR: $GENE_DATA does not have 'symbol' column"
    echo "Header: $(head -1 "$GENE_DATA")"
    exit 1
fi

module load singularity

singularity exec \
    --nv \
    ${BIND_MOUNTS} \
    --bind "${PYTHON_LIBS_SCGPT_GPU_2025Q3}:/home/jovyan/python_user:ro" \
    --env "PYTHONUSERBASE=/home/jovyan/python_user" \
    --env "PYTHONUNBUFFERED=1" \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
    "${CONTAINER_PYTHON_SCVI_GPU_2025Q3}" \
    python "$SCVI_SCRIPT" \
        --counts-npz "$COUNTS_NPZ" \
        --gene-data "$GENE_DATA" \
        --metadata "$METADATA" \
        --compartment "$COMP" \
        --n-latent "$N_LATENT" \
        --n-hvg 4000 \
        --seed 42 \
        --output-dir "$OUTPUT_DIR"

echo ""
echo "=== Complete: $(date) ==="
ls -lh "$OUTPUT_DIR/"
