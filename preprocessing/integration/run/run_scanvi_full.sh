#!/bin/bash
#SBATCH --job-name=scanvi_full
#SBATCH --account=dalawson_lab
#SBATCH --partition=free-gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --gres=gpu:1
#SBATCH --time=08:00:00
#SBATCH --array=0-2
#SBATCH --output=${HPC_LOGS_DIR}/%x_%A_%a.out
#SBATCH --error=${HPC_LOGS_DIR}/%x_%A_%a.err

# =============================================================================
# Per-compartment scANVI fine-tune on top of pre-trained scVI models.
# =============================================================================
# Loads the scVI model from outputs/scvi/{comp}/n_latent_{N}/scvi_model/ and
# fine-tunes it as a semi-supervised scANVI model using L1 cell-type labels
# (level1_annotation). Produces the canonical X_scANVI latent representation,
# UMAP, predicted labels, and Leiden clusterings consumed by Track A's L2
# annotation pipeline.
#
# Resources: Tier 5+GPU (128G / 8 CPU / 8h / 1 A100) — fine-tune is fast vs.
# scVI's full training (max_epochs=30 vs. 300).
#
# Requires N_LATENT to match a completed scVI run (default 50).
#
# Submit pattern:
#   source publication/config/load_paths.sh
#   sbatch publication/preprocessing/integration/run/run_scanvi_full.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../config/load_paths.sh"

INTEGRATION_ROOT="${PROJECT_PUBLICATION}/preprocessing/integration"
SCANVI_SCRIPT="${INTEGRATION_ROOT}/scripts/02_compartment_scanvi.py"

COMPONENTS="${SOURCE_IHBCAV1_ASSEMBLY_COMPONENTS}"
GENE_DATA="${COMPONENTS}gene_data.csv"

COMPARTMENTS=("Immune" "Epithelial" "Stromal")
SHORTS=("imm" "epi" "str")
COMP="${COMPARTMENTS[$SLURM_ARRAY_TASK_ID]}"
SHORT="${SHORTS[$SLURM_ARRAY_TASK_ID]}"

N_LATENT="${N_LATENT:-50}"
COUNTS_NPZ="${COMPONENTS}${SHORT}_counts.npz"
METADATA="${COMPONENTS}${SHORT}_metadata_enriched.csv"
SCVI_MODEL_DIR="${SOURCE_IHBCAV1_ASSEMBLY_ROOT}scvi/${SHORT}/n_latent_${N_LATENT}/scvi_model"
OUTPUT_DIR="${SOURCE_IHBCAV1_ASSEMBLY_SCANVI}${SHORT}/n_latent_${N_LATENT}"

mkdir -p "$OUTPUT_DIR"

echo "=== scANVI: $COMP n_latent=$N_LATENT (task $SLURM_ARRAY_TASK_ID) ==="
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "Counts: $COUNTS_NPZ"
echo "Metadata: $METADATA"
echo "scVI model: $SCVI_MODEL_DIR"
echo "Output: $OUTPUT_DIR"
echo ""

if [ ! -d "$SCVI_MODEL_DIR" ]; then
    echo "ERROR: scVI model dir not found: $SCVI_MODEL_DIR"
    echo "scANVI requires a completed scVI run at this path."
    exit 1
fi

for f in "$COUNTS_NPZ" "$METADATA" "$GENE_DATA"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: $f not found"
        exit 1
    fi
done

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
    python "$SCANVI_SCRIPT" \
        --counts-npz "$COUNTS_NPZ" \
        --gene-data "$GENE_DATA" \
        --metadata "$METADATA" \
        --compartment "$COMP" \
        --n-latent "$N_LATENT" \
        --scvi-model-dir "$SCVI_MODEL_DIR" \
        --output-dir "$OUTPUT_DIR" \
        --labels-key level1_annotation \
        --seed 42

echo ""
echo "=== Complete ==="
ls -lh "$OUTPUT_DIR/"
