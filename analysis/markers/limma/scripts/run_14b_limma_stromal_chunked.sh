#!/bin/bash
#SBATCH --job-name=14b_limma_str_c
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --array=0-41
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --output=slurm-%x_%j_%a.out
#SBATCH --error=slurm-%x_%j_%a.err

# =============================================================================
# Pseudobulk limma-voom one-vs-rest markers for iHBCA v1 stromal compartment
# across 7 Leiden resolutions, chunked for SLURM array parallelism.
# See run_14b_limma_epithelial_chunked.sh header for full pipeline description.
#
# Canonical resolution for the T2 supp table: leiden_1.0 (26 clusters → 12 L2).
# =============================================================================

set -euo pipefail

if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$PWD"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"

RESOLUTIONS=(0.1 0.2 0.3 0.5 0.8 1.0 1.5)
MAX_CHUNKS=6
CHUNK_SIZE=8

RES_IDX=$((SLURM_ARRAY_TASK_ID / MAX_CHUNKS))
CHUNK_IDX=$((SLURM_ARRAY_TASK_ID % MAX_CHUNKS))
RES="${RESOLUTIONS[$RES_IDX]}"
GROUPBY_COL="leiden_${RES}"
COMPARTMENT="stromal"
TAG="${COMPARTMENT}_${GROUPBY_COL}"

INDIR="${SOURCE_IHBCA_V1_PSEUDOBULK_ROOT}/${COMPARTMENT}/lineage_markers/deseq2/${GROUPBY_COL}"
OUTDIR="${PROJECT_PUBLICATION}/publication/analysis/markers/limma/outputs/${COMPARTMENT}/${GROUPBY_COL}"

module load singularity

R_SING_EXEC="singularity exec --no-home \
  ${BIND_MOUNTS} \
  --bind ${R_LIBS_R_SPATIAL_4_3_3}:/home/jovyan/R/userlib:ro \
  --env R_LIBS_USER=/home/jovyan/R/userlib \
  ${CONTAINER_R_SPATIAL_4_3_3}"

echo "=== 14b limma-voom ${COMPARTMENT} chunked: task=${SLURM_ARRAY_TASK_ID} res=${RES} chunk=${CHUNK_IDX}/${MAX_CHUNKS} groupby=${GROUPBY_COL} ==="
echo "Start: $(date)"

mkdir -p "$OUTDIR"

$R_SING_EXEC Rscript "${SCRIPT_DIR}/14b_marker_limma.R" \
  --counts  "$INDIR/${TAG}_pseudobulk_counts.csv" \
  --meta    "$INDIR/${TAG}_pseudobulk_meta.csv" \
  --output-dir "$OUTDIR" \
  --level-tag  "$GROUPBY_COL" \
  --n-top 20 \
  --chunk-idx  "$CHUNK_IDX" \
  --chunk-size "$CHUNK_SIZE"

echo "Done: $(date)"
