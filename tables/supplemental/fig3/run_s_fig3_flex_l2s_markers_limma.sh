#!/bin/bash
#SBATCH --job-name=flex_l2s_limma_ann
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/share/crsp/lab/dalawson/nwechter/iHBCA_publication/coordination/handoff/kai_20260523_l1_markers/logs/%x_%j.out
#SBATCH --error=/share/crsp/lab/dalawson/nwechter/iHBCA_publication/coordination/handoff/kai_20260523_l1_markers/logs/%x_%j.err

# limma-voom on FLEX L2S using ANNOTATED vocabulary (flex_l2s_labels.csv).
# Excludes ARTIFACT_* groups. Prior run used raw cluster IDs.

set -euo pipefail

SCRIPT_DIR="/share/crsp/lab/dalawson/nwechter/iHBCA_publication/coordination/handoff/kai_20260523_l1_markers"
source "/share/crsp/lab/dalawson/nwechter/iHBCA_publication/publication/config/load_paths.sh"

SCRIPT="${SCRIPT_DIR}/run_flex_l2s_limma.R"
SEURAT="${PROJECT_SPATIAL_HBCA}/project/01_Preprocessing/outputs/18_IntegrationPreview/seurat_objects/integrated_seurat.rds"
L2S_CSV="${PROJECT_PUBLICATION}/publication/analysis/annotation/flex/outputs/flex_l2s_labels.csv"
OUT_DIR="${SCRIPT_DIR}/outputs"

mkdir -p "${OUT_DIR}"

for f in "$SCRIPT" "$SEURAT" "$L2S_CSV"; do
    [ -f "$f" ] || { echo "ERROR: missing input: $f"; exit 1; }
done

module load singularity/3.11.3

singularity exec \
    --cleanenv --containall --no-home \
    --bind "${R_LIBS_USER_DEFAULT}:/home/jovyan/R/library:ro" \
    ${BIND_MOUNTS} \
    --env "HOME=/home/jovyan" \
    --env "R_LIBS_USER=/home/jovyan/R/library" \
    --env "LANG=en_US.UTF-8" \
    --env "LC_ALL=en_US.UTF-8" \
    "${CONTAINER_R_SPATIAL_4_3_3}" \
    Rscript "$SCRIPT" \
        --seurat-path "$SEURAT" \
        --l2s-csv "$L2S_CSV" \
        --donor-col "orig.ident" \
        --min-cells-per-pb 10 \
        --min-donors-per-L2S 3 \
        --outdir "$OUT_DIR"

echo "=== Complete ==="
ls -la "${OUT_DIR}/flex_l2s_limma_de_annotated.csv"
