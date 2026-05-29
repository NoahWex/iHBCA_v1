#!/bin/bash
#SBATCH --job-name=gsea_heatmap
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519/logs/16_gsea_heatmap_%j.out
#SBATCH --error=/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519/logs/16_gsea_heatmap_%j.err

# Render GSEA Hallmark heatmap for one track.
# Env vars: GSEA_CSV (input pathways_long.csv), OUT_PDF (output path),
# optional CELL_TYPES (comma-separated subset).

set -euo pipefail
PROJECT_ROOT="/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519"
PUBLICATION_ROOT="/share/crsp/lab/dalawson/nwechter/iHBCA_publication"

if [[ -z "${GSEA_CSV:-}" || -z "${OUT_PDF:-}" ]]; then
    echo "ERROR: GSEA_CSV and OUT_PDF env vars required" >&2
    exit 1
fi

module load singularity
CONTAINER="/dfs7/singularity_containers/rcic/JHUB3/Rocky8_jupyter_base_R4.3.3_Spatial.sif"

singularity exec \
    --bind /share/crsp/lab/dalawson:/share/crsp/lab/dalawson:rw \
    --bind /dfs7:/dfs7:ro \
    --bind /pub/nwechter:/pub/nwechter:ro \
    --env "R_LIBS_USER=/pub/nwechter/biojhub4_dir/Rocky8_jupyter_base_R4.3.3_Spatial.sif/R/library" \
    "${CONTAINER}" \
    Rscript "${PROJECT_ROOT}/scripts/16_gsea_heatmap.R" \
        --gsea-csv          "${GSEA_CSV}" \
        --out-pdf           "${OUT_PDF}" \
        --publication-root  "${PUBLICATION_ROOT}" \
        --cell-types        "${CELL_TYPES:-}" \
        --cell-types-file   "${CELL_TYPES_FILE:-}" \
        --pathways-file     "${PATHWAYS_FILE:-}"

ls -la "${OUT_PDF}"
