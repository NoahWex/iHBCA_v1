#!/bin/bash
#SBATCH --job-name=aucell_l1_fibro
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519/logs/13_aucell_l1_%j.out
#SBATCH --error=/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519/logs/13_aucell_l1_%j.err

set -euo pipefail

PROJECT_ROOT="/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519"
PUBLICATION_ROOT="/share/crsp/lab/dalawson/nwechter/iHBCA_publication"

SCORE_DIR="${SCORE_DIR:-${PROJECT_ROOT}/outputs/plots/aucell_preneoplastic}"
OUT_DIR="${OUT_DIR:-${PROJECT_ROOT}/outputs/plots/aucell_preneoplastic}"

module load singularity
CONTAINER="/dfs7/singularity_containers/rcic/JHUB3/Rocky8_jupyter_base_R4.3.3_Spatial.sif"

singularity exec \
    --bind /share/crsp/lab/dalawson:/share/crsp/lab/dalawson:rw \
    --bind /dfs7:/dfs7:ro \
    --bind /pub/nwechter:/pub/nwechter:ro \
    --env "R_LIBS_USER=/pub/nwechter/biojhub4_dir/Rocky8_jupyter_base_R4.3.3_Spatial.sif/R/library" \
    "${CONTAINER}" \
    Rscript "${PROJECT_ROOT}/scripts/13_aucell_l1_fibroblast.R" \
        --score-dir         "${SCORE_DIR}" \
        --publication-root  "${PUBLICATION_ROOT}" \
        --out-dir           "${OUT_DIR}" \
        --panels-file       "${PANELS_FILE:-}" \
        --out-stem          "${OUT_STEM:-fibro_l1_beeswarm}"

ls -la "${OUT_DIR}"/fibro_l1_beeswarm.pdf "${OUT_DIR}"/l1_*.csv
