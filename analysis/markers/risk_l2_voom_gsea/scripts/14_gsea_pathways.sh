#!/bin/bash
#SBATCH --job-name=gsea_track
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519/logs/14_gsea_%j.out
#SBATCH --error=/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519/logs/14_gsea_%j.err

# fgsea pre-ranked GSEA on one cohort × resolution track.
# Parameterize via env vars (SLURM doesn't pass positional args cleanly):
#   DE_DIR  — voom DE output directory for this track (one CSV per cell type)
#   OUT_DIR — pathway output directory
# Submit via: hpc submit ... --sbatch "--export=ALL,DE_DIR=<>,OUT_DIR=<>"
#
# Tier 2 sizing (16G, 4 CPUs, 1h): fgsea is fast (~1-2 min per cell type),
# but msigdbr load is ~30s and we iterate over up to 42 L2 cell types.

set -euo pipefail
PROJECT_ROOT="/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519"

if [[ -z "${DE_DIR:-}" || -z "${OUT_DIR:-}" ]]; then
    echo "ERROR: DE_DIR and OUT_DIR env vars required" >&2
    echo "Example: --sbatch '--export=ALL,DE_DIR=outputs/de_results/A_full/A2_full,OUT_DIR=outputs/gsea/A_full_L1'" >&2
    exit 1
fi

mkdir -p "${OUT_DIR}"

module load singularity
CONTAINER="/dfs7/singularity_containers/rcic/JHUB3/Rocky8_jupyter_base_R4.3.3_Spatial.sif"

echo "=== GSEA track ==="
echo "DE_DIR:  ${DE_DIR}"
echo "OUT_DIR: ${OUT_DIR}"
echo "Start:   $(date)"

singularity exec \
    --bind /share/crsp/lab/dalawson:/share/crsp/lab/dalawson:rw \
    --bind /dfs7:/dfs7:ro \
    --bind /pub/nwechter:/pub/nwechter:ro \
    --env "R_LIBS_USER=/pub/nwechter/biojhub4_dir/Rocky8_jupyter_base_R4.3.3_Spatial.sif/R/library" \
    "${CONTAINER}" \
    Rscript "${PROJECT_ROOT}/scripts/14_gsea_pathways.R" \
        --de-dir      "${DE_DIR}" \
        --out-dir     "${OUT_DIR}" \
        --collections "H"

echo "Done: $(date)"
ls -la "${OUT_DIR}/pathways_long.csv" "${OUT_DIR}/audit.csv" 2>/dev/null || true
