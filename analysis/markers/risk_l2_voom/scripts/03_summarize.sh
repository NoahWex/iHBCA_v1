#!/bin/bash
#SBATCH --job-name=risk_l2_summary
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519/logs/03_summary_%j.out
#SBATCH --error=/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519/logs/03_summary_%j.err

# Step 3: cross-formula comparison + recommended-formula report.
# Tier 1 sizing per ~/.claude/rules/hpc-resource-rules.md.

set -euo pipefail

PROJECT_ROOT="/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519"

module load singularity
CONTAINER="/dfs7/singularity_containers/rcic/JHUB3/Rocky8_jupyter_base_R4.3.3_Spatial.sif"
R_LIBS_USER_HOST="/pub/nwechter/biojhub4_dir/Rocky8_jupyter_base_R4.3.3_Spatial.sif/R/library"

echo "=== risk_l2_voom_20260519 / Step 3 summary ==="
echo "Start: $(date)"

singularity exec \
    --bind /share/crsp/lab/dalawson:/share/crsp/lab/dalawson:rw \
    --bind /dfs7:/dfs7:ro \
    --bind /dfs8:/dfs8:ro \
    --bind "${R_LIBS_USER_HOST}:/home/jovyan/R/library:ro" \
    --env "R_LIBS_USER=/home/jovyan/R/library" \
    "${CONTAINER}" \
    Rscript "${PROJECT_ROOT}/scripts/03_summarize.R" \
        --project-root "${PROJECT_ROOT}" \
        --inquiry-yaml "${PROJECT_ROOT}/config/inquiry.yaml" \
        --min-sig 20

echo "Done: $(date)"
echo
echo "Comparison outputs:"
ls -la "${PROJECT_ROOT}/outputs/comparisons/" 2>/dev/null
echo
echo "Report:"
ls -la "${PROJECT_ROOT}/reports/CP3_formula_choice.md" 2>/dev/null
