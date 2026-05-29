#!/bin/bash
#SBATCH --job-name=geno_tested_de
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519/logs/07_geno_tested_%j.out
#SBATCH --error=/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519/logs/07_geno_tested_%j.err

set -euo pipefail
PROJECT_ROOT="/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519"
PB_DIR="${PROJECT_ROOT}/outputs/pseudobulk"
GENE_MAPPING="/share/crsp/lab/dalawson/nwechter/iHBCAv1_upload/iHBCA_upload/publication/mappings/gene_symbol_to_ensembl.tsv"
# Use derived_donor_metadata.csv (287-donor full atlas with risk_class +
# brca_genotype string-coded). NB: parity_binary in this file is numeric 0/1;
# the cohort_BR1_vs_AR_design.csv has it string-coded. For this restricted
# query we need brca_genotype which IS string-coded in both — but parity_binary
# differs. Use the cohort design file (already filtered to risk_class AR+BR1)
# and join the brca_genotype back from the derived file.
DONOR_META_PRIMARY="/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/da_pipeline/inquiries/risk_main_effects_20260507/outputs/stageA_cohorts/cohort_BR1_vs_AR_design.csv"

module load singularity
CONTAINER="/dfs7/singularity_containers/rcic/JHUB3/Rocky8_jupyter_base_R4.3.3_Spatial.sif"
R_LIBS_USER_HOST="/pub/nwechter/biojhub4_dir/Rocky8_jupyter_base_R4.3.3_Spatial.sif/R/library"

singularity exec \
    --bind /share/crsp/lab/dalawson:/share/crsp/lab/dalawson:rw \
    --bind /dfs7:/dfs7:ro \
    --bind "${R_LIBS_USER_HOST}:/home/jovyan/R/library:ro" \
    --env "R_LIBS_USER=/home/jovyan/R/library" \
    "${CONTAINER}" \
    Rscript "${PROJECT_ROOT}/scripts/07_genotyped_cohort_fibros.R" \
        --project-root "${PROJECT_ROOT}" \
        --counts       "${PB_DIR}/risk_l2_voom_pseudobulk_counts.csv" \
        --pb-meta      "${PB_DIR}/risk_l2_voom_pseudobulk_meta.csv" \
        --donor-meta   "${DONOR_META_PRIMARY}" \
        --gene-mapping "${GENE_MAPPING}"
