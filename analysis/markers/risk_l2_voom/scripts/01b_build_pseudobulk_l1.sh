#!/bin/bash
#SBATCH --job-name=risk_l1_pseudobulk
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519/logs/01b_pseudobulk_l1_%j.out
#SBATCH --error=/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519/logs/01b_pseudobulk_l1_%j.err

# L1-resolution pseudobulk. Mirror of 01_build_pseudobulk.sh — same aggregator
# (pseudobulk_aggregate_full_object.py), same NPZ inputs; swap cell_to_l2.csv
# for cell_to_l1.csv and --groupby-col L2_label for L1_label. Output tag
# risk_l1_voom so artifacts don't collide with the L2 build.
# Pattern source: 01_build_pseudobulk.sh (this stage).
# Tier 4 sizing matches the L2 build — same cell count, just 11 groups instead
# of 42, so memory profile is similar (matrix sizes per group are larger).

set -euo pipefail

PROJECT_ROOT="/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519"
COMPONENTS="/share/crsp/lab/dalawson/nwechter/iHBCAv1_upload/.dev/assembly_rebuild_20260405/outputs/components"
AGGREGATOR_PY="/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Annotation/dev/full_object_markers_20260421/scripts/pseudobulk_aggregate_full_object.py"

CELL_TO_L1="${PROJECT_ROOT}/outputs/pseudobulk/cell_to_l1.csv"
OUT_DIR="${PROJECT_ROOT}/outputs/pseudobulk"

TAG="risk_l1_voom"

mkdir -p "${OUT_DIR}"

module load singularity
CONTAINER="/dfs7/singularity_containers/rcic/JHUB3/Rocky8_jupyter_base_R4.3.3_Spatial.sif"

echo "=== risk_l1_voom / Step 1b pseudobulk (L1 resolution) ==="
echo "Host: $(hostname); Start: $(date)"

for f in "${AGGREGATOR_PY}" "${CELL_TO_L1}" \
         "${COMPONENTS}/epi_counts.npz" "${COMPONENTS}/imm_counts.npz" \
         "${COMPONENTS}/str_counts.npz" "${COMPONENTS}/gene_data.csv"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: missing input: ${f}" >&2
        exit 1
    fi
done

singularity exec \
    --bind /share/crsp/lab/dalawson:/share/crsp/lab/dalawson:rw \
    --bind /dfs7:/dfs7:ro \
    --bind /dfs8:/dfs8:ro \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
    "${CONTAINER}" \
    python3 "${AGGREGATOR_PY}" \
        --epi-npz   "${COMPONENTS}/epi_counts.npz" \
        --epi-meta  "${COMPONENTS}/epi_metadata_enriched.csv" \
        --imm-npz   "${COMPONENTS}/imm_counts.npz" \
        --imm-meta  "${COMPONENTS}/imm_metadata_enriched.csv" \
        --str-npz   "${COMPONENTS}/str_counts.npz" \
        --str-meta  "${COMPONENTS}/str_metadata_enriched.csv" \
        --gene-data-csv    "${COMPONENTS}/gene_data.csv" \
        --leiden-aligned   "${CELL_TO_L1}" \
        --groupby-col      L1_label \
        --donor-key        patientID \
        --study-key        dataset \
        --min-cells        5 \
        --output-dir       "${OUT_DIR}" \
        --tag              "${TAG}"

echo "Done: $(date)"
echo
echo "Outputs:"
ls -la "${OUT_DIR}/${TAG}_pseudobulk_counts.csv" "${OUT_DIR}/${TAG}_pseudobulk_meta.csv" || true
