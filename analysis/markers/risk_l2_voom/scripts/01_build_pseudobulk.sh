#!/bin/bash
#SBATCH --job-name=risk_l2_pseudobulk
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519/logs/01_pseudobulk_%j.out
#SBATCH --error=/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519/logs/01_pseudobulk_%j.err

# Step 1: aggregate (L2 x donor) pseudobulk counts from per-compartment NPZ.
# Reuses pseudobulk_aggregate_full_object.py unchanged — only the
# --leiden-aligned (here cell_to_l2.csv) and --groupby-col (level2) inputs
# change. Pattern source: this script is identical in shape to
# full_object_markers_20260421/run/*.sh that produced the leiden-resolution
# pseudobulks.
#
# Tier 4 sizing (64G, 4 CPUs, 2h) per ~/.claude/rules/hpc-resource-rules.md
# for full-cohort pseudobulk on 2.12M cells.
# Subset-testable: --n-cells N flag exists in the aggregator for smoke tests.

set -euo pipefail

PROJECT_ROOT="/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Abundance/DE_Testing/risk_l2_voom_20260519"
COMPONENTS="/share/crsp/lab/dalawson/nwechter/iHBCAv1_upload/.dev/assembly_rebuild_20260405/outputs/components"
AGGREGATOR_PY="/share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Annotation/dev/full_object_markers_20260421/scripts/pseudobulk_aggregate_full_object.py"

CELL_TO_L2="${PROJECT_ROOT}/outputs/pseudobulk/cell_to_l2.csv"
OUT_DIR="${PROJECT_ROOT}/outputs/pseudobulk"

# Optional smoke-test mode via SUBSET_N env var:
#   hpc submit ... --sbatch "--export=ALL,SUBSET_N=50000"
SUBSET_N="${SUBSET_N:-}"
SUBSET_FLAGS=""
TAG="risk_l2_voom"
if [[ -n "${SUBSET_N}" ]]; then
    SUBSET_FLAGS="--n-cells ${SUBSET_N}"
    TAG="risk_l2_voom_smoke${SUBSET_N}"
    echo "*** SMOKE TEST MODE: subsampling each compartment to ${SUBSET_N} cells ***"
fi

mkdir -p "${OUT_DIR}"

module load singularity

CONTAINER="/dfs7/singularity_containers/rcic/JHUB3/Rocky8_jupyter_base_R4.3.3_Spatial.sif"

echo "=== risk_l2_voom_20260519 / Step 1 pseudobulk ==="
echo "Host: $(hostname)"
echo "Start: $(date)"
echo "Tag: ${TAG}"

# Sanity-check inputs exist before kicking off the heavy job
for f in "${AGGREGATOR_PY}" "${CELL_TO_L2}" \
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
        --leiden-aligned   "${CELL_TO_L2}" \
        --groupby-col      L2_label \
        --donor-key        patientID \
        --study-key        dataset \
        --min-cells        5 \
        --output-dir       "${OUT_DIR}" \
        --tag              "${TAG}" \
        ${SUBSET_FLAGS}

echo "Done: $(date)"
echo
echo "Outputs:"
ls -la "${OUT_DIR}/${TAG}_pseudobulk_counts.csv" "${OUT_DIR}/${TAG}_pseudobulk_meta.csv" || true
