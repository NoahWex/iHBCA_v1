#!/bin/bash
#SBATCH --job-name=aucell_full
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=96G
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Full-cohort AUCell scoring (sensitivity mirror of the tested-cohort canonical
# run). Same engine and signature panels as build_ihbca_fig2_aucell_scores.sh,
# with --cohort full (n = 159 donors: 125 AR + 34 BR1). Memory is bumped to 96G
# because the full cohort has ~3x more stromal cells than the tested cohort.

set -euo pipefail

# Hardcoded HPC root — sbatch doesn't propagate caller env. Path locked by paths.yaml.
PUBLICATION_ROOT="${PUBLICATION_ROOT:-${PROJECT_PUBLICATION:-/share/crsp/lab/dalawson/nwechter/iHBCA_publication}}"
if [ ! -f "${PUBLICATION_ROOT}/publication/config/load_paths.sh" ]; then
  echo "ERROR: PUBLICATION_ROOT=${PUBLICATION_ROOT} missing publication/config/load_paths.sh." >&2
  exit 1
fi
echo "[init] PUBLICATION_ROOT=${PUBLICATION_ROOT}"
# shellcheck disable=SC1091
source "${PUBLICATION_ROOT}/publication/config/load_paths.sh"

COMPONENTS="${SOURCE_IHBCAV1_ASSEMBLY_COMPONENTS}"
LABELS="${PUBLICATION_ROOT}/publication/analysis/annotation/labels_full.csv"
DONOR_META="${SOURCE_IHBCA_V1_ABUNDANCE}/da_pipeline/inquiries/risk_main_effects_20260507/config/derived_donor_metadata.csv"
GENE_MAPPING="${USER_ROOT}/iHBCAv1_upload/iHBCA_upload/publication/mappings/gene_symbol_to_ensembl.tsv"

OUT_DIR="${OUT_DIR:-${PUBLICATION_ROOT}/publication/analysis/markers/aucell_methods_qc/outputs/plots/aucell_preneoplastic_full}"
WORKDIR="${OUT_DIR}/_workdir"
mkdir -p "${OUT_DIR}" "${WORKDIR}"

module load singularity
CONTAINER="${CONTAINER_R_SPATIAL_4_3_3}"
R_LIBS_BIND="${R_LIBS_R_SPATIAL_4_3_3}"

PY_SCRIPT="${PUBLICATION_ROOT}/publication/analysis/markers/aucell_methods_qc/scripts/build_ihbca_fig2_aucell_scores.py"
R_SCRIPT="${PUBLICATION_ROOT}/publication/analysis/markers/aucell_methods_qc/scripts/build_ihbca_fig2_aucell_scores.R"

singularity exec \
    ${BIND_MOUNTS} \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
    --env "R_LIBS_USER=${R_LIBS_BIND}" \
    "${CONTAINER}" \
    python3 "${PY_SCRIPT}" \
        --str-npz       "${COMPONENTS}/str_counts.npz" \
        --str-meta      "${COMPONENTS}/str_metadata_enriched.csv" \
        --labels-full   "${LABELS}" \
        --donor-meta    "${DONOR_META}" \
        --gene-data     "${COMPONENTS}/gene_data.csv" \
        --gene-mapping  "${GENE_MAPPING}" \
        --output-dir    "${OUT_DIR}" \
        --r-script      "${R_SCRIPT}" \
        --workdir       "${WORKDIR}" \
        --cohort        full

ls -la "${OUT_DIR}"
