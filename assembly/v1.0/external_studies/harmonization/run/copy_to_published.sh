#!/bin/bash
#SBATCH --job-name=copy_pub
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1:00:00

# =============================================================================
# copy_to_published.sh
# =============================================================================
# Phase 2.2 → published/: Copy rebuilt study outputs into the standardized
# published output directories.
#
# Copies per study:
#   harmonization/outputs/study_objects/{study}/  →  outputs/{study}/published/
#     seurat.rds, metadata.csv, embedding_native.csv, embedding_joint.csv
#   harmonization/studies/{study}/outputs/{study}_cells.csv  →  cell_id_mapping.csv
#
# Then validates cell counts against expected values from studies.yaml.
#
# Why did the file copy go to published/?
# Because it finally got its act together.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../../config/load_paths.sh"

BASE="${SOURCE_IHBCAV1_EXTERNAL_STUDIES%/}"
SRC_DIR="${BASE}/harmonization/outputs/study_objects"
MAPPING_DIR="${BASE}/harmonization/studies"
DST_DIR="${BASE}/outputs"

STUDIES=(gray kumar murrow nee twigger pal_norm_epi pal_norm_total pal_norm_b1 reed)
EXPECTED_CELLS=(52681 714331 86136 230100 110744 53716 54332 59766 803283)

FILES=(seurat.rds metadata.csv embedding_native.csv embedding_joint.csv)

echo "=============================================="
echo "Copy to Published: Phase 2.2"
echo "=============================================="
echo "Start time: $(date)"
echo ""

FAIL=0

for i in "${!STUDIES[@]}"; do
  STUDY="${STUDIES[$i]}"
  EXPECTED="${EXPECTED_CELLS[$i]}"
  SRC="${SRC_DIR}/${STUDY}"
  PUB="${DST_DIR}/${STUDY}/published"

  echo "--- ${STUDY} ---"

  # Check source exists
  if [ ! -d "$SRC" ]; then
    echo "  ERROR: source dir missing: $SRC"
    FAIL=1
    continue
  fi

  # Create published dir
  mkdir -p "$PUB"

  # Copy study outputs
  for F in "${FILES[@]}"; do
    if [ -f "${SRC}/${F}" ]; then
      cp "${SRC}/${F}" "${PUB}/${F}"
      echo "  Copied: ${F} ($(du -h "${PUB}/${F}" | cut -f1))"
    else
      echo "  WARNING: ${F} not found in source"
      FAIL=1
    fi
  done

  # Copy cell ID mapping
  MAPPING_FILE="${MAPPING_DIR}/${STUDY}/outputs/${STUDY}_cells.csv"
  if [ -f "$MAPPING_FILE" ]; then
    cp "$MAPPING_FILE" "${PUB}/cell_id_mapping.csv"
    echo "  Copied: cell_id_mapping.csv ($(du -h "${PUB}/cell_id_mapping.csv" | cut -f1))"
  else
    echo "  WARNING: cell mapping not found: $MAPPING_FILE"
    FAIL=1
  fi

  # Copy prepare_summary.yaml if exists
  if [ -f "${SRC}/prepare_summary.yaml" ]; then
    cp "${SRC}/prepare_summary.yaml" "${PUB}/prepare_summary.yaml"
    echo "  Copied: prepare_summary.yaml"
  fi

  # Validate cell count from metadata.csv header
  if [ -f "${PUB}/metadata.csv" ]; then
    # Count data rows (total lines minus header)
    ACTUAL=$(( $(wc -l < "${PUB}/metadata.csv") - 1 ))
    if [ "$ACTUAL" -eq "$EXPECTED" ]; then
      echo "  Cell count: ${ACTUAL} ✓ (matches expected)"
    else
      echo "  ERROR: Cell count ${ACTUAL} != expected ${EXPECTED}"
      FAIL=1
    fi
  fi

  echo ""
done

# =============================================================================
# Harmonized metadata: copy to consumer-facing location
# =============================================================================
HARM_SRC="${BASE}/harmonization/outputs/harmonized_metadata"
HARM_DST="${DST_DIR}/harmonized_metadata"

echo "--- harmonized_metadata ---"
if [ -d "$HARM_SRC" ]; then
  mkdir -p "$HARM_DST"
  for F in harmonized_donor_metadata.csv harmonized_cell_metadata.csv coverage_matrix.csv no_contrast_notes.yaml; do
    if [ -f "${HARM_SRC}/${F}" ]; then
      cp "${HARM_SRC}/${F}" "${HARM_DST}/${F}"
      echo "  Copied: ${F} ($(du -h "${HARM_DST}/${F}" | cut -f1))"
    else
      echo "  Skipped: ${F} (not found)"
    fi
  done
  # Copy HTML report if present
  for F in "${HARM_SRC}"/*.html; do
    [ -f "$F" ] && cp "$F" "$HARM_DST/" && echo "  Copied: $(basename $F)"
  done
else
  echo "  ERROR: harmonized metadata source not found: $HARM_SRC"
  FAIL=1
fi
echo ""

echo "=============================================="
if [ $FAIL -eq 0 ]; then
  echo "SUCCESS: All studies + metadata copied and validated"
else
  echo "WARNING: Some checks failed — review above"
fi
echo "End time: $(date)"
echo "=============================================="

exit $FAIL
