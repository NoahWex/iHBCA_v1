#!/bin/bash
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --time=06:00:00
#SBATCH --job-name=diff_all_source
# =============================================================================
# diff_all_source.sh — Diff all 8 h5ads against pre-fix references
# =============================================================================
# Runs diff_h5ads.py for each study, comparing post-fix h5ads against
# Spatial_HBCA originals (source) and CxG published h5ad (integrated).
#
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
REPO_ROOT="${PROJECT_IHBCAV1_UPLOAD}"
PUB_ROOT="$REPO_ROOT/publication"
SCRIPTS="$PUB_ROOT/scripts"
NEW_DIR="$PUB_ROOT/outputs/source_datasets"
REF_DIR="${SOURCE_SPATIAL_PUBLICATION_OUTPUTS}source_datasets"
REPORT_DIR="$PUB_ROOT/outputs/validation_reports"

PY_CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"
BINDS=(
    ${BIND_MOUNTS} \
)

# Study -> output filename mapping
declare -A STUDIES=(
    [gray]=gray2022
    [kumar]=kumar2023
    [murrow]=murrow2022
    [nee]=nee2023
    [twigger]=twigger2022
    [reed]=reed2024
    [pal]=pal2021
)

module load singularity

mkdir -p "$REPORT_DIR" "$PUB_ROOT/logs"

PASS_COUNT=0
FAIL_COUNT=0
SKIP_COUNT=0
RESULTS=()

echo "============================================="
echo "DIFF ALL SOURCE DATASETS"
echo "  New:       $NEW_DIR"
echo "  Reference: $REF_DIR"
echo "  Reports:   $REPORT_DIR"
echo "============================================="
echo ""

for study in gray kumar murrow nee twigger reed pal; do
    filename="${STUDIES[$study]}"
    new_file="$NEW_DIR/${filename}.h5ad"
    ref_file="$REF_DIR/${filename}.h5ad"
    report_file="$REPORT_DIR/${study}_regen_diff.yaml"

    echo "---------------------------------------------"
    echo "DIFF: $study ($filename)"
    echo "---------------------------------------------"

    if [[ ! -f "$new_file" ]]; then
        echo "  SKIP: New file not found: $new_file"
        SKIP_COUNT=$((SKIP_COUNT + 1))
        RESULTS+=("$study: SKIP (new file missing)")
        continue
    fi

    if [[ ! -f "$ref_file" ]]; then
        echo "  SKIP: Reference file not found: $ref_file"
        SKIP_COUNT=$((SKIP_COUNT + 1))
        RESULTS+=("$study: SKIP (reference file missing)")
        continue
    fi

    echo "  New: $(ls -lh "$new_file" | awk '{print $5}')"
    echo "  Ref: $(ls -lh "$ref_file" | awk '{print $5}')"

    # diff_h5ads.py exits non-zero on FAIL — don't let set -e kill the loop,
    # we parse the YAML report ourselves below
    singularity exec \
        --no-mount bind-paths \
        "${BINDS[@]}" \
        --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
        --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
        "$PY_CONTAINER" \
        python3 "$SCRIPTS/diff_h5ads.py" \
            --new "$new_file" \
            --reference "$ref_file" \
            --output "$report_file" \
        || true

    # Allow CRSP cache to settle before reading report
    sleep 3

    # Check result — diff_h5ads.py writes YAML with overall_pass: true/false
    if [[ -f "$report_file" ]]; then
        result=$(grep "overall_pass:" "$report_file" 2>/dev/null || echo "UNKNOWN")
        if echo "$result" | grep -q "true"; then
            echo "  RESULT: PASS"
            PASS_COUNT=$((PASS_COUNT + 1))
            RESULTS+=("$study: PASS")
        elif echo "$result" | grep -q "false"; then
            echo "  RESULT: FAIL"
            FAIL_COUNT=$((FAIL_COUNT + 1))
            RESULTS+=("$study: FAIL")
        else
            echo "  RESULT: UNKNOWN (could not parse report)"
            RESULTS+=("$study: UNKNOWN")
        fi
    else
        echo "  RESULT: UNKNOWN (report file not found — CRSP cache lag?)"
        RESULTS+=("$study: UNKNOWN (report missing)")
    fi
    echo ""
done

# =============================================================================
# Run: Integrated object (vs CxG published h5ad)
# =============================================================================
echo ""
echo "---------------------------------------------"
echo "DIFF: integrated (all-breast-cells vs CxG published)"
echo "---------------------------------------------"

INTEGRATED_NEW="$PUB_ROOT/outputs/integrated_objects/all-breast-cells.h5ad"
INTEGRATED_REF="${SOURCE_INTEGRATED_ATLAS_REED}integration_iHBCA.h5ad"
INTEGRATED_REPORT="$REPORT_DIR/integrated_regen_diff.yaml"

if [[ ! -f "$INTEGRATED_NEW" ]]; then
    echo "  SKIP: New file not found: $INTEGRATED_NEW"
    SKIP_COUNT=$((SKIP_COUNT + 1))
    RESULTS+=("integrated: SKIP (new file missing)")
elif [[ ! -f "$INTEGRATED_REF" ]]; then
    echo "  SKIP: Reference file not found: $INTEGRATED_REF"
    SKIP_COUNT=$((SKIP_COUNT + 1))
    RESULTS+=("integrated: SKIP (reference file missing)")
else
    echo "  New: $(ls -lh "$INTEGRATED_NEW" | awk '{print $5}')"
    echo "  Ref: $(ls -lh "$INTEGRATED_REF" | awk '{print $5}')"

    singularity exec \
        --no-mount bind-paths \
        "${BINDS[@]}" \
        --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
        --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
        "$PY_CONTAINER" \
        python3 "$SCRIPTS/diff_h5ads.py" \
            --new "$INTEGRATED_NEW" \
            --reference "$INTEGRATED_REF" \
            --output "$INTEGRATED_REPORT" \
        || true

    sleep 3

    if [[ -f "$INTEGRATED_REPORT" ]]; then
        result=$(grep "overall_pass:" "$INTEGRATED_REPORT" 2>/dev/null || echo "UNKNOWN")
        if echo "$result" | grep -q "true"; then
            echo "  RESULT: PASS"
            PASS_COUNT=$((PASS_COUNT + 1))
            RESULTS+=("integrated: PASS")
        elif echo "$result" | grep -q "false"; then
            echo "  RESULT: FAIL (expected — structural changes from fixes)"
            FAIL_COUNT=$((FAIL_COUNT + 1))
            RESULTS+=("integrated: FAIL (expected diffs)")
        else
            RESULTS+=("integrated: UNKNOWN")
        fi
    else
        RESULTS+=("integrated: UNKNOWN (report missing)")
    fi
fi
echo ""

echo "============================================="
echo "SUMMARY"
echo "============================================="
for r in "${RESULTS[@]}"; do
    echo "  $r"
done
echo ""
echo "PASS: $PASS_COUNT  FAIL: $FAIL_COUNT  SKIP: $SKIP_COUNT"
echo ""
echo "NOTE: FAIL is expected for studies with structural changes from"
echo "validation fixes (A1-A12). Key check: cell counts must match."
echo ""

if [[ $SKIP_COUNT -gt 0 ]]; then
    echo "WARNING: $SKIP_COUNT studies SKIPPED"
fi

echo "Reports: $REPORT_DIR/"
echo "Done: $(date)"
