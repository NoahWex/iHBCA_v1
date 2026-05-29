#!/bin/bash
#SBATCH --job-name=validate_sketch
#SBATCH --account=dalawson_lab
#SBATCH --partition=free
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=02:00:00

# =============================================================================
# Validate sketch object with CxG, CAP, and HCA validators
# =============================================================================

set -uo pipefail
# Note: not using -e because validators return non-zero on validation failures

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
REPO_ROOT="${PROJECT_IHBCAV1_UPLOAD}"
SKETCH="${REPO_ROOT}/publication/outputs/integrated_objects/all-breast-cells-sketch.h5ad"
REPORTS="${REPO_ROOT}/publication/outputs/validation_reports"
HCA_WRAPPER="${REPO_ROOT}/publication/scripts/run_hca_validator.py"

echo "============================================="
echo "Validate sketch object"
echo "  Date: $(date)"
echo "  Node: $(hostname)"
echo "  Job ID: ${SLURM_JOB_ID}"
echo "  Input: ${SKETCH}"
echo "============================================="

if [ ! -f "${SKETCH}" ]; then
    echo "ERROR: Sketch h5ad not found: ${SKETCH}"
    exit 1
fi

echo "Sketch size: $(ls -lh "${SKETCH}" | awk '{print $5}')"
mkdir -p "${REPORTS}"

module load mamba/24.3.0
source activate hca_validators

PASS_COUNT=0
FAIL_COUNT=0

# --- CxG Schema Validator ---
echo ""
echo "=== CxG Schema Validator ==="
CXG_LOG="${REPORTS}/all-breast-cells-sketch_cxg.log"
cellxgene-schema validate "${SKETCH}" 2>&1 | tee "${CXG_LOG}" || true
if grep -q "is_valid=True" "${CXG_LOG}"; then
    echo "RESULT: PASS"
    PASS_COUNT=$((PASS_COUNT + 1))
elif grep -q "is_valid=False" "${CXG_LOG}"; then
    echo "RESULT: FAIL"
    FAIL_COUNT=$((FAIL_COUNT + 1))
else
    echo "RESULT: ERR (could not determine)"
    FAIL_COUNT=$((FAIL_COUNT + 1))
fi

# --- CAP Validator ---
echo ""
echo "=== CAP Validator ==="
CAP_LOG="${REPORTS}/all-breast-cells-sketch_cap.log"
python3 -m cap_upload_validator "${SKETCH}" 2>&1 | tee "${CAP_LOG}" || true
if grep -q "Validation passed" "${CAP_LOG}"; then
    echo "RESULT: PASS"
    PASS_COUNT=$((PASS_COUNT + 1))
else
    echo "RESULT: FAIL"
    FAIL_COUNT=$((FAIL_COUNT + 1))
fi

# --- HCA Validator ---
echo ""
echo "=== HCA Validator ==="
HCA_LOG="${REPORTS}/all-breast-cells-sketch_hca.log"
python3 "${HCA_WRAPPER}" "${SKETCH}" 2>&1 | tee "${HCA_LOG}" || true
if grep -q "is_valid: True" "${HCA_LOG}"; then
    echo "RESULT: PASS"
    PASS_COUNT=$((PASS_COUNT + 1))
elif grep -q "is_valid: False" "${HCA_LOG}"; then
    echo "RESULT: FAIL"
    FAIL_COUNT=$((FAIL_COUNT + 1))
else
    echo "RESULT: ERR (could not determine)"
    FAIL_COUNT=$((FAIL_COUNT + 1))
fi

# --- Summary ---
echo ""
echo "============================================="
echo "Validation Summary"
echo "  PASS: ${PASS_COUNT}/3"
echo "  FAIL: ${FAIL_COUNT}/3"
echo "  Date: $(date)"
echo "============================================="
