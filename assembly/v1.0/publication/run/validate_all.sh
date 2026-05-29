#!/bin/bash
#SBATCH --job-name=validate_all
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=180G
#SBATCH --time=06:00:00

# =============================================================================
# Run all 3 validators (CxG, CAP, HCA) on all 8 h5ads (7 source + 1 integrated)
# =============================================================================
# Logs per file per validator to publication/outputs/validation_reports/
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
REPO_ROOT="${PROJECT_IHBCAV1_UPLOAD}"
source "${REPO_ROOT}/publication/run/load_config.sh"
PUB_ROOT="${REPO_ROOT}/publication"

# Env-var overrides for alternate builds (e.g., SOURCE_DIR=.../cxg_build/source_datasets)
SOURCE_DIR="${SOURCE_DIR:-${OUTPUT_SOURCE_DATASETS}}"
INTEGRATED_DIR="${INTEGRATED_DIR:-${OUTPUT_INTEGRATED_OBJECTS}}"
REPORT_DIR="${REPORT_DIR:-${OUTPUT_VALIDATION_REPORTS}}"

# --- File lists ---

SOURCE_FILES=(
    "gray2022.h5ad"
    "kumar2023.h5ad"
    "murrow2022.h5ad"
    "nee2023.h5ad"
    "twigger2022.h5ad"
    "reed2024.h5ad"
    "pal2021.h5ad"
)

INTEGRATED_FILES=(
    "all-breast-cells.h5ad"
)

echo "=== Full Validation Run: 8 h5ads x 3 validators ==="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Memory: ${SLURM_MEM_PER_NODE}MB"
echo ""

mkdir -p "${REPORT_DIR}"

# --- Activate validator environment ---

module load mamba/24.3.0
source activate hca_validators

echo "Python: $(which python3)"
echo "Python version: $(python3 --version)"
echo "anndata version: $(python3 -c 'import anndata; print(anndata.__version__)' 2>/dev/null || echo 'not found')"
echo "cellxgene-schema version: $(python3 -c 'import cellxgene_schema; print(cellxgene_schema.__version__)' 2>/dev/null || echo 'not found')"
echo ""

# --- Determine CAP command (module invocation preferred over CLI entry point) ---

if python3 -c "import cap_upload_validator" 2>/dev/null; then
    CAP_CMD="python3 -m cap_upload_validator"
elif command -v capval &> /dev/null; then
    CAP_CMD="capval"
else
    echo "WARNING: No CAP validator found. CAP checks will be skipped."
    CAP_CMD=""
fi
echo "CAP command: ${CAP_CMD:-NONE}"
echo ""

# --- Tracking arrays ---

declare -a RESULT_FILE=()
declare -a RESULT_CXG=()
declare -a RESULT_CAP=()
declare -a RESULT_HCA=()

# =============================================================================
# Validate one file with all 3 validators
# Args: $1 = full path, $2 = basename (without .h5ad)
# =============================================================================
validate_file() {
    local INPUT_PATH="$1"
    local BASENAME="$2"

    echo "========================================================================"
    echo "Validating: ${BASENAME}.h5ad"
    echo "  Path: ${INPUT_PATH}"
    echo "========================================================================"

    if [ ! -f "${INPUT_PATH}" ]; then
        echo "  SKIPPED: File not found"
        RESULT_FILE+=("${BASENAME}")
        RESULT_CXG+=("SKIP")
        RESULT_CAP+=("SKIP")
        RESULT_HCA+=("SKIP")
        return
    fi

    echo "  Size: $(ls -lh "${INPUT_PATH}" | awk '{print $5}')"
    echo "  Started: $(date)"

    # --- CxG validator ---
    echo ""
    echo "  --- CxG Validator ---"
    cellxgene-schema validate "${INPUT_PATH}" \
        2>&1 | tee "${REPORT_DIR}/${BASENAME}_cxg.log" || true
    CXG_EXIT=${PIPESTATUS[0]}
    echo "  CxG exit code: ${CXG_EXIT}"

    # --- CAP validator ---
    echo ""
    echo "  --- CAP Validator ---"
    if [ -n "${CAP_CMD}" ]; then
        ${CAP_CMD} "${INPUT_PATH}" \
            2>&1 | tee "${REPORT_DIR}/${BASENAME}_cap.log" || true
        CAP_EXIT=${PIPESTATUS[0]}
    else
        echo "  SKIPPED: CAP validator not available"
        CAP_EXIT=-1
    fi
    echo "  CAP exit code: ${CAP_EXIT}"

    # --- HCA validator (via wrapper — hca-schema-validator lacks __main__.py) ---
    echo ""
    echo "  --- HCA Validator ---"
    python3 "${PUB_ROOT}/scripts/run_hca_validator.py" "${INPUT_PATH}" \
        2>&1 | tee "${REPORT_DIR}/${BASENAME}_hca.log" || true
    HCA_EXIT=${PIPESTATUS[0]}
    echo "  HCA exit code: ${HCA_EXIT}"

    echo "  Finished: $(date)"
    echo ""

    # --- Record results (parse log content, not exit codes) ---
    RESULT_FILE+=("${BASENAME}")

    # CxG: parse is_valid from log
    CXG_LOG="${REPORT_DIR}/${BASENAME}_cxg.log"
    if grep -q "is_valid=False" "$CXG_LOG" 2>/dev/null; then
        RESULT_CXG+=("FAIL")
    elif grep -q "is_valid=True" "$CXG_LOG" 2>/dev/null; then
        RESULT_CXG+=("PASS")
    else
        RESULT_CXG+=("ERR")
    fi

    # CAP: "Validation passed!" = pass, anything else = fail
    CAP_LOG="${REPORT_DIR}/${BASENAME}_cap.log"
    if [ -z "${CAP_CMD}" ]; then
        RESULT_CAP+=("SKIP")
    elif grep -q "Validation passed" "$CAP_LOG" 2>/dev/null; then
        RESULT_CAP+=("PASS")
    else
        RESULT_CAP+=("FAIL")
    fi

    # HCA: "is_valid: True" vs "is_valid: False"
    HCA_LOG="${REPORT_DIR}/${BASENAME}_hca.log"
    if grep -q "is_valid: False" "$HCA_LOG" 2>/dev/null; then
        RESULT_HCA+=("FAIL")
    elif grep -q "is_valid: True" "$HCA_LOG" 2>/dev/null; then
        RESULT_HCA+=("PASS")
    else
        RESULT_HCA+=("ERR")
    fi
}

# =============================================================================
# Run: Source datasets
# =============================================================================
echo ""
echo "================================================================"
echo "  PHASE 1: Source Datasets (${#SOURCE_FILES[@]} files)"
echo "================================================================"
echo ""

for F in "${SOURCE_FILES[@]}"; do
    BASENAME="${F%.h5ad}"
    validate_file "${SOURCE_DIR}/${F}" "${BASENAME}"
done

# =============================================================================
# Run: Integrated objects
# =============================================================================
echo ""
echo "================================================================"
echo "  PHASE 2: Integrated Objects (${#INTEGRATED_FILES[@]} files)"
echo "================================================================"
echo ""

for F in "${INTEGRATED_FILES[@]}"; do
    BASENAME="${F%.h5ad}"
    validate_file "${INTEGRATED_DIR}/${F}" "${BASENAME}"
done

# =============================================================================
# Summary table
# =============================================================================
echo ""
echo "========================================================================"
echo "VALIDATION SUMMARY"
echo "========================================================================"
echo ""

TOTAL_FILES=${#RESULT_FILE[@]}
PASS_COUNT=0
FAIL_COUNT=0

printf "%-35s  %-6s  %-6s  %-6s\n" "FILE" "CxG" "CAP" "HCA"
printf "%-35s  %-6s  %-6s  %-6s\n" "-----------------------------------" "------" "------" "------"

for i in $(seq 0 $((TOTAL_FILES - 1))); do
    printf "%-35s  %-6s  %-6s  %-6s\n" \
        "${RESULT_FILE[$i]}" "${RESULT_CXG[$i]}" "${RESULT_CAP[$i]}" "${RESULT_HCA[$i]}"

    if [[ "${RESULT_CXG[$i]}" == "PASS" && "${RESULT_CAP[$i]}" != "FAIL" && "${RESULT_HCA[$i]}" == "PASS" ]]; then
        PASS_COUNT=$((PASS_COUNT + 1))
    else
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
done

echo ""
echo "Files validated: ${TOTAL_FILES}"
echo "All-pass: ${PASS_COUNT}"
echo "Has failures: ${FAIL_COUNT}"
echo "Reports: ${REPORT_DIR}/"
echo ""
ls -la "${REPORT_DIR}/"*.log 2>/dev/null || echo "No log files found"
echo ""
echo "Done: $(date)"
