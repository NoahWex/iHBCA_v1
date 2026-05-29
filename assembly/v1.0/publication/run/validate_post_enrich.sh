#!/bin/bash
#SBATCH --job-name=validate_post_enrich
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --time=06:00:00

# Post-enrichment validation: CAP validator on all 11 h5ads
# Uses validate_all.sh logic but with 256GB for enlarged integrated objects

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
REPO_ROOT="${PROJECT_IHBCAV1_UPLOAD}"
PUB_ROOT="${REPO_ROOT}/publication"
SOURCE_DIR="${PUB_ROOT}/outputs/source_datasets"
INTEGRATED_DIR="${PUB_ROOT}/outputs/integrated_objects"
REPORT_DIR="${PUB_ROOT}/outputs/validation_reports"

SOURCE_FILES=(gray2022.h5ad kumar2023.h5ad murrow2022.h5ad nee2023.h5ad twigger2022.h5ad reed2024.h5ad pal2021.h5ad)
INTEGRATED_FILES=(all-breast-cells.h5ad breast-epithelial-lineage.h5ad breast-stromal-lineage.h5ad breast-immune-lineage.h5ad)

echo "=== Post-Enrichment Validation ==="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo ""

mkdir -p "${REPORT_DIR}"

module load mamba/24.3.0
source activate hca_validators

echo "Python: $(which python3)"
echo ""

# Determine CAP command
if python3 -c "import cap_upload_validator" 2>/dev/null; then
    CAP_CMD="python3 -m cap_upload_validator"
elif command -v capval &> /dev/null; then
    CAP_CMD="capval"
else
    echo "ERROR: No CAP validator found"
    exit 1
fi

declare -a RESULT_FILE=()
declare -a RESULT_CAP=()

validate_cap() {
    local INPUT_PATH="$1"
    local BASENAME="$2"

    echo "--- Validating: ${BASENAME} ---"
    echo "  Size: $(ls -lh "${INPUT_PATH}" | awk '{print $5}')"

    ${CAP_CMD} "${INPUT_PATH}" \
        2>&1 | tee "${REPORT_DIR}/${BASENAME}_cap_post_enrich.log" || true

    CAP_LOG="${REPORT_DIR}/${BASENAME}_cap_post_enrich.log"
    RESULT_FILE+=("${BASENAME}")
    if grep -q "Validation passed" "$CAP_LOG" 2>/dev/null; then
        RESULT_CAP+=("PASS")
        echo "  -> PASS"
    else
        RESULT_CAP+=("FAIL")
        echo "  -> FAIL"
    fi
    echo ""
}

echo "=== Source Datasets ==="
for F in "${SOURCE_FILES[@]}"; do
    validate_cap "${SOURCE_DIR}/${F}" "${F%.h5ad}"
done

echo "=== Integrated Objects ==="
for F in "${INTEGRATED_FILES[@]}"; do
    validate_cap "${INTEGRATED_DIR}/${F}" "${F%.h5ad}"
done

echo ""
echo "========================================"
echo "POST-ENRICHMENT VALIDATION SUMMARY"
echo "========================================"
echo ""

TOTAL=${#RESULT_FILE[@]}
PASS_COUNT=0

printf "%-35s  %-6s\n" "FILE" "CAP"
printf "%-35s  %-6s\n" "---" "---"

for i in $(seq 0 $((TOTAL - 1))); do
    printf "%-35s  %-6s\n" "${RESULT_FILE[$i]}" "${RESULT_CAP[$i]}"
    if [[ "${RESULT_CAP[$i]}" == "PASS" ]]; then
        PASS_COUNT=$((PASS_COUNT + 1))
    fi
done

echo ""
echo "Total: ${TOTAL}, Passed: ${PASS_COUNT}, Failed: $((TOTAL - PASS_COUNT))"
echo "Done: $(date)"
