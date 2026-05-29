#!/bin/bash
# =============================================================================
# Dual Build Validation: CxG build + triple-validate both HCA and CxG sets
# =============================================================================
#
# Orchestrates the full CxG build pipeline via SLURM dependency chains:
#   Phase 1: Assemble CxG h5ads (source × 7 + integrated × 1)
#   Phase 2: Enrich CxG source embeddings (needs integrated from phase 1)
#   Phase 3: Enrich CxG source metadata + integrated
#   Phase 4: Triple-validate both HCA and CxG builds
#
# Usage:
#   bash publication/run/dual_build_validation.sh [--dry-run]
#   PARTITION=free bash publication/run/dual_build_validation.sh
#
# All sbatch commands use env-var overrides — no clone scripts needed.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
REPO_ROOT="${PROJECT_IHBCAV1_UPLOAD}"
RUN_DIR="${REPO_ROOT}/publication/run"
CXG_SOURCE="${REPO_ROOT}/publication/outputs/cxg_build/source_datasets"
CXG_INTEGRATED="${REPO_ROOT}/publication/outputs/cxg_build/integrated_objects"
CXG_REPORTS="${REPO_ROOT}/publication/outputs/validation_reports/cxg_build"

# Partition override (free = no billing limits but preemptible)
PARTITION="${PARTITION:-standard}"
PARTITION_FLAG="--partition=${PARTITION}"

DRY_RUN=false
if [[ "${1:-}" == "--dry-run" ]]; then
    DRY_RUN=true
    echo "=== DRY RUN — commands will be printed but not executed ==="
    echo ""
fi

submit() {
    local desc="$1"
    shift
    echo "--- ${desc} ---" >&2
    echo "  $*" >&2
    if $DRY_RUN; then
        echo "  [dry-run: skipped]" >&2
        echo "" >&2
        echo "DRY_$$_${desc// /_}"
        return 0
    fi
    local JOB_ID
    JOB_ID=$("$@")
    echo "  Job ID: ${JOB_ID}" >&2
    echo "" >&2
    echo "$JOB_ID"
}

echo "============================================================"
echo "DUAL BUILD VALIDATION PIPELINE"
echo "Date: $(date)"
echo "============================================================"
echo ""
echo "Partition:           ${PARTITION}"
echo "CxG source dir:     ${CXG_SOURCE}"
echo "CxG integrated dir: ${CXG_INTEGRATED}"
echo "CxG reports dir:    ${CXG_REPORTS}"
echo ""

# Create output directories
mkdir -p "$CXG_SOURCE" "$CXG_INTEGRATED" "$CXG_REPORTS"

# =========================================================================
# Phase 1: Assemble CxG h5ads
# =========================================================================
echo "======== PHASE 1: CxG Assembly ========"

P1_SRC=$(submit "CxG source assembly (array 0-6)" \
    sbatch --parsable ${PARTITION_FLAG} \
    --export="ALL,TARGET=cxg,OUTPUT_DIR=${CXG_SOURCE}" \
    "${RUN_DIR}/reassemble_all_source.sh")

P1_INT=$(submit "CxG integrated assembly" \
    sbatch --parsable ${PARTITION_FLAG} \
    --export="ALL,TARGET=cxg,OUTPUT_DIR=${CXG_INTEGRATED}" \
    "${RUN_DIR}/assemble_integrated.sh")

# =========================================================================
# Phase 2: Enrich CxG source embeddings (depends on integrated + source)
# =========================================================================
echo "======== PHASE 2: CxG Source Embedding Enrichment ========"

P2_EMB=$(submit "CxG embedding enrichment (array 0-6, after phase 1)" \
    sbatch --parsable ${PARTITION_FLAG} \
    --dependency="afterok:${P1_SRC},afterok:${P1_INT}" \
    --export="ALL,SOURCE_DIR=${CXG_SOURCE}" \
    "${RUN_DIR}/enrich_source_embeddings.sh")

# =========================================================================
# Phase 3: Enrich CxG metadata (depends on embedding enrichment)
# =========================================================================
echo "======== PHASE 3: CxG Metadata Enrichment ========"

P3_SRC=$(submit "CxG source enrichment (array 0-6, after phase 2)" \
    sbatch --parsable ${PARTITION_FLAG} \
    --dependency="afterok:${P2_EMB}" \
    --export="ALL,SOURCE_DIR=${CXG_SOURCE}" \
    "${RUN_DIR}/enrich_source.sh")

P3_INT=$(submit "CxG integrated enrichment (after phase 1)" \
    sbatch --parsable ${PARTITION_FLAG} \
    --dependency="afterok:${P1_INT}" \
    --export="ALL,INTEGRATED_DIR=${CXG_INTEGRATED}" \
    "${RUN_DIR}/enrich_integrated.sh")

# =========================================================================
# Phase 4: Triple-validate both builds (depends on all enrichment)
# =========================================================================
echo "======== PHASE 4: Triple Validation ========"

P4_HCA=$(submit "HCA build validation (re-run)" \
    sbatch --parsable ${PARTITION_FLAG} \
    "${RUN_DIR}/validate_all.sh")

P4_CXG=$(submit "CxG build validation (after all CxG enrichment)" \
    sbatch --parsable ${PARTITION_FLAG} \
    --dependency="afterok:${P3_SRC},afterok:${P3_INT}" \
    --export="ALL,SOURCE_DIR=${CXG_SOURCE},INTEGRATED_DIR=${CXG_INTEGRATED},REPORT_DIR=${CXG_REPORTS}" \
    "${RUN_DIR}/validate_all.sh")

# =========================================================================
# Summary
# =========================================================================
echo "======== SUBMISSION SUMMARY ========"
echo ""
echo "Phase 1 — CxG Assembly:"
echo "  Source (array):     ${P1_SRC}"
echo "  Integrated:        ${P1_INT}"
echo ""
echo "Phase 2 — Embedding Enrichment:"
echo "  Source embeddings:  ${P2_EMB}  (after ${P1_SRC}, ${P1_INT})"
echo ""
echo "Phase 3 — Metadata Enrichment:"
echo "  Source metadata:    ${P3_SRC}  (after ${P2_EMB})"
echo "  Integrated:        ${P3_INT}  (after ${P1_INT})"
echo ""
echo "Phase 4 — Validation:"
echo "  HCA build:         ${P4_HCA}  (no deps)"
echo "  CxG build:         ${P4_CXG}  (after ${P3_SRC}, ${P3_INT})"
echo ""
echo "Monitor: squeue -u \$USER"
echo "After completion: parse logs in"
echo "  HCA: ${REPO_ROOT}/publication/outputs/validation_reports/"
echo "  CxG: ${CXG_REPORTS}/"
echo ""
echo "Total: 6 submissions, ~22 SLURM tasks"
