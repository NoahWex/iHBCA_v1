#!/bin/bash
# =============================================================================
# submit_build_pipeline.sh — Orchestrator for parallel Milo build pipeline
# =============================================================================
# Submits Phase 1 (prepare) then Phase 2 (Milo build) with SLURM dependency.
#
# Usage:
#   ./submit_build_pipeline.sh                    # All studies, full run
#   ./submit_build_pipeline.sh --dry-run          # All studies, dry-run Phase 2
#   ./submit_build_pipeline.sh --study 0          # Single study (gray)
#   ./submit_build_pipeline.sh --study 1,8        # XL studies (kumar, reed)
#   ./submit_build_pipeline.sh --study 8 --xl     # Reed with XL resources
#   ./submit_build_pipeline.sh --phase2-only 0    # Skip Phase 1, run Phase 2 for gray
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../../config/load_paths.sh"
BASE_PATH="${SOURCE_IHBCAV1_EXTERNAL_STUDIES%/}"
RUN_DIR="${BASE_PATH}/harmonization/run"

# Defaults
DRY_RUN=0
STUDY_INDICES=""
XL_MODE=0
PHASE2_ONLY=""

# Parse args
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        --study)
            STUDY_INDICES="$2"
            shift 2
            ;;
        --xl)
            XL_MODE=1
            shift
            ;;
        --phase2-only)
            PHASE2_ONLY="$2"
            shift 2
            ;;
        *)
            echo "Unknown arg: $1"
            echo "Usage: $0 [--dry-run] [--study N[,M]] [--xl] [--phase2-only N[,M]]"
            exit 1
            ;;
    esac
done

# Map study indices to Milo array indices (study N -> milo 2N, 2N+1)
compute_milo_indices() {
    local study_indices="$1"
    local milo_indices=""
    IFS=',' read -ra INDICES <<< "$study_indices"
    for idx in "${INDICES[@]}"; do
        local native=$((idx * 2))
        local joint=$((idx * 2 + 1))
        if [ -n "$milo_indices" ]; then
            milo_indices="${milo_indices},${native},${joint}"
        else
            milo_indices="${native},${joint}"
        fi
    done
    echo "$milo_indices"
}

# Resource overrides for XL mode
XL_PREP_ARGS=""
XL_MILO_ARGS=""
if [ "$XL_MODE" = "1" ]; then
    XL_PREP_ARGS="--mem=280G --time=4:00:00"
    XL_MILO_ARGS="--mem=280G --time=3-00:00:00 --cpus-per-task=16"
fi

# --- Phase 2 only mode ---
if [ -n "$PHASE2_ONLY" ]; then
    MILO_ARRAY=$(compute_milo_indices "$PHASE2_ONLY")
    echo "=== Phase 2 Only ==="
    echo "Study indices: $PHASE2_ONLY"
    echo "Milo array: $MILO_ARRAY"

    EXPORT_ARGS=""
    [ "$DRY_RUN" = "1" ] && EXPORT_ARGS="--export=ALL,DRY_RUN=1"

    JOB2=$(sbatch --parsable \
        --array="$MILO_ARRAY" \
        $XL_MILO_ARGS \
        $EXPORT_ARGS \
        "${RUN_DIR}/build_milo_parallel.sh")
    echo "Phase 2 job: $JOB2 (array=$MILO_ARRAY)"
    echo ""
    echo "Monitor: squeue -j $JOB2"
    exit 0
fi

# --- Full pipeline ---
echo "=== Parallel Milo Build Pipeline ==="
echo "Dry run: $DRY_RUN"
echo "XL mode: $XL_MODE"

# Phase 1: Prepare
PREP_ARRAY_ARGS=""
if [ -n "$STUDY_INDICES" ]; then
    PREP_ARRAY_ARGS="--array=$STUDY_INDICES"
fi

JOB1=$(sbatch --parsable \
    $PREP_ARRAY_ARGS \
    $XL_PREP_ARGS \
    "${RUN_DIR}/prepare_study.sh")
echo "Phase 1 (prepare): job $JOB1"

# Phase 2: Milo build (depends on Phase 1)
MILO_ARRAY_ARGS=""
if [ -n "$STUDY_INDICES" ]; then
    MILO_INDICES=$(compute_milo_indices "$STUDY_INDICES")
    MILO_ARRAY_ARGS="--array=$MILO_INDICES"
fi

EXPORT_ARGS=""
[ "$DRY_RUN" = "1" ] && EXPORT_ARGS="--export=ALL,DRY_RUN=1"

JOB2=$(sbatch --parsable \
    --dependency=afterok:${JOB1} \
    $MILO_ARRAY_ARGS \
    $XL_MILO_ARGS \
    $EXPORT_ARGS \
    "${RUN_DIR}/build_milo_parallel.sh")
echo "Phase 2 (milo):    job $JOB2 (depends on $JOB1)"

echo ""
echo "=== Pipeline Submitted ==="
echo "Phase 1: $JOB1"
echo "Phase 2: $JOB2 (starts after Phase 1 completes)"
echo ""
echo "Monitor:"
echo "  squeue -j ${JOB1},${JOB2}"
echo "  # Or tail logs:"
echo "  tail -f ${BASE_PATH}/harmonization/logs/prep_study_${JOB1}_*.out"
