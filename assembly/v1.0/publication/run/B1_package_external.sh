#!/bin/bash
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
# =============================================================================
# B1_package_external.sh — Package external studies as CxG-compliant h5ad
# =============================================================================
# Two-phase pipeline:
#   Phase 1 (extract): R extracts raw counts from seurat.rds → intermediates/
#   Phase 2 (assemble): Python assembles intermediates + metadata → h5ad
#
# Usage:
#   B1_package_external.sh submit       Submit both phases (run from login node)
#   B1_package_external.sh extract      Phase 1 (run as SLURM array task)
#   B1_package_external.sh assemble     Phase 2 (run as SLURM array task)
#
# Phase 1 array: 0-8 (9 extraction tasks: 6 standalone + 3 Pal sub-studies)
# Phase 2 array: 0-6 (7 assembly tasks: 6 standalone + 1 Pal consolidated)
#
# =============================================================================

set -euo pipefail

# === Paths ===
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
REPO_ROOT="${PROJECT_IHBCAV1_UPLOAD}"
PUB_ROOT="$REPO_ROOT/publication"
EXT_ROOT="$REPO_ROOT/external_studies/outputs"
SCRIPTS="$PUB_ROOT/scripts"
OUTPUTS="$PUB_ROOT/outputs/source_datasets"
LOG_DIR="$PUB_ROOT/logs"
INTERMEDIATES="$OUTPUTS/intermediates"

# === Containers ===
R_CONTAINER="${CONTAINER_R_SPATIAL_4_3_3}"
PY_CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"
R_LIBS_USER="${R_LIBS_R_SPATIAL_4_3_3}"

# === Common bind mounts (no /scratch, no /dfs6 — they don't exist everywhere) ===
BINDS=(
    ${BIND_MOUNTS} \
)

# === Study arrays ===
# Phase 1: All studies to extract (including 3 Pal sub-studies)
EXTRACT_STUDIES=(gray kumar murrow nee twigger reed pal_norm_epi pal_norm_total pal_norm_b1)
# Phase 2: Output studies to assemble (Pal consolidated into 1)
ASSEMBLE_STUDIES=(gray kumar murrow nee twigger reed pal)

# =============================================================================
# Accept mode from $1 (direct invocation) or B1_MODE env var (sbatch --export).
# sbatch sometimes drops positional args — env var is more reliable.
# NOTE: Cannot use {braces} in :? error message — bash treats } as closing the
# parameter expansion, corrupting MODE with a trailing }.
if [[ -n "${1:-}" ]]; then
    MODE="$1"
elif [[ -n "${B1_MODE:-}" ]]; then
    MODE="$B1_MODE"
else
    echo "Usage: $0 <submit|extract|assemble>" >&2
    echo "  OR:  B1_MODE=extract sbatch $0" >&2
    exit 1
fi
# =============================================================================

case "$MODE" in

# ---- SUBMIT: Launch both phases from login node ----------------------------
submit)
    mkdir -p "$LOG_DIR" "$INTERMEDIATES"

    echo "==========================================="
    echo "B1: Package External Studies as CxG h5ad"
    echo "==========================================="
    echo ""
    echo "Phase 1: Extract counts (${#EXTRACT_STUDIES[@]} studies)"
    echo "  Studies: ${EXTRACT_STUDIES[*]}"
    echo ""
    echo "Phase 2: Assemble h5ad (${#ASSEMBLE_STUDIES[@]} outputs)"
    echo "  Studies: ${ASSEMBLE_STUDIES[*]}"
    echo ""

    SCRIPT_PATH="$(readlink -f "$0")"

    # Phase 1: R extraction array job
    # Use --wrap to set B1_MODE inline — avoids sbatch --export parsing issues
    # that cause B1_MODE=extract to arrive as $1 instead of an env var.
    # 4h/96G: reed (803K cells) produces a 24 GB mtx that needs time for
    # writeMM + gzip + copy to CRSP, and peaks at ~60-80 GB during extraction.
    PHASE1=$(sbatch --parsable \
        --account=dalawson_lab \
        --partition=standard \
        --job-name=B1_extract \
        --array=0-$((${#EXTRACT_STUDIES[@]} - 1)) \
        --cpus-per-task=4 \
        --mem=96G \
        --time=04:00:00 \
        --output="$LOG_DIR/B1_extract_%a_%j.out" \
        --error="$LOG_DIR/B1_extract_%a_%j.err" \
        --wrap "export B1_MODE=extract; bash $SCRIPT_PATH")
    echo "Phase 1 submitted: Job $PHASE1 (array 0-$((${#EXTRACT_STUDIES[@]} - 1)))"

    # Phase 2: Python assembly array job (depends on all Phase 1 tasks)
    PHASE2=$(sbatch --parsable \
        --account=dalawson_lab \
        --partition=standard \
        --job-name=B1_assemble \
        --dependency=afterok:${PHASE1} \
        --array=0-$((${#ASSEMBLE_STUDIES[@]} - 1)) \
        --cpus-per-task=4 \
        --mem=64G \
        --time=02:00:00 \
        --output="$LOG_DIR/B1_assemble_%a_%j.out" \
        --error="$LOG_DIR/B1_assemble_%a_%j.err" \
        --wrap "export B1_MODE=assemble; bash $SCRIPT_PATH")
    echo "Phase 2 submitted: Job $PHASE2 (array 0-$((${#ASSEMBLE_STUDIES[@]} - 1)), depends on $PHASE1)"

    echo ""
    echo "Monitor:"
    echo "  squeue -u \$USER -n B1_extract,B1_assemble"
    echo "  ls -lt $LOG_DIR/B1_*"
    echo ""
    echo "Tail a specific task (e.g., extract task 0):"
    echo "  tail -f $LOG_DIR/B1_extract_0_*.out"
    ;;

# ---- EXTRACT: R phase — export count matrix from seurat.rds ----------------
extract)
    INDEX="${SLURM_ARRAY_TASK_ID:?Must run as SLURM array job}"
    STUDY="${EXTRACT_STUDIES[$INDEX]}"

    INPUT_RDS="$EXT_ROOT/$STUDY/published/seurat.rds"
    OUT_DIR="$INTERMEDIATES/$STUDY"

    echo "============================================="
    echo "Phase 1: Extract counts"
    echo "  Study:  $STUDY (index $INDEX)"
    echo "  Input:  $INPUT_RDS"
    echo "  Output: $OUT_DIR"
    echo "============================================="

    if [[ ! -f "$INPUT_RDS" ]]; then
        echo "FATAL: Input RDS not found: $INPUT_RDS"
        exit 1
    fi

    mkdir -p "$OUT_DIR"

    # Kumar's published metadata.csv is 0 bytes — extract from RDS
    EXTRA_FLAGS=""
    if [[ "$STUDY" == "kumar" ]]; then
        EXTRA_FLAGS="--extract-metadata"
    fi

    module load singularity

    # Use local scratch for temp files. SLURM may set $TMPDIR to a per-job
    # directory; fall back to /tmp. Bind-mount it explicitly since we use
    # --no-mount bind-paths.
    LOCAL_TMP="${TMPDIR:-/tmp}"
    singularity exec \
        --no-mount bind-paths \
        "${BINDS[@]}" \
        --bind "$R_LIBS_USER:/home/jovyan/R/library:ro" \
        --bind "$LOCAL_TMP:$LOCAL_TMP:rw" \
        --env "R_LIBS_USER=/home/jovyan/R/library" \
        --env "TMPDIR=$LOCAL_TMP" \
        "$R_CONTAINER" \
        Rscript "$SCRIPTS/extract_counts_for_h5ad.R" \
            --study "$STUDY" \
            --input-rds "$INPUT_RDS" \
            --output-dir "$OUT_DIR" \
            $EXTRA_FLAGS

    echo ""
    echo "Phase 1 COMPLETE: $STUDY"
    echo "  Files:"
    ls -lh "$OUT_DIR/"
    ;;

# ---- ASSEMBLE: Python phase — build CxG-compliant h5ad ---------------------
assemble)
    INDEX="${SLURM_ARRAY_TASK_ID:?Must run as SLURM array job}"
    STUDY="${ASSEMBLE_STUDIES[$INDEX]}"

    echo "============================================="
    echo "Phase 2: Assemble h5ad"
    echo "  Study: $STUDY (index $INDEX)"
    echo "============================================="

    module load singularity

    singularity exec \
        --no-mount bind-paths \
        "${BINDS[@]}" \
        --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
        --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
        "$PY_CONTAINER" \
        python3 "$SCRIPTS/assemble_h5ad.py" \
            --study "$STUDY" \
            --repo-root "$REPO_ROOT"

    echo ""
    echo "Phase 2 COMPLETE: $STUDY"
    if [[ "$STUDY" == "pal" ]]; then
        ls -lh "$OUTPUTS/pal2021.h5ad" 2>/dev/null || true
    else
        FILENAME=$(python3 -c "
import yaml
with open('$PUB_ROOT/config/source_dataset_registry.yaml') as f:
    reg = yaml.safe_load(f)
print(reg['datasets']['$STUDY']['output_filename'])
" 2>/dev/null || echo "${STUDY}.h5ad")
        ls -lh "$OUTPUTS/$FILENAME" 2>/dev/null || true
    fi
    ;;

*)
    echo "Usage: $0 {submit|extract|assemble}"
    echo ""
    echo "  submit   — Submit both phases to SLURM (run from login node)"
    echo "  extract  — Phase 1: R extraction (run as SLURM array task)"
    echo "  assemble — Phase 2: Python assembly (run as SLURM array task)"
    exit 1
    ;;

esac
