#!/bin/bash
#SBATCH --job-name=s17_post
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=01:30:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/s17_post_%A_%a.out
#SBATCH --error=.slurm_stubs/s17_post_%A_%a.err
#SBATCH --array=0-3

# ============================================================================
# step_17_integration_sweep — winner postprocess + aggregate + SVD
#
# Runs the post-scoring pipeline for one target:
#   1. scib_aggregate.py  -> scib_comparison.csv + winner_selection.json
#   2. svd_diagnostic.py  -> stdout report + svd.json, fails loud on PC1>0.90
#   3. postprocess_winner.py -> umap.csv + leiden_multi.csv + metadata_with_umap.csv
#
# Array layout:
#   0  full
#   1  epi
#   2  str
#   3  imm
#
# Resources (Tier 3): 4 CPU / 48 GB / 1.5h
#   postprocess_winner on ~265K cells: neighbors+UMAP peak ~30 GB,
#   7-resolution Leiden sweep ~15 min. Compartment runs (~46-120K cells)
#   are proportionally smaller but share the same wall-time tier for
#   safety (UMAP on scVI n200 compartment embeddings is the slowest step).
#
# CPU-only (no GPU). Chains aggregate + SVD before postprocess so the
# winner is selected here rather than hardcoded in the SLURM script.
# Runs as a per-target array.
#
# Env vars:
#   MANUAL_OVERRIDE : optional config label to force as winner. Passed to
#                     scib_aggregate.py --manual-override. Use when the SVD
#                     diagnostic flags PC1 dominance and you've decided to
#                     pick a specific config by hand.
#   DRY_RUN=1       : validate inputs, skip compute
# ============================================================================

if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CONTAINER_TYPE=scgpt_gpu
source "${PIPELINE_ROOT}/run/_common.sh"


# Redirect stdout/stderr to canonical publication logs location.
# SBATCH directives above land tiny stubs in .slurm_stubs/; real output goes
# to the publication-tier ${HPC_LOGS_DIR}/ which is resolved by _common.sh.
mkdir -p "${HPC_LOGS_DIR}"
exec >> "${HPC_LOGS_DIR}/s17_post_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.out" 2>> "${HPC_LOGS_DIR}/s17_post_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.err"
TARGETS=(full epi str imm)
TARGET="${TARGETS[$SLURM_ARRAY_TASK_ID]}"

case "$TARGET" in
    full)
        SCIB_DIR="${CFG_CANONICAL_SWEEP_SCIB}/full_object"
        SWEEP_DIR="${CFG_CANONICAL_SWEEP_FULL_OBJECT}"
        WINNER_OUT="${CFG_CANONICAL_SWEEP_WINNER}/full_object"
        ;;
    epi)
        SCIB_DIR="${CFG_CANONICAL_SWEEP_SCIB}/compartments/Epithelial"
        SWEEP_DIR="${CFG_CANONICAL_SWEEP_COMPARTMENTS}/Epithelial"
        WINNER_OUT="${CFG_CANONICAL_SWEEP_WINNER}/compartments/Epithelial"
        ;;
    str)
        SCIB_DIR="${CFG_CANONICAL_SWEEP_SCIB}/compartments/Stromal"
        SWEEP_DIR="${CFG_CANONICAL_SWEEP_COMPARTMENTS}/Stromal"
        WINNER_OUT="${CFG_CANONICAL_SWEEP_WINNER}/compartments/Stromal"
        ;;
    imm)
        SCIB_DIR="${CFG_CANONICAL_SWEEP_SCIB}/compartments/Immune"
        SWEEP_DIR="${CFG_CANONICAL_SWEEP_COMPARTMENTS}/Immune"
        WINNER_OUT="${CFG_CANONICAL_SWEEP_WINNER}/compartments/Immune"
        ;;
    *)
        echo "FATAL: unknown target $TARGET" >&2
        exit 2
        ;;
esac

SVD_DIR="${CFG_CANONICAL_SWEEP_SVD}/${TARGET}"
mkdir -p "$WINNER_OUT" "$SVD_DIR"

echo "Task:       $SLURM_ARRAY_TASK_ID"
echo "Target:     $TARGET"
echo "SCIB dir:   $SCIB_DIR"
echo "Sweep dir:  $SWEEP_DIR"
echo "Winner out: $WINNER_OUT"
echo "SVD out:    $SVD_DIR"

AGG_SCRIPT="${CFG_SCRIPTS_DIR}/step_17_integration_sweep/source/scib_aggregate.py"
SVD_SCRIPT="${CFG_SCRIPTS_DIR}/step_17_integration_sweep/source/svd_diagnostic.py"
POST_SCRIPT="${CFG_SCRIPTS_DIR}/step_17_integration_sweep/source/postprocess_winner.py"

for s in "$AGG_SCRIPT" "$SVD_SCRIPT" "$POST_SCRIPT"; do
    if [ ! -f "$s" ]; then
        echo "FATAL: script not found: $s" >&2
        exit 2
    fi
done

# Optional manual override passed through to the aggregator
OVERRIDE_FLAG=""
if [ -n "${MANUAL_OVERRIDE:-}" ]; then
    OVERRIDE_FLAG="--manual-override ${MANUAL_OVERRIDE}"
fi

# Allow NaN bio metrics when no labels CSV is available (e.g. pre-L1-delivery).
# Falls back to batch-only scoring. Set ALLOW_NAN_BIO=0 to restore strict mode.
NAN_BIO_FLAG="--allow-nan-bio"
if [ "${ALLOW_NAN_BIO:-1}" = "0" ]; then
    NAN_BIO_FLAG=""
fi

# -----------------------------------------------------------------
# Stage 1: aggregate scIB JSONs into wide CSV + winner_selection.json
# -----------------------------------------------------------------
echo
echo "=== Stage 1: scib_aggregate.py ==="
run_python_singularity \
    "$AGG_SCRIPT" \
    --scib-dir "$SCIB_DIR" \
    --output-name "scib_comparison.csv" \
    $OVERRIDE_FLAG \
    $NAN_BIO_FLAG \
    ${DRY_RUN_FLAG}

# -----------------------------------------------------------------
# Stage 2: SVD diagnostic. fail_on_pc1 is pinned in sweep_configs.yaml,
# so this stage may exit non-zero and abort the pipeline before stage 3.
# That's intentional — the coordinator must inspect a failing SVD
# manually and decide whether to re-run with MANUAL_OVERRIDE.
# -----------------------------------------------------------------
echo
echo "=== Stage 2: svd_diagnostic.py ==="
run_python_singularity \
    "$SVD_SCRIPT" \
    --csv "${SCIB_DIR}/scib_comparison.csv" \
    --output "${SVD_DIR}/svd.json" \
    ${DRY_RUN_FLAG}

# -----------------------------------------------------------------
# Stage 3: winner postprocess. Reads winner_selection.json to find
# the selected config, then runs UMAP + multi-res Leiden on its h5ad.
# -----------------------------------------------------------------
echo
echo "=== Stage 3: postprocess_winner.py ==="
if [ "${DRY_RUN:-0}" = "1" ]; then
    echo "DRY RUN: skipping winner resolution; would read winner_selection.json"
    exit 0
fi

WINNER_CFG=$(python3 -c "import json,sys; d=json.load(open('${SCIB_DIR}/winner_selection.json')); print(d['selected_config'])")
if [ -z "$WINNER_CFG" ]; then
    echo "FATAL: could not read selected_config from ${SCIB_DIR}/winner_selection.json" >&2
    exit 2
fi
WINNER_H5AD="${SWEEP_DIR}/${WINNER_CFG}/integrated.h5ad"
if [ ! -f "$WINNER_H5AD" ]; then
    echo "FATAL: winner H5AD not found: $WINNER_H5AD" >&2
    exit 2
fi

echo "Winner config: $WINNER_CFG"
echo "Winner h5ad:   $WINNER_H5AD"

run_python_singularity \
    "$POST_SCRIPT" \
    --h5ad "$WINNER_H5AD" \
    --output-dir "$WINNER_OUT"

echo "Exit: $?, End: $(date)"
