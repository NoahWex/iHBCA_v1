#!/bin/bash
#SBATCH --job-name=s17_scib
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/s17_scib_%A_%a.out
#SBATCH --error=.slurm_stubs/s17_scib_%A_%a.err
# Default array: 17 configs x 8 metrics = 136 tasks (for SWEEP_TARGET=full).
# For SWEEP_TARGET=all submit --array=0-543 (17*8*4).
#SBATCH --array=0-135

# ============================================================================
# step_17_integration_sweep — scIB metrics stage
#
# Computes ONE scIB metric for ONE config (+ target) per array task.
# Each task writes {scib_dir}/{target_subdir}/{config}/{metric}.json.
# The aggregation step (scib_aggregate.py) runs AFTER this array completes.
#
# Reads the 8-metric grid (kBET included as a sidecar). Per-target subdir
# under the scib_results root. CPU-only (no GPU needed for scIB metrics).
#
# Resources (Tier 3): 4 CPU / 64 GB / 2h
#   Peak usage: PCR on 265K cells holds two AnnData copies in memory (~50 GB).
#   isolated_labels/graph_conn peak ~15 GB. kBET ~20 GB via pynndescent.
#
# Env vars:
#   SWEEP_TARGET : full|epi|str|imm (default: full). One target per array.
#                  For SWEEP_TARGET=all, submit 4 arrays separately or expand
#                  the index math below.
#   DRY_RUN=1    : validate inputs, skip compute
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
exec >> "${HPC_LOGS_DIR}/s17_scib_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.out" 2>> "${HPC_LOGS_DIR}/s17_scib_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.err"
SWEEP_TARGET="${SWEEP_TARGET:-full}"

# Must match sweep_configs.yaml:configs ordering and metrics union
CONFIGS=(scvi_n20 scvi_n50 scvi_n100 scvi_n200 \
         scvi_austin_n20 scvi_austin_n50 scvi_austin_n100 scvi_austin_n200 \
         harmony_n20 harmony_n50 harmony_n100 harmony_n200 \
         harmony_theta1 harmony_theta5 \
         concord_n50 concord_n100 pca_n50)
# 8 metrics = 5 scored + 3 sidecar (per sweep_configs.yaml:metrics)
METRICS=(ASW_batch ASW_label NMI ARI graph_conn isolated_labels PCR kBET)

N_CONFIGS=${#CONFIGS[@]}   # 17
N_METRICS=${#METRICS[@]}  # 8

CONFIG_IDX=$((SLURM_ARRAY_TASK_ID / N_METRICS))
METRIC_IDX=$((SLURM_ARRAY_TASK_ID % N_METRICS))

if [ "$CONFIG_IDX" -ge "$N_CONFIGS" ]; then
    echo "FATAL: CONFIG_IDX=$CONFIG_IDX >= N_CONFIGS=$N_CONFIGS" >&2
    echo "  Check --array range vs N_CONFIGS*N_METRICS" >&2
    exit 2
fi

CONFIG=${CONFIGS[$CONFIG_IDX]}
METRIC=${METRICS[$METRIC_IDX]}

# Resolve target subdir and input H5AD path
case "$SWEEP_TARGET" in
    full)
        H5AD="${CFG_CANONICAL_SWEEP_FULL_OBJECT}/${CONFIG}/integrated.h5ad"
        SCIB_OUT_BASE="${CFG_CANONICAL_SWEEP_SCIB}/full_object"
        ;;
    epi)
        H5AD="${CFG_CANONICAL_SWEEP_COMPARTMENTS}/Epithelial/${CONFIG}/integrated.h5ad"
        SCIB_OUT_BASE="${CFG_CANONICAL_SWEEP_SCIB}/compartments/Epithelial"
        ;;
    str)
        H5AD="${CFG_CANONICAL_SWEEP_COMPARTMENTS}/Stromal/${CONFIG}/integrated.h5ad"
        SCIB_OUT_BASE="${CFG_CANONICAL_SWEEP_SCIB}/compartments/Stromal"
        ;;
    imm)
        H5AD="${CFG_CANONICAL_SWEEP_COMPARTMENTS}/Immune/${CONFIG}/integrated.h5ad"
        SCIB_OUT_BASE="${CFG_CANONICAL_SWEEP_SCIB}/compartments/Immune"
        ;;
    *)
        echo "FATAL: unknown SWEEP_TARGET=$SWEEP_TARGET" >&2
        exit 2
        ;;
esac

OUTPUT_DIR="${SCIB_OUT_BASE}/${CONFIG}"
mkdir -p "$OUTPUT_DIR"

# Resolve the bio-labels CSV.
# Primary: CFG_CANONICAL_L1_LABELS (iHBCA_V1 L2.0c-derived, post-step_20).
# Fallback: per-compartment Kumar proxy labels (label_kumar_2023 normalized),
#   used when the primary doesn't exist yet. These give a genuine bio signal
#   today; re-run with primary labels when L2.0c delivery is complete.
LABELS_CSV="${CFG_CANONICAL_L1_LABELS}"
if [ ! -f "${LABELS_CSV}" ]; then
    KUMAR_LABELS_DIR="${CFG_PREPROCESSING_ROOT}/preprocessing_wrapup/labels_for_scib"
    case "$SWEEP_TARGET" in
        epi) LABELS_CSV="${KUMAR_LABELS_DIR}/labels_epithelial.csv" ;;
        str) LABELS_CSV="${KUMAR_LABELS_DIR}/labels_stromal.csv" ;;
        imm) LABELS_CSV="${KUMAR_LABELS_DIR}/labels_immune.csv" ;;
    esac
    echo "NOTE: Using Kumar proxy labels: $LABELS_CSV"
fi

echo "Task:       $SLURM_ARRAY_TASK_ID"
echo "Target:     $SWEEP_TARGET"
echo "Config:     $CONFIG"
echo "Metric:     $METRIC"
echo "H5AD:       $H5AD"
echo "Output dir: $OUTPUT_DIR"
echo "Labels:     $LABELS_CSV"

SCRIPT="${CFG_SCRIPTS_DIR}/step_17_integration_sweep/source/scib_metric.py"
if [ ! -f "$SCRIPT" ]; then
    echo "FATAL: $SCRIPT not found" >&2
    exit 2
fi

# Build optional --labels-csv arg only if it exists on disk (otherwise
# scib_metric.py will log the absence and return NaN, per fail-loud docs).
LABELS_ARG=""
if [ -n "$LABELS_CSV" ] && [ -f "$LABELS_CSV" ]; then
    LABELS_ARG="--labels-csv $LABELS_CSV"
fi

run_python_singularity \
    "$SCRIPT" \
    --h5ad "$H5AD" \
    --metric "$METRIC" \
    --batch-key patient_id \
    --label-key consensus_label \
    $LABELS_ARG \
    --output-dir "$OUTPUT_DIR" \
    ${DRY_RUN_FLAG}

echo "Exit: $?, End: $(date)"
