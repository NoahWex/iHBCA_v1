#!/bin/bash
#SBATCH --job-name=s17_sweep
#SBATCH --partition=free-gpu
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --gres=gpu:1
#SBATCH --time=08:00:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/s17_sweep_%A_%a.out
#SBATCH --error=.slurm_stubs/s17_sweep_%A_%a.err
# Default array: 0-16 (17 configs x 1 target = 17 tasks). For target=all,
# submit with --array=0-67 (17 configs x 4 targets = 68 tasks) and set
# SWEEP_TARGET=all.
#SBATCH --array=0-16

# ============================================================================
# step_17_integration_sweep — sweep stage
#
# Runs one (config, target) pair of the integration sweep per array task.
# Delegates to scripts/step_17_integration_sweep/wrapper.py which enforces
# the target x config grid and resolves all paths from config/paths.yaml.
#
# Full-object and per-compartment runs are unified into a single wrapper
# selected by the SWEEP_TARGET env var.
#
# Resources (Tier 5): 8 CPU / 128 GB / 1 GPU / 8h
#   scVI training on ~265K cells with n_latent=200 peaks ~80 GB and runs ~2-4h.
#   Harmony/PCA stages peak ~40 GB and run in 20-40 min — this profile covers
#   the worst case (scvi_austin_n200 with batch_size=1024).
#
# SLURM array layout options (controlled by SWEEP_TARGET):
#   SWEEP_TARGET=full  (default): --array=0-11  (12 configs x 1 target)
#   SWEEP_TARGET=epi              --array=0-11
#   SWEEP_TARGET=str              --array=0-16
#   SWEEP_TARGET=imm              --array=0-16
#   SWEEP_TARGET=all              --array=0-67 (17 configs x 4 targets)
#                                  Task i -> target = TARGETS[i / 17],
#                                            config = i % 17
# Config order (0-16): scvi_n20, scvi_n50, scvi_n100, scvi_n200,
#   scvi_austin_n20..n200, harmony_n20, harmony_n50, harmony_n100, harmony_n200,
#   harmony_theta1, harmony_theta5, concord_n50, concord_n100, pca_n50
#
# Required env vars (pass via sbatch --export=ALL,SWEEP_TARGET=full):
#   SWEEP_TARGET: one of full|epi|str|imm|all (default: full)
# Optional:
#   DRY_RUN=1     : validate inputs, skip training
#   TEST_MODE=1   : pass --test to wrapper (5000 cells, 10 epochs)
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
exec >> "${HPC_LOGS_DIR}/s17_sweep_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.out" 2>> "${HPC_LOGS_DIR}/s17_sweep_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.err"
SWEEP_TARGET="${SWEEP_TARGET:-full}"
N_CONFIGS=17   # must match len(sweep_configs.yaml:configs)

# Resolve (target, config_idx) for this task
TARGETS_ALL=(full epi str imm)
if [ "$SWEEP_TARGET" = "all" ]; then
    TARGET_IDX=$(( SLURM_ARRAY_TASK_ID / N_CONFIGS ))
    CONFIG_IDX=$(( SLURM_ARRAY_TASK_ID % N_CONFIGS ))
    if [ "$TARGET_IDX" -ge 4 ]; then
        echo "FATAL: TARGET_IDX=$TARGET_IDX out of range for --array 0-$(( 4 * N_CONFIGS - 1 ))" >&2
        exit 2
    fi
    TARGET="${TARGETS_ALL[$TARGET_IDX]}"
else
    TARGET="$SWEEP_TARGET"
    CONFIG_IDX="$SLURM_ARRAY_TASK_ID"
fi

if [ "$CONFIG_IDX" -ge "$N_CONFIGS" ]; then
    echo "FATAL: CONFIG_IDX=$CONFIG_IDX exceeds N_CONFIGS=$N_CONFIGS" >&2
    echo "  Check that --array matches sweep_configs.yaml length." >&2
    exit 2
fi

echo "Task:    $SLURM_ARRAY_TASK_ID"
echo "Target:  $TARGET"
echo "Config:  index $CONFIG_IDX"
echo "GPU:     ${CUDA_VISIBLE_DEVICES:-none}"

# Optional test / dry-run
EXTRA_FLAGS=""
if [ "${TEST_MODE:-0}" = "1" ]; then
    EXTRA_FLAGS="$EXTRA_FLAGS --test"
fi

# Path validation
STEP_17_WRAPPER="${CFG_SCRIPTS_DIR}/step_17_integration_sweep/wrapper.py"
if [ ! -f "$STEP_17_WRAPPER" ]; then
    echo "FATAL: wrapper not found: $STEP_17_WRAPPER" >&2
    exit 2
fi

run_python_singularity \
    "$STEP_17_WRAPPER" \
    --target "$TARGET" \
    --config-idx "$CONFIG_IDX" \
    --project-root "$CFG_PROJECT_ROOT" \
    --temp-dir "$JOB_TEMP_DIR" \
    $EXTRA_FLAGS \
    ${DRY_RUN_FLAG}

echo "Exit: $?, End: $(date)"
