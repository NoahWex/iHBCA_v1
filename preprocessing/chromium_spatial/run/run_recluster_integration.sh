#!/bin/bash
#SBATCH --job-name=recluster
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:45:00
#SBATCH --array=0-2
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/recluster_%A_%a.out
#SBATCH --error=.slurm_stubs/recluster_%A_%a.err

# Cluster patch job. Loads the existing compartment integrated.h5ad,
# runs sc.pp.neighbors + sc.tl.leiden at user-specified params, and writes
# a sidecar CSV under outputs/integration_intermediate/.../patches/. The
# Rmd's load_patches chunk picks it up automatically on the next render.
#
# Container: scgpt_gpu (Python + scanpy + anndata)
# Tier:      4cpu/24G/45min — neighbor graph on 100K cells in 100D fits in <8G
#
# Required env vars (export via --sbatch "--export=ALL,VAR=value"):
#   PATCH_ID         identifier for this patch (also column prefix unless overridden)
#   RESOLUTIONS      COLON-separated leiden resolutions, e.g. "0.4:0.6:1.5:2.0"
#                    (colons because sbatch --export uses commas as key separators
#                    and shell-escaping the commas does not survive parsing).
# Optional:
#   N_NEIGHBORS      neighbor count (default 30)
#   USE_REP          obsm key (default X_emb)
#   COL_PREFIX       column name prefix (default = PATCH_ID)
#   RANDOM_SEED      leiden seed (default 0)
#   DRY_RUN=1        validate inputs and exit
#
# Tasks:
#   0  Epithelial    102,952 cells
#   1  Stromal       119,141 cells
#   2  Immune         47,942 cells

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
exec >> "${HPC_LOGS_DIR}/recluster_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.out" 2>> "${HPC_LOGS_DIR}/recluster_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.err"
if [ -z "${PATCH_ID:-}" ]; then
    echo "FATAL: PATCH_ID env var not set" >&2
    exit 2
fi
if [ -z "${RESOLUTIONS:-}" ]; then
    echo "FATAL: RESOLUTIONS env var not set" >&2
    exit 2
fi

N_NEIGHBORS="${N_NEIGHBORS:-30}"
USE_REP="${USE_REP:-X_emb}"
COL_PREFIX="${COL_PREFIX:-${PATCH_ID}}"
RANDOM_SEED="${RANDOM_SEED:-0}"

# Translate colon-separated RESOLUTIONS into the comma-separated form
# 06c_recluster_integration.py expects.
RESOLUTIONS_PY=$(echo "${RESOLUTIONS}" | tr ':' ',')

COMPARTMENTS=(Epithelial Stromal Immune)
COMPARTMENT=${COMPARTMENTS[$SLURM_ARRAY_TASK_ID]}
CONFIG_LABEL=scvi_n100

INPUT_H5AD="${CFG_UPSTREAM_ROOT}/compartments_step15/${COMPARTMENT}/integration/${CONFIG_LABEL}/integrated.h5ad"
PATCHES_DIR="${CFG_PROJECT_ROOT}/outputs/integration_intermediate/${COMPARTMENT}/${CONFIG_LABEL}/patches"

mkdir -p "${PATCHES_DIR}"

echo "Compartment:  ${COMPARTMENT}"
echo "Config:       ${CONFIG_LABEL}"
echo "Input h5ad:   ${INPUT_H5AD}"
echo "Patches dir:  ${PATCHES_DIR}"
echo "Patch id:     ${PATCH_ID}"
echo "use_rep:      ${USE_REP}"
echo "n_neighbors:  ${N_NEIGHBORS}"
echo "resolutions:  ${RESOLUTIONS_PY}  (raw env: ${RESOLUTIONS})"
echo "col_prefix:   ${COL_PREFIX}"
echo "random_seed:  ${RANDOM_SEED}"

run_python_singularity \
    "${CFG_SCRIPTS_DIR}/06c_recluster_integration.py" \
    --h5ad         "${INPUT_H5AD}" \
    --patches-dir  "${PATCHES_DIR}" \
    --patch-id     "${PATCH_ID}" \
    --use-rep      "${USE_REP}" \
    --n-neighbors  "${N_NEIGHBORS}" \
    --resolutions  "${RESOLUTIONS_PY}" \
    --col-prefix   "${COL_PREFIX}" \
    --random-seed  "${RANDOM_SEED}" \
    ${DRY_RUN_FLAG}

echo "Exit: $?, End: $(date)"
