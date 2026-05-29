#!/bin/bash
#SBATCH --job-name=s20b_promote
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/s20b_promote_%j.out
#SBATCH --error=.slurm_stubs/s20b_promote_%j.err

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
exec >> "${HPC_LOGS_DIR}/s20b_promote_${SLURM_JOB_ID}.out" 2>> "${HPC_LOGS_DIR}/s20b_promote_${SLURM_JOB_ID}.err"
SWEEP_H5AD="${CFG_CANONICAL_SWEEP_FULL_OBJECT}/scvi_n100/integrated.h5ad"
WINNER_DIR="${CFG_PREPROCESSING_ROOT}/17_ScviIntegration/winner/full_object"
OUT_DIR="${PIPELINE_ROOT}/outputs/integration_intermediate/full/scvi_n100"

echo "Sweep h5ad:  $SWEEP_H5AD"
echo "Winner dir:  $WINNER_DIR"
echo "Out dir:     $OUT_DIR"

SCRIPT="${CFG_SCRIPTS_DIR}/step_20_finalize/20b_promote_full_scvi_n100.py"
run_python_singularity \
    "$SCRIPT" \
    --sweep-h5ad "$SWEEP_H5AD" \
    --winner-dir "$WINNER_DIR" \
    --out-dir "$OUT_DIR"

echo "Exit: $?, End: $(date)"
