#!/bin/bash
#SBATCH --job-name=svd_only
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:05:00
#SBATCH --output=.slurm_stubs/svd_only_%j.out
#SBATCH --error=.slurm_stubs/svd_only_%j.err

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
exec >> "${HPC_LOGS_DIR}/svd_only_${SLURM_JOB_ID}.out" 2>> "${HPC_LOGS_DIR}/svd_only_${SLURM_JOB_ID}.err"
SVD_SCRIPT="${CFG_SCRIPTS_DIR}/step_17_integration_sweep/source/svd_diagnostic.py"

for comp in Epithelial Stromal Immune; do
    target=$(echo "$comp" | tr '[:upper:]' '[:lower:]' | cut -c1-3)
    csv="${CFG_CANONICAL_SWEEP_SCIB}/compartments/${comp}/scib_comparison.csv"
    out="${CFG_CANONICAL_SWEEP_SVD}/${target}/svd.json"
    mkdir -p "$(dirname "$out")"
    echo "=== SVD: $comp ==="
    run_python_singularity "$SVD_SCRIPT" --csv "$csv" --output "$out"
done
