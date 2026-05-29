#!/bin/bash
#
# Step 20 — Canonical manifest (finalization).
#
# Pattern src: run_step_18_deseq.sh (single-task, Tier 1)
#
# Walks integration_intermediate/ and deseq/, validates all expected files
# exist, records cell counts, writes outputs/manifest.yaml.
# Run after step 18 DESeq completes.
#
# Tier 1: validation + YAML write only, no heavy computation.

#SBATCH --job-name=step_20_finalize
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=.slurm_stubs/step_20_finalize_%j.out
#SBATCH --error=.slurm_stubs/step_20_finalize_%j.err

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
exec >> "${HPC_LOGS_DIR}/step_20_finalize_${SLURM_JOB_ID}.out" 2>> "${HPC_LOGS_DIR}/step_20_finalize_${SLURM_JOB_ID}.err"
SCRIPT="${CFG_SCRIPTS_DIR}/step_20_finalize/20_build_manifest.py"

echo "--- Step 20 Canonical Manifest ---"
echo "Outputs root: ${CFG_OUTPUTS_ROOT}"
echo "Manifest out: ${CFG_CANONICAL_MANIFEST}"
echo "----------------------------------"

run_python_singularity python "${SCRIPT}" \
    --outputs-root "${CFG_OUTPUTS_ROOT}" \
    --out          "${CFG_CANONICAL_MANIFEST}"

echo "Exit: $?, End: $(date)"
