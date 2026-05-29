#!/bin/bash
#
# 4a-counts step 1b: re-ingest UCI604 per-position count bundles.
#
# UCI604's xenium-ranger transcripts.parquet has cell_id stored as int32 while
# other patients use string hashes. The legacy producer cast cells side to str
# but left transcripts side as int, producing zero-count bundles for all UCI604
# positions. This script is the legacy 03a producer with both sides coerced
# to str. Outputs land in xenium_counts_v2_uci604/<xenium_id>/.
#
# Tier 2: 4 CPUs / 16G / 1h per task.

#SBATCH --job-name=xen_reingest_uci604
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=.slurm_stubs/reingest_uci604_%j.out
#SBATCH --error=.slurm_stubs/reingest_uci604_%j.err

if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CONTAINER_TYPE=python_scvi_gpu_2025Q3
source "${PIPELINE_ROOT}/run/_common.sh"

mkdir -p "${CFG_HPC_LOGS_DIR}" "${CFG_XENIUM_COUNTS_V2_UCI604}"
exec >> "${CFG_HPC_LOGS_DIR}/reingest_uci604_${SLURM_JOB_ID}.out" \
    2>> "${CFG_HPC_LOGS_DIR}/reingest_uci604_${SLURM_JOB_ID}.err"

run_python_singularity \
    "${PIPELINE_ROOT}/4a-counts/scripts/02_reingest_uci604.py" \
        --manifest "${CFG_XENIUM_MANIFEST}" \
        --patient  UCI604 \
        --out-dir  "${CFG_XENIUM_COUNTS_V2_UCI604}"

echo "Done reingest_uci604: $(date)"
