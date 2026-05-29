#!/bin/bash
#
# 4a-counts step 3: merge UCI604 v2 bundles into cohort xenium_qc.csv.
#
# Replaces the zero-count UCI604 rows in xenium_qc.csv with the correctly-typed
# v2 obs from xenium_counts_v2_uci604/, applies the Phase 1 pass_qc_whole
# threshold, and rewrites xenium_qc.csv in place.
#
# Tier 1: 2 CPUs / 4G / 15min.

#SBATCH --job-name=xen_merge_uci604
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=.slurm_stubs/merge_uci604_%j.out
#SBATCH --error=.slurm_stubs/merge_uci604_%j.err

if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CONTAINER_TYPE=python_spatial_2025Q4E2
source "${PIPELINE_ROOT}/run/_common.sh"

mkdir -p "${CFG_HPC_LOGS_DIR}"
exec >> "${CFG_HPC_LOGS_DIR}/merge_uci604_${SLURM_JOB_ID}.out" \
    2>> "${CFG_HPC_LOGS_DIR}/merge_uci604_${SLURM_JOB_ID}.err"

run_python_singularity \
    "${PIPELINE_ROOT}/4a-counts/scripts/03_merge_uci604_and_rescan.py" \
        --v2-root "${CFG_XENIUM_COUNTS_V2_UCI604}" \
        --xqc     "${CFG_XENIUM_QC}" \
        --fqc     "${CFG_FLEX_QC}" \
        --audit   "${CFG_COHORT_QC_AUDIT}"

echo "Done merge_uci604: $(date)"
