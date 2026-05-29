#!/bin/bash
#
# 4a-counts step 5: drop UCI604 Tumor + Ipsilateral cells from cohort.
#
# Excludes cells whose ID begins with UCI604_Tumor_xenium_1 or
# UCI604_Ipsilateral_xenium_1 (cancer + contralateral samples not part of the
# clean 65-sample reference cohort). Writes xenium_qc_clean.csv.
#
# Tier 1: 2 CPUs / 2G / 5min.

#SBATCH --job-name=xen_filter_clean
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=.slurm_stubs/filter_clean_%j.out
#SBATCH --error=.slurm_stubs/filter_clean_%j.err

if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CONTAINER_TYPE=python_spatial_2025Q4E2
source "${PIPELINE_ROOT}/run/_common.sh"

mkdir -p "${CFG_HPC_LOGS_DIR}"
exec >> "${CFG_HPC_LOGS_DIR}/filter_clean_${SLURM_JOB_ID}.out" \
    2>> "${CFG_HPC_LOGS_DIR}/filter_clean_${SLURM_JOB_ID}.err"

run_python_singularity \
    "${PIPELINE_ROOT}/4a-counts/scripts/filter_clean_cohort.py" \
        --input  "${CFG_XENIUM_QC}" \
        --output "${CFG_XENIUM_QC_CLEAN}"

echo "Done filter_clean_cohort: $(date)"
