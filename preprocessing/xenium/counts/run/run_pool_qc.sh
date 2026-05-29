#!/bin/bash
#
# 4a-counts step 2: noise-floor QC consolidation across the cohort.
#
# Aggregates per-position obs.csv files from xenium_counts_per_position/ +
# xenium_counts_v2_uci604/ into cohort xenium_qc.csv with the soft pass_qc_whole
# threshold (nFeature_Whole > 5 AND tx_per_gene_Whole > 1). FLEX QC consolidated
# from per-compartment CSVs.
#
# Tier 1: 2 CPUs / 4G / 30min.

#SBATCH --job-name=xen_pool_qc
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=.slurm_stubs/pool_qc_%j.out
#SBATCH --error=.slurm_stubs/pool_qc_%j.err

if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CONTAINER_TYPE=python_spatial_2025Q4E2
source "${PIPELINE_ROOT}/run/_common.sh"

mkdir -p "${CFG_HPC_LOGS_DIR}" "$(dirname "${CFG_XENIUM_QC}")"
exec >> "${CFG_HPC_LOGS_DIR}/pool_qc_${SLURM_JOB_ID}.out" \
    2>> "${CFG_HPC_LOGS_DIR}/pool_qc_${SLURM_JOB_ID}.err"

# 01_pool_qc.py expects a pre-pooled xenium_obs.csv; build it on-the-fly by
# concatenating per-position obs.csv files from the 3-seg producer.
TMP_OBS="${JOB_TEMP_DIR}/xenium_obs_pooled.csv"
python3 -c "
import os, pandas as pd
parts = []
for root in ('${CFG_XENIUM_COUNTS_PER_POSITION}', '${CFG_XENIUM_COUNTS_V2_UCI604}'):
    if not os.path.isdir(root): continue
    for xid in sorted(os.listdir(root)):
        p = os.path.join(root, xid, 'obs.csv')
        if os.path.isfile(p):
            parts.append(pd.read_csv(p))
df = pd.concat(parts, ignore_index=True)
df.to_csv('${TMP_OBS}', index=False)
print(f'Pooled {len(df)} cells from {len(parts)} positions to ${TMP_OBS}')
"

run_python_singularity \
    "${PIPELINE_ROOT}/4a-counts/scripts/01_pool_qc.py" \
        --xenium-obs "${TMP_OBS}" \
        --flex-qc-dir "${CFG_FLEX_CONCORD_DIR}" \
        --out-dir "$(dirname "${CFG_XENIUM_QC}")"

echo "Done pool_qc: $(date)"
