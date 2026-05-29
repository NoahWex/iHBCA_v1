#!/bin/bash
#SBATCH --job-name=inspect_obs_cols
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --time=00:15:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1

set -euo pipefail

SIF="${CONTAINER_PYTHON_SPATIAL_2025Q2}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
H5AD="${SOURCE_IHBCA_INTEGRATED}"

module load singularity

singularity exec \
    ${BIND_MOUNTS} \
    "$SIF" python3 -c "
import anndata as ad
import numpy as np
import pandas as pd

print('Loading h5ad (backed mode)...')
adata = ad.read_h5ad('$H5AD', backed='r')
obs = adata.obs

print()
print('=== FACS_status vs facs_status ===')
print()

for col in ['facs_status', 'FACS_status']:
    print(f'--- {col} ---')
    print(f'  dtype: {obs[col].dtype}')
    vc = obs[col].value_counts(dropna=False)
    print(f'  value_counts:')
    for val, cnt in vc.items():
        print(f'    {repr(val)}: {cnt:,}')
    print()

print('=== Cross-tab facs_status vs FACS_status ===')
ct = pd.crosstab(obs['facs_status'].astype(str), obs['FACS_status'].astype(str), margins=True)
print(ct.to_string())

print()
print('=== Null coverage comparison ===')
for col in ['facs_status', 'FACS_status']:
    null_mask = obs[col].isna() | (obs[col].astype(str) == 'nan') | (obs[col].astype(str) == '')
    print(f'  {col}: {null_mask.sum():,} nulls / {len(obs):,} total ({100*null_mask.sum()/len(obs):.1f}%)')

print()
print('=== Numeric cols dtype check ===')
for col in ['n_genes', 'n_counts', 'percent_mito']:
    print(f'  {col}: dtype={obs[col].dtype}, sample={list(obs[col].head(3))}')
"
