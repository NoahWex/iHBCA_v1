#!/bin/bash
#SBATCH --job-name=inspect_facs
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
import pandas as pd

print('Loading...')
adata = ad.read_h5ad('$H5AD', backed='r')
obs = adata.obs

disagree = obs[obs['facs_status'].astype(str) != obs['FACS_status'].astype(str)][
    ['dataset', 'facs_status', 'FACS_status']
]

print(f'Total disagreements: {len(disagree):,} / {len(obs):,} cells')
print()

print('=== By dataset ===')
print(disagree.groupby('dataset').size().sort_values(ascending=False).to_string())
print()

print('=== Disagreement patterns by dataset ===')
for ds, grp in disagree.groupby('dataset'):
    print(f'--- {ds} ({len(grp):,} cells) ---')
    ct = pd.crosstab(grp['facs_status'].astype(str), grp['FACS_status'].astype(str))
    print(ct.to_string())
    print()
"
