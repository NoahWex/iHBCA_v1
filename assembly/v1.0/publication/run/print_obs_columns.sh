#!/bin/bash
#SBATCH --job-name=print_obs_cols
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --time=00:10:00
#SBATCH --mem=8G
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
import h5py
import numpy as np

path = '$H5AD'
with h5py.File(path, 'r') as f:
    col_order = list(f['obs'].attrs.get('column-order', []))
    all_keys = sorted(f['obs'].keys())

    # Sample values for context
    print(f'=== obs columns ({len(col_order)} in column-order, {len(all_keys)} total keys) ===')
    print()
    for i, col in enumerate(col_order):
        if col in f['obs']:
            obj = f['obs'][col]
            # Try to get a few sample values
            try:
                if hasattr(obj, 'keys'):  # categorical (has codes/categories)
                    if 'categories' in obj:
                        cats = list(obj['categories'][:5])
                        n_cats = len(obj['categories'])
                        print(f'  {i+1:3d}. {col}  [categorical, {n_cats} levels, e.g. {cats}]')
                    else:
                        print(f'  {i+1:3d}. {col}  [group]')
                else:
                    vals = obj[:3]
                    print(f'  {i+1:3d}. {col}  [e.g. {list(vals)}]')
            except Exception as e:
                print(f'  {i+1:3d}. {col}  [err: {e}]')
        else:
            print(f'  {i+1:3d}. {col}  [MISSING from keys]')

    extra = [k for k in all_keys if k not in col_order and not k.startswith('_')]
    if extra:
        print()
        print(f'Keys not in column-order: {extra}')
"
