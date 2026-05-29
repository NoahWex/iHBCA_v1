#!/bin/bash
#SBATCH --job-name=inspect_raw_var
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --time=00:10:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
REPO="${PROJECT_IHBCAV1_UPLOAD}"
SIF="${CONTAINER_PYTHON_SPATIAL_2025Q2}"

INTEGRATED="$REPO/publication/outputs/integrated_objects/all-breast-cells.h5ad"
SKETCH="$REPO/publication/outputs/integrated_objects/sketch.h5ad"

echo "=== Inspecting raw/var column-order ==="

for H5AD in "$INTEGRATED" "$SKETCH"; do
    echo ""
    echo "--- $H5AD ---"
    module load singularity
    singularity exec \
        ${BIND_MOUNTS} \
        "$SIF" python3 -c "
import h5py
path = '$H5AD'
with h5py.File(path, 'r') as f:
    has_raw = 'raw' in f
    has_raw_var = has_raw and 'var' in f['raw']
    print('has_raw:', has_raw)
    print('has_raw_var:', has_raw_var)
    if has_raw_var:
        col_order = list(f['raw/var'].attrs.get('column-order', []))
        keys = sorted(f['raw/var'].keys())
        print('raw/var keys:', keys[:5], '...' if len(keys) > 5 else '')
        print('feature_is_filtered in keys:', 'feature_is_filtered' in keys)
        print('feature_is_filtered in column-order:', 'feature_is_filtered' in col_order)
        print('column-order length:', len(col_order))
        if 'feature_is_filtered' in col_order:
            idx = col_order.index('feature_is_filtered')
            print('position in column-order:', idx)
"
done

echo ""
echo "=== Done ==="
