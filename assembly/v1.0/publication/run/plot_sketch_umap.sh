#!/bin/bash
#SBATCH --job-name=plot_sketch_umap
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --time=00:15:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
PROJECT="${PROJECT_IHBCAV1_UPLOAD}"
SKETCH=$PROJECT/publication/outputs/integrated_objects/all-breast-cells-sketch.h5ad
OUTDIR=$PROJECT/publication/outputs/figures
CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"

mkdir -p "$OUTDIR"

module load singularity

singularity exec \
    ${BIND_MOUNTS} \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
    "$CONTAINER" python3 -c "
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

print('Loading sketch...')
adata = ad.read_h5ad('$SKETCH', backed='r')
print(f'  {adata.n_obs:,} cells, obsm keys: {list(adata.obsm.keys())}')

# Pick UMAP key
umap_key = None
for k in ['X_umap', 'X_umap_joint', 'X_umap_scVI']:
    if k in adata.obsm:
        umap_key = k
        break
if umap_key is None:
    print('ERROR: no UMAP found in obsm')
    exit(1)
print(f'  Using: {umap_key}')

coords = adata.obsm[umap_key][:, :2]

COLOR_COLS = [
    ('dataset',             'tab10',  False),
    ('level0_annotation',   'Set1',   False),
    ('level1_annotation',   'tab20',  False),
    ('cell_type_ontology_term_id', 'tab20b', False),
]

fig, axes = plt.subplots(2, 2, figsize=(18, 14))
axes = axes.flatten()

for ax, (col, cmap, _) in zip(axes, COLOR_COLS):
    if col not in adata.obs.columns:
        ax.set_visible(False)
        continue
    labels = adata.obs[col].astype(str).values
    unique = sorted(set(labels))
    cmap_obj = plt.get_cmap(cmap, len(unique))
    label_to_color = {l: cmap_obj(i) for i, l in enumerate(unique)}
    colors = [label_to_color[l] for l in labels]

    ax.scatter(coords[:, 0], coords[:, 1], c=colors, s=0.3, linewidths=0, rasterized=True)
    ax.set_title(col, fontsize=11)
    ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
    ax.set_xticks([]); ax.set_yticks([])

    # legend (outside right for large sets, inside for small)
    handles = [plt.Line2D([0],[0], marker='o', color='w',
               markerfacecolor=label_to_color[l], markersize=6, label=l)
               for l in unique]
    ncol = max(1, len(unique) // 20)
    ax.legend(handles=handles, fontsize=5, ncol=ncol,
              loc='lower right' if len(unique) <= 10 else 'lower left',
              markerscale=1.5, framealpha=0.7)

plt.suptitle(f'iHBCA sketch ({adata.n_obs:,} cells) — {umap_key}', fontsize=13)
plt.tight_layout()

out = '$OUTDIR/sketch_umap.png'
plt.savefig(out, dpi=150, bbox_inches='tight')
print(f'Saved: {out}')
"

echo "Done: $(date)"
