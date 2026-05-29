#!/bin/bash
#SBATCH --job-name=sketch_abund
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --time=00:15:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=1

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
PROJECT="${PROJECT_IHBCAV1_UPLOAD}"
INTEGRATED=$PROJECT/publication/outputs/integrated_objects/all-breast-cells.h5ad
SKETCH=$PROJECT/publication/outputs/integrated_objects/all-breast-cells-sketch.h5ad
CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"

module load singularity

singularity exec \
    ${BIND_MOUNTS} \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
    "$CONTAINER" python3 -c "
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

print('Loading...')
full = ad.read_h5ad('$INTEGRATED', backed='r')
sk   = ad.read_h5ad('$SKETCH',     backed='r')

print(f'Full: {full.n_obs:,} cells')
print(f'Sketch: {sk.n_obs:,} cells ({100*sk.n_obs/full.n_obs:.1f}%)')

COLS = ['dataset', 'level0_annotation', 'level1_annotation', 'cell_type_ontology_term_id',
        'sex_ontology_term_id', 'assay_ontology_term_id', 'self_reported_ethnicity_ontology_term_id']

results = []
for col in COLS:
    if col not in full.obs.columns or col not in sk.obs.columns:
        continue
    full_pct = full.obs[col].astype(str).value_counts(normalize=True).sort_index()
    sk_pct   = sk.obs[col].astype(str).value_counts(normalize=True).sort_index()
    combined = pd.DataFrame({'full': full_pct, 'sketch': sk_pct}).fillna(0)
    # Max absolute deviation across categories
    mad = (combined['full'] - combined['sketch']).abs().max()
    # Pearson correlation of proportions
    corr = combined.corr().loc['full', 'sketch'] if len(combined) > 1 else 1.0
    results.append({'column': col, 'n_categories': len(combined),
                    'max_abs_dev': mad, 'pearson_r': corr})

df = pd.DataFrame(results)
print()
print('=== Sketch fidelity summary ===')
print(df.to_string(index=False, float_format='{:.4f}'.format))

# Detailed breakdown for key columns
print()
for col in ['dataset', 'level0_annotation', 'level1_annotation']:
    if col not in full.obs.columns:
        continue
    full_pct = full.obs[col].astype(str).value_counts(normalize=True).sort_index() * 100
    sk_pct   = sk.obs[col].astype(str).value_counts(normalize=True).sort_index() * 100
    combined = pd.DataFrame({'full_%': full_pct, 'sketch_%': sk_pct}).fillna(0)
    combined['diff_%'] = (combined['sketch_%'] - combined['full_%']).round(2)
    combined = combined.sort_values('full_%', ascending=False)
    print(f'--- {col} ---')
    print(combined.round(2).to_string())
    print()

# Save figure: scatter of full vs sketch proportions for level1_annotation
col = 'level1_annotation'
full_pct = full.obs[col].astype(str).value_counts(normalize=True).sort_index() * 100
sk_pct   = sk.obs[col].astype(str).value_counts(normalize=True).sort_index() * 100
combined = pd.DataFrame({'full': full_pct, 'sketch': sk_pct}).fillna(0)

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

ax = axes[0]
ax.scatter(combined['full'], combined['sketch'], s=60, alpha=0.8)
lim = max(combined.max()) * 1.05
ax.plot([0, lim], [0, lim], 'k--', lw=1, alpha=0.5)
for _, row in combined.iterrows():
    if abs(row['full'] - row['sketch']) > 1.0:
        ax.annotate(row.name, (row['full'], row['sketch']),
                    fontsize=7, ha='left', va='bottom')
ax.set_xlabel('Full integrated (%)')
ax.set_ylabel('Sketch (%)')
ax.set_title(f'level1_annotation proportions\n(full vs sketch, r={combined.corr().loc[\"full\",\"sketch\"]:.4f})')

ax = axes[1]
col2 = 'dataset'
full_pct2 = full.obs[col2].astype(str).value_counts(normalize=True).sort_index() * 100
sk_pct2   = sk.obs[col2].astype(str).value_counts(normalize=True).sort_index() * 100
combined2 = pd.DataFrame({'full': full_pct2, 'sketch': sk_pct2}).fillna(0)
x = np.arange(len(combined2))
w = 0.35
ax.bar(x - w/2, combined2['full'],  w, label='full',   alpha=0.8)
ax.bar(x + w/2, combined2['sketch'], w, label='sketch', alpha=0.8)
ax.set_xticks(x); ax.set_xticklabels(combined2.index, rotation=30, ha='right', fontsize=9)
ax.set_ylabel('% of cells'); ax.set_title('Dataset proportions')
ax.legend()

plt.tight_layout()
out = '$PROJECT/publication/outputs/figures/sketch_abundance_fidelity.png'
plt.savefig(out, dpi=150, bbox_inches='tight')
print(f'Saved: {out}')
"
echo "Done: $(date)"
