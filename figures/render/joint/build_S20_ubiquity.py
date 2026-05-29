"""
§2 Supp S20 — "Is M1+M2 (baseline architecture) present in every sample?"

Single question. Per-sample matrix: fraction of non-epi cells assigned to
each motif. Rows = samples grouped by patient, columns = 9 motifs.
M1 and M2 columns should be filled at every row (universal baseline);
M4 and M6 columns should be sporadic (selective niche).
"""
import warnings, numpy as np, pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Rectangle
from pathlib import Path
warnings.filterwarnings('ignore')

mpl.rcParams.update({
    'font.family': 'Helvetica', 'font.size': 6,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
    'axes.linewidth': 0.4, 'axes.labelsize': 7,
    'xtick.labelsize': 5.5, 'ytick.labelsize': 5.5,
    'xtick.major.width': 0.4, 'ytick.major.width': 0.4,
    'xtick.major.size': 1.5, 'ytick.major.size': 1.5,
    'legend.fontsize': 6, 'legend.frameon': False,
})

OUT = Path(__file__).resolve().parent.parent / 'drafts'
MOTIFS = [f'M{i}' for i in range(9)]
PATIENTS = ['Pat1','Pat2','UCI220228','UCI604']
PAT_COLORS = {'Pat1':'#9b9b9b','Pat2':'#7B2CBF','UCI220228':'#cc4444','UCI604':'#BA68C8'}

mc = pd.read_parquet('/tmp/motif_cells_post.parquet')

# Use ALL non-epi cells (in OR out of unit) — argmax_P is the motif assignment
# Per-sample total + per-motif counts
per_sample_total = mc.groupby('sample_id').size().rename('n_total')
per_sample_motif = (mc[mc['unit_id'] != '']
                       .groupby(['sample_id','argmax_P']).size()
                       .unstack(fill_value=0)
                       .reindex(columns=MOTIFS, fill_value=0))
# Fraction = (cells in motif unit) / (total non-epi cells in sample)
frac = per_sample_motif.div(per_sample_total, axis=0)

# Attach patient + sample ordering: group by patient, then sample order
sample_to_pat = (mc[['sample_id','patient_id']].drop_duplicates()
                  .set_index('sample_id')['patient_id'])
sample_order = []
for p in PATIENTS:
    s_pat = sorted(sample_to_pat[sample_to_pat == p].index.tolist())
    sample_order.extend(s_pat)
frac = frac.reindex(sample_order)

print(f'Samples × motifs: {frac.shape}')
print(f'Fraction range: {frac.values.min():.4f} – {frac.values.max():.4f}')

# Compute M1+M2 union per sample (the "baseline"):
m1m2 = frac['M1'] + frac['M2']
print('\nM1+M2 baseline fraction:')
print(f'  min: {m1m2.min():.3f}, max: {m1m2.max():.3f}, median: {m1m2.median():.3f}')
print(f'  samples with M1+M2 == 0:  {(m1m2 == 0).sum()} / {len(m1m2)}')
print(f'  samples with M1+M2 < 0.05: {(m1m2 < 0.05).sum()} / {len(m1m2)}')
print('\nM6 niche fraction:')
m6 = frac['M6']
print(f'  min: {m6.min():.3f}, max: {m6.max():.3f}, median: {m6.median():.3f}')
print(f'  samples with M6 == 0:  {(m6 == 0).sum()} / {len(m6)}')

# Compose figure: motif heatmap + patient strip on left
n_samples = len(sample_order)
fig_h = max(3.0, 0.10 * n_samples + 1.0)
fig = plt.figure(figsize=(5.5, fig_h))
# layout: patient strip | heatmap | colorbar
gs = fig.add_gridspec(1, 3, width_ratios=[0.6, 18, 0.5], wspace=0.04,
                       left=0.18, right=0.94, top=0.92, bottom=0.10)
ax_pat = fig.add_subplot(gs[0])
ax = fig.add_subplot(gs[1], sharey=ax_pat)
ax_cb = fig.add_subplot(gs[2])

# Heatmap
mat = frac.values  # n_samples × 9
im = ax.imshow(mat, aspect='auto', cmap='viridis',
                vmin=0, vmax=np.percentile(mat, 99))
ax.set_xticks(range(len(MOTIFS)))
ax.set_xticklabels(MOTIFS, fontsize=6.5)
ax.set_yticks(range(n_samples))
ax.set_yticklabels([s.replace('_xenium_1','') for s in sample_order], fontsize=4)
ax.set_xlabel('Motif', fontsize=7)
ax.set_title('Sample × motif fraction (in-unit / total non-epi)', fontsize=7, pad=4)

# Patient strip on left
ax_pat.set_xlim(0, 1); ax_pat.set_ylim(-0.5, n_samples - 0.5)
for i, s in enumerate(sample_order):
    p = sample_to_pat[s]
    ax_pat.add_patch(Rectangle((0, i-0.5), 1, 1,
                                  facecolor=PAT_COLORS[p],
                                  edgecolor='none'))
ax_pat.set_xticks([]); ax_pat.set_yticks([])
ax_pat.invert_yaxis()
for s in ax_pat.spines.values():
    s.set_visible(False)
# Patient labels (one per group, centered)
for p in PATIENTS:
    rows = [i for i, s in enumerate(sample_order) if sample_to_pat[s] == p]
    if not rows: continue
    mid = (min(rows) + max(rows)) / 2
    ax_pat.text(-0.6, mid, p, ha='right', va='center', fontsize=6,
                 color=PAT_COLORS[p], fontweight='bold')

# Colorbar
cb = plt.colorbar(im, cax=ax_cb)
cb.outline.set_linewidth(0.3)
cb.ax.tick_params(labelsize=5.5, length=2, width=0.3)
cb.set_label('Fraction', fontsize=6, labelpad=3)

ax.invert_yaxis()

fig.savefig(OUT/'S20_ubiquity.pdf', bbox_inches='tight')
fig.savefig(OUT/'S20_ubiquity.png', dpi=600, bbox_inches='tight')
print(f'\n  → drafts/S20_ubiquity.pdf  (samples={n_samples})')

# Audit table
print('\nM1+M2 union vs M6 across samples (head):')
audit = pd.DataFrame({'M1+M2': m1m2.round(3), 'M6': m6.round(3),
                       'patient': sample_to_pat[sample_order].values})
print(audit.head(15).to_string())
