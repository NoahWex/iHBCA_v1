"""
§2 Supp S25 — "What L1.5 cell types define each motif's NMF factor?"

Loading heatmap: 9 motifs × 21 L1.5 non-epi cell types, color = basis weight.
Rows ordered M0..M8 (motif palette in y-tick color). Cols ordered by lineage
(Stromal → Vascular → T → DC → B/Plas → Myeloid).
Cells where basis ≥ 0.05 are flagged with a thin black outline (the "signal
carrier" threshold).
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
    'xtick.labelsize': 6, 'ytick.labelsize': 7,
    'xtick.major.width': 0.4, 'ytick.major.width': 0.4,
    'xtick.major.size': 2, 'ytick.major.size': 2,
    'legend.fontsize': 6, 'legend.frameon': False,
})

OUT = Path(__file__).resolve().parent.parent / 'drafts'
MOTIFS = [f'M{i}' for i in range(9)]
PROG_PAL = {
    'M0':'#E6194B','M1':'#3CB44B','M2':'#4363D8','M3':'#F58231',
    'M4':'#911EB4','M5':'#42D4F4','M6':'#F032E6','M7':'#BFEF45','M8':'#FABED4',
}
# Lineage ordering of non-epi L1.5 columns (left → right)
L1P5_ORDER = [
    'Fb','Fb_Activated','Fb_SFRP4','Adipo',
    'EC','PV','LEC','Vas-cap',
    'CD4T','CD8T','Treg','T-NK',
    'cDC','cDC1','cDC2','pDC',
    'B','Plas',
    'Mac','Mast','Neu',
]
SIG_THRESHOLD = 0.05

basis = pd.read_csv('/Users/NoahWechter/Library/Application Support/CRSP Desktop/Volumes.noindex/CRSP Lab - dalawson.localized/nwechter/Spatial_HBCA_Xenium/NicheFramework/fig4_lane_c_nmf_20260515/data/nmf/basis.csv').set_index('program')
basis.index = [p.replace('P','M') for p in basis.index]
basis = basis.reindex(MOTIFS)
# Reorder columns by L1P5_ORDER
cols_present = [c for c in L1P5_ORDER if c in basis.columns]
basis = basis[cols_present]

# Figure: wider than tall for the matrix layout
fig, ax = plt.subplots(figsize=(6.8, 3.4))
mat = basis.values
im = ax.imshow(mat, aspect='auto', cmap='viridis',
                vmin=0, vmax=max(0.5, mat.max()))
# Outline signal carriers (basis ≥ threshold)
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        if mat[i, j] >= SIG_THRESHOLD:
            ax.add_patch(Rectangle((j-0.5, i-0.5), 1, 1,
                                      facecolor='none', edgecolor='black',
                                      linewidth=0.5, zorder=3))
# Annotate weights ≥ 0.10 (top carriers) inline
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        if mat[i, j] >= 0.10:
            color = 'white' if mat[i, j] >= 0.50 else 'black'
            ax.text(j, i, f'{mat[i,j]:.2f}', ha='center', va='center',
                     fontsize=4.5, color=color)
# Axis ticks
ax.set_xticks(range(len(cols_present)))
ax.set_xticklabels(cols_present, fontsize=5.5, rotation=45, ha='right')
ax.set_yticks(range(len(MOTIFS)))
ax.set_yticklabels(MOTIFS, fontsize=7)
# Color y-tick labels by motif palette
for i, m in enumerate(MOTIFS):
    ax.get_yticklabels()[i].set_color(PROG_PAL[m])
    ax.get_yticklabels()[i].set_fontweight('bold')
ax.set_xlabel('L1.5 cell type', fontsize=7)
ax.set_ylabel('Motif', fontsize=7)

# Colorbar (right margin)
cax = fig.add_axes([0.97, 0.30, 0.012, 0.45])
cb = plt.colorbar(im, cax=cax)
cb.outline.set_linewidth(0.3)
cb.ax.tick_params(labelsize=5.5, length=2, width=0.3)
cb.set_label('Basis weight (H)', fontsize=6.5, labelpad=3)

# Lineage block separators (vertical dashed lines)
LINEAGE_BREAKS = [4, 8, 12, 16, 18]
for b in LINEAGE_BREAKS:
    ax.axvline(b - 0.5, color='black', lw=0.4, ls=':', alpha=0.5)

# Compact legend for the outline meaning
ax.text(1.02, -0.16, f'■ basis ≥ {SIG_THRESHOLD} (signal carrier)',
         transform=ax.transAxes, fontsize=5, ha='right', va='top')

plt.tight_layout(rect=[0, 0, 0.94, 1])
fig.savefig(OUT/'S25_loading_heatmap.pdf', bbox_inches='tight')
fig.savefig(OUT/'S25_loading_heatmap.png', dpi=600, bbox_inches='tight')
print(f'  → drafts/S25_loading_heatmap.pdf')
