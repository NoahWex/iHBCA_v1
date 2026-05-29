"""
§2 Supp S6 v2 — "Does each motif carry the same L1.5 composition (incl epi)
across patients?"

Per-motif × per-patient stacked bar: fraction of cells associated with the
motif unit, including:
  - non-epi cells inside the unit (motif_cells_post)
  - epi cells within R=20µm of the unit boundary (per_unit_nearby_epi_counts)

Epi bars sit at the top of each stack (warm reds — visually distinct from
the non-epi cool tones). Stacked stable across patients = motif is cellular
reproducible. Different patient stacks = motif is patient-specific.
"""
import argparse, warnings, numpy as np, pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path
warnings.filterwarnings('ignore')

mpl.rcParams.update({
    'font.family': 'Helvetica', 'font.size': 6,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
    'axes.linewidth': 0.4, 'axes.labelsize': 7, 'axes.titlesize': 7,
    'xtick.labelsize': 5.5, 'ytick.labelsize': 5.5,
    'xtick.major.width': 0.4, 'ytick.major.width': 0.4,
    'xtick.major.size': 2, 'ytick.major.size': 2,
    'legend.fontsize': 5.5, 'legend.frameon': False,
})

ap = argparse.ArgumentParser()
ap.add_argument('--basis-csv', required=True, help='NMF basis matrix (program × L1.5)')
ap.add_argument('--motif-cells', required=True, help='motif_cells_post.parquet')
ap.add_argument('--epi-counts', required=True, help='per_unit_nearby_epi_counts.parquet')
ap.add_argument('--out-dir', required=True, help='Output directory')
args = ap.parse_args()
OUT = Path(args.out_dir)
MOTIFS = [f'M{i}' for i in range(9)]
PATIENTS = ['Pat1','Pat2','UCI220228','UCI604']
PROG_PAL = {
    'M0':'#E6194B','M1':'#3CB44B','M2':'#4363D8','M3':'#F58231',
    'M4':'#911EB4','M5':'#42D4F4','M6':'#F032E6','M7':'#BFEF45','M8':'#FABED4',
}
L1P5_PAL = {
    'BMYO-myo':'#8B1A1A','LASP-basal':'#CC3333','LASP':'#E55555','LHS':'#F08A8A',
    'Fb':'#4A6FA5','Fb_Activated':'#2C7FB8','Fb_SFRP4':'#6BAED6','Adipo':'#D4BC4F',
    'EC':'#2E8B57','PV':'#1A5276','LEC':'#76D7C4','Vas-cap':'#16A085',
    'CD4T':'#C71585','CD8T':'#DB7093','Treg':'#7B2CBF','T-NK':'#DA70D6',
    'cDC':'#B8860B','cDC1':'#DAA520','cDC2':'#FFD700','pDC':'#FFA500',
    'B':'#6A1B9A','Plas':'#BA68C8','Mac':'#6D4C41','Mast':'#8D6E63','Neu':'#A1887F',
}
EPI = ['LASP-basal','LASP','LHS','BMYO-myo']  # top-of-stack
basis = pd.read_csv(args.basis_csv).set_index('program')

mc = pd.read_parquet(args.motif_cells)
mc_in = mc[mc['unit_id']!=''].copy()
ne = pd.read_parquet(args.epi_counts)

# Non-epi composition per (motif, patient)
ne_comp = (mc_in.groupby(['argmax_P','patient_id','l1p5_short']).size()
              .unstack(fill_value=0))
# Epi composition per (motif, patient) — aggregated across units
epi_comp = (ne.groupby(['motif','patient_id','epi_type'])['n_epi'].sum()
              .unstack(fill_value=0))

# Build per-motif data
fig, axes = plt.subplots(3, 3, figsize=(8.5, 7.5), sharey=True)
for idx, m in enumerate(MOTIFS):
    ax = axes[idx//3, idx%3]
    # Pull non-epi composition for this motif
    if m in ne_comp.index.get_level_values('argmax_P'):
        nep = ne_comp.loc[m]
    else:
        nep = pd.DataFrame(0, index=PATIENTS, columns=[])
    nep = nep.reindex(PATIENTS, fill_value=0)
    # Pull epi composition for this motif
    if m in epi_comp.index.get_level_values('motif'):
        ep = epi_comp.loc[m]
    else:
        ep = pd.DataFrame(0, index=PATIENTS, columns=[])
    ep = ep.reindex(PATIENTS, fill_value=0)
    # Combine — make sure all columns exist
    nep = nep.reindex(columns=[c for c in nep.columns if c in L1P5_PAL], fill_value=0)
    ep = ep.reindex(columns=EPI, fill_value=0)
    combined = pd.concat([nep, ep], axis=1)
    # Fraction
    row_tot = combined.sum(axis=1).replace(0, 1)
    frac = combined.div(row_tot, axis=0)
    # Order columns: non-epi by basis weight desc, then epi (top of stack)
    if m.replace('M','P') in basis.index:
        nep_order = basis.loc[m.replace('M','P')].sort_values(ascending=False).index.tolist()
        nep_order = [c for c in nep_order if c in nep.columns]
    else:
        nep_order = nep.columns.tolist()
    col_order = nep_order + EPI
    frac = frac.reindex(columns=col_order, fill_value=0)
    # Stack
    bottom = np.zeros(len(PATIENTS))
    for col in frac.columns:
        if frac[col].sum() == 0: continue
        color = L1P5_PAL.get(col, '#cccccc')
        # Separator line between non-epi and epi: thicker edge when entering epi
        is_epi = col in EPI
        ax.bar(range(len(PATIENTS)), frac[col].values, bottom=bottom,
                color=color, edgecolor='white', linewidth=0.3 if not is_epi else 0.6,
                label=col, width=0.7)
        bottom += frac[col].values
    ax.set_title(m, fontsize=8, color=PROG_PAL[m], fontweight='bold')
    ax.set_xticks(range(len(PATIENTS)))
    ax.set_xticklabels(PATIENTS, fontsize=5, rotation=30, ha='right')
    ax.set_ylim(0, 1.02)
    if idx%3 == 0:
        ax.set_ylabel('L1.5 fraction\n(non-epi + nearby epi)', fontsize=6.5)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    # n unit-cells + n nearby-epi annotation
    for j, p in enumerate(PATIENTS):
        n_unit = int(nep.loc[p].sum()) if p in nep.index else 0
        n_epi  = int(ep.loc[p].sum()) if p in ep.index else 0
        ax.text(j, 1.03, f'{n_unit}+{n_epi}', ha='center', va='bottom', fontsize=4.5)

# Build a single combined legend with both non-epi + epi types present
present = set()
for m in MOTIFS:
    nep = ne_comp.loc[m].reindex(PATIENTS, fill_value=0) if m in ne_comp.index.get_level_values('argmax_P') else pd.DataFrame()
    if not nep.empty:
        for c in nep.columns:
            if (nep[c].sum() > 0) and (c in L1P5_PAL):
                present.add(c)
    ep = epi_comp.loc[m].reindex(PATIENTS, fill_value=0) if m in epi_comp.index.get_level_values('motif') else pd.DataFrame()
    if not ep.empty:
        for c in ep.columns:
            if ep[c].sum() > 0:
                present.add(c)
# Order: non-epi by lineage, then epi at bottom
LINEAGE_ORDER = ['Fb','Fb_Activated','Fb_SFRP4','Adipo',
                  'EC','PV','LEC','Vas-cap',
                  'CD4T','CD8T','Treg','T-NK',
                  'cDC','cDC1','cDC2','pDC',
                  'B','Plas','Mac','Mast','Neu',
                  'LASP-basal','LASP','LHS','BMYO-myo']
ord_present = [c for c in LINEAGE_ORDER if c in present]
import matplotlib.patches as mpatches
handles = [mpatches.Patch(color=L1P5_PAL[c], label=c) for c in ord_present]
fig.legend(handles=handles, loc='center right', bbox_to_anchor=(1.04, 0.5),
            fontsize=5, handlelength=0.8, handletextpad=0.4,
            labelspacing=0.25, ncol=1, frameon=False)
plt.tight_layout(rect=[0, 0, 0.93, 0.97])
OUT.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT/'S6_xpatient_composition.pdf', bbox_inches='tight')
fig.savefig(OUT/'S6_xpatient_composition.png', dpi=600, bbox_inches='tight')
print(f'  → {OUT}/S6_xpatient_composition.pdf')
