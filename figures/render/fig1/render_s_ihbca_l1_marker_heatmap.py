#!/usr/bin/env python3
"""
Expression heatmap keyed on working labels from build_label_map.py.

Columns = working labels (aggregated from res 1.0 clusters, weighted by cell count).
Rows    = top N marker genes per label by pooled LOR (ENSG excluded).
Column order = lineage groups: CD8 T → CD4 T → NK/ILC → B → Plasma → Myeloid → artifacts.

Output: synthesis/figures/label_heatmap_{mode}.png

Modes:
  markers     top N per label by pooled LOR (default)
  lineage     curated lineage gene lists
  compartment curated compartment gene lists
  artifact    curated artifact gene lists
"""

import argparse
import csv
import os
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import yaml
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist


# Lineage display order for column grouping (artifact labels at end)
LINEAGE_ORDER = [
    'CD8_T', 'T_ambiguous', 'CD4_T', 'Treg',
    'NK', 'ILC',
    'B', 'plasma',
    'macrophage', 'DC', 'monocyte', 'mast', 'neutrophil', 'pDC',
    'artifact',
]

LINEAGE_COLORS = {
    'CD8_T':      '#aed6f1',
    'T_ambiguous':'#85c1e9',
    'CD4_T':      '#d6eaf8',
    'Treg':       '#d5d8dc',
    'NK':         '#a9dfbf',
    'ILC':        '#82e0aa',
    'B':          '#d7bde2',
    'plasma':     '#c9b1d9',
    'macrophage': '#fad7a0',
    'DC':         '#f9e4b7',
    'monocyte':   '#f5cba7',
    'mast':       '#fadbd8',
    'neutrophil': '#fdebd0',
    'pDC':        '#e8daef',
    'artifact':   '#eaeded',
}

FLAG_COLORS = {
    'IEG':        '#e74c3c',
    'heat_shock': '#e67e22',
    'proliferation': '#9b59b6',
    'MT_high':    '#7f8c8d',
}


def load_yaml(path):
    with open(path) as f:
        return yaml.safe_load(f)


def load_expression(expr_csv):
    df = pd.read_csv(expr_csv, index_col='gene_id')
    sym_map = {}
    if 'symbol' in df.columns:
        sym_map = df['symbol'].dropna().to_dict()
        df = df.drop(columns=['symbol'])
    df.columns = [str(c) for c in df.columns]
    return df, sym_map


def build_label_expression(label_map_csv, expr_df):
    """Weighted-mean expression per label across contributing res 1.0 clusters."""
    lmap = pd.read_csv(label_map_csv, index_col='cell_id')
    counts = lmap.groupby(['label', 'cluster_1.0']).size().reset_index(name='n')

    label_expr = {}
    for label, grp in counts.groupby('label'):
        weighted = np.zeros(len(expr_df))
        total = 0
        for _, row in grp.iterrows():
            cid = str(int(row['cluster_1.0']))
            if cid in expr_df.columns:
                weighted += expr_df[cid].values * row['n']
                total += row['n']
        if total > 0:
            label_expr[label] = weighted / total
    return label_expr


def ordered_labels(label_map_csv, ann_yaml):
    """Return labels sorted by lineage group order, with cell counts."""
    lmap = pd.read_csv(label_map_csv, index_col='cell_id')
    label_counts = lmap['label'].value_counts().to_dict()

    clusters = ann_yaml['clusters']
    # Build label → lineage from YAML (stripped _1 from working_label)
    label_lineage = {}
    for info in clusters.values():
        wl = info.get('working_label', '')
        lineage = info.get('lineage', 'artifact')
        label_lineage[wl] = lineage

    result = []
    for lineage in LINEAGE_ORDER:
        group = [(lbl, n) for lbl, n in label_counts.items()
                 if label_lineage.get(lbl, 'artifact') == lineage]
        group.sort(key=lambda x: -x[1])
        result.extend(group)

    # Any labels not captured by lineage order (fallback)
    seen = {lbl for lbl, _ in result}
    for lbl, n in sorted(label_counts.items(), key=lambda x: -x[1]):
        if lbl not in seen:
            result.append((lbl, n))

    return result  # [(label, n_cells), ...]


def load_pooled_logreg(logreg_csv, sym_map, exclude_ensg=True):
    """Returns dict: cluster -> list of {gene, symbol, lor}"""
    data = defaultdict(list)
    with open(logreg_csv) as f:
        for row in csv.DictReader(f):
            lor_s = row.get('pooled_log_or', '')
            if not lor_s:
                continue
            lor = float(lor_s)
            if lor <= 0:
                continue
            gene = row['gene']
            sym = sym_map.get(gene, gene)
            if exclude_ensg and sym.startswith('ENSG'):
                continue
            data[row['cluster']].append({'gene': gene, 'symbol': sym, 'lor': lor})
    return data


def select_marker_genes(logreg, label_map_csv, ann_yaml, top_per_label, sym_map,
                        exclude_ensg=True):
    """Top N genes per label by pooled LOR, ordered by cross-label breadth."""
    lmap = pd.read_csv(label_map_csv, index_col='cell_id')

    # Map cluster → label(s)
    cluster_to_labels = defaultdict(set)
    for _, row in lmap[['label', 'cluster_1.0']].drop_duplicates().iterrows():
        cluster_to_labels[str(int(row['cluster_1.0']))].add(row['label'])

    # Build label → top genes
    label_genes = defaultdict(list)
    for cluster, rows in logreg.items():
        labels = cluster_to_labels.get(cluster, {cluster})
        top = sorted(rows, key=lambda r: -r['lor'])[:top_per_label]
        for lbl in labels:
            label_genes[lbl].extend(top)

    # Collect candidates; track max LOR and cross-label count
    gene_max_lor = defaultdict(float)
    gene_label_count = defaultdict(int)
    candidate_syms = set()

    for lbl, rows in label_genes.items():
        seen_in_label = set()
        top = sorted(rows, key=lambda r: -r['lor'])[:top_per_label]
        for r in top:
            sym = r['symbol']
            if sym not in seen_in_label:
                candidate_syms.add(sym)
                seen_in_label.add(sym)
                gene_label_count[sym] += 1
                gene_max_lor[sym] = max(gene_max_lor[sym], r['lor'])

    ordered = sorted(candidate_syms,
                     key=lambda s: (-gene_label_count[s], -gene_max_lor[s]))
    return ordered


def select_curated_genes(gene_lists, mode, sym_map, logreg):
    if mode == 'lineage':
        categories = gene_lists.get('lineage_markers', {})
    elif mode == 'compartment':
        categories = gene_lists.get('compartment_markers', {})
    elif mode == 'artifact':
        categories = gene_lists.get('artifact_markers', {})
    else:
        raise ValueError(mode)

    sym_to_max_lor = defaultdict(float)
    for rows in logreg.values():
        for r in rows:
            sym_to_max_lor[r['symbol']] = max(sym_to_max_lor[r['symbol']], r['lor'])

    result = []
    seen = set()
    for cat, genes in categories.items():
        for sym in genes:
            if sym in seen or sym not in sym_to_max_lor:
                continue
            seen.add(sym)
            result.append(sym)
    return result


def select_canonical_genes(ann_yaml, label_map_csv):
    """One entry per canonical marker per label, ordered by label similarity position."""
    clusters = ann_yaml['clusters']
    lmap = pd.read_csv(label_map_csv, index_col='cell_id')
    cluster_to_label = {}
    for _, row in lmap[['label', 'cluster_1.0']].drop_duplicates().iterrows():
        cluster_to_label[str(int(row['cluster_1.0']))] = row['label']

    # Collect (label, symbol) pairs — one row per symbol per label
    seen = set()
    result = []  # list of (label, symbol)
    for key, info in clusters.items():
        cid = key.lstrip('c')
        lbl = cluster_to_label.get(cid)
        if lbl is None:
            continue
        for sym in info.get('canonical_markers', []):
            if isinstance(sym, str) and not sym.startswith('ENSG'):
                if (lbl, sym) not in seen:
                    seen.add((lbl, sym))
                    result.append((lbl, sym))
    return result  # [(label, symbol), ...]


def get_artifact_flags(sym, gene_lists):
    flags = []
    for cat in ['IEG', 'heat_shock', 'proliferation', 'MT_high']:
        if sym in gene_lists.get('artifact_markers', {}).get(cat, []):
            flags.append(cat)
    return flags


def render_label_heatmap(label_map_csv, ann_yaml_path, logreg_csv, expr_csv,
                         gene_lists_path, out_path, mode='markers',
                         top_per_label=8, exclude_artifacts=True):
    ann = load_yaml(ann_yaml_path)
    gene_lists = load_gene_lists_safe(gene_lists_path)
    expr_df, sym_map = load_expression(expr_csv)
    logreg = load_pooled_logreg(logreg_csv, sym_map)

    label_order = ordered_labels(label_map_csv, ann)
    if exclude_artifacts:
        # Get artifact lineage labels
        lmap = load_yaml(ann_yaml_path)
        artifact_labels = set()
        for info in lmap['clusters'].values():
            if info.get('lineage') == 'artifact':
                wl = info.get('working_label', '')
                artifact_labels.add(wl)
        label_order = [(lbl, n) for lbl, n in label_order
                       if lbl not in artifact_labels]

    labels = [lbl for lbl, _ in label_order]
    label_ncells = {lbl: n for lbl, n in label_order}

    # Build gene × label expression matrix
    label_expr = build_label_expression(label_map_csv, expr_df)

    if mode == 'top_markers':
        symbols = select_marker_genes(logreg, label_map_csv, ann,
                                      top_per_label, sym_map)
    elif mode == 'canonical_markers':
        canonical_pairs = select_canonical_genes(ann, label_map_csv)
        symbols = list(dict.fromkeys(sym for _, sym in canonical_pairs))
    else:
        symbols = select_curated_genes(gene_lists, mode, sym_map, logreg)

    # Filter to symbols with expression data
    sym_to_ensg = {v: k for k, v in sym_map.items()}
    def has_expr(sym):
        ensg = sym_to_ensg.get(sym, sym)
        return ensg in expr_df.index or sym in expr_df.index
    symbols = [s for s in symbols if has_expr(s)]

    if not symbols:
        print(f'No genes found for mode {mode}')
        return

    # Build matrix
    def get_gene_vector(sym):
        ensg = sym_to_ensg.get(sym, sym)
        if ensg in expr_df.index:
            idx = expr_df.index.get_loc(ensg)
        elif sym in expr_df.index:
            idx = expr_df.index.get_loc(sym)
        else:
            return None
        return idx

    n_genes = len(symbols)
    n_labels = len(labels)
    matrix = np.zeros((n_genes, n_labels))
    for i, sym in enumerate(symbols):
        idx = get_gene_vector(sym)
        if idx is None:
            continue
        for j, lbl in enumerate(labels):
            if lbl in label_expr:
                matrix[i, j] = label_expr[lbl][idx]

    # Drop all-zero rows
    nonzero = np.any(matrix != 0, axis=1)
    symbols = [symbols[i] for i in range(len(symbols)) if nonzero[i]]
    matrix = matrix[nonzero, :]

    if len(symbols) == 0:
        print('No expression data found')
        return

    # Cluster columns by expression similarity, then diagonalize rows
    if n_labels > 1:
        col_dist = pdist(matrix.T, metric='correlation')
        col_dist = np.nan_to_num(col_dist, nan=1.0)
        col_link = linkage(col_dist, method='average')
        col_order = leaves_list(col_link)
    else:
        col_order = list(range(n_labels))

    labels = [labels[i] for i in col_order]
    label_ncells = {lbl: label_ncells[lbl] for lbl in labels}
    matrix = matrix[:, col_order]

    # Re-diagonalize rows against clustered column order
    row_order = np.argsort(np.argmax(matrix, axis=1), kind='stable')
    symbols = [symbols[i] for i in row_order]
    matrix = matrix[row_order, :]

    # Lineage color bar for columns
    clusters_ann = ann['clusters']
    label_lin = {}
    for info in clusters_ann.values():
        wl = info.get('working_label', '')
        lin = info.get('lineage', 'artifact')
        label_lin[wl] = lin

    n_genes = len(symbols)
    fig_h = max(6, n_genes * 0.22 + 2.5)
    fig_w = max(10, n_labels * 0.5 + 4)
    vmax = min(abs(matrix).max(), 3.0)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    cmap = plt.get_cmap('RdBu_r').copy()
    im = ax.imshow(matrix, aspect='auto', cmap=cmap,
                   vmin=-vmax, vmax=vmax, interpolation='nearest')

    # Lineage color strip above columns
    for j, lbl in enumerate(labels):
        lin = label_lin.get(lbl, 'artifact')
        color = LINEAGE_COLORS.get(lin, '#f0f0f0')
        ax.add_patch(mpatches.Rectangle(
            (j - 0.5, -1.5), 1, 1,
            color=color, clip_on=False, zorder=3))

    # Lineage group separators
    prev_lin = None
    for j, lbl in enumerate(labels):
        lin = label_lin.get(lbl, 'artifact')
        if lin != prev_lin and j > 0:
            ax.axvline(j - 0.5, color='white', lw=1.5, zorder=4)
        prev_lin = lin

    ax.set_xlim(-0.5, n_labels - 0.5)
    ax.set_xticks(range(n_labels))
    ax.set_xticklabels(
        [f'{lbl}\n({label_ncells.get(lbl,0):,})' for lbl in labels],
        fontsize=5.5, rotation=45, ha='right')
    ax.set_yticks(range(n_genes))
    ax.set_yticklabels(symbols, fontsize=6)

    # Flag badges (top_markers mode)
    if mode == 'top_markers':
        ax.set_xlim(-0.5, n_labels - 0.5 + 1.5)
        for i, sym in enumerate(symbols):
            flags = get_artifact_flags(sym, gene_lists)
            for k, flag in enumerate(flags[:2]):
                color = FLAG_COLORS.get(flag, '#bdc3c7')
                rect = mpatches.FancyBboxPatch(
                    (n_labels + 0.1 + k * 0.45, i - 0.35), 0.4, 0.7,
                    boxstyle='round,pad=0.05',
                    facecolor=color, edgecolor='none', alpha=0.85,
                    transform=ax.transData, clip_on=False)
                fig.add_artist(rect)
                ax.text(n_labels + 0.3 + k * 0.45, i,
                        flag[:4], ha='center', va='center',
                        fontsize=4, color='white', fontweight='bold',
                        transform=ax.transData, clip_on=False)

    # Lineage legend
    seen_lin = []
    for lbl in labels:
        lin = label_lin.get(lbl, 'artifact')
        if lin not in seen_lin:
            seen_lin.append(lin)
    handles = [mpatches.Patch(color=LINEAGE_COLORS.get(l, '#f0f0f0'), label=l)
               for l in seen_lin]
    fig.legend(handles=handles, title='Lineage', fontsize=5, title_fontsize=5.5,
               loc='lower right', bbox_to_anchor=(1.0, 0.0), ncol=2)

    mode_label = {'top_markers': f'top {top_per_label}/label',
                  'canonical_markers': 'canonical markers',
                  'lineage': 'lineage markers',
                  'compartment': 'compartment markers',
                  'artifact': 'artifact markers'}.get(mode, mode)
    ax.set_title(f'Immune — working labels  ({mode_label}, z-scored)',
                 fontsize=9, pad=16)
    plt.colorbar(im, ax=ax, label='z-scored mean expression', shrink=0.5, pad=0.01)
    plt.tight_layout()

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out_path}')


def load_gene_lists_safe(path):
    if path and os.path.exists(path):
        with open(path) as f:
            return yaml.safe_load(f)
    return {}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--label-map', required=True)
    parser.add_argument('--annotation-yaml', required=True)
    parser.add_argument('--data-dir', required=True)
    parser.add_argument('--out-dir', required=True)
    parser.add_argument('--resolution', default='1.0')
    parser.add_argument('--gene-lists', default=None)
    parser.add_argument('--top-per-label', type=int, default=6)
    parser.add_argument('--mode', default='top_markers',
                        choices=['top_markers', 'canonical_markers',
                                 'lineage', 'compartment', 'artifact'])
    parser.add_argument('--include-artifacts', action='store_true')
    args = parser.parse_args()

    logreg_csv = os.path.join(args.data_dir, 'logreg',
                              f'logreg_markers_leiden_{args.resolution}.csv')
    expr_csv = os.path.join(args.data_dir, 'mean_expression',
                            f'scaled_mean_expr_leiden_{args.resolution}.csv')
    out_path = os.path.join(args.out_dir, 'figures',
                            f'label_heatmap_{args.mode}.png')

    render_label_heatmap(
        label_map_csv=args.label_map,
        ann_yaml_path=args.annotation_yaml,
        logreg_csv=logreg_csv,
        expr_csv=expr_csv,
        gene_lists_path=args.gene_lists,
        out_path=out_path,
        mode=args.mode,
        top_per_label=args.top_per_label,
        exclude_artifacts=not args.include_artifacts,
    )


if __name__ == '__main__':
    main()
