"""Xenium x FLEX 2x2 marker comparison + cytoplasm artifact screen + FLEX-context clusters.

2x2 layout:
  (Xenium nuclear markers, Xenium nuclear expression)   | (Xenium nuclear markers, FLEX expression)
  (FLEX markers, FLEX expression)                       | (FLEX markers, Xenium nuclear expression)

Expression plots (2x2 and cyto screen) are generated at every requested resolution.
Feature plots are generated at the primary --res only.

FLEX-context section: FLEX cells clustered in their own latent subspace at
multiple resolutions, with coverage/compression stats vs joint Leiden clusters.

Usage:
  python3 05_2x2_marker_comparison.py \\
      --compartment Immune \\
      --latent outputs/Immune_clean/joint_latent.csv \\
      --obs    outputs/Immune_clean/joint_obs.csv \\
      --nuc-bundle  <pooled_nuclear_clean/> \\
      --cyto-bundle <pooled_cytoplasmic_clean/> \\
      --flex-dir    <integration_intermediate/> \\
      --flex-qc     <flex_qc.csv> \\
      --out-dir outputs/Immune_clean/panels2x2 \\
      --html    outputs/Immune_clean/immune_2x2.html \\
      -k 30 --res 0.3 --n-top 4
"""
import argparse, base64, gzip, io, os
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as _mp
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scipy.io, scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ── parallel Wilcoxon worker ───────────────────────────────────────────────────

def _wilcoxon_worker(task):
    """Top-level picklable worker: reconstruct sparse AnnData, run Wilcoxon, return DataFrame."""
    X_data, X_indices, X_indptr, shape, var_names, cluster_labels, n_genes = task
    import scipy.sparse as _sp, anndata as _ad, pandas as _pd, scanpy as _sc
    X = _sp.csr_matrix((X_data, X_indices, X_indptr), shape=shape)
    obs = _pd.DataFrame({"leiden": _pd.Categorical(cluster_labels)})
    a = _ad.AnnData(X=X, obs=obs)
    a.var_names = _pd.Index(var_names)
    _sc.tl.rank_genes_groups(a, "leiden", method="wilcoxon", n_genes=n_genes, pts=True)
    return _sc.get.rank_genes_groups_df(a, group=None)


def _adata_to_task(adata, cluster_labels, n_genes):
    """Package an AnnData into a picklable task tuple for _wilcoxon_worker."""
    X = sp.csr_matrix(adata.X)
    return (X.data, X.indices, X.indptr, X.shape,
            list(adata.var_names), list(cluster_labels), n_genes)


def run_markers_parallel(tasks_dict, n_workers):
    """Run multiple Wilcoxon DE tests in parallel.

    tasks_dict: {key: task_tuple}  (key = arbitrary label)
    Returns:    {key: markers_DataFrame}
    """
    results = {}
    ctx = _mp.get_context("spawn")
    with ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as pool:
        futures = {pool.submit(_wilcoxon_worker, task): key
                   for key, task in tasks_dict.items()}
        for fut in as_completed(futures):
            key = futures[fut]
            results[key] = fut.result()
    return results


# ── helpers ────────────────────────────────────────────────────────────────────

def load_mtx(mtx_gz, cells_tsv, genes_tsv):
    with gzip.open(mtx_gz, "rb") as gz:
        mat = scipy.io.mmread(io.BytesIO(gz.read()))
    cells = pd.read_csv(cells_tsv, header=None)[0].tolist()
    genes = pd.read_csv(genes_tsv, header=None)[0].tolist()
    a = ad.AnnData(X=sp.csr_matrix(mat), obs=pd.DataFrame(index=cells))
    a.var_names = genes
    return a


def load_flex(flex_dir, comp, qc_path):
    d = os.path.join(flex_dir, comp, "scvi_n100")
    a = load_mtx(
        os.path.join(d, "counts.mtx.gz"),
        os.path.join(d, "cells.tsv"),
        os.path.join(d, "genes.tsv"),
    )
    qc = pd.read_csv(qc_path, index_col="cell_id")
    keep = qc.index[qc["pass_qc"]].intersection(a.obs_names)
    return a[keep].copy()


def top_genes(markers_df, n):
    """Top-n per cluster, deduplicated, preserving score order."""
    top = (markers_df.groupby("group")
           .apply(lambda g: g.nlargest(n, "scores"))
           .reset_index(drop=True))
    seen, ordered = set(), []
    for g in top["names"]:
        if g not in seen:
            seen.add(g)
            ordered.append(g)
    return ordered


def run_markers(a, groupby, n_genes=25):
    sc.tl.rank_genes_groups(a, groupby=groupby, method="wilcoxon",
                            n_genes=n_genes, pts=True)
    return sc.get.rank_genes_groups_df(a, group=None)


def dotplot(a, genes, groupby, title, outpath):
    genes_in = [g for g in genes if g in a.var_names]
    if not genes_in:
        print(f"  SKIP {os.path.basename(outpath)}: no matching genes", flush=True)
        return False
    sc.pl.dotplot(a, var_names=genes_in, groupby=groupby, show=False,
                  figsize=(max(6, len(genes_in) * 0.45 + 2),
                           max(4, a.obs[groupby].nunique() * 0.35 + 2)))
    fig = plt.gcf()
    fig.suptitle(title, fontsize=8, y=1.01)
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  {os.path.basename(outpath)}", flush=True)
    return True


def pct_nz_per_cluster(adata, groupby, genes):
    rows = {}
    sub = adata[:, genes]
    for cl, grp in sub.obs.groupby(groupby):
        X = sub[grp.index].X
        if sp.issparse(X):
            rows[cl] = np.array((X > 0).mean(axis=0)).ravel()
        else:
            rows[cl] = (X > 0).mean(axis=0)
    return pd.DataFrame(rows, index=genes).T


def parse_position(cell_id, patient_id):
    rest = cell_id[len(patient_id) + 1:]
    parts = rest.split("_")
    pos_parts = []
    for part in parts:
        if part.lower() in ("xenium", "flex") or len(part) > 12:
            break
        pos_parts.append(part)
    return "_".join(pos_parts) if pos_parts else "unknown"


def patient_position_grid(x, y, obs, out_dir):
    patients = sorted(obs["patient_id"].dropna().unique())
    positions = sorted(obs["position"].dropna().unique())
    cmap = plt.get_cmap("tab20", len(positions))
    pos_palette = {pos: cmap(i) for i, pos in enumerate(positions)}

    ncols = min(4, len(patients))
    nrows = (len(patients) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 3.5), squeeze=False)
    axes_flat = axes.ravel()

    for idx, pat in enumerate(patients):
        ax = axes_flat[idx]
        mask = (obs["patient_id"] == pat).values
        xp, yp = x[mask], y[mask]
        pos_vals = obs.loc[obs["patient_id"] == pat, "position"].fillna("NA").values
        colors = [pos_palette.get(pos, (0.7, 0.7, 0.7, 1.0)) for pos in pos_vals]
        order = np.random.default_rng(0).permutation(mask.sum())
        ax.scatter(xp[order], yp[order], c=[colors[i] for i in order],
                   s=0.8, alpha=0.3, rasterized=True, linewidths=0)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(pat, fontsize=8)
        for spine in ax.spines.values():
            spine.set_linewidth(0.25)
        for pos in pd.unique(pos_vals):
            pmask = pos_vals == pos
            if pmask.sum() < 5:
                continue
            ax.text(xp[pmask].mean(), yp[pmask].mean(), pos,
                    fontsize=4.5, ha="center", va="center", fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.1", facecolor="white",
                              alpha=0.6, linewidth=0))

    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], marker="o", color="w",
                      markerfacecolor=pos_palette.get(pos, (0.7, 0.7, 0.7)),
                      markersize=5, label=pos)
               for pos in positions]
    fig.legend(handles=handles, loc="lower center", ncol=min(8, len(positions)),
               fontsize=6, frameon=False, bbox_to_anchor=(0.5, -0.02))

    for idx in range(len(patients), len(axes_flat)):
        axes_flat[idx].axis("off")

    fig.suptitle("UMAP split by patient — colored by position", fontsize=9)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    outpath = os.path.join(out_dir, "umap_patient_x_position.png")
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  umap_patient_x_position.png", flush=True)
    return outpath


def scatter_cat(ax, x, y, cat, title, s=1.0, alpha=0.3, seed=0, label_centroids=False):
    from matplotlib.lines import Line2D
    cats = pd.Categorical(cat)
    codes, names = cats.codes, cats.categories
    cmap = plt.get_cmap("tab20", max(len(names), 1))
    palette = [cmap(i) for i in range(len(names))]
    order = np.random.default_rng(seed).permutation(len(x))
    ax.scatter(x[order], y[order], c=[palette[c] for c in codes[order]],
               s=s, alpha=alpha, rasterized=True, linewidths=0, edgecolors="none")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, fontsize=8)
    for spine in ax.spines.values():
        spine.set_linewidth(0.25)
    if label_centroids:
        x_arr, y_arr = np.array(x, dtype=float), np.array(y, dtype=float)
        for i, name in enumerate(names):
            mask = codes == i
            if mask.sum() < 5:
                continue
            cx, cy = x_arr[mask].mean(), y_arr[mask].mean()
            ax.text(cx, cy, str(name), fontsize=5, ha="center", va="center",
                    fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.1", facecolor="white",
                              alpha=0.65, linewidth=0))
    handles = [Line2D([0], [0], marker="o", color="w",
                      markerfacecolor=palette[i], markersize=4, label=names[i])
               for i in range(len(names))]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5),
              fontsize=5, frameon=False, markerscale=1.5)


def scatter_cont(ax, x, y, vals, title, cmap="viridis", vmin=None, vmax=None,
                 s=0.5, alpha=0.25, seed=0):
    m = np.isfinite(vals)
    order = np.random.default_rng(seed).permutation(len(x))
    sc_obj = ax.scatter(x[order][m[order]], y[order][m[order]],
                        c=vals[order][m[order]], cmap=cmap,
                        vmin=vmin, vmax=vmax, s=s, alpha=alpha,
                        rasterized=True, linewidths=0)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, fontsize=8)
    for spine in ax.spines.values():
        spine.set_linewidth(0.25)
    plt.colorbar(sc_obj, ax=ax, shrink=0.6, pad=0.02)


def umap_cat_panel(x, y, obs, key, title, out_dir, s=0.5, alpha=0.25, label_centroids=False):
    outpath = os.path.join(out_dir, f"umap_{key}.png")
    fig, ax = plt.subplots(figsize=(5, 4))
    scatter_cat(ax, x, y, obs[key].fillna("NA").values, title,
                s=s, alpha=alpha, label_centroids=label_centroids)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  umap_{key}.png", flush=True)
    return outpath


def umap_cont_panel(x, y, vals, title, out_dir, stem, cmap="viridis",
                    vmin=None, vmax=None, s=0.5, alpha=0.25):
    outpath = os.path.join(out_dir, f"umap_{stem}.png")
    fig, ax = plt.subplots(figsize=(5, 4))
    scatter_cont(ax, x, y, vals, title, cmap=cmap, vmin=vmin, vmax=vmax, s=s, alpha=alpha)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  umap_{stem}.png", flush=True)
    return outpath


def feature_plot_grid(x, y, adata, genes, leiden_labels, out_dir, stem, s=0.5, alpha=0.3):
    genes_in = [g for g in genes if g in adata.var_names]
    if not genes_in:
        return None
    gene_idx = {g: adata.var_names.get_loc(g) for g in genes_in}
    ncols = min(6, len(genes_in))
    nrows = (len(genes_in) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3, nrows * 2.8), squeeze=False)
    axes_flat = axes.ravel()

    centroids = {}
    if leiden_labels is not None:
        for cl in pd.unique(leiden_labels):
            mask = leiden_labels == cl
            if mask.sum() >= 5:
                centroids[cl] = (x[mask].mean(), y[mask].mean())

    for i, g in enumerate(genes_in):
        col = gene_idx[g]
        if sp.issparse(adata.X):
            vals = np.array(adata.X[:, col].todense()).ravel().astype(float)
        else:
            vals = adata.X[:, col].astype(float)
        scatter_cont(axes_flat[i], x, y, vals, g, cmap="Reds", s=s, alpha=alpha)
        for cl, (cx, cy) in centroids.items():
            axes_flat[i].text(cx, cy, str(cl), fontsize=4.5, ha="center", va="center",
                              fontweight="bold",
                              bbox=dict(boxstyle="round,pad=0.1", facecolor="white",
                                        alpha=0.55, linewidth=0))
    for i in range(len(genes_in), len(axes_flat)):
        axes_flat[i].axis("off")
    fig.tight_layout()
    outpath = os.path.join(out_dir, f"featureplot_{stem}.png")
    fig.savefig(outpath, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  featureplot_{stem}.png", flush=True)
    return outpath


def flex_cluster_umap(ux, uy, flex_mask, fl_cats, fl_codes, title, out_dir, stem):
    """UMAP with Xenium cells gray, FLEX cells colored by FLEX-only leiden."""
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(ux, uy, s=0.3, alpha=0.1, color="#cccccc", rasterized=True, linewidths=0)
    cmap = plt.get_cmap("tab20", max(len(fl_cats), 1))
    palette = [cmap(i) for i in range(len(fl_cats))]
    ux_f, uy_f = ux[flex_mask], uy[flex_mask]
    order = np.random.default_rng(0).permutation(flex_mask.sum())
    ax.scatter(ux_f[order], uy_f[order], c=[palette[fl_codes[i]] for i in order],
               s=1.0, alpha=0.5, rasterized=True, linewidths=0)
    for i, name in enumerate(fl_cats):
        mask = fl_codes == i
        if mask.sum() < 5:
            continue
        cx, cy = ux_f[mask].mean(), uy_f[mask].mean()
        ax.text(cx, cy, str(name), fontsize=5, ha="center", va="center",
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.1", facecolor="white", alpha=0.7, linewidth=0))
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, fontsize=8)
    fig.tight_layout()
    outpath = os.path.join(out_dir, f"umap_flex_only_{stem}.png")
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  umap_flex_only_{stem}.png", flush=True)
    return outpath


def flex_coverage_heatmap(fl_obs, fl_leiden_key, joint_leiden_key, out_dir, stem):
    """Heatmap: FLEX-only cluster x joint cluster, row-normalized (coverage)."""
    ct = pd.crosstab(fl_obs[fl_leiden_key], fl_obs[joint_leiden_key], normalize="index")
    fig, ax = plt.subplots(figsize=(max(6, len(ct.columns) * 0.45 + 2),
                                    max(4, len(ct) * 0.4 + 2)))
    im = ax.imshow(ct.values, aspect="auto", cmap="YlOrRd", vmin=0, vmax=1)
    ax.set_xticks(range(len(ct.columns)))
    ax.set_xticklabels(ct.columns, fontsize=6, rotation=90)
    ax.set_yticks(range(len(ct)))
    ax.set_yticklabels(ct.index, fontsize=6)
    ax.set_xlabel(f"Joint Leiden ({joint_leiden_key})", fontsize=8)
    ax.set_ylabel(f"FLEX-only cluster ({fl_leiden_key})", fontsize=8)
    ax.set_title("Coverage: fraction of FLEX-only cluster cells per joint cluster (row-normalized)", fontsize=8)
    plt.colorbar(im, ax=ax, shrink=0.8, label="Fraction")
    fig.tight_layout()
    outpath = os.path.join(out_dir, f"flex_coverage_{stem}.png")
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  flex_coverage_{stem}.png", flush=True)
    return outpath


def flex_compression_bar(fl_obs, fl_leiden_key, joint_leiden_key, out_dir, stem):
    """Bar chart: number of joint clusters covered (>=5% of cells) per FLEX cluster."""
    ct = pd.crosstab(fl_obs[fl_leiden_key], fl_obs[joint_leiden_key], normalize="index")
    covered = (ct >= 0.05).sum(axis=1)
    fig, ax = plt.subplots(figsize=(max(5, len(covered) * 0.6 + 1), 3))
    ax.bar(range(len(covered)), covered.values, color="#4c72b0")
    ax.set_xticks(range(len(covered)))
    ax.set_xticklabels(covered.index, fontsize=7)
    ax.set_xlabel("FLEX-only cluster", fontsize=8)
    ax.set_ylabel("# joint clusters (>=5%)", fontsize=8)
    ax.set_title(f"Coverage breadth: FLEX cluster compression in joint embedding ({fl_leiden_key} x {joint_leiden_key})", fontsize=8)
    ax.axhline(covered.mean(), color="k", linestyle="--", lw=0.8, label=f"Mean={covered.mean():.1f}")
    ax.legend(fontsize=7)
    fig.tight_layout()
    outpath = os.path.join(out_dir, f"flex_compression_{stem}.png")
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  flex_compression_{stem}.png", flush=True)
    return outpath


def b64(path):
    if not path or not os.path.exists(path):
        name = os.path.basename(path) if path else "unknown"
        return f'<p class="missing">{name} not found</p>'
    with open(path, "rb") as f:
        data = base64.b64encode(f.read()).decode()
    return f'<img src="data:image/png;base64,{data}" style="max-width:100%;">'


def section(title, body, note=None):
    note_html = f'<p class="note">{note}</p>' if note else ""
    return f'<div class="section"><h3>{title}</h3>{note_html}{body}</div>'


def row2(*imgs):
    cells = "".join(f'<div class="cell2">{i}</div>' for i in imgs)
    return f'<div class="row2">{cells}</div>'


def row3(*imgs):
    cells = "".join(f'<div class="cell3">{i}</div>' for i in imgs)
    return f'<div class="row3">{cells}</div>'


def build_tab_section(ns, title, note, tab_content_map, all_res, primary_res):
    """Build a tabbed HTML section. tab_content_map: {res_float -> html_string}."""
    btns, panes = [], []
    for r in all_res:
        stem = str(r).replace(".", "p")
        active_cls = " active" if r == primary_res else ""
        btns.append(
            '<button class="tab-btn' + active_cls + '" data-ns="' + ns + '" '
            'id="' + ns + '-btn-' + stem + '" '
            'onclick="showTab(\'' + ns + '\',\'' + stem + '\')">' + str(r) + '</button>'
        )
        content = tab_content_map.get(r, "<p class='missing'>not generated</p>")
        panes.append(
            '<div class="tab-pane' + active_cls + '" data-ns="' + ns + '" '
            'id="' + ns + '-pane-' + stem + '">' + content + '</div>'
        )
    bar = '<div class="tab-bar">' + "".join(btns) + "</div>"
    note_html = f'<p class="note">{note}</p>' if note else ""
    return '<div class="section"><h3>' + title + '</h3>' + note_html + bar + "".join(panes) + '</div>'


# ── main ───────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--compartment", required=True)
    ap.add_argument("--latent", required=True)
    ap.add_argument("--obs",    required=True)
    ap.add_argument("--nuc-bundle",  required=True)
    ap.add_argument("--cyto-bundle", required=True)
    ap.add_argument("--flex-dir",    required=True)
    ap.add_argument("--flex-qc",     required=True)
    ap.add_argument("--out-dir",     required=True)
    ap.add_argument("--html",        required=True)
    ap.add_argument("--umap",  required=True)
    ap.add_argument("--l0p5", default=None)
    ap.add_argument("--flex-cluster-summary", default=None)
    ap.add_argument("--flex-leiden-joint", default=None)
    ap.add_argument("--nmp",  default=None)
    ap.add_argument("--leiden-assignments", default=None,
                    help="Pre-computed leiden_assignments.csv from 02_compartment_analysis.py. "
                         "If provided, skips neighbor graph + leiden and loads cluster columns "
                         "directly. Required to keep cascade YAML cluster numbers consistent.")
    ap.add_argument("-k", type=int, default=30)
    ap.add_argument("--res", type=float, default=0.3)
    ap.add_argument("--resolutions", default="0.1,0.3,0.5,0.7,1.0,1.5",
                    help="Leiden resolutions for joint expression plots")
    ap.add_argument("--flex-only-resolutions", default="0.1,0.3,0.5",
                    help="FLEX-only Leiden resolutions for FLEX-context section")
    ap.add_argument("--n-top", type=int, default=4)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    rk = f"leiden_{args.res}"

    all_res = sorted(set([float(r) for r in args.resolutions.split(",")]))
    if args.res not in all_res:
        all_res = sorted(all_res + [args.res])

    flex_only_res = sorted(set([float(r) for r in args.flex_only_resolutions.split(",")]))

    def p(name): return os.path.join(args.out_dir, f"{name}.png")
    def rstem(r): return str(r).replace(".", "p")

    # 1. Joint Leiden at all resolutions ─────────────────────────────────────
    print("[1/8] Loading latent + joint Leiden...", flush=True)
    obs     = pd.read_csv(args.obs,    index_col=0)
    latent  = pd.read_csv(args.latent, index_col=0).loc[obs.index]

    a = ad.AnnData(obs=obs.copy())
    a.obsm["X_concord"] = latent.values.astype(np.float32)

    if args.leiden_assignments:
        # Load pre-computed leiden from 02_compartment_analysis.py so that
        # cluster numbers here match leiden_assignments.csv exactly.
        asgn = pd.read_csv(args.leiden_assignments, index_col=0, dtype=str)
        for r in all_res:
            key = f"leiden_{r}"
            if key in asgn.columns:
                a.obs[key] = asgn[key].reindex(a.obs_names).fillna("-1")
                print(f"   joint res={r}: {a.obs[key].nunique()} clusters (pre-computed)", flush=True)
            else:
                print(f"   WARNING: {key} not in leiden_assignments.csv, skipping", flush=True)
        # Restrict all_res to what was actually loaded — avoids KeyError downstream.
        all_res = [r for r in all_res if f"leiden_{r}" in a.obs.columns]
        if args.res not in all_res:
            args.res = all_res[-1]
        rk = f"leiden_{args.res}"
    else:
        sc.pp.neighbors(a, n_neighbors=args.k, use_rep="X_concord", metric="cosine")
        for r in all_res:
            key = f"leiden_{r}"
            sc.tl.leiden(a, resolution=r, key_added=key, random_state=42)
            a.obs[key] = a.obs[key].astype(str)
            print(f"   joint res={r}: {a.obs[key].nunique()} clusters", flush=True)

    primary_clusters = sorted(a.obs[rk].unique(), key=lambda x: int(x))

    # L0.5 labels
    if args.l0p5:
        l0p5_df = pd.read_csv(args.l0p5, usecols=["cell_id", "l0p5_final"]).set_index("cell_id")
        a.obs["l0p5"] = l0p5_df["l0p5_final"].reindex(a.obs_names)

    if args.flex_cluster_summary and args.flex_leiden_joint:
        fl_j = pd.read_csv(args.flex_leiden_joint, index_col="cell_id")
        fl_sum_df = pd.read_csv(args.flex_cluster_summary)
        fl_sum_raw = fl_sum_df.set_index("cluster")["top_l0p5"]
        fl_sum = fl_sum_raw.squeeze() if isinstance(fl_sum_raw, pd.DataFrame) else fl_sum_raw
        leiden_col = fl_j.columns[0]
        flex_l0p5 = fl_j[leiden_col].map(fl_sum)
        flex_mask_idx = a.obs_names[a.obs["platform"] == "flex"]
        if "l0p5" not in a.obs.columns:
            a.obs["l0p5"] = pd.NA
        a.obs.loc[flex_mask_idx, "l0p5"] = flex_l0p5.reindex(flex_mask_idx).values
        print(f"   FLEX l0p5: {flex_l0p5.notna().sum():,}/{len(flex_l0p5):,}", flush=True)

    if args.nmp:
        nmp_df = pd.read_csv(args.nmp, usecols=["cell_id", "nmp_nuclear", "nmp_cyto"]).set_index("cell_id")
        a.obs["nmp_nuclear"] = pd.to_numeric(nmp_df["nmp_nuclear"].reindex(a.obs_names), errors="coerce")
        a.obs["nmp_cyto"]    = pd.to_numeric(nmp_df["nmp_cyto"].reindex(a.obs_names), errors="coerce")

    a.obs["position"] = [parse_position(cid, pid)
                         for cid, pid in zip(a.obs_names, a.obs["patient_id"])]

    # 2. FLEX-only Leiden (on joint latent, FLEX cells only) ──────────────────
    print("[2/8] FLEX-only Leiden clustering...", flush=True)
    flex_mask_bool = (a.obs["platform"] == "flex").values
    flex_ids  = a.obs_names[flex_mask_bool]
    xen_ids   = a.obs_names[~flex_mask_bool]

    fl_a = ad.AnnData(obs=a.obs.loc[flex_ids].copy())
    fl_a.obsm["X_concord"] = latent.loc[flex_ids].values.astype(np.float32)
    sc.pp.neighbors(fl_a, n_neighbors=min(args.k, len(flex_ids) - 1),
                    use_rep="X_concord", metric="cosine")

    for r in flex_only_res:
        key = f"fl_leiden_{r}"
        sc.tl.leiden(fl_a, resolution=r, key_added=key, random_state=42)
        fl_a.obs[key] = fl_a.obs[key].astype(str)
        a.obs.loc[flex_ids, key] = fl_a.obs[key].values
        a.obs[key] = a.obs[key].fillna("Xenium")
        print(f"   FLEX-only res={r}: {fl_a.obs[key].nunique()} clusters", flush=True)

    # 3. UMAP panels ──────────────────────────────────────────────────────────
    print("[3/8] Rendering UMAPs...", flush=True)
    umap_df = pd.read_csv(args.umap, index_col=0).reindex(a.obs_names)
    ux, uy = umap_df["UMAP1"].values, umap_df["UMAP2"].values

    has_l0p5 = "l0p5" in a.obs.columns
    meta_paths = []
    for key, title, centroids in [
            ("platform", "Platform", False),
            *([("l0p5", "L0.5", True)] if has_l0p5 else []),
            ("patient_id", "Patient", False),
            ("position",   "Position", True),
    ]:
        if key in a.obs.columns:
            meta_paths.append(umap_cat_panel(ux, uy, a.obs, key, title, args.out_dir,
                                             label_centroids=centroids))

    pat_pos_path = patient_position_grid(ux, uy, a.obs, args.out_dir)

    leiden_paths = []
    for r in all_res:
        key = f"leiden_{r}"
        leiden_paths.append(umap_cat_panel(ux, uy, a.obs, key, f"Leiden res={r}",
                                           args.out_dir, label_centroids=True))

    nmp_paths = []
    if args.nmp:
        for col, title, cmap in [("nmp_nuclear", "NMP nuclear", "Reds"),
                                  ("nmp_cyto",    "NMP cyto",    "Blues")]:
            nmp_paths.append(
                umap_cont_panel(ux, uy, a.obs[col].values.astype(float),
                                title, args.out_dir, col, cmap=cmap, vmin=0))

    # FLEX-only cluster UMAPs (Xenium grayed)
    flex_cl_umap_paths = {}
    for r in flex_only_res:
        fl_key = f"fl_leiden_{r}"
        fl_obs_sub = fl_a.obs.copy()
        cats = pd.Categorical(fl_obs_sub[fl_key])
        stem = rstem(r)
        path = flex_cluster_umap(ux, uy, flex_mask_bool, cats.categories,
                                 cats.codes, f"FLEX-only Leiden res={r}",
                                 args.out_dir, stem)
        flex_cl_umap_paths[r] = path

    # 4–6. Load count matrices ─────────────────────────────────────────────────
    print("[4/8] Loading count matrices...", flush=True)
    nuc_raw = load_mtx(
        os.path.join(args.nuc_bundle, "xenium_nuclear_counts.mtx.gz"),
        os.path.join(args.nuc_bundle, "xenium_cells.tsv"),
        os.path.join(args.nuc_bundle, "xenium_genes.tsv"),
    )
    common_nuc = nuc_raw.obs_names.intersection(xen_ids)
    nuc = nuc_raw[common_nuc].copy()
    sc.pp.normalize_total(nuc, target_sum=1e4); sc.pp.log1p(nuc)

    cyto_raw = load_mtx(
        os.path.join(args.cyto_bundle, "xenium_cytoplasmic_counts.mtx.gz"),
        os.path.join(args.cyto_bundle, "xenium_cells.tsv"),
        os.path.join(args.cyto_bundle, "xenium_genes.tsv"),
    )
    common_cyto = cyto_raw.obs_names.intersection(xen_ids)
    cyto = cyto_raw[common_cyto].copy()
    sc.pp.normalize_total(cyto, target_sum=1e4); sc.pp.log1p(cyto)

    flex_raw = load_flex(args.flex_dir, args.compartment, args.flex_qc)
    common_flex = flex_raw.obs_names.intersection(flex_ids)
    flex = flex_raw[common_flex].copy()
    sc.pp.normalize_total(flex, target_sum=1e4); sc.pp.log1p(flex)

    print(f"   nuc={nuc.shape}  cyto={cyto.shape}  flex={flex.shape}", flush=True)

    panel_set = set(nuc.var_names)
    nuc_umap  = umap_df.reindex(nuc.obs_names)
    flex_umap = umap_df.reindex(flex.obs_names)
    cyto_umap = umap_df.reindex(cyto.obs_names)

    # 7. Joint per-resolution markers — parallel Wilcoxon, serial dotplots ────
    n_workers = max(1, min(len(all_res) + len(flex_only_res), os.cpu_count() or 4) - 1)
    print(f"[5/8] Joint markers ({len(all_res)} resolutions, {n_workers} parallel workers)...",
          flush=True)

    # Assign cluster labels onto count objects for each resolution
    for r in all_res:
        key = f"leiden_{r}"
        nuc.obs[key]  = a.obs.loc[common_nuc, key].values
        cyto.obs[key] = a.obs.loc[common_cyto, key].values
        flex.obs[key] = a.obs.loc[common_flex, key].values

    for r in flex_only_res:
        fl_key = f"fl_leiden_{r}"
        if fl_key in a.obs.columns:
            flex.obs[fl_key] = a.obs.loc[common_flex, fl_key].values

    # Build all DE tasks: {(r, assay): task_tuple}
    de_tasks = {}
    valid_res = []
    for r in all_res:
        key = f"leiden_{r}"
        if a.obs[key].nunique() < 2:
            continue
        valid_res.append(r)
        de_tasks[(r, "nuc")]  = _adata_to_task(nuc,  nuc.obs[key].values,  args.n_top * 4)
        de_tasks[(r, "cyto")] = _adata_to_task(cyto, cyto.obs[key].values, args.n_top * 4)
        de_tasks[(r, "flex")] = _adata_to_task(flex, flex.obs[key].values, args.n_top * 4)

    for r in flex_only_res:
        fl_key = f"fl_leiden_{r}"
        if fl_key in flex.obs.columns and flex.obs[fl_key].nunique() >= 2:
            de_tasks[(r, "flex_only")] = _adata_to_task(flex, flex.obs[fl_key].values, args.n_top * 4)

    print(f"   Submitting {len(de_tasks)} DE tasks...", flush=True)
    de_results = run_markers_parallel(de_tasks, n_workers)
    print(f"   All DE tasks complete.", flush=True)

    # Save marker CSVs
    per_res = {}
    for r in all_res:
        key  = f"leiden_{r}"
        stem = rstem(r)
        if a.obs[key].nunique() < 2:
            per_res[r] = {"skipped": True}
            continue

        nuc_mdf  = de_results[(r, "nuc")]
        cyto_mdf = de_results[(r, "cyto")]
        flex_mdf = de_results[(r, "flex")]

        nuc_mdf.to_csv( os.path.join(args.out_dir, f"markers_nuclear_res{stem}.csv"),  index=False)
        cyto_mdf.to_csv(os.path.join(args.out_dir, f"markers_cyto_res{stem}.csv"),     index=False)
        flex_mdf.to_csv(os.path.join(args.out_dir, f"markers_flex_res{stem}.csv"),     index=False)

        nuc_genes  = top_genes(nuc_mdf,  args.n_top)
        cyto_genes = top_genes(cyto_mdf, args.n_top)
        flex_genes = top_genes(flex_mdf, args.n_top)

        fp_nuc = fp_flex = fp_cyto = None
        if r == args.res:
            fp_nuc  = feature_plot_grid(nuc_umap["UMAP1"].values,  nuc_umap["UMAP2"].values,
                                        nuc,  nuc_genes,  nuc.obs[key].values,
                                        args.out_dir, f"nuc_markers_{stem}")
            fp_flex = feature_plot_grid(flex_umap["UMAP1"].values, flex_umap["UMAP2"].values,
                                        flex, flex_genes, flex.obs[key].values,
                                        args.out_dir, f"flex_markers_{stem}")
            fp_cyto = feature_plot_grid(cyto_umap["UMAP1"].values, cyto_umap["UMAP2"].values,
                                        cyto, cyto_genes, cyto.obs[key].values,
                                        args.out_dir, f"cyto_markers_{stem}")

        flex_in_panel = [g for g in flex_genes if g in panel_set]
        dp_nn  = p(f"2x2_{stem}_nuc_x_nuc")
        dp_nf  = p(f"2x2_{stem}_nuc_x_flex")
        dp_ff  = p(f"2x2_{stem}_flex_x_flex")
        dp_fn  = p(f"2x2_{stem}_flex_x_nuc")
        dp_cc  = p(f"cyto_{stem}_x_cyto")
        dp_cn  = p(f"cyto_{stem}_x_nuc")
        dp_cf  = p(f"cyto_{stem}_x_flex")

        print(f"   dotplots res={r}...", flush=True)
        dotplot(nuc,  nuc_genes,     key, f"Nuclear markers - Xenium nuclear (res={r})", dp_nn)
        dotplot(flex, nuc_genes,     key, f"Nuclear markers - FLEX full-assay (res={r})", dp_nf)
        dotplot(flex, flex_genes,    key, f"FLEX markers - FLEX full-assay (res={r})", dp_ff)
        dotplot(nuc,  flex_in_panel, key,
                f"FLEX markers - Xenium nuclear ({len(flex_in_panel)}/{len(flex_genes)} in panel, res={r})", dp_fn)
        dotplot(cyto, cyto_genes, key, f"Cyto markers - cytoplasmic (res={r})", dp_cc)
        dotplot(nuc,  cyto_genes, key, f"Cyto markers - nuclear (res={r})", dp_cn)
        dotplot(flex, cyto_genes, key, f"Cyto markers - FLEX (res={r})", dp_cf)

        per_res[r] = {
            "key": key, "stem": stem, "skipped": False,
            "fp_nuc": fp_nuc, "fp_flex": fp_flex, "fp_cyto": fp_cyto,
            "dp_nn": dp_nn, "dp_nf": dp_nf, "dp_ff": dp_ff, "dp_fn": dp_fn,
            "dp_cc": dp_cc, "dp_cn": dp_cn, "dp_cf": dp_cf,
        }

    # 8. FLEX-context markers + coverage ──────────────────────────────────────
    print(f"[6/8] FLEX-context clusters ({len(flex_only_res)} resolutions)...", flush=True)
    per_flex_res = {}

    for r in flex_only_res:
        fl_key = f"fl_leiden_{r}"
        stem   = rstem(r)
        print(f"   FLEX res={r}...", flush=True)

        n_fl_cl = fl_a.obs[fl_key].nunique()
        if n_fl_cl < 2:
            per_flex_res[r] = {"skipped": True}
            continue

        flex_mdf_fl = de_results.get((r, "flex_only"))
        if flex_mdf_fl is None:
            per_flex_res[r] = {"skipped": True}
            continue
        flex_mdf_fl.to_csv(os.path.join(args.out_dir, f"markers_flex_only_res{stem}.csv"), index=False)
        flex_genes_fl = top_genes(flex_mdf_fl, args.n_top)

        fp_fl = feature_plot_grid(flex_umap["UMAP1"].values, flex_umap["UMAP2"].values,
                                  flex, flex_genes_fl, flex.obs[fl_key].values,
                                  args.out_dir, f"flex_only_markers_{stem}")

        dp_fl = p(f"flex_only_{stem}_dotplot")
        dotplot(flex, flex_genes_fl, fl_key, f"FLEX-only cluster markers (res={r})", dp_fl)

        # Coverage and compression vs joint primary resolution
        fl_obs_sub = a.obs.loc[flex_ids, [fl_key, rk]].dropna()
        cov_path  = flex_coverage_heatmap(fl_obs_sub, fl_key, rk, args.out_dir, stem)
        comp_path = flex_compression_bar(fl_obs_sub, fl_key, rk, args.out_dir, stem)

        per_flex_res[r] = {
            "skipped": False, "fl_key": fl_key, "stem": stem,
            "fp_fl": fp_fl, "dp_fl": dp_fl,
            "umap": flex_cl_umap_paths.get(r),
            "coverage": cov_path, "compression": comp_path,
        }

    # 9. Artifact screen (primary resolution) ─────────────────────────────────
    print("[7/8] Artifact screen...", flush=True)

    bal = a.obs.groupby([rk, "platform"]).size().unstack(fill_value=0)
    for col in ["xenium", "flex"]:
        if col not in bal.columns:
            bal[col] = 0
    bal = bal.loc[primary_clusters]
    bal["total"] = bal["xenium"] + bal["flex"]
    bal["xen_frac"] = bal["xenium"] / bal["total"]
    overall_xen = flex_mask_bool.mean()

    fig, axes = plt.subplots(1, 2, figsize=(16, 4))
    xv = np.arange(len(primary_clusters))
    ax = axes[0]
    ax.bar(xv, bal["xen_frac"], color="#4c72b0", label="Xenium")
    ax.bar(xv, 1 - bal["xen_frac"], bottom=bal["xen_frac"], color="#dd8452", label="FLEX")
    ax.axhline(1 - overall_xen, color="k", linestyle="--", lw=0.8,
               label=f"Expected Xenium ({1-overall_xen:.0%})")
    ax.set_xticks(xv); ax.set_xticklabels(primary_clusters, fontsize=7)
    ax.set_xlabel("Cluster"); ax.set_ylabel("Cell fraction")
    ax.set_title(f"Platform balance (res={args.res})"); ax.legend(fontsize=7)

    ax = axes[1]
    ax.bar(xv, bal["xenium"], color="#4c72b0", label="Xenium", alpha=0.8)
    ax.bar(xv, bal["flex"], bottom=bal["xenium"], color="#dd8452", label="FLEX", alpha=0.8)
    ax.set_xticks(xv); ax.set_xticklabels(primary_clusters, fontsize=7)
    ax.set_xlabel("Cluster"); ax.set_ylabel("N cells (log)")
    ax.set_title(f"Cell counts (res={args.res})"); ax.legend(fontsize=7); ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(p("artifact_platform_balance"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    shared = [g for g in nuc.var_names if g in flex.var_names]
    nuc_pct  = pct_nz_per_cluster(nuc,  rk, shared)
    flex_pct = pct_nz_per_cluster(flex, rk, shared)
    common_cl = nuc_pct.index.intersection(flex_pct.index)
    conc = pd.DataFrame({
        "cluster": np.repeat(common_cl, len(shared)),
        "gene":    np.tile(shared, len(common_cl)),
        "xen_pct": nuc_pct.loc[common_cl].values.ravel(),
        "flex_pct": flex_pct.loc[common_cl].values.ravel(),
    })
    conc["disc"] = np.abs(conc["xen_pct"] - conc["flex_pct"])
    conc.to_csv(os.path.join(args.out_dir, "artifact_pct_concordance.csv"), index=False)

    ncols = min(5, len(common_cl))
    nrows = (len(common_cl) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3, nrows * 2.8), squeeze=False)
    axes_flat = axes.ravel()
    for idx, cl in enumerate(common_cl):
        ax = axes_flat[idx]
        sub = conc[conc["cluster"] == cl]
        ax.scatter(sub["xen_pct"], sub["flex_pct"], s=10, alpha=0.5, linewidths=0)
        for _, row in sub.sort_values(by="disc", ascending=False).head(5).iterrows():
            ax.annotate(row["gene"], (row["xen_pct"], row["flex_pct"]),
                        fontsize=4.5, ha="center", va="bottom")
        ax.plot([0, 1], [0, 1], "k--", lw=0.5, alpha=0.4)
        ax.set_title(f"Cluster {cl}", fontsize=7)
        ax.set_xlabel("Xenium nuc pct_nz", fontsize=6)
        ax.set_ylabel("FLEX pct_nz", fontsize=6)
        ax.tick_params(labelsize=5)
    for idx in range(len(common_cl), len(axes_flat)):
        axes_flat[idx].axis("off")
    fig.suptitle(f"pct_nz concordance: Xenium nuclear vs FLEX (res={args.res})", fontsize=8)
    fig.tight_layout()
    fig.savefig(p("artifact_pct_concordance"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    nuc_pct2  = pct_nz_per_cluster(nuc,  rk, list(nuc.var_names))
    cyto_pct2 = pct_nz_per_cluster(cyto, rk, list(cyto.var_names))
    ncols2 = min(5, len(primary_clusters))
    nrows2 = (len(primary_clusters) + ncols2 - 1) // ncols2
    fig, axes = plt.subplots(nrows2, ncols2, figsize=(ncols2 * 3, nrows2 * 2.8), squeeze=False)
    axes_flat = axes.ravel()
    for idx, cl in enumerate(primary_clusters):
        ax = axes_flat[idx]
        if cl not in nuc_pct2.index or cl not in cyto_pct2.index:
            ax.axis("off"); continue
        xv2, yv2 = nuc_pct2.loc[cl], cyto_pct2.loc[cl]
        ax.scatter(xv2, yv2, s=10, alpha=0.5, linewidths=0)
        disc_s = pd.Series(np.abs(xv2.values - yv2.values), index=xv2.index).nlargest(5)
        for g in disc_s.index:
            ax.annotate(g, (xv2[g], yv2[g]), fontsize=4.5, ha="center", va="bottom")
        ax.plot([0, 1], [0, 1], "k--", lw=0.5, alpha=0.4)
        ax.set_title(f"Cluster {cl}", fontsize=7)
        ax.set_xlabel("Nuclear pct_nz", fontsize=6)
        ax.set_ylabel("Cyto pct_nz", fontsize=6)
        ax.tick_params(labelsize=5)
    for idx in range(len(primary_clusters), len(axes_flat)):
        axes_flat[idx].axis("off")
    fig.suptitle(f"Nuclear vs cytoplasmic pct_nz (res={args.res})", fontsize=8)
    fig.tight_layout()
    fig.savefig(p("artifact_nuc_vs_cyto"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # 10. Render HTML ──────────────────────────────────────────────────────────
    print("[8/8] Rendering HTML...", flush=True)

    CSS = """
body { font-family: Arial, sans-serif; font-size: 13px; margin: 20px; background: #f5f5f5; }
h1   { color: #2c3e50; margin-bottom: 4px; }
h3   { color: #2c3e50; border-bottom: 1px solid #ddd; padding-bottom: 4px; }
.section { background: white; border: 1px solid #ddd; border-radius: 4px;
           padding: 16px; margin-bottom: 24px; }
.note { color: #666; font-style: italic; font-size: 11px; margin: 0 0 8px; }
.row2 { display: flex; gap: 8px; flex-wrap: wrap; }
.cell2 { flex: 1 1 45%; min-width: 280px; }
.row3 { display: flex; gap: 8px; flex-wrap: wrap; }
.cell3 { flex: 1 1 30%; min-width: 240px; }
.missing { color: #aaa; font-style: italic; font-size: 11px; }
img { max-width: 100%; }
.tab-bar { display: flex; background: #34495e; border-radius: 4px 4px 0 0;
           flex-wrap: wrap; margin-bottom: 0; }
.tab-btn { padding: 7px 18px; cursor: pointer; color: #ccc; border: none;
           background: transparent; font-size: 13px; }
.tab-btn.active { background: white; color: #2c3e50; font-weight: bold; }
.tab-pane { display: none; padding: 12px 0 0 0; }
.tab-pane.active { display: block; }
    """

    JS = """
function showTab(ns, id) {
  document.querySelectorAll('[data-ns="' + ns + '"].tab-btn').forEach(function(b) {
    b.classList.remove('active');
  });
  document.querySelectorAll('[data-ns="' + ns + '"].tab-pane').forEach(function(c) {
    c.classList.remove('active');
  });
  document.getElementById(ns + '-btn-' + id).classList.add('active');
  document.getElementById(ns + '-pane-' + id).classList.add('active');
}
    """

    meta_imgs   = "".join('<div class="cell3">' + b64(pp) + '</div>' for pp in meta_paths)
    leiden_imgs = "".join('<div class="cell3">' + b64(pp) + '</div>' for pp in leiden_paths)
    nmp_imgs    = "".join('<div class="cell3">' + b64(pp) + '</div>' for pp in nmp_paths)

    nmp_btn_str = (
        '<button class="tab-btn" data-ns="emb" id="emb-btn-nmp" '
        'onclick="showTab(\'emb\',\'nmp\')">NMP</button>'
    ) if nmp_paths else ""

    nmp_pane_str = (
        '<div class="tab-pane" data-ns="emb" id="emb-pane-nmp">'
        '<p class="note">NMP nuclear (Reds) and NMP cyto (Blues). FLEX cells are NaN.</p>'
        '<div class="row3">' + nmp_imgs + '</div></div>'
    ) if nmp_paths else ""

    resolutions_str = args.resolutions
    k_str = str(args.k)
    res_str = str(args.res)

    embed_tabs = (
        '<div class="tab-bar">'
        '<button class="tab-btn active" data-ns="emb" id="emb-btn-meta" '
        'onclick="showTab(\'emb\',\'meta\')">Metadata</button>'
        '<button class="tab-btn" data-ns="emb" id="emb-btn-leiden" '
        'onclick="showTab(\'emb\',\'leiden\')">Leiden</button>'
        '<button class="tab-btn" data-ns="emb" id="emb-btn-patient" '
        'onclick="showTab(\'emb\',\'patient\')">Patient \u00d7 Position</button>'
        + nmp_btn_str +
        '</div>'
        '<div class="tab-pane active" data-ns="emb" id="emb-pane-meta">'
        '<p class="note">Platform, L0.5, patient, and position coloring.</p>'
        '<div class="row3">' + meta_imgs + '</div></div>'
        '<div class="tab-pane" data-ns="emb" id="emb-pane-leiden">'
        '<p class="note">Leiden at resolutions: ' + resolutions_str + '. k=' + k_str + ', cosine.</p>'
        '<div class="row3">' + leiden_imgs + '</div></div>'
        '<div class="tab-pane" data-ns="emb" id="emb-pane-patient">'
        '<p class="note">One panel per patient. Cells colored by anatomical position.</p>'
        + b64(pat_pos_path) +
        '</div>'
        + nmp_pane_str
    )
    embed_section = '<div class="section"><h3>Embedding</h3>' + embed_tabs + '</div>'

    # Feature plots (primary resolution)
    pr = per_res.get(args.res, {})
    fp_section = ""
    if not pr.get("skipped"):
        fp_note_nuc = "Top " + str(args.n_top) + " nuclear markers per cluster, UMAP (Xenium, res=" + res_str + ")."
        fp_note_fl  = "Top " + str(args.n_top) + " FLEX markers per cluster, UMAP (FLEX cells, res=" + res_str + ")."
        fp_note_cy  = "Top " + str(args.n_top) + " cyto markers per cluster, UMAP (Xenium, res=" + res_str + ")."
        fp_section = "".join([
            section("Feature plots \u2014 Xenium nuclear top markers",
                    b64(pr.get("fp_nuc")) if pr.get("fp_nuc") else "<p class='missing'>not generated</p>",
                    note=fp_note_nuc),
            section("Feature plots \u2014 FLEX top markers",
                    b64(pr.get("fp_flex")) if pr.get("fp_flex") else "<p class='missing'>not generated</p>",
                    note=fp_note_fl),
            section("Feature plots \u2014 Cytoplasm top markers",
                    b64(pr.get("fp_cyto")) if pr.get("fp_cyto") else "<p class='missing'>not generated</p>",
                    note=fp_note_cy),
        ])

    # 2x2 tabbed section
    tab_2x2 = {}
    for r in all_res:
        rr = per_res.get(r, {})
        if rr.get("skipped"):
            tab_2x2[r] = "<p class='missing'>skipped</p>"
        else:
            tab_2x2[r] = (row2(b64(rr.get("dp_nn")), b64(rr.get("dp_nf"))) +
                          row2(b64(rr.get("dp_ff")), b64(rr.get("dp_fn"))))

    dotplot_2x2_section = build_tab_section(
        "d2x2",
        "2\u00d72: Marker source \u00d7 Expression platform",
        ("Row 1: nuclear-derived markers (left: Xenium nuclear; right: FLEX). "
         "Row 2: FLEX-derived markers (left: FLEX; right: Xenium nuclear). Tab = resolution."),
        tab_2x2, all_res, args.res
    )

    # Cyto tabbed section
    tab_cyto = {}
    for r in all_res:
        rr = per_res.get(r, {})
        if rr.get("skipped"):
            tab_cyto[r] = "<p class='missing'>skipped</p>"
        else:
            tab_cyto[r] = row3(b64(rr.get("dp_cc")), b64(rr.get("dp_cn")), b64(rr.get("dp_cf")))

    cyto_section = build_tab_section(
        "dcyto", "Cytoplasm marker screen",
        "Markers from cytoplasmic expression. Left: cyto; center: nuclear; right: FLEX.",
        tab_cyto, all_res, args.res
    )

    # FLEX-context section
    tab_flex = {}
    primary_flex_res = flex_only_res[0] if flex_only_res else None
    for r in flex_only_res:
        rfr = per_flex_res.get(r, {})
        if rfr.get("skipped"):
            tab_flex[r] = "<p class='missing'>skipped</p>"
            continue
        umap_img    = b64(rfr.get("umap"))
        cov_img     = b64(rfr.get("coverage"))
        comp_img    = b64(rfr.get("compression"))
        dp_img      = b64(rfr.get("dp_fl"))
        fp_img      = b64(rfr.get("fp_fl")) if rfr.get("fp_fl") else "<p class='missing'>not generated</p>"
        tab_flex[r] = (
            '<p class="note">FLEX-only Leiden clusters (joint latent, FLEX subset) at res=' + str(r) + '.</p>'
            '<div class="row3">' + umap_img + cov_img + comp_img + '</div>'
            '<h4 style="color:#555;margin:12px 0 4px">FLEX cluster markers</h4>'
            + fp_img + dp_img
        )

    if flex_only_res and tab_flex:
        flex_context_section = build_tab_section(
            "dflx",
            "FLEX-context clusters",
            ("FLEX cells clustered independently in their joint-latent subspace. "
             "Coverage heatmap: fraction of each FLEX cluster's cells per joint Leiden cluster (row-normalized). "
             "Compression bar: how many joint clusters each FLEX cluster spans (>=5%). "
             "Markers from FLEX counts."),
            tab_flex, flex_only_res, primary_flex_res
        )
    else:
        flex_context_section = ""

    artifact_html = "".join([
        section("Artifact: platform balance per cluster",
                b64(p("artifact_platform_balance")),
                note="Left: Xenium vs FLEX fraction (primary res). Right: absolute cell counts."),
        section("Artifact: pct_nz concordance \u2014 Xenium nuclear vs FLEX",
                b64(p("artifact_pct_concordance")),
                note="Each point = shared panel gene. Above diagonal: higher in FLEX."),
        section("Artifact: nuclear vs cytoplasmic pct_nz",
                b64(p("artifact_nuc_vs_cyto")),
                note="Above diagonal: more cytoplasmic than nuclear."),
    ])

    body = (embed_section + fp_section + dotplot_2x2_section +
            cyto_section + flex_context_section + artifact_html)

    html = (
        '<!DOCTYPE html>\n<html><head><meta charset="utf-8">\n'
        '<title>' + args.compartment + ' \u2014 Xenium \u00d7 FLEX Marker Comparison</title>\n'
        '<style>' + CSS + '</style>\n'
        '</head><body>\n'
        '<h1>' + args.compartment + ' \u2014 Xenium \u00d7 FLEX Marker Comparison</h1>\n'
        '<p style="color:#666;font-size:11px;">Primary res=' + res_str + ' \u00b7 k=' + k_str +
        ' \u00b7 top ' + str(args.n_top) + ' markers/cluster'
        ' \u00b7 joint resolutions: ' + resolutions_str +
        ' \u00b7 FLEX-only resolutions: ' + args.flex_only_resolutions + '</p>\n'
        + body +
        '\n<script>' + JS + '</script>\n'
        '</body></html>'
    )

    os.makedirs(os.path.dirname(os.path.abspath(args.html)), exist_ok=True)
    with open(args.html, "w") as f:
        f.write(html)
    print(f"Written: {args.html}", flush=True)


if __name__ == "__main__":
    main()
