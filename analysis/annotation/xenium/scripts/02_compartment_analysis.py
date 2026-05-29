"""Per-compartment analysis: Leiden clustering, marker calling, diagnostic UMAPs.

Pattern refs:
  Leiden + neighbors: 03_joint_umap.py (RAPIDS neighbors, scanpy CPU UMAP)
  Marker calling: scanpy rank_genes_groups (Wilcoxon)
  Panels: 04_render_umap_panels_v5.py (scatter_cat / scatter_cont / save)

Outputs in --out-dir (panels/):
  umap_platform, umap_patient, umap_l0p5
  umap_leiden_{res} x6 (0.1, 0.3, 0.5, 0.7, 1.0, 1.5)
  umap_n_counts, umap_n_genes
  umap_nmp_nuclear, umap_nmp_cyto, umap_nmp_exclude
  umap_knn_confidence
  umap_flex_leiden_0.3  (FLEX cells colored by Phase 3 Leiden cluster)
  dotplot_markers_{primary_res}
  markers_{primary_res}.csv
"""
import argparse, gzip, io, os
import numpy as np
import pandas as pd
import anndata as ad
import rapids_singlecell as rsc
import scanpy as sc
import scipy.io, scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

RESOLUTIONS = [0.1, 0.3, 0.5, 0.7, 1.0, 1.5]

COMP_L0P5 = {
    "Epithelial": ["LASP", "LHS", "BMYO"],
    "Stromal":    ["FB", "PV", "LE", "VE", "Adipocyte"],
    "Immune":     ["T_cell", "B_cell", "Plasma", "NK", "Myeloid", "Mast", "Neutro"],
}


# ── Rendering helpers ──────────────────────────────────────────────────────────

def scatter_cat(ax, x, y, cat, title, palette=None, s=1.0, alpha=0.3, seed=0):
    cats = pd.Categorical(cat)
    codes = cats.codes
    names = cats.categories
    if palette is None:
        cmap = plt.get_cmap("tab20", len(names))
        palette = [cmap(i) for i in range(len(names))]
    order = np.random.default_rng(seed).permutation(len(x))
    ax.scatter(x[order], y[order], c=[palette[c] for c in codes[order]],
               s=s, alpha=alpha, rasterized=True, linewidths=0, edgecolors="none")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, fontsize=8)
    for sp in ax.spines.values(): sp.set_linewidth(0.25)
    handles = [Line2D([0], [0], marker="o", color="w",
                      markerfacecolor=palette[i], markersize=4, label=names[i])
               for i in range(len(names))]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5),
              fontsize=5, frameon=False, markerscale=1.5)


def scatter_cont(ax, x, y, vals, title, cmap="viridis", vmin=None, vmax=None,
                 s=1.0, alpha=0.3, seed=0):
    m = np.isfinite(vals)
    order = np.random.default_rng(seed).permutation(len(x))
    sc_obj = ax.scatter(x[order][m[order]], y[order][m[order]],
                        c=vals[order][m[order]], cmap=cmap,
                        vmin=vmin, vmax=vmax,
                        s=s, alpha=alpha, rasterized=True, linewidths=0)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, fontsize=8)
    for sp in ax.spines.values(): sp.set_linewidth(0.25)
    plt.colorbar(sc_obj, ax=ax, shrink=0.6, pad=0.02)


def save(fig, path_base):
    fig.tight_layout()
    fig.savefig(path_base + ".pdf", dpi=200, bbox_inches="tight")
    fig.savefig(path_base + ".png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  {os.path.basename(path_base)}.png", flush=True)


# ── Data loading ───────────────────────────────────────────────────────────────

def load_counts_subset(bundle_dir, keep_ids):
    """Load pooled nuclear MTX, return subset AnnData (normalized)."""
    mtx_path = os.path.join(bundle_dir, "xenium_nuclear_counts.mtx.gz")
    with gzip.open(mtx_path, "rb") as gz:
        mat = scipy.io.mmread(io.BytesIO(gz.read()))
    X = sp.csr_matrix(mat)
    genes = pd.read_csv(os.path.join(bundle_dir, "xenium_genes.tsv"), header=None)[0].tolist()
    cells = pd.read_csv(os.path.join(bundle_dir, "xenium_cells.tsv"), header=None)[0].tolist()
    cell_idx = pd.Index(cells)
    keep = cell_idx.isin(keep_ids)
    X_sub = X[keep]
    obs_sub = pd.DataFrame(index=cell_idx[keep])
    a = ad.AnnData(X=X_sub, obs=obs_sub)
    a.var_names = genes
    a.obs["n_counts"] = np.array(X_sub.sum(axis=1)).flatten()
    a.obs["n_genes"]  = np.array((X_sub > 0).sum(axis=1)).flatten()
    return a


def load_flex_counts(flex_dir: str, comp: str, flex_qc_path: str) -> ad.AnnData:
    """Load full FLEX count matrix for one compartment (pass_qc cells only)."""
    cdir = os.path.join(flex_dir, comp, "scvi_n100")
    with gzip.open(os.path.join(cdir, "counts.mtx.gz"), "rb") as gz:
        mat = scipy.io.mmread(io.BytesIO(gz.read()))
    X = sp.csr_matrix(mat)
    genes = pd.read_csv(os.path.join(cdir, "genes.tsv"), header=None)[0].tolist()
    cells = pd.read_csv(os.path.join(cdir, "cells.tsv"), header=None)[0].tolist()
    obs = pd.read_csv(os.path.join(cdir, "obs.csv"), index_col=0)
    flex_qc = pd.read_csv(flex_qc_path)
    qc = flex_qc[flex_qc["compartment"] == comp].set_index("cell_id")
    keep_ids = set(qc.index[qc["pass_qc"]])
    keep = pd.Index(cells).isin(keep_ids)
    a = ad.AnnData(X=X[keep], obs=obs.loc[pd.Index(cells)[keep]].copy())
    a.var_names = genes
    print(f"  FLEX/{comp}: {a.shape[0]:,} pass_qc cells × {a.shape[1]:,} genes", flush=True)
    return a


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--compartment", required=True, choices=list(COMP_L0P5.keys()))
    ap.add_argument("--latent",       required=True)
    ap.add_argument("--obs",          required=True)
    ap.add_argument("--umap",         required=True)
    ap.add_argument("--xenium-bundle", required=True)
    ap.add_argument("--l0p5",         required=True)
    ap.add_argument("--flex-leiden",  required=True)
    ap.add_argument("--flex-dir",     required=True,
                    help="integration_intermediate dir with {comp}/scvi_n100/ subdirs")
    ap.add_argument("--flex-qc",      required=True)
    ap.add_argument("--nmp",          required=True)
    ap.add_argument("--out-dir",      required=True)
    ap.add_argument("-k", type=int, default=30)
    args = ap.parse_args()

    comp = args.compartment
    os.makedirs(args.out_dir, exist_ok=True)

    # Load latent + obs
    print(f"[{comp}] Loading latent + obs...", flush=True)
    lat = pd.read_csv(args.latent, index_col=0)
    obs = pd.read_csv(args.obs).set_index("cell_id").reindex(lat.index)
    umap_df = pd.read_csv(args.umap).set_index("cell_id").reindex(lat.index)
    print(f"  {len(lat):,} cells (platforms: {obs['platform'].value_counts().to_dict()})",
          flush=True)

    # Join metadata
    l0 = pd.read_csv(args.l0p5, usecols=["cell_id", "l0p5_final", "knn_confidence"]).set_index("cell_id")
    obs["l0p5"] = l0["l0p5_final"].reindex(obs.index)
    obs["knn_confidence"] = pd.to_numeric(l0["knn_confidence"].reindex(obs.index), errors="coerce")

    # FLEX leiden clusters — use Int64 before string to avoid "0.0" labels
    fl = pd.read_csv(args.flex_leiden).set_index("cell_id")
    for col in [f"leiden_{r}" for r in ["0.1", "0.3", "0.5", "1.0"]]:
        if col in fl.columns:
            obs[f"flex_{col}"] = fl[col].reindex(obs.index).astype("Int64").astype("string")

    nmp = pd.read_csv(args.nmp, usecols=["cell_id", "nmp_nuclear", "nmp_cyto"]).set_index("cell_id")
    obs["nmp_nuclear"] = pd.to_numeric(nmp["nmp_nuclear"].reindex(obs.index), errors="coerce")
    obs["nmp_cyto"]    = pd.to_numeric(nmp["nmp_cyto"].reindex(obs.index),    errors="coerce")
    obs["nmp_exclude"] = (nmp["nmp_nuclear"] > nmp["nmp_cyto"]).reindex(obs.index)

    # QC from count matrix (Xenium only)
    xen_ids = obs.index[obs["platform"] == "xenium"]
    print(f"  Loading counts for {len(xen_ids):,} Xenium cells...", flush=True)
    counts_a = load_counts_subset(args.xenium_bundle, set(xen_ids))
    obs.loc[counts_a.obs_names, "n_counts"] = counts_a.obs["n_counts"].values
    obs.loc[counts_a.obs_names, "n_genes"]  = counts_a.obs["n_genes"].values

    # Build AnnData for Leiden
    a = ad.AnnData(X=np.zeros((lat.shape[0], 1), dtype=np.float32), obs=obs.copy())
    a.obsm["X_concord"] = lat.values.astype(np.float32)
    a.obsm["X_umap"]    = umap_df[["UMAP1", "UMAP2"]].to_numpy(dtype=np.float32)

    # RAPIDS neighbors + scanpy CPU Leiden at each resolution
    print(f"  RAPIDS neighbors k={args.k}...", flush=True)
    rsc.pp.neighbors(a, n_neighbors=args.k, use_rep="X_concord", metric="cosine")

    for res in RESOLUTIONS:
        key = f"leiden_{res}"
        print(f"  Leiden res={res}...", flush=True)
        sc.tl.leiden(a, resolution=res, key_added=key, random_state=42)
        a.obs[key] = a.obs[key].astype(str)
        n_cl = a.obs[key].nunique()
        print(f"    {n_cl} clusters", flush=True)

    # Save leiden assignments for L1.5 label mapping
    leiden_cols = [f"leiden_{r}" for r in RESOLUTIONS]
    leiden_df = a.obs[leiden_cols].copy()
    leiden_df.index.name = "cell_id"
    leiden_df.to_csv(os.path.join(args.out_dir, "..", "leiden_assignments.csv"))
    print(f"  Saved leiden_assignments.csv", flush=True)

    primary_res = "leiden_0.3"

    xen_mask = a.obs["platform"] == "xenium"
    # Build marker AnnData from counts_a (proper var_names) + leiden labels.
    # Cannot assign counts_a.X to a_xen.X because a was built with 1 dummy variable.
    xen_obs_names = a.obs_names[xen_mask.to_numpy(dtype=bool)]
    common = counts_a.obs_names.intersection(xen_obs_names)
    counts_sub = counts_a[common].copy()
    a_xen = ad.AnnData(X=counts_sub.X, obs=counts_sub.obs.copy())
    a_xen.var_names = counts_sub.var_names
    leiden_labels = a.obs.loc[common, primary_res]
    a_xen.obs[primary_res] = leiden_labels.values
    sc.pp.normalize_total(a_xen, target_sum=1e4)
    sc.pp.log1p(a_xen)
    print(f"  {len(common):,} Xenium cells for marker calling", flush=True)
    print(f"  rank_genes_groups (Wilcoxon, {primary_res})...", flush=True)
    sc.tl.rank_genes_groups(a_xen, groupby=primary_res, method="wilcoxon",
                             n_genes=20, pts=True)
    markers_df = sc.get.rank_genes_groups_df(a_xen, group=None)
    markers_stem = f"markers_{primary_res}"
    markers_df.to_csv(os.path.join(args.out_dir, f"{markers_stem}.csv"), index=False)
    print(f"  Saved {markers_stem}.csv ({len(markers_df)} rows)", flush=True)

    # FLEX markers — full assay (~6000 genes), joint-defined clusters
    print(f"  Loading FLEX counts (full assay)...", flush=True)
    flex_counts = load_flex_counts(args.flex_dir, comp, args.flex_qc)
    flex_obs_names = a.obs_names[a.obs["platform"] == "flex"]
    flex_common = flex_counts.obs_names.intersection(flex_obs_names)
    flex_sub = flex_counts[flex_common].copy()
    a_flex = ad.AnnData(X=flex_sub.X, obs=flex_sub.obs.copy())
    a_flex.var_names = flex_sub.var_names
    flex_leiden = a.obs.loc[flex_common, primary_res]
    a_flex.obs[primary_res] = flex_leiden.values
    sc.pp.normalize_total(a_flex, target_sum=1e4)
    sc.pp.log1p(a_flex)
    print(f"  {len(flex_common):,} FLEX cells for marker calling", flush=True)
    print(f"  rank_genes_groups FLEX (Wilcoxon, {primary_res})...", flush=True)
    sc.tl.rank_genes_groups(a_flex, groupby=primary_res, method="wilcoxon",
                             n_genes=50, pts=True)
    flex_markers_df = sc.get.rank_genes_groups_df(a_flex, group=None)
    flex_stem = f"markers_flex_{primary_res}"
    flex_markers_df.to_csv(os.path.join(args.out_dir, f"{flex_stem}.csv"), index=False)
    print(f"  Saved {flex_stem}.csv ({len(flex_markers_df)} rows)", flush=True)

    # FLEX marker dotplot (top 5 per cluster)
    flex_top5 = (flex_markers_df.groupby("group")
                                .apply(lambda g: g.nlargest(5, "scores"))
                                .reset_index(drop=True))
    flex_genes = flex_top5["names"].unique().tolist()[:60]
    if flex_genes:
        sc.pl.dotplot(a_flex, var_names=flex_genes, groupby=primary_res, show=False)
        fig = plt.gcf()
        fig.savefig(os.path.join(args.out_dir, f"dotplot_{flex_stem}.pdf"),
                    dpi=200, bbox_inches="tight")
        fig.savefig(os.path.join(args.out_dir, f"dotplot_{flex_stem}.png"),
                    dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  dotplot_{flex_stem}.png", flush=True)

    # ── Render panels ──────────────────────────────────────────────────────────
    x = umap_df["UMAP1"].to_numpy(dtype=float)
    y = umap_df["UMAP2"].to_numpy(dtype=float)
    xlim_raw = (float(np.percentile(x, 0.1)), float(np.percentile(x, 99.9)))
    ylim_raw = (float(np.percentile(y, 0.1)), float(np.percentile(y, 99.9)))
    px = 0.03 * (xlim_raw[1] - xlim_raw[0])
    py = 0.03 * (ylim_raw[1] - ylim_raw[0])
    xlim = (xlim_raw[0] - px, xlim_raw[1] + px)
    ylim = (ylim_raw[0] - py, ylim_raw[1] + py)

    print("  Rendering panels...", flush=True)

    # Platform
    fig, ax = plt.subplots(figsize=(7, 5.5))
    scatter_cat(ax, x, y, obs["platform"].fillna("na"),
                f"{comp} — platform (n={len(x):,})",
                palette=["#1f77b4", "#ff7f0e"])
    ax.set_xlim(xlim); ax.set_ylim(ylim)
    save(fig, os.path.join(args.out_dir, "umap_platform"))

    # Patient
    fig, ax = plt.subplots(figsize=(7, 5.5))
    scatter_cat(ax, x, y, obs["patient_id"].fillna("na"), f"{comp} — patient")
    ax.set_xlim(xlim); ax.set_ylim(ylim)
    save(fig, os.path.join(args.out_dir, "umap_patient"))

    # L0.5
    has_l0 = obs["l0p5"].notna().to_numpy(dtype=bool)
    if has_l0.sum() > 0:
        fig, ax = plt.subplots(figsize=(8, 5.5))
        scatter_cat(ax, x[has_l0], y[has_l0],
                    obs["l0p5"].to_numpy(dtype=str)[has_l0],
                    f"{comp} — L0.5 (n={has_l0.sum():,})")
        ax.set_xlim(xlim); ax.set_ylim(ylim)
        save(fig, os.path.join(args.out_dir, "umap_l0p5"))

    # Leiden at each resolution
    for res in RESOLUTIONS:
        key = f"leiden_{res}"
        n_cl = a.obs[key].nunique()
        fig, ax = plt.subplots(figsize=(8, 5.5))
        scatter_cat(ax, x, y, a.obs[key], f"{comp} — Leiden res={res} ({n_cl} clusters)")
        ax.set_xlim(xlim); ax.set_ylim(ylim)
        save(fig, os.path.join(args.out_dir, f"umap_{key}"))

    # QC: n_counts, n_genes (Xenium only)
    xen_vals = xen_mask.to_numpy(dtype=bool)
    for col, label, cm in [("n_counts", "n_counts (nuclear)", "plasma"),
                            ("n_genes",  "n_genes (nuclear)",  "viridis")]:
        vals = obs[col].to_numpy(dtype=float)
        fig, ax = plt.subplots(figsize=(7, 5.5))
        scatter_cont(ax, x[xen_vals], y[xen_vals], vals[xen_vals],
                     f"{comp} — {label} (Xenium)", cmap=cm,
                     vmin=float(np.nanpercentile(vals[xen_vals], 5)),
                     vmax=float(np.nanpercentile(vals[xen_vals], 95)))
        ax.set_xlim(xlim); ax.set_ylim(ylim)
        save(fig, os.path.join(args.out_dir, f"umap_{col}"))

    # NMP nuclear + cytoplasmic scores (Xenium only)
    for col, label, cm in [("nmp_nuclear", "NMP nuclear score", "Reds"),
                            ("nmp_cyto",   "NMP cyto score",    "Blues")]:
        vals = obs[col].to_numpy(dtype=float)
        if np.isfinite(vals[xen_vals]).sum() > 0:
            fig, ax = plt.subplots(figsize=(7, 5.5))
            scatter_cont(ax, x[xen_vals], y[xen_vals], vals[xen_vals],
                         f"{comp} — {label} (Xenium)", cmap=cm,
                         vmin=float(np.nanpercentile(vals[xen_vals], 5)),
                         vmax=float(np.nanpercentile(vals[xen_vals], 95)))
            ax.set_xlim(xlim); ax.set_ylim(ylim)
            save(fig, os.path.join(args.out_dir, f"umap_{col}"))

    # NMP exclude flag (Xenium)
    nmp_flag = obs.loc[xen_mask, "nmp_exclude"].astype("boolean").fillna(False).astype(str).values
    fig, ax = plt.subplots(figsize=(7, 5.5))
    scatter_cat(ax, x[xen_vals], y[xen_vals], nmp_flag,
                f"{comp} — NMP exclude (Xenium)",
                palette=["#bbbbbb", "#d62728"])
    ax.set_xlim(xlim); ax.set_ylim(ylim)
    save(fig, os.path.join(args.out_dir, "umap_nmp_exclude"))

    # KNN confidence (Xenium only — FLEX cells are NaN)
    knn_vals = obs["knn_confidence"].to_numpy(dtype=float)
    if np.isfinite(knn_vals[xen_vals]).sum() > 0:
        fig, ax = plt.subplots(figsize=(7, 5.5))
        scatter_cont(ax, x[xen_vals], y[xen_vals], knn_vals[xen_vals],
                     f"{comp} — KNN confidence (Xenium)", cmap="magma",
                     vmin=0.0, vmax=1.0)
        ax.set_xlim(xlim); ax.set_ylim(ylim)
        save(fig, os.path.join(args.out_dir, "umap_knn_confidence"))

    # FLEX Leiden 0.3 (integer labels via Int64 to avoid "0.0" strings)
    flex_vals = obs["platform"] == "flex"
    if flex_vals.sum() > 0 and "flex_leiden_0.3" in obs.columns:
        flex_mask = flex_vals.to_numpy(dtype=bool)
        fig, ax = plt.subplots(figsize=(8, 5.5))
        scatter_cat(ax, x[flex_mask], y[flex_mask],
                    obs.loc[flex_vals, "flex_leiden_0.3"].fillna("na"),
                    f"{comp} — FLEX Leiden 0.3 (n={flex_vals.sum():,})")
        ax.set_xlim(xlim); ax.set_ylim(ylim)
        save(fig, os.path.join(args.out_dir, "umap_flex_leiden_0.3"))

    # Marker dot plot (top 5 genes per cluster)
    top5 = (markers_df.groupby("group")
                      .apply(lambda g: g.nlargest(5, "scores"))
                      .reset_index(drop=True))
    genes_to_plot = top5["names"].unique().tolist()[:50]
    if genes_to_plot:
        sc.pl.dotplot(a_xen, var_names=genes_to_plot, groupby=primary_res, show=False)
        fig = plt.gcf()
        fig.savefig(os.path.join(args.out_dir, f"dotplot_{markers_stem}.pdf"),
                    dpi=200, bbox_inches="tight")
        fig.savefig(os.path.join(args.out_dir, f"dotplot_{markers_stem}.png"),
                    dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  dotplot_{markers_stem}.png", flush=True)

    print(f"[{comp}] Done: all panels written to {args.out_dir}", flush=True)


if __name__ == "__main__":
    main()
