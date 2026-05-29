"""Render P3 UOQ vs. rest UMAP panel on joint embedding.

Excludes Tumor and Ipsilateral cells from both categories.
Produces two subpanels: all cells, Xenium-only.
Pattern: scatter_cat / save from 04_render_umap_panels_v5.py
"""
import argparse, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

EXCLUDE = ("UCI604_Tumor_xenium_1", "UCI604_Ipsilateral_xenium_1")


def is_excluded(cell_id):
    return any(cell_id.startswith(p) for p in EXCLUDE)


def scatter_cat(ax, x, y, cat, title, palette, s=0.8, alpha=0.3,
                xlim=None, ylim=None, seed=0):
    cats = pd.Categorical(cat)
    codes = cats.codes
    names = cats.categories
    order = np.random.default_rng(seed).permutation(len(x))
    ax.scatter(x[order], y[order], c=[palette[c] for c in codes[order]],
               s=s, alpha=alpha, rasterized=True, linewidths=0, edgecolors="none")
    if xlim: ax.set_xlim(xlim)
    if ylim: ax.set_ylim(ylim)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, fontsize=8)
    for sp in ax.spines.values(): sp.set_linewidth(0.25)
    handles = [Line2D([0], [0], marker="o", color="w",
                      markerfacecolor=palette[i], markersize=4, label=names[i])
               for i in range(len(names))]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5),
              fontsize=6, frameon=False, markerscale=1.5)


def save(fig, path_base):
    fig.tight_layout()
    fig.savefig(path_base + ".pdf", dpi=300, bbox_inches="tight")
    fig.savefig(path_base + ".png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {os.path.basename(path_base)}.png", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--umap", required=True)
    ap.add_argument("--obs", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print("Loading...", flush=True)
    u = pd.read_csv(args.umap).set_index("cell_id")
    obs = pd.read_csv(args.obs).set_index("cell_id").reindex(u.index)

    excluded = pd.Series([is_excluded(c) for c in u.index], index=u.index)
    mask = ~excluded

    x = u.loc[mask, "UMAP1"].values
    y = u.loc[mask, "UMAP2"].values
    cell_ids = u.index[mask]
    platform = obs.loc[mask, "platform"].values

    is_p3_uoq = np.array(
        [("_P3_" in c) and ("_UOQ" in c) for c in cell_ids]
    )
    label = np.where(is_p3_uoq, "P3 UOQ", "rest")

    n_total = mask.sum()
    n_uoq = is_p3_uoq.sum()
    n_excluded = excluded.sum()
    print(f"  {n_excluded:,} excluded (Tumor/Ipsilateral)", flush=True)
    print(f"  {n_total:,} cells plotted  ({n_uoq:,} P3 UOQ, {n_total-n_uoq:,} rest)", flush=True)

    x0, x1 = float(np.percentile(x, 0.1)), float(np.percentile(x, 99.9))
    y0, y1 = float(np.percentile(y, 0.1)), float(np.percentile(y, 99.9))
    pad = 0.03
    xlim = (x0 - pad*(x1-x0), x1 + pad*(x1-x0))
    ylim = (y0 - pad*(y1-y0), y1 + pad*(y1-y0))

    palette = {"rest": "#cccccc", "P3 UOQ": "#ff7f0e"}
    pal_list = [palette["P3 UOQ"], palette["rest"]]  # alphabetical: P3 UOQ, rest

    # Panel 1: all cells
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    scatter_cat(axes[0], x, y, label,
                f"P3 UOQ vs rest — all cells (n={n_total:,})",
                palette=pal_list, xlim=xlim, ylim=ylim)

    # Panel 2: Xenium only
    xen = platform == "xenium"
    scatter_cat(axes[1], x[xen], y[xen], label[xen],
                f"P3 UOQ vs rest — Xenium only (n={xen.sum():,})",
                palette=pal_list, xlim=xlim, ylim=ylim)

    save(fig, os.path.join(args.out_dir, "umap_p3_uoq"))
    print("Done.", flush=True)


if __name__ == "__main__":
    main()
