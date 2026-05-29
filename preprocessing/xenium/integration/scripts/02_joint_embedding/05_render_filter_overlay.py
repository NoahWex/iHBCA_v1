"""Render filter-overlay UMAPs: show which cells each axis catches.

Panels:
  1. passes_all (gray=pass, red=any-fail)
  2. filter status: all_pass / prox_fail / nmp_fail / both_fail
  3. just NMP fail
  4. just proximity fail
"""
import argparse, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def scatter(ax, x, y, cat, title, palette, s=1.0, alpha=0.3, xlim=None, ylim=None,
            legend=True):
    cats = pd.Categorical(cat)
    codes = cats.codes
    names = list(cats.categories)
    # palette: dict keyed by name → build list aligned to codes
    pal_list = [palette[n] for n in names]
    order = np.random.default_rng(0).permutation(len(x))
    ax.scatter(x[order], y[order], c=[pal_list[c] for c in codes[order]],
               s=s, alpha=alpha, rasterized=True, linewidths=0, edgecolors="none")
    if xlim: ax.set_xlim(xlim)
    if ylim: ax.set_ylim(ylim)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, fontsize=8)
    for sp in ax.spines.values(): sp.set_linewidth(0.25)
    if legend:
        handles = [plt.Line2D([0], [0], marker="o", color="w",
                              markerfacecolor=pal_list[i], markersize=5, label=names[i])
                   for i in range(len(names))]
        ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5),
                  fontsize=6, frameon=False, markerscale=1.5)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--umap", required=True)
    ap.add_argument("--filter-csv", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    u = pd.read_csv(args.umap).set_index("cell_id")
    f = pd.read_csv(args.filter_csv).set_index("cell_id")
    df = u.join(f, how="inner")
    print(f"xenium cells w/ UMAP: {len(df)}", flush=True)

    x = df["UMAP1"].values; y = df["UMAP2"].values
    xlim = (float(np.percentile(x, 0.1)), float(np.percentile(x, 99.9)))
    ylim = (float(np.percentile(y, 0.1)), float(np.percentile(y, 99.9)))
    px = 0.03 * (xlim[1] - xlim[0]); py = 0.03 * (ylim[1] - ylim[0])
    xlim = (xlim[0]-px, xlim[1]+px); ylim = (ylim[0]-py, ylim[1]+py)

    # Panel 1: passes_all
    pa = df["passes_all"].astype(bool)
    fig, ax = plt.subplots(figsize=(7, 5.5), dpi=300)
    scatter(ax, x, y, np.where(pa, "pass", "fail"),
            f"passes_all ({int(pa.sum()):,}/{len(pa):,} = {100*pa.mean():.1f}%)",
            palette={"pass": "#bbbbbb", "fail": "#d62728"},
            xlim=xlim, ylim=ylim)
    plt.tight_layout()
    plt.savefig(os.path.join(args.out_dir, "umap_filter_passes_all.pdf"), dpi=300)
    plt.close()

    # Panel 2: axis combo status
    prox = df["passes_proximity"].astype(bool)
    nmp = df["passes_nmp"].astype(bool)
    status = np.full(len(df), "all_pass", dtype=object)
    status[(~prox) & (~nmp)] = "both_fail"
    status[(~prox) & nmp] = "prox_fail_only"
    status[prox & (~nmp)] = "nmp_fail_only"
    fig, ax = plt.subplots(figsize=(8, 5.5), dpi=300)
    scatter(ax, x, y, status, "Per-axis filter status",
            palette={"all_pass": "#cccccc",
                     "prox_fail_only": "#1f77b4",
                     "nmp_fail_only": "#ff7f0e",
                     "both_fail": "#d62728"},
            xlim=xlim, ylim=ylim)
    plt.tight_layout()
    plt.savefig(os.path.join(args.out_dir, "umap_filter_status.pdf"), dpi=300)
    plt.close()

    # Panel 3: only NMP fail highlighted
    fig, ax = plt.subplots(figsize=(7, 5.5), dpi=300)
    scatter(ax, x, y, np.where(nmp, "nmp_pass", "nmp_fail"),
            f"Axis 2 NMP (fail = nuclear > cyto): {int((~nmp).sum()):,} cells",
            palette={"nmp_pass": "#cccccc", "nmp_fail": "#d62728"},
            xlim=xlim, ylim=ylim)
    plt.tight_layout()
    plt.savefig(os.path.join(args.out_dir, "umap_filter_nmp_only.pdf"), dpi=300)
    plt.close()

    # Panel 4: only proximity fail
    fig, ax = plt.subplots(figsize=(7, 5.5), dpi=300)
    scatter(ax, x, y, np.where(prox, "prox_pass", "prox_fail"),
            f"Axis 1 FLEX proximity (fail): {int((~prox).sum()):,} cells",
            palette={"prox_pass": "#cccccc", "prox_fail": "#1f77b4"},
            xlim=xlim, ylim=ylim)
    plt.tight_layout()
    plt.savefig(os.path.join(args.out_dir, "umap_filter_prox_only.pdf"), dpi=300)
    plt.close()
    print("Done.", flush=True)


if __name__ == "__main__":
    main()
