"""Per-compartment UMAP panels colored by L0.5 type.

Produces one panel per compartment (Epithelial, Stromal, Immune) plus a
3-panel overview, all colored by l0p5 within the joint_v6 UMAP space.
Cells from other compartments are shown as light gray background.

Panels produced (PDF + PNG):
  umap_epithelial_l0p5   — Epithelial cells colored by L0.5; rest gray
  umap_stromal_l0p5      — Stromal cells colored by L0.5; rest gray
  umap_immune_l0p5       — Immune cells colored by L0.5; rest gray
  umap_compartments_overview — 3-panel row summary
"""
import argparse, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

GRAY = "#cccccc"

COMPARTMENT_ORDER = ["Epithelial", "Stromal", "Immune"]


def scatter_bg_fg(ax, x_all, y_all, x_fg, y_fg, cat_fg, title, palette,
                  s=0.8, alpha_bg=0.08, alpha_fg=0.35, xlim=None, ylim=None, seed=0):
    """Plot all cells in gray, then foreground cells colored by category."""
    bg_order = np.random.default_rng(seed).permutation(len(x_all))
    ax.scatter(x_all[bg_order], y_all[bg_order], c=GRAY, s=s,
               alpha=alpha_bg, rasterized=True, linewidths=0, edgecolors="none")

    cats = pd.Categorical(cat_fg)
    codes = cats.codes
    names = cats.categories
    fg_order = np.random.default_rng(seed + 1).permutation(len(x_fg))
    ax.scatter(x_fg[fg_order], y_fg[fg_order],
               c=[palette[c] for c in codes[fg_order]],
               s=s, alpha=alpha_fg, rasterized=True, linewidths=0, edgecolors="none")

    if xlim: ax.set_xlim(xlim)
    if ylim: ax.set_ylim(ylim)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, fontsize=8)
    for sp in ax.spines.values(): sp.set_linewidth(0.25)

    handles = [Line2D([0], [0], marker="o", color="w",
                      markerfacecolor=palette[i], markersize=4, label=names[i])
               for i in range(len(names))]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5),
              fontsize=5, frameon=False, markerscale=1.5)


def save(fig, path_base):
    fig.tight_layout()
    fig.savefig(path_base + ".pdf", dpi=300, bbox_inches="tight")
    fig.savefig(path_base + ".png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  {os.path.basename(path_base)}.pdf", flush=True)


def build_palette(categories):
    cmap = plt.get_cmap("tab20", len(categories))
    return [cmap(i) for i in range(len(categories))]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--umap", required=True)
    ap.add_argument("--obs", required=True)
    ap.add_argument("--l0p5", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print("Loading...", flush=True)
    u = pd.read_csv(args.umap).set_index("cell_id")
    obs = pd.read_csv(args.obs).set_index("cell_id").reindex(u.index)
    l0 = pd.read_csv(args.l0p5, usecols=["cell_id", "l0p5_final"]).set_index("cell_id")
    obs["l0p5"] = l0["l0p5_final"].reindex(obs.index)

    x = u["UMAP1"].to_numpy(dtype=float)
    y = u["UMAP2"].to_numpy(dtype=float)
    print(f"  {len(x):,} cells total", flush=True)

    xlim = tuple(np.percentile(x, [0.1, 99.9]))
    ylim = tuple(np.percentile(y, [0.1, 99.9]))
    pad = lambda lim, r=0.03: (lim[0] - r*(lim[1]-lim[0]), lim[1] + r*(lim[1]-lim[0]))
    xlim = pad(xlim); ylim = pad(ylim)
    # Build a shared palette across all L0.5 types so colors are consistent across panels
    all_types = sorted(obs["l0p5"].dropna().unique())
    palette_map = {t: c for t, c in zip(all_types, build_palette(all_types))}

    # Per-compartment panels
    for comp in COMPARTMENT_ORDER:
        mask = (obs["compartment"] == comp).values & obs["l0p5"].notna().values
        x_fg = x[mask]; y_fg = y[mask]
        cat_fg = obs["l0p5"].values[mask]
        types_in_comp = sorted(pd.Categorical(cat_fg).categories)
        palette = [palette_map[t] for t in types_in_comp]

        n_fg = mask.sum()
        print(f"  {comp}: {n_fg:,} cells ({len(types_in_comp)} L0.5 types)", flush=True)

        fig, ax = plt.subplots(figsize=(8, 5.5))
        scatter_bg_fg(ax, x, y, x_fg, y_fg,
                      pd.Categorical(cat_fg, categories=types_in_comp),
                      f"{comp} (n={n_fg:,})", palette,
                      xlim=xlim, ylim=ylim)
        save(fig, os.path.join(args.out_dir, f"umap_{comp.lower()}_l0p5"))

    # 3-panel overview
    fig, axes = plt.subplots(1, 3, figsize=(22, 5.5))
    for ax, comp in zip(axes, COMPARTMENT_ORDER):
        mask = (obs["compartment"] == comp).values & obs["l0p5"].notna().values
        x_fg = x[mask]; y_fg = y[mask]
        cat_fg = obs["l0p5"].values[mask]
        types_in_comp = sorted(pd.Categorical(cat_fg).categories)
        palette = [palette_map[t] for t in types_in_comp]
        scatter_bg_fg(ax, x, y, x_fg, y_fg,
                      pd.Categorical(cat_fg, categories=types_in_comp),
                      f"{comp} (n={mask.sum():,})", palette,
                      xlim=xlim, ylim=ylim)
    save(fig, os.path.join(args.out_dir, "umap_compartments_overview"))

    print("Done.", flush=True)


if __name__ == "__main__":
    main()
