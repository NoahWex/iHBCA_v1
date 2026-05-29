"""Render diagnostic UMAP panels for joint_nuc_n100 winner.

Panels:
  1. platform     (flex vs xenium)
  2. compartment  (Epi / Stroma / Immune)
  3. patient
  4. l0p5_final   (13 types, xenium cells)
  5. nmp_exclude  (xenium: nmp_nuclear > nmp_cyto)
"""
import argparse, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def scatter_panel(ax, x, y, cat, title, palette=None, s=1.0, alpha=0.3, legend=True,
                  xlim=None, ylim=None):
    cats = pd.Categorical(cat)
    codes = cats.codes
    names = cats.categories
    if palette is None:
        cmap = plt.get_cmap("tab20", len(names))
        palette = [cmap(i) for i in range(len(names))]
    # random draw order so no layer dominates
    order = np.random.default_rng(0).permutation(len(x))
    ax.scatter(x[order], y[order], c=[palette[c] for c in codes[order]],
               s=s, alpha=alpha, rasterized=True, linewidths=0, edgecolors="none")
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, fontsize=8)
    for sp in ax.spines.values():
        sp.set_linewidth(0.25)
    if legend:
        handles = [plt.Line2D([0], [0], marker="o", color="w",
                              markerfacecolor=palette[i], markersize=4, label=names[i])
                   for i in range(len(names))]
        ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5),
                  fontsize=5, frameon=False, markerscale=1.5)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--umap", required=True)
    ap.add_argument("--obs", required=True)
    ap.add_argument("--l0p5", required=True)
    ap.add_argument("--nmp", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print("Loading...", flush=True)
    u = pd.read_csv(args.umap).set_index("cell_id")
    obs = pd.read_csv(args.obs).set_index("cell_id").reindex(u.index)
    l0 = pd.read_csv(args.l0p5, usecols=["cell_id", "l0p5_final"]).set_index("cell_id")
    nmp = pd.read_csv(args.nmp, usecols=["cell_id", "nmp_nuclear", "nmp_cyto"]).set_index("cell_id")
    obs["l0p5"] = l0["l0p5_final"].reindex(obs.index)
    obs["nmp_exclude"] = (nmp["nmp_nuclear"] > nmp["nmp_cyto"]).reindex(obs.index)

    x = u["UMAP1"].values; y = u["UMAP2"].values
    print(f"  {len(x):,} cells", flush=True)
    # Clip axes to bulk (ignore outliers that squash the embedding)
    xlim = (float(np.percentile(x, 0.1)), float(np.percentile(x, 99.9)))
    ylim = (float(np.percentile(y, 0.1)), float(np.percentile(y, 99.9)))
    pad_x = 0.03 * (xlim[1] - xlim[0])
    pad_y = 0.03 * (ylim[1] - ylim[0])
    xlim = (xlim[0] - pad_x, xlim[1] + pad_x)
    ylim = (ylim[0] - pad_y, ylim[1] + pad_y)
    print(f"  xlim={xlim}, ylim={ylim}", flush=True)

    # --- Panel 1: platform ---
    fig, ax = plt.subplots(figsize=(7, 5.5), dpi=300)
    scatter_panel(ax, x, y, obs["platform"].fillna("na"),
                  f"Platform (joint_nuc_n100, n={len(x):,})",
                  palette=["#1f77b4", "#ff7f0e"], xlim=xlim, ylim=ylim)
    plt.tight_layout()
    plt.savefig(os.path.join(args.out_dir, "umap_platform.pdf"), dpi=300)
    plt.savefig(os.path.join(args.out_dir, "umap_platform.png"), dpi=300)
    plt.close()
    print("  platform.pdf", flush=True)

    # --- Panel 2: compartment ---
    fig, ax = plt.subplots(figsize=(7, 5.5), dpi=300)
    scatter_panel(ax, x, y, obs["compartment"].fillna("na"),
                  "Compartment",
                  palette=["#d62728", "#2ca02c", "#9467bd", "#8c564b"],
                  xlim=xlim, ylim=ylim)
    plt.tight_layout()
    plt.savefig(os.path.join(args.out_dir, "umap_compartment.pdf"), dpi=300)
    plt.savefig(os.path.join(args.out_dir, "umap_compartment.png"), dpi=300)
    plt.close()
    print("  compartment.pdf", flush=True)

    # --- Panel 3: patient ---
    fig, ax = plt.subplots(figsize=(7, 5.5), dpi=300)
    scatter_panel(ax, x, y, obs["patient_id"].fillna("na"), "Patient",
                  xlim=xlim, ylim=ylim)
    plt.tight_layout()
    plt.savefig(os.path.join(args.out_dir, "umap_patient.pdf"), dpi=300)
    plt.savefig(os.path.join(args.out_dir, "umap_patient.png"), dpi=300)
    plt.close()
    print("  patient.pdf", flush=True)

    # --- Panel 4: l0p5 (all cells that have a label, both platforms) ---
    has_l0 = obs["l0p5"].notna().values
    fig, ax = plt.subplots(figsize=(8, 5.5), dpi=300)
    scatter_panel(ax, x[has_l0], y[has_l0], obs["l0p5"][has_l0],
                  f"L0.5 (labeled cells, n={int(has_l0.sum()):,})",
                  xlim=xlim, ylim=ylim)
    plt.tight_layout()
    plt.savefig(os.path.join(args.out_dir, "umap_l0p5.pdf"), dpi=300)
    plt.savefig(os.path.join(args.out_dir, "umap_l0p5.png"), dpi=300)
    plt.close()
    print("  l0p5.pdf", flush=True)

    # --- Panel 5: NMP exclude (xenium only) ---
    xen_mask = (obs["platform"] == "xenium").values
    flag = obs.loc[xen_mask, "nmp_exclude"].astype("boolean").fillna(False).astype(str).values
    fig, ax = plt.subplots(figsize=(7, 5.5), dpi=300)
    scatter_panel(ax, x[xen_mask], y[xen_mask], flag,
                  "NMP EXCLUDE (nuclear > cyto) — Xenium",
                  palette=["#bbbbbb", "#d62728"], xlim=xlim, ylim=ylim)
    plt.tight_layout()
    plt.savefig(os.path.join(args.out_dir, "umap_nmp_exclude.pdf"), dpi=300)
    plt.savefig(os.path.join(args.out_dir, "umap_nmp_exclude.png"), dpi=300)
    plt.close()
    print("  nmp_exclude.pdf", flush=True)

    print("Done.", flush=True)


if __name__ == "__main__":
    main()
