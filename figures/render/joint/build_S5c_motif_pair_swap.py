"""S5c (true) — Motif-pair NMF refit swap-survival heatmap.

For each motif pair (i, j): swap positions of M_i and M_j cells, recompute
the K-NN composition substrate, refit NMF k=9, measure per-motif
top-cluster retention. Plot the mean retention of (M_i, M_j) under their
own pairwise swap: lower = more compositionally distinct.
"""
from __future__ import annotations
import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

mpl.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 6, "axes.labelsize": 7, "axes.titlesize": 7,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
    "pdf.fonttype": 42,
})

PROG_PAL = {
    "M0": "#E6194B", "M1": "#3CB44B", "M2": "#4363D8", "M3": "#F58231",
    "M4": "#911EB4", "M5": "#42D4F4", "M6": "#F032E6", "M7": "#BFEF45",
    "M8": "#FABED4",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.tsv, sep="\t")
    df = df.drop_duplicates(subset=["mode", "seed", "radius", "pair"])
    p = df[df["mode"] == "motif-pair"].copy()
    print(f"S5c: {len(p)} motif-pair tasks")

    def parse_pair(s):
        return tuple(int(x) for x in str(s).replace(",", "-").split("-"))
    p["i"] = p["pair"].apply(lambda s: parse_pair(s)[0])
    p["j"] = p["pair"].apply(lambda s: parse_pair(s)[1])

    mat = np.full((9, 9), np.nan)
    n_seeds = np.zeros((9, 9), dtype=int)
    for _, r in p.iterrows():
        i, j = int(r["i"]), int(r["j"])
        ri, rj = r[f"ret_M{i}"], r[f"ret_M{j}"]
        val = (ri + rj) / 2
        if np.isnan(mat[i, j]):
            mat[i, j] = val
            n_seeds[i, j] = 1
        else:
            n = n_seeds[i, j]
            mat[i, j] = (mat[i, j] * n + val) / (n + 1)
            n_seeds[i, j] = n + 1
        mat[j, i] = mat[i, j]
        n_seeds[j, i] = n_seeds[i, j]

    motifs = [f"M{i}" for i in range(9)]

    # 90mm × 90mm Nature single-column square heatmap
    fig, ax = plt.subplots(figsize=(95 / 25.4, 90 / 25.4))

    cmap = mpl.colormaps["RdBu_r"]
    vmin, vmax = 0.3, 0.9
    im = ax.imshow(mat, cmap=cmap, vmin=vmin, vmax=vmax, aspect="equal")

    # Diagonal mask
    for i in range(9):
        ax.add_patch(plt.Rectangle((i - 0.5, i - 0.5), 1, 1,
                                     facecolor="#222222", edgecolor="none"))

    # Annotate cells
    for i in range(9):
        for j in range(9):
            if i == j:
                continue
            v = mat[i, j]
            if np.isnan(v):
                continue
            txt_color = "white" if (v < 0.45 or v > 0.75) else "black"
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                    color=txt_color, fontsize=5)

    ax.set_xticks(np.arange(9))
    ax.set_yticks(np.arange(9))
    ax.set_xticklabels(motifs)
    ax.set_yticklabels(motifs)
    for tick, m in zip(ax.get_xticklabels(), motifs):
        tick.set_color(PROG_PAL[m]); tick.set_fontweight("bold")
    for tick, m in zip(ax.get_yticklabels(), motifs):
        tick.set_color(PROG_PAL[m]); tick.set_fontweight("bold")
    ax.set_xlabel("M_j")
    ax.set_ylabel("M_i")

    cbar = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02, aspect=18,
                         orientation="vertical")
    cbar.set_label("swap survival", fontsize=6)
    cbar.ax.tick_params(width=0.5, length=2, labelsize=6)
    cbar.outline.set_linewidth(0.4)

    for spine in ["top", "right", "left", "bottom"]:
        ax.spines[spine].set_linewidth(0.5)
    ax.tick_params(width=0.5, length=2)

    plt.tight_layout()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out, dpi=600, bbox_inches="tight")
    plt.savefig(str(out).replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    print(f"saved {out}")

    # Print summary
    print("\nMost distinct pairs (lowest swap survival):")
    pairs = [(i, j, mat[i, j]) for i in range(9) for j in range(i + 1, 9) if not np.isnan(mat[i, j])]
    pairs.sort(key=lambda x: x[2])
    for i, j, v in pairs[:6]:
        print(f"  M{i} <-> M{j}: {v:.2f}")
    print("\nMost interchangeable pairs (highest swap survival):")
    for i, j, v in pairs[-5:]:
        print(f"  M{i} <-> M{j}: {v:.2f}")


if __name__ == "__main__":
    main()
