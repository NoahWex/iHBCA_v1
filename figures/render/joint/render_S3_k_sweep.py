#!/usr/bin/env python3
"""Render Fig 3 supp — NMF rank selection (K-sweep).

Two-panel: (a) reconstruction error vs K, with chosen K marked; (b) second
difference of recon_err (elbow detector), peak marks the chosen K.

Reads canonical K-sweep diagnostics produced by
fig4_lane_c_nmf_20260515/scripts/01_cell_nmf.py.
"""
from __future__ import annotations
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--diagnostics", required=True,
                    help="diagnostics.csv (k, recon_err, second_diff)")
    ap.add_argument("--k-chosen", type=int, default=9)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.diagnostics)

    fig, axes = plt.subplots(1, 2, figsize=(4.5, 1.8), dpi=300,
                              gridspec_kw={"wspace": 0.35})

    # (a) recon_err vs K
    ax = axes[0]
    ax.plot(df["k"], df["recon_err"], color="#222", linewidth=0.8,
             marker="o", markersize=3, markerfacecolor="#222")
    ax.axvline(args.k_chosen, color="#cc0000", linewidth=0.6,
                linestyle="--", alpha=0.7)
    ax.text(args.k_chosen + 0.15, df["recon_err"].max() * 0.95,
             f"k={args.k_chosen}", color="#cc0000", fontsize=6,
             ha="left", va="top")
    ax.set_xlabel("rank k", fontsize=6)
    ax.set_ylabel("reconstruction error", fontsize=6)
    ax.tick_params(labelsize=5, length=1.5)
    for spine in ax.spines.values():
        spine.set_linewidth(0.3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # (b) second difference
    ax = axes[1]
    sub = df.dropna(subset=["second_diff"])
    ax.bar(sub["k"], sub["second_diff"], width=0.7,
            color=["#cc0000" if k == args.k_chosen else "#999"
                    for k in sub["k"]], linewidth=0)
    ax.axhline(0, color="black", linewidth=0.3)
    ax.set_xlabel("rank k", fontsize=6)
    ax.set_ylabel("second difference\nof reconstruction error", fontsize=6)
    ax.tick_params(labelsize=5, length=1.5)
    for spine in ax.spines.values():
        spine.set_linewidth(0.3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
