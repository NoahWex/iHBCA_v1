#!/usr/bin/env python3
"""s_S5d — Local-radius shuffle null: per-motif retention as a function of
shuffle radius. For each radius r ∈ {25, 50, 100, 200, 500, 1000 µm}, cell
positions are shuffled within local neighbourhoods of radius r (cell-type
identities preserved) and the dominant motif assignment is re-derived under
projection through the canonical basis. The fraction of cells whose
dominant-motif assignment is unchanged is reported per motif.

Reads pre-computed retention data from null_results.tsv (mode=local, 3 seeds
per radius), aggregates per-motif mean ± range across seeds, and renders a
line plot with one line per motif. Companion to s_S5c (paired positional
swap) — together they form the distribution side of motif validation.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

SHARE = ("/Users/NoahWechter/Library/Application Support/CRSP Desktop/"
         "Volumes.noindex/CRSP Lab - dalawson.localized/nwechter")
SUB = (f"{SHARE}/Spatial_HBCA_Xenium/NicheFramework/"
       "fig4_lane_c_nmf_20260515/reports/validation/null_results.tsv")
OUT_DIR = "/tmp/s_S5d_local_radius_shuffle"

MOTIFS = [f"M{i}" for i in range(9)]
MOTIF_LABELS = {
    "M0": "M0  Fb_Activated", "M1": "M1  EC", "M2": "M2  Fb",
    "M3": "M3  PV", "M4": "M4  T+DC", "M5": "M5  LEC",
    "M6": "M6  B", "M7": "M7  Mac", "M8": "M8  Plas",
}
# Diverging palette: baseline (M0/M1/M2) in stromal-grey/blue tones,
# specialized (M3-M8) in distinct hues
MOTIF_COLORS = {
    "M0": "#888888",  # stromal grey
    "M1": "#4477AA",  # endothelial blue
    "M2": "#888844",  # fibroblast olive
    "M3": "#CC6677",  # perivascular pink
    "M4": "#117733",  # T+DC green
    "M5": "#66CCEE",  # lymphatic light-blue
    "M6": "#AA3377",  # B-cell magenta
    "M7": "#EE8866",  # macrophage orange
    "M8": "#882255",  # plasma dark-burgundy
}


def main():
    Path(OUT_DIR).mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(SUB, sep="\t")
    loc = df[df["mode"] == "local"].copy()
    loc["radius"] = loc["radius"].astype(float)
    print(f"[load] {len(loc)} local rows across "
          f"{loc['radius'].nunique()} radii × {loc['seed'].nunique()} seeds",
          flush=True)

    # Aggregate: per-motif mean across seeds at each radius
    ret_cols = [f"ret_{m}" for m in MOTIFS]
    agg = loc.groupby("radius")[ret_cols].agg(["mean", "min", "max"])
    radii = sorted(loc["radius"].unique())
    print(f"[radii] {radii}", flush=True)

    # === Render ===
    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    for m in MOTIFS:
        means = [agg.loc[r, (f"ret_{m}", "mean")] for r in radii]
        mins  = [agg.loc[r, (f"ret_{m}", "min")]  for r in radii]
        maxs  = [agg.loc[r, (f"ret_{m}", "max")]  for r in radii]
        ax.plot(radii, means, "-o", color=MOTIF_COLORS[m],
                 label=MOTIF_LABELS[m], linewidth=1.4, markersize=4)
        ax.fill_between(radii, mins, maxs, color=MOTIF_COLORS[m], alpha=0.15)

    ax.set_xscale("log")
    ax.set_xlim(20, 1200)
    ax.set_ylim(0, 1.0)
    ax.set_xlabel("Local shuffle radius (µm)", fontsize=7)
    ax.set_ylabel("Dominant-motif retention", fontsize=7)
    ax.set_xticks(radii)
    ax.set_xticklabels([str(int(r)) for r in radii], fontsize=6)
    ax.tick_params(axis="y", labelsize=6)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5),
               fontsize=5.5, frameon=False, ncol=1)

    plt.tight_layout()
    out_pdf = f"{OUT_DIR}/s_S5d_local_radius_shuffle.pdf"
    fig.savefig(out_pdf, format="pdf", bbox_inches="tight")
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[write] {out_pdf}", flush=True)

    # Long-form CSV for substrate
    long_rows = []
    for r in radii:
        for m in MOTIFS:
            long_rows.append({
                "radius_um": int(r),
                "motif": m,
                "retention_mean": float(agg.loc[r, (f"ret_{m}", "mean")]),
                "retention_min":  float(agg.loc[r, (f"ret_{m}", "min")]),
                "retention_max":  float(agg.loc[r, (f"ret_{m}", "max")]),
                "n_seeds": int(loc[loc["radius"] == r].shape[0]),
            })
    pd.DataFrame(long_rows).to_csv(
        f"{OUT_DIR}/s_S5d_long_table.csv", index=False)
    print(f"[write] long table", flush=True)
    print("[done]", flush=True)


if __name__ == "__main__":
    main()
