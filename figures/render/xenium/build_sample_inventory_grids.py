#!/usr/bin/env python3
"""Sample-inventory grid renders — three plot types as separate PDFs.

For each plot type, produces a single PDF with a grid of small panels:
  rows = patient_id (4), columns = canonical region (17 max).  Empty
  cell when a patient × region combo isn't sampled.

Plot types:
  1. xenium_spatial:  per-sample x/y centroids, colored by L1.5
  2. flex_umap:       joint UMAP filtered to FLEX cells of that region,
                      with the full FLEX cohort as light-grey backdrop
  3. joint_umap:      joint UMAP filtered to all cells of that region
                      (both platforms), with the full joint cohort
                      as light-grey backdrop

Coloring is shared across all three grids (same L1.5 palette).
"""
import pandas as pd
import numpy as np
import os
import sys
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SHARE = ("/Users/NoahWechter/Library/Application Support/CRSP Desktop/"
         "Volumes.noindex/CRSP Lab - dalawson.localized/nwechter")
PREVIEW = (f"{SHARE}/iHBCA_publication/coordination/staging/"
           "fig3_promotion_preview_20260520")
JL = "/tmp/joint_l1p5.csv"
XO = "/tmp/xenium_obs.csv"
UMAP = (f"{SHARE}/Spatial_HBCA_Xenium/pipeline/outputs/embedding/"
        "joint_umap.csv")
UMAP_TMP = "/tmp/joint_umap.csv"
OUT_DIR = (f"{PREVIEW}/figures/output/fig3/supplemental/"
           "sample_inventory")

PATIENTS = ["Pat1", "Pat2", "UCI604", "UCI220228"]
# Canonical region key derived from library_id (strip patient + xenium_N
# suffix).  Sorted canonical order spans 17 regions per breast.
REGION_ORDER = [
    "P1",
    "P2_Upper", "P2_Lower",
    "P2_UOQ", "P2_UIQ", "P2_LOQ", "P2_LIQ",
    "P3_A_UOQ", "P3_A_UIQ", "P3_A_LOQ", "P3_A_LIQ",
    "P3_M_UOQ", "P3_M_UIQ", "P3_M_LOQ", "P3_M_LIQ",
    "P3_P_UOQ", "P3_P_UIQ", "P3_P_LOQ", "P3_P_LIQ",
]


def library_to_region(library_id):
    """Strip patient prefix to get canonical region key."""
    if not isinstance(library_id, str):
        return None
    s = library_id
    # Strip patient prefix (Pat1_, Pat2_, UCI604_, UCI220228_)
    for p in PATIENTS:
        if s.startswith(p + "_"):
            s = s[len(p) + 1:]
            break
    # Strip xenium_N suffix
    s = s.replace("_xenium_1", "").replace("_xenium_2", "")
    return s


# Qualitative palette for L1.5 cell types (Tableau20-style; consistent
# across grids).  Stable assignment regardless of which subset is
# rendered.
L1P5_PALETTE = {
    "Fb": "#1f77b4", "Fb_Activated": "#aec7e8", "Fb_SFRP4": "#5599c9",
    "Adipo": "#ff9896",
    "EC": "#2ca02c", "Vas-cap": "#98df8a", "PV": "#bcbd22",
    "LEC": "#17becf",
    "BMYO-myo": "#8c564b", "LASP-basal": "#c49c94", "LASP": "#d62728",
    "LHS": "#e377c2",
    "CD4T": "#9467bd", "CD8T": "#c5b0d5", "Treg": "#dbdb8d",
    "T-NK": "#7f7f7f",
    "cDC": "#ff7f0e", "cDC1": "#ffbb78", "cDC2": "#f7b6d2",
    "pDC": "#393b79",
    "B": "#637939", "Plas": "#8c6d31",
    "Mac": "#843c39", "Mac_art": "#bd9e39",
    "Mast": "#7b4173", "Neu": "#ad494a",
    "Stromal_art": "#cccccc", "Immune_art": "#dddddd",
    "BMYO_art": "#bbbbbb",
}


def _ensure_umap():
    """Ensure joint umap CSV is cached locally (CRSP read can stale)."""
    if os.path.exists(UMAP_TMP) and os.path.getsize(UMAP_TMP) > 1_000_000:
        return UMAP_TMP
    if os.path.exists(UMAP) and os.path.getsize(UMAP) > 1_000_000:
        return UMAP
    print("[fetch] joint_umap via hpc file cat ...")
    os.system(f"~/.claude/hpc-toolkit/bin/hpc file cat "
              f"{XENIUM_ROOT}/"
              f"pipeline/outputs/embedding/joint_umap.csv "
              f"> {UMAP_TMP}")
    return UMAP_TMP


def _color_array(l1p5_series):
    return np.array([L1P5_PALETTE.get(c, "#999999") for c in l1p5_series])


def render_grid(plot_type, jl, umap_df=None, xo_df=None):
    out_path = f"{OUT_DIR}/inventory_{plot_type}.pdf"

    n_rows = len(PATIENTS)
    n_cols = len(REGION_ORDER)
    fig, axes = plt.subplots(n_rows, n_cols,
                              figsize=(n_cols * 1.3, n_rows * 1.4),
                              squeeze=False)

    if plot_type in ("flex_umap", "joint_umap") and umap_df is not None:
        # Backdrop = all cells in the relevant platform set
        if plot_type == "flex_umap":
            backdrop = umap_df[umap_df["platform"] == "flex"]
        else:
            backdrop = umap_df
        backdrop_xy = backdrop[["UMAP1", "UMAP2"]].values

    for i, pat in enumerate(PATIENTS):
        for j, region in enumerate(REGION_ORDER):
            ax = axes[i, j]
            # Subset cells: patient + region
            sub_jl = jl[(jl["patient_id"] == pat) &
                        (jl["region_key"] == region)]
            if len(sub_jl) == 0:
                ax.set_axis_off()
                continue

            if plot_type == "xenium_spatial":
                xen = sub_jl[sub_jl["platform"] == "xenium"]
                if len(xen) == 0 or xo_df is None:
                    ax.set_axis_off()
                    continue
                # Merge centroids
                m = xen.merge(xo_df, on="cell_id", how="inner")
                if len(m) == 0:
                    ax.set_axis_off()
                    continue
                ax.scatter(m["x_centroid"], m["y_centroid"],
                           c=_color_array(m["l1p5_short"]), s=0.2,
                           alpha=0.6, linewidths=0)
                ax.set_aspect("equal")

            elif plot_type == "flex_umap":
                # Background
                ax.scatter(backdrop_xy[:, 0], backdrop_xy[:, 1],
                           c="#EEEEEE", s=0.15, alpha=0.4, linewidths=0)
                flex = sub_jl[sub_jl["platform"] == "flex"]
                if len(flex) == 0 or umap_df is None:
                    ax.set_xticks([]); ax.set_yticks([])
                    for s in ax.spines.values():
                        s.set_visible(False)
                    continue
                m = flex.merge(umap_df[["cell_id", "UMAP1", "UMAP2"]],
                               on="cell_id", how="inner")
                if len(m) > 0:
                    ax.scatter(m["UMAP1"], m["UMAP2"],
                               c=_color_array(m["l1p5_short"]), s=0.3,
                               alpha=0.75, linewidths=0)

            elif plot_type == "joint_umap":
                ax.scatter(backdrop_xy[:, 0], backdrop_xy[:, 1],
                           c="#EEEEEE", s=0.15, alpha=0.4, linewidths=0)
                m = sub_jl.merge(
                    umap_df[["cell_id", "UMAP1", "UMAP2"]],
                    on="cell_id", how="inner")
                if len(m) > 0:
                    ax.scatter(m["UMAP1"], m["UMAP2"],
                               c=_color_array(m["l1p5_short"]), s=0.3,
                               alpha=0.75, linewidths=0)

            ax.set_xticks([]); ax.set_yticks([])
            for s in ax.spines.values():
                s.set_visible(False)
            if i == 0:
                ax.set_title(region, fontsize=5, pad=2)
            if j == 0:
                ax.set_ylabel(pat, fontsize=6, rotation=90,
                              labelpad=4)

    fig.suptitle(f"Sample inventory — {plot_type}",
                 fontsize=8, y=0.995)
    plt.tight_layout(rect=(0.01, 0.01, 0.99, 0.97), h_pad=0.3, w_pad=0.3)
    fig.savefig(out_path, format="pdf", bbox_inches="tight")
    fig.savefig(out_path.replace(".pdf", ".png"), dpi=200,
                 bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plot-type",
                    choices=["xenium_spatial", "flex_umap", "joint_umap",
                             "all"],
                    default="all")
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)

    print("[load] joint_l1p5 ...")
    jl = pd.read_csv(JL, usecols=[
        "cell_id", "platform", "patient_id", "library_id",
        "p_level", "quadrant", "l1p5_short"])
    jl["region_key"] = jl["library_id"].map(library_to_region)
    print(f"   {len(jl)} cells | regions present: "
          f"{sorted(jl['region_key'].dropna().unique())}")

    umap_df = None
    xo_df = None
    if args.plot_type in ("flex_umap", "joint_umap", "all"):
        umap_path = _ensure_umap()
        print(f"[load] joint_umap from {umap_path} ...")
        umap_df = pd.read_csv(umap_path)
        umap_df = umap_df.merge(jl[["cell_id", "platform"]],
                                 on="cell_id", how="inner")
        print(f"   {len(umap_df)} umap cells")
    if args.plot_type in ("xenium_spatial", "all"):
        print(f"[load] xenium_obs ...")
        xo_df = pd.read_csv(XO, usecols=[
            "cell_id", "x_centroid", "y_centroid", "qc_pass_nuclear"])
        xo_df = xo_df[xo_df["qc_pass_nuclear"] == True]  # noqa: E712
        print(f"   {len(xo_df)} xenium cells with positions")

    plot_types = (["xenium_spatial", "flex_umap", "joint_umap"]
                  if args.plot_type == "all" else [args.plot_type])
    for pt in plot_types:
        print(f"\n[render] {pt} ...")
        render_grid(pt, jl, umap_df=umap_df, xo_df=xo_df)
    return 0


if __name__ == "__main__":
    sys.exit(main())
