"""
render_supp_xenium_three_axis_filter.py — Fig 3 Supp, Xenium three-axis filter.

Phase X4 cohort-filter receipt: three-panel strip showing each filter axis
in its native space, with pass/fail colouring.

    (a) proximity hist        — mean_knn_dist_to_flex, q95 FLEX-self threshold
    (b) NMP scatter           — nmp_nuclear vs nmp_cyto, y=x diagonal
    (c) per-library QC bars   — pass_qc_whole rate per library, colored by patient

Resolves the M-2 §"cohort filter and joint embedding" [VERIFY:per-axis-counts]
flag by surfacing per-axis fail counts as text annotations on each panel.

Inputs:
    three_axis_filter.csv  (1,076,449 rows × 9 cols)
        cell_id, mean_knn_dist_to_flex, nmp_nuclear, nmp_cyto,
        pass_qc_whole, passes_proximity, passes_nmp, passes_qc_whole, passes_all

Outputs:
    s_supp_xenium_three_axis_filter.pdf
    s_supp_xenium_three_axis_filter_stats.csv  (per-axis pass/fail counts)
"""
from __future__ import annotations

import argparse
import logging
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
# Honor PROJECT_ROOT env var (set by SLURM wrapper); fall back to inferring
# from script location (parents[4] = iHBCA_publication root from staging/render/).
import os as _os
_PROJECT_ROOT_ENV = _os.environ.get("PROJECT_ROOT", "")
PROJECT_ROOT_DEFAULT = (Path(_PROJECT_ROOT_ENV) if _PROJECT_ROOT_ENV
                        else SCRIPT_DIR.parents[4])
sys.path.insert(0, str(PROJECT_ROOT_DEFAULT / "publication" / "config"))
from load_aesthetics import (  # noqa: E402
    get_matplotlib_theme,
    get_palette,
    load_aesthetics,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("supp_three_axis")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--three-axis-csv", required=True,
                   help="three_axis_filter.csv (per-cell axis pass/fail, post-X1)")
    p.add_argument("--qc-clean-csv", required=True,
                   help="xenium_qc_clean.csv (pre-X1 cohort, real pass_qc_whole "
                        "variation across libraries — used for panel c)")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--panel-id", default="s_supp_xenium_three_axis_filter")
    p.add_argument("--test", action="store_true",
                   help="Subsample to 10K cells for quick local validation")
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def parse_library(cell_id: str) -> str:
    """`Pat1_P1_xenium_1__aaacmemg-1` -> `Pat1_P1_xenium_1`."""
    return cell_id.split("__", 1)[0]


def parse_patient(cell_id: str) -> str:
    """`Pat1_P1_xenium_1__aaacmemg-1` -> `Pat1`."""
    return cell_id.split("_", 1)[0]


_P3_DEPTH_FULL = {"A": "Anterior", "M": "Middle", "P": "Posterior"}


def parse_position_short(library_id: str) -> str:
    """`Pat1_P3_A_LIQ_xenium_1` -> `P3-A-LIQ`. Compact label for x-axis grouping.

    Handles all observed formats:
      Pat1_P1_xenium_1            -> P1
      Pat1_P2_Upper_xenium_1      -> P2-Upper
      Pat2_P2_LIQ_xenium_1        -> P2-LIQ
      Pat1_P3_M_LOQ_xenium_1      -> P3-M-LOQ
    """
    # Strip `_xenium_N` suffix
    base = re.sub(r"_xenium_\d+$", "", library_id)
    # Strip patient prefix
    parts = base.split("_", 1)
    return parts[1].replace("_", "-") if len(parts) > 1 else "?"


def parse_p_level(library_id: str) -> str:
    """`Pat1_P3_A_LIQ_xenium_1` -> `P3`."""
    m = re.search(r"_(P\d)_", library_id + "_")
    return m.group(1) if m else "?"


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    log.info("Loading %s", args.three_axis_csv)
    df = pd.read_csv(args.three_axis_csv)
    log.info("Loaded %d cells, columns: %s", len(df), list(df.columns))

    if args.test:
        rng = np.random.default_rng(args.seed)
        idx = rng.choice(len(df), size=min(10_000, len(df)), replace=False)
        df = df.iloc[idx].reset_index(drop=True)
        log.info("Test mode: subsampled to %d cells", len(df))

    df["library_id"] = df["cell_id"].map(parse_library)
    df["patient_id"] = df["cell_id"].map(parse_patient)

    aesthetics = load_aesthetics()
    get_matplotlib_theme(aesthetics)

    qc_pal = get_palette("qc_status", aesthetics)
    color_pass = qc_pal["kept"]
    color_fail = qc_pal["removed_initial"]

    try:
        patient_pal = get_palette("patient", aesthetics)
    except KeyError:
        patient_pal = None
        log.warning("No 'patient' palette; falling back to tab10")

    # ---- Stats ----
    n = len(df)
    stats_rows = []
    for axis in ("passes_proximity", "passes_nmp", "passes_qc_whole", "passes_all"):
        n_pass = int(df[axis].sum())
        stats_rows.append({"axis": axis, "n_pass": n_pass, "n_fail": n - n_pass,
                           "pass_rate": n_pass / n})
    stats = pd.DataFrame(stats_rows)
    stats_path = out_dir / f"{args.panel_id}_stats.csv"
    stats.to_csv(stats_path, index=False)
    log.info("Stats written: %s", stats_path)
    log.info("\n%s", stats.to_string(index=False))

    # ---- Figure ----
    # Double-column 2-row layout (Nature 183mm = 7.2 in).
    #   Top row: panels (a) and (b) side by side (continuous-axis filters).
    #   Bottom row: panel (c) full width — per-library attrition bars.
    fig = plt.figure(figsize=(7.2, 4.8))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.75],
                          hspace=0.65, wspace=0.35,
                          left=0.07, right=0.99,
                          bottom=0.14, top=0.93)

    # ---- (a) Proximity hist ----
    ax_a = fig.add_subplot(gs[0, 0])
    x = df["mean_knn_dist_to_flex"].to_numpy(dtype=float)
    # Threshold = max(x[passes_proximity]) — this is the q95 cut from the filter step
    thr = float(df.loc[df["passes_proximity"], "mean_knn_dist_to_flex"].max())
    bins = np.linspace(np.nanmin(x), np.nanquantile(x, 0.999), 60)
    ax_a.hist(x[df["passes_proximity"]], bins=bins, color=color_pass,
              alpha=0.85, label="pass")
    ax_a.hist(x[~df["passes_proximity"]], bins=bins, color=color_fail,
              alpha=0.85, label="fail")
    ax_a.axvline(thr, color="black", lw=0.6, ls="--")
    ax_a.text(thr, ax_a.get_ylim()[1] * 0.95, f" q95 = {thr:.3f}",
              fontsize=5, va="top", ha="left")
    ax_a.set_xlabel("mean kNN dist to FLEX (joint latent)", fontsize=6)
    ax_a.set_ylabel("cells", fontsize=6)
    fail_n = int((~df["passes_proximity"]).sum())
    ax_a.text(0.98, 0.98, f"fail: {fail_n:,}", transform=ax_a.transAxes,
              fontsize=5, ha="right", va="top")
    ax_a.legend(fontsize=5, frameon=False, loc="upper left",
                handlelength=1.0, handletextpad=0.4)

    # ---- (b) NMP scatter ----
    ax_b = fig.add_subplot(gs[0, 1])
    nuc = df["nmp_nuclear"].to_numpy(dtype=float)
    cyt = df["nmp_cyto"].to_numpy(dtype=float)
    fail_mask = ~df["passes_nmp"].to_numpy()
    # Shuffle to avoid one class painting over the other
    order = np.random.default_rng(args.seed).permutation(len(df))
    nuc_s, cyt_s, fail_s = nuc[order], cyt[order], fail_mask[order]
    colors = np.where(fail_s, color_fail, color_pass)
    ax_b.scatter(nuc_s, cyt_s, s=0.6, c=colors, alpha=0.25,
                 linewidths=0, rasterized=True)
    lim_hi = float(np.nanquantile(np.concatenate([nuc, cyt]), 0.99))
    ax_b.plot([0, lim_hi], [0, lim_hi], color="black", lw=0.6, ls="--")
    ax_b.set_xlim(0, lim_hi)
    ax_b.set_ylim(0, lim_hi)
    ax_b.set_xlabel("NMP nuclear (%)", fontsize=6)
    ax_b.set_ylabel("NMP cytoplasmic (%)", fontsize=6)
    ax_b.set_aspect("equal", adjustable="box")
    fail_b = int(fail_mask.sum())
    ax_b.text(0.98, 0.02, f"fail (nuc > cyto): {fail_b:,}",
              transform=ax_b.transAxes, fontsize=5, ha="right", va="bottom")

    # ---- (c) Per-library QC pass rate ----
    # Substrate: xenium_qc_clean.csv (pre-X1 cohort with real pass_qc_whole
    # variation). df (three_axis_filter) is post-X1 and has pass_qc_whole == True
    # for every row, which would render flat bars at 1.0.
    log.info("Loading QC-clean substrate for panel c: %s", args.qc_clean_csv)
    qc = pd.read_csv(args.qc_clean_csv,
                     usecols=["cell_id", "pass_qc_whole"])
    qc["library_id"] = qc["cell_id"].map(parse_library)
    qc["patient_id"] = qc["cell_id"].map(parse_patient)
    log.info("QC-clean: %d cells, %d libraries; pass_qc_whole rate = %.3f",
             len(qc), qc["library_id"].nunique(),
             qc["pass_qc_whole"].astype(bool).mean())

    ax_c = fig.add_subplot(gs[1, :])

    # Build per-cell fate categories by joining qc_clean (pre-X1 cohort) with
    # three_axis_filter (post-X1, with per-axis status):
    #   fail_qc        : in qc but pass_qc_whole=False           (X1 noise-floor drop)
    #   fail_prox_only : in three_axis, fails proximity only
    #   fail_nmp_only  : in three_axis, fails NMP only
    #   fail_multi     : in three_axis, fails 2+ axes (prox + NMP)
    #   kept           : in three_axis, passes_all
    fate = qc[["cell_id", "pass_qc_whole", "library_id", "patient_id"]].copy()
    fate["pass_qc_whole"] = fate["pass_qc_whole"].astype(bool)
    tax = df[["cell_id", "passes_proximity", "passes_nmp", "passes_all"]].copy()
    fate = fate.merge(tax, on="cell_id", how="left")
    n_fail_axes = (
        (~fate["passes_proximity"].fillna(True).astype(bool)).astype(int) +
        (~fate["passes_nmp"].fillna(True).astype(bool)).astype(int)
    )
    fate["category"] = "kept"
    fate.loc[~fate["pass_qc_whole"], "category"] = "fail_qc"
    in_tax = fate["passes_all"].notna()
    fate.loc[in_tax & ~fate["passes_all"].fillna(False).astype(bool) &
             (n_fail_axes == 1) & ~fate["passes_proximity"].fillna(True).astype(bool),
             "category"] = "fail_prox_only"
    fate.loc[in_tax & ~fate["passes_all"].fillna(False).astype(bool) &
             (n_fail_axes == 1) & ~fate["passes_nmp"].fillna(True).astype(bool),
             "category"] = "fail_nmp_only"
    fate.loc[in_tax & ~fate["passes_all"].fillna(False).astype(bool) &
             (n_fail_axes >= 2),
             "category"] = "fail_multi"

    cat_order = ["kept", "fail_qc", "fail_prox_only", "fail_nmp_only", "fail_multi"]
    cat_colors = {
        "kept":           qc_pal["kept"],
        "fail_qc":        qc_pal["removed_initial"],
        "fail_prox_only": qc_pal["removed_manifold"],
        "fail_nmp_only":  qc_pal["removed_contam"],
        "fail_multi":     qc_pal["removed_doublet"],
    }
    cat_labels = {
        "kept":           "kept (pass all 3 axes)",
        "fail_qc":        "fail X1 QC (noise floor)",
        "fail_prox_only": "fail proximity only",
        "fail_nmp_only":  "fail NMP only",
        "fail_multi":     "fail 2+ axes",
    }

    by_lib_cat = (
        fate.groupby(["library_id", "patient_id", "category"]).size()
            .unstack(fill_value=0)
    )
    for c in cat_order:
        if c not in by_lib_cat.columns:
            by_lib_cat[c] = 0
    by_lib_cat = by_lib_cat[cat_order]
    totals = by_lib_cat.sum(axis=1).replace(0, 1)
    props = by_lib_cat.div(totals, axis=0)

    libs = props.index.get_level_values("library_id").tolist()
    pats = props.index.get_level_values("patient_id").tolist()
    positions = [parse_position_short(l) for l in libs]
    order_df = pd.DataFrame({"library_id": libs, "patient_id": pats,
                             "position": positions, "total": totals.values})
    order_df = order_df.sort_values(["patient_id", "position", "library_id"]).reset_index(drop=True)
    order_df["x"] = np.arange(len(order_df))

    props = props.loc[
        pd.MultiIndex.from_arrays([order_df["library_id"], order_df["patient_id"]],
                                  names=["library_id", "patient_id"])
    ]

    # Stacked bars
    bottom = np.zeros(len(order_df))
    for c in cat_order:
        ax_c.bar(order_df["x"], props[c].values, bottom=bottom,
                 color=cat_colors[c], width=0.92, edgecolor="none",
                 label=cat_labels[c])
        bottom = bottom + props[c].values

    ax_c.set_ylim(0, 1.0)
    ax_c.set_xlim(-0.6, len(order_df) - 0.4)
    ax_c.set_ylabel("cell fate (proportion)", fontsize=6)
    ax_c.set_xlabel(f"library (n={len(order_df)}, ordered patient → position)", fontsize=6)
    ax_c.set_xticks(order_df["x"].values)
    ax_c.set_xticklabels(order_df["position"].values, rotation=45,
                         ha="right", rotation_mode="anchor", fontsize=4)
    ax_c.tick_params(axis="x", length=2, pad=1)

    # Patient bracket labels + separator lines
    uniq_pats = list(dict.fromkeys(order_df["patient_id"]))
    if patient_pal is not None:
        pat_to_color = {p: patient_pal.get(p, "#888888") for p in uniq_pats}
    else:
        cmap = plt.get_cmap("tab10")
        pat_to_color = {p: cmap(i) for i, p in enumerate(uniq_pats)}
    pat_starts = order_df.groupby("patient_id")["x"].agg(["min", "max"])
    for pat in uniq_pats:
        lo = pat_starts.loc[pat, "min"]
        hi = pat_starts.loc[pat, "max"]
        ax_c.text((lo + hi) / 2, 1.05, pat, transform=ax_c.transData,
                  fontsize=6, ha="center", va="bottom",
                  fontweight="bold", color=pat_to_color[pat])
        if pat != uniq_pats[-1]:
            ax_c.axvline(hi + 0.5, color="#cccccc", lw=0.4, zorder=0)

    # Bottom legend for fate categories
    legend_handles = [
        Rectangle((0, 0), 1, 1, fc=cat_colors[c], ec="none", label=cat_labels[c])
        for c in cat_order
    ]
    fig.legend(handles=legend_handles, loc="lower center",
               ncol=len(cat_order), fontsize=5.5, frameon=False,
               handlelength=1.0, handletextpad=0.4, columnspacing=1.2,
               bbox_to_anchor=(0.5, -0.01))

    # ---- Per-axis "intersection" footnote (top, above the layout) ----
    n_all = int(df["passes_all"].sum())
    footnote = (
        f"n = {n:,} post-X1 cells | pass all 3 axes: {n_all:,} ({100*n_all/n:.1f}%)"
    )
    fig.text(0.5, 0.97, footnote, fontsize=5, ha="center", va="top")

    # ---- Save ----
    pdf_path = out_dir / f"{args.panel_id}.pdf"
    png_path = out_dir / f"{args.panel_id}.png"
    fig.savefig(pdf_path, format="pdf", dpi=600, bbox_inches="tight")
    fig.savefig(png_path, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved %s and %s", pdf_path, png_path)


if __name__ == "__main__":
    main()
