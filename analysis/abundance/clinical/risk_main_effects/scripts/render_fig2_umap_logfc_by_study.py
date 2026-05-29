"""
render_fig2_umap_logfc_by_study.py — Per-study faceted logFC UMAP.

Layout: 4 columns x 2 rows.
  [0,0] = main effect (all studies combined)
  [0,1]-[0,3] = studies 1-3 (alphabetical)
  [1,0]-[1,3] = studies 4-7

Per-study panel rendering:
  - Out-of-study cells: medium grey (#C0C0C0)
  - Within-study, NS (no significant neighborhoods): white
  - Within-study, significant: blue-white-red logFC diverging scale

Color scale is GLOBAL across all subpanels (1st-99th percentile of full
main-effect logFC) so panels are directly comparable.

One PDF per context. Stack 3 contexts vertically in Illustrator for the
final composite (BR1_vs_AR_tested, HRS_vs_AR, parity_in_AR).

Inputs:
  --substrate-dir containing cell_projection.csv (UMAP1, UMAP2, <context>__cell_logfc)
  --labels-csv with cell_id + study column (for study assignment)

Output:
  <out-dir>/s_ihbca_umap_<context>_by_study.{pdf,png}
"""
from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib_config")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("fig2_umap_by_study")

CONTEXT_DISPLAY = {
    "parity_in_AR": "Parity in AR",
    "BR1_vs_AR": "BRCA1 vs AR",
    "BR1_vs_AR_tested": "BRCA1 vs AR-tested",
    "HRS_vs_AR": "HRS vs AR",
    "parity_x_HR_BRCA1": "Parity x BRCA1",
}

STUDY_ORDER = ["gray", "kumar", "murrow", "nee", "pal", "reed", "twigger"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument("--context", required=True)
    p.add_argument("--substrate-dir", required=True, type=Path)
    p.add_argument("--labels-csv", default=None, type=Path,
                   help="CSV with cell_id + index_study columns. "
                        "Default: <project-root>/publication/analysis/annotation/labels_full.csv")
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--color-cap", type=float, default=None)
    p.add_argument("--dpi", type=int, default=600)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    project_root = args.project_root.resolve()

    aes_path = project_root / "publication" / "config" / "aesthetics.yaml"
    import yaml
    with open(aes_path) as f:
        aes = yaml.safe_load(f)

    pt_size = float(aes["dimensions"]["panel_types"]["umap"]["pt_size"])
    scales = aes["scales"]["expression_diverging"]
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "expression_diverging",
        [scales["low"], scales["mid"], scales["high"]],
        N=256,
    )

    log.info("loading cell_projection.csv")
    cell_proj = pd.read_csv(args.substrate_dir / "cell_projection.csv")
    log.info("  cells: %d", len(cell_proj))

    color_col = f"{args.context}__cell_logfc"
    if color_col not in cell_proj.columns:
        available = [c for c in cell_proj.columns if c.endswith("__cell_logfc")]
        raise KeyError(f"Missing {color_col}; available: {available}")

    if "study" not in cell_proj.columns and "index_study" not in cell_proj.columns:
        labels_csv = args.labels_csv or (
            project_root / "publication" / "analysis" / "annotation" / "labels_full.csv"
        )
        log.info("joining study from %s", labels_csv)
        labels = pd.read_csv(labels_csv)
        if "index_study" in labels.columns:
            labels = labels[["cell_id", "index_study"]].rename(columns={"index_study": "study"})
        elif "study" in labels.columns:
            labels = labels[["cell_id", "study"]]
        else:
            raise KeyError(f"No study column in {labels_csv}; cols: {list(labels.columns)}")
        cell_proj = cell_proj.merge(labels, on="cell_id", how="left")
        log.info("  study coverage: %d / %d", cell_proj["study"].notna().sum(), len(cell_proj))
    elif "index_study" in cell_proj.columns:
        cell_proj = cell_proj.rename(columns={"index_study": "study"})

    cell_proj["study"] = cell_proj["study"].str.lower()
    studies_present = sorted(cell_proj["study"].dropna().unique())
    studies = [s for s in STUDY_ORDER if s in studies_present]
    log.info("  studies: %s", studies)

    have_vals = cell_proj[color_col].dropna().to_numpy()
    if args.color_cap is not None:
        lfc_max = float(args.color_cap)
    else:
        q01, q99 = np.nanpercentile(have_vals, [1.0, 99.0])
        lfc_max = float(max(abs(q01), abs(q99), 0.1))
    log.info("  color cap: +-%.3f", lfc_max)
    norm = mcolors.TwoSlopeNorm(vmin=-lfc_max, vcenter=0.0, vmax=lfc_max)

    nrows, ncols = 2, 4
    panel_w, panel_h = 2.2, 2.2
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(ncols * panel_w, nrows * panel_h),
                             squeeze=False)

    umap1 = cell_proj["UMAP1"].to_numpy()
    umap2 = cell_proj["UMAP2"].to_numpy()
    logfc = cell_proj[color_col].to_numpy()
    study_arr = cell_proj["study"].to_numpy()

    def render_panel(ax, title, mask_in_study, logfc_vals):
        out_of_study = ~mask_in_study
        in_study_ns = mask_in_study & np.isnan(logfc_vals)
        in_study_sig = mask_in_study & ~np.isnan(logfc_vals)

        if out_of_study.any():
            ax.scatter(umap1[out_of_study], umap2[out_of_study],
                       s=pt_size * 0.3, c="#C0C0C0", alpha=0.15,
                       linewidths=0, rasterized=True)
        if in_study_ns.any():
            ax.scatter(umap1[in_study_ns], umap2[in_study_ns],
                       s=pt_size * 0.6, c="white", edgecolors="#E0E0E0",
                       linewidths=0.1, alpha=0.4, rasterized=True)
        if in_study_sig.any():
            vals = np.clip(logfc_vals[in_study_sig], -lfc_max, lfc_max)
            ax.scatter(umap1[in_study_sig], umap2[in_study_sig],
                       c=vals, cmap=cmap, norm=norm,
                       s=pt_size * 1.2, alpha=0.75,
                       linewidths=0, rasterized=True)

        ax.set_title(title, fontsize=6, fontweight="bold", pad=2)
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_linewidth(0.3)
        ax.set_aspect("equal", adjustable="datalim")

    panels = ["All studies"] + studies
    for idx, panel_label in enumerate(panels):
        row, col = divmod(idx, ncols)
        if row >= nrows:
            break
        ax = axes[row][col]

        if panel_label == "All studies":
            mask = np.ones(len(cell_proj), dtype=bool)
        else:
            mask = study_arr == panel_label

        display_name = panel_label.capitalize() if panel_label != "All studies" else "All studies"
        render_panel(ax, display_name, mask, logfc)

    for idx in range(len(panels), nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row][col].set_visible(False)

    cbar_ax = fig.add_axes([0.92, 0.25, 0.015, 0.5])
    cbar = fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=cmap),
        cax=cbar_ax,
    )
    cbar.set_label(f"logFC ({CONTEXT_DISPLAY.get(args.context, args.context)})",
                   fontsize=5)
    cbar.ax.tick_params(labelsize=4)

    fig.subplots_adjust(left=0.02, right=0.90, top=0.95, bottom=0.02,
                        wspace=0.08, hspace=0.12)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    panel_id = f"s_ihbca_umap_{args.context}_by_study"
    pdf_path = args.out_dir / f"{panel_id}.pdf"
    png_path = args.out_dir / f"{panel_id}.png"
    fig.savefig(pdf_path, dpi=args.dpi, bbox_inches="tight")
    fig.savefig(png_path, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    log.info("wrote %s", pdf_path)
    log.info("wrote %s", png_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
