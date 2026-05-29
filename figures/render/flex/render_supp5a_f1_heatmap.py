"""render_supp5a_f1_heatmap.py — Fig 3 Supp 5a.

Per-compartment F1 heatmap: V1 L2 labels (rows, UCell-positive at threshold)
× Leiden clusters (columns) at the anchor resolution 1.0. Validates whether
the resolution-1.0 partition resolves each V1 label cleanly.

The F1 matrix is the deterministic output of `compute_v1_label_f1_stack` in
publication/analysis/annotation/flex/scripts/05_render_structural.py:96-129.
This render reads a pre-computed CSV with that same shape rather than
re-running scripts 01-02 on the FLEX intermediates (substrate decision logged
in coordination/reports/CP_supp5a_substrate.md).

Pattern reference:
- Plot styling (sns.heatmap with vmin=0, vmax=1, Blues cmap): publication/
  analysis/annotation/flex/scripts/05_render_structural.py:286-296
  (emit_f1_heatmaps).
- Diagonal ordering of rows + columns (place each row's argmax column near
  the diagonal): same script lines 62-81 (diagonal_ordering).

No on-plot title or panel letter — Illustrator composition adds those.

Outputs:
    supp5a_flex_f1_{epi,imm,str}.pdf
    supp5a_flex_f1_{epi,imm,str}_data.csv (provenance copy of F1 matrix)
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import List, Tuple

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib_config")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("supp5a")


COMPARTMENT_FULL = {"epi": "Epithelial", "imm": "Immune", "str": "Stromal"}


def _setup_aesthetics(project_root: Path) -> dict:
    sys.path.insert(0, str(project_root / "publication" / "config"))
    from load_aesthetics import load_aesthetics, get_matplotlib_theme  # type: ignore
    cfg = load_aesthetics()
    plt.rcParams.update(get_matplotlib_theme(cfg))
    return cfg


def diagonal_ordering(mat: pd.DataFrame) -> Tuple[List[str], List[str]]:
    """Order rows by argmax-column then descending max; order columns by row of
    first appearance. Direct port of 05_render_structural.py:62-81."""
    col_names = list(mat.columns)
    col_to_pos = {c: i for i, c in enumerate(col_names)}
    row_argmax = mat.idxmax(axis=1).map(col_to_pos)
    row_max = mat.max(axis=1)
    row_order = (
        pd.DataFrame({"argmax": row_argmax, "neg_max": -row_max}, index=mat.index)
        .sort_values(by=["argmax", "neg_max"])
        .index.tolist()
    )
    mat_r = mat.loc[row_order]
    col_first: dict = {}
    for ridx, cname in enumerate(mat_r.idxmax(axis=1).values):
        col_first.setdefault(cname, ridx)
    col_max = mat.max(axis=0)
    cols_used = sorted(
        [c for c in col_names if c in col_first], key=lambda c: col_first[c]
    )
    cols_unused = sorted(
        [c for c in col_names if c not in col_first], key=lambda c: -float(col_max[c])
    )
    return row_order, cols_used + cols_unused


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project-root", type=Path, required=True)
    p.add_argument("--compartment", required=True, choices=["epi", "imm", "str"])
    p.add_argument("--f1-csv", type=Path, required=True,
                   help="Pre-computed F1 matrix CSV (rows=V1 labels, cols=clusters at "
                        "Leiden res 1.0). Produced by 05_render_structural.py "
                        "(compute_v1_label_f1_stack at threshold 0.5).")
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--test", action="store_true",
                   help="Trim columns to 5 for fast iteration")
    p.add_argument("--reorder", action="store_true",
                   help="Re-apply diagonal_ordering. Off by default — dev CSVs are "
                        "already diagonal-ordered by 05_render_structural emit_f1_heatmaps.")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    project_root = args.project_root.resolve()
    cfg = _setup_aesthetics(project_root)

    if not args.f1_csv.exists():
        sys.exit(f"FATAL: F1 CSV missing: {args.f1_csv}")

    f1 = pd.read_csv(args.f1_csv, index_col=0)
    log.info("Loaded F1 matrix: %d V1 labels × %d clusters", *f1.shape)

    # Cast cluster columns to string for stable label rendering.
    f1.columns = [str(c) for c in f1.columns]

    if args.reorder:
        row_order, col_order = diagonal_ordering(f1)
        f1 = f1.loc[row_order, col_order]
        log.info("Re-applied diagonal_ordering")

    if args.test:
        f1 = f1.iloc[:, :5]
        log.info("--test: trimmed to %d clusters", f1.shape[1])

    # Per Noah v3 review (2026-05-07): bump panel widths from 60mm (v2) to 75mm
    # so F1 confusion cells have more legible per-cluster ticks. Compose with
    # 3mm gutters at v3 layout (3×75 + 2×3 = 231mm > 183mm page; composer
    # tightens gutters and reduces inter-panel margin to fit, see compose_*.py
    # ROW_HEIGHT_MM and ROW_GUTTER_MM_BY_ROW). Height scales with native aspect.
    n_rows, n_cols = f1.shape
    width_in = 75.0 / 25.4   # 75mm
    aspect = max(0.4, min(1.6, n_rows / max(n_cols, 1)))
    height_in = width_in * aspect + 0.4   # +0.4in for axis labels + colorbar

    fig, ax = plt.subplots(figsize=(width_in, height_in))
    sns.heatmap(
        f1, cmap="Blues", vmin=0, vmax=1,
        cbar_kws={"label": "F1", "shrink": 0.4},
        linewidths=0, ax=ax,
        xticklabels=True, yticklabels=True,
    )
    ax.set_xlabel("Leiden cluster (res 1.0)", fontsize=7)
    ax.set_ylabel("V1 L2 label (UCell-positive)", fontsize=7)
    ax.tick_params(axis="x", labelsize=5, rotation=90)
    ax.tick_params(axis="y", labelsize=5)
    plt.tight_layout()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_pdf = args.out_dir / f"supp5a_flex_f1_{args.compartment}.pdf"
    out_csv = args.out_dir / f"supp5a_flex_f1_{args.compartment}_data.csv"
    fig.savefig(out_pdf, bbox_inches="tight",
                dpi=int(cfg.get("rendering", {}).get("dpi", 600)))
    plt.close(fig)
    f1.to_csv(out_csv, index_label="v1_label")
    log.info("Wrote %s + %s", out_pdf, out_csv)
    return 0


if __name__ == "__main__":
    sys.exit(main())
