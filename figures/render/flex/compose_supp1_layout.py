"""compose_supp1_layout.py — assemble Supp 1 panels into supp1_layout_v1.pdf.

V1 layout (locked 2026-05-08):
- Row 1 (top): 1a pair plot (left) + 1c per-patient UMAP grid (right). Both
  approximately square, ~85mm tall.
- Row 2 (bottom): 1b cell-attrition waterfall, full double-column width,
  ~110mm tall.

Reviewer narrative:
- 1a: pre-integration QC in feature space (kept / MAD / doublet)
- 1c: same attrition story in patient-batch (pre-integration) latent space
- 1b: cohort-wide cell-attrition cascade through all filters (full 5-way)

Total page: 183mm × 247mm (Nature double-column, supp may use full page).

pikepdf preserves vector content of each panel as form XObjects.

Outputs:
    supp1_layout_v1.pdf

Pattern source:
    publication/figures/render/flex/compose_supp5_layout.py
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import List

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("compose_supp1")


MM_PER_INCH = 25.4
PT_PER_INCH = 72.0
PAGE_WIDTH_MM = 183.0
PAGE_HEIGHT_MM = 247.0
GUTTER_MM = 5.0
BORDER_MM = 5.0

ROW_HEIGHT_MM = {
    "row1": 95.0,
    "row2": 110.0,
}

ROWS = [
    ("row1", [
        "supp1_flex_cell_qc_pairplot.pdf",
        "supp1_flex_filter_umap.pdf",
    ]),
    ("row2", [
        "supp1_flex_cell_attrition.pdf",
    ]),
]


def mm_to_pt(mm: float) -> float:
    return mm / MM_PER_INCH * PT_PER_INCH


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--panel-dir", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--page-width-mm", type=float, default=PAGE_WIDTH_MM)
    p.add_argument("--page-height-mm", type=float, default=PAGE_HEIGHT_MM)
    p.add_argument("--gutter-mm", type=float, default=GUTTER_MM)
    p.add_argument("--border-mm", type=float, default=BORDER_MM)
    return p.parse_args()


def collect_panels(panel_dir: Path) -> List[List[Path]]:
    grid: List[List[Path]] = []
    for _row_id, fnames in ROWS:
        row_paths: List[Path] = []
        for fn in fnames:
            p = panel_dir / fn
            if not p.exists():
                sys.exit(f"FATAL: missing panel {p}")
            row_paths.append(p)
        grid.append(row_paths)
    return grid


def compose(
    grid: List[List[Path]], out: Path,
    page_w_mm: float, page_h_mm: float, gutter_mm: float, border_mm: float,
) -> None:
    from pikepdf import Pdf, Rectangle  # type: ignore

    page_w_pt = mm_to_pt(page_w_mm)
    page_h_pt = mm_to_pt(page_h_mm)
    border_pt = mm_to_pt(border_mm)
    gutter_pt = mm_to_pt(gutter_mm)
    inner_w = page_w_pt - 2 * border_pt

    out_pdf = Pdf.new()
    out_pdf.add_blank_page(page_size=(page_w_pt, page_h_pt))
    page = out_pdf.pages[0]

    cur_top = page_h_pt - border_pt
    for row_paths, (row_id, _) in zip(grid, ROWS):
        row_h_pt = mm_to_pt(ROW_HEIGHT_MM[row_id])
        n_cols = len(row_paths)
        cell_w = (inner_w - (n_cols - 1) * gutter_pt) / n_cols
        log.info(
            "row %s: %d cells, cell %.1f×%.1f pt (top=%.1f)",
            row_id, n_cols, cell_w, row_h_pt, cur_top,
        )

        for c, panel_path in enumerate(row_paths):
            src = Pdf.open(panel_path)
            src_page = src.pages[0]
            mb = src_page.mediabox
            src_w = float(mb[2]) - float(mb[0])
            src_h = float(mb[3]) - float(mb[1])

            x0 = border_pt + c * (cell_w + gutter_pt)
            y0 = cur_top - row_h_pt

            scale = min(cell_w / src_w, row_h_pt / src_h)
            placed_w = src_w * scale
            placed_h = src_h * scale
            tx = x0 + (cell_w - placed_w) / 2.0
            ty = y0 + (row_h_pt - placed_h) / 2.0

            page.add_overlay(
                src_page,
                Rectangle(tx, ty, tx + placed_w, ty + placed_h),
            )
            log.info(
                "  placed %s (%.1f×%.1f → %.1f×%.1f at %.1f,%.1f)",
                panel_path.name, src_w, src_h, placed_w, placed_h, tx, ty,
            )

        cur_top = cur_top - row_h_pt - gutter_pt

    out.parent.mkdir(parents=True, exist_ok=True)
    out_pdf.save(out)
    log.info("Wrote %s", out)


def main() -> int:
    args = parse_args()
    grid = collect_panels(args.panel_dir)
    log.info("Collected panels: %d rows", len(grid))
    compose(
        grid, args.out,
        args.page_width_mm, args.page_height_mm,
        args.gutter_mm, args.border_mm,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
