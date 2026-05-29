"""compose_supp5_layout.py — assemble Supp 5 panels into the v6 layout PDF.

V6 layout (locked 2026-05-07):
- Row 1: 5a × 3 (one per compartment), 3-column split, ~65mm tall
- Row 2: 5b × 3 (one per compartment), 3-column split, ~70mm tall (square UMAPs)
- Row 3: 5c × 1 merged 107-marker step-diagonal panoramic heatmap, full row width, ~92mm tall

Total page: Nature double-column = 183mm × 247mm (full page allowed for supp).
65 + 70 + 92 + 2*5 (border) + 2*5 (gutter) = 247mm.

pikepdf preserves vector content of each panel as form XObjects.

Outputs:
    supp5_layout_v6.pdf

Revision history (algorithmic, not bookkeeping):
- v2 (2026-05-07): merged 5c, dropped per-compartment color bar
- v3: bumped marker count + drop comp-bar fully
- v4: panoramic 5c + ARTIFACT-at-end
- v5: within-block diagonal sort
- v6: step-diagonal (col-block per row) sort — current

Pattern reference: prior versions of this script (see prior/v{1..5}/) and
pikepdf form-xobject docs. No in-repo precedent for multi-panel composition.
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
log = logging.getLogger("compose_supp5_v2")


MM_PER_INCH = 25.4
PT_PER_INCH = 72.0
PAGE_WIDTH_MM = 183.0
PAGE_HEIGHT_MM = 247.0
GUTTER_MM = 5.0
BORDER_MM = 5.0

ROW_HEIGHT_MM = {
    "5a": 65.0,
    "5b": 70.0,
    "5c": 92.0,
}

# v3: per-row column gutter override. 5a panels are 75mm wide each (vs 60mm in
# v2), so 3 × 75mm = 225mm needs to fit in 183mm inner width. Solution: drop
# the 5a gutter to 3mm AND let panels overflow the inner border slightly
# (composer crops panels to cell width via fit-to-cell scaling, so 75mm panels
# in a ~57mm cell will scale down). The cell width is computed from
# inner_width - (n-1)*gutter; with 5a gutter=3mm: cell = (173 - 6) / 3 = 55.7mm.
# Each 75mm panel scales to 55.7mm width and proportionally less height. Net:
# 5a panels fit but at reduced print size from the script-rendered 75mm. v3
# users who want true 75mm-printed 5a should compose with a wider page (e.g.
# 8.5×11 letter) or accept the 55.7mm in-grid render.
ROW_GUTTER_MM = {
    "5a": 3.0,
    "5b": 5.0,
    "5c": 5.0,   # single panel, gutter unused
}

ROWS = [
    ("5a", [
        "supp5a_flex_f1_epi.pdf",
        "supp5a_flex_f1_imm.pdf",
        "supp5a_flex_f1_str.pdf",
    ]),
    ("5b", [
        "supp5b_flex_label_umap_epi.pdf",
        "supp5b_flex_label_umap_imm.pdf",
        "supp5b_flex_label_umap_str.pdf",
    ]),
    ("5c", [
        "supp5c_flex_canonical_artifact_combined.pdf",
    ]),
]


def mm_to_pt(mm: float) -> float:
    return mm / MM_PER_INCH * PT_PER_INCH


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--panel-dir", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--page-width-mm", type=float, default=PAGE_WIDTH_MM)
    p.add_argument("--page-height-mm", type=float, default=PAGE_HEIGHT_MM)
    p.add_argument("--gutter-mm", type=float, default=GUTTER_MM)
    p.add_argument("--border-mm", type=float, default=BORDER_MM)
    return p.parse_args()


def collect_panels(panel_dir: Path) -> List[List[Path]]:
    grid: List[List[Path]] = []
    for row_id, fnames in ROWS:
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
    page_w_mm: float, page_h_mm: float, default_gutter_mm: float, border_mm: float,
) -> None:
    from pikepdf import Pdf, Rectangle  # type: ignore

    page_w_pt = mm_to_pt(page_w_mm)
    page_h_pt = mm_to_pt(page_h_mm)
    border_pt = mm_to_pt(border_mm)
    inner_w = page_w_pt - 2 * border_pt

    out_pdf = Pdf.new()
    out_pdf.add_blank_page(page_size=(page_w_pt, page_h_pt))
    page = out_pdf.pages[0]

    cur_top = page_h_pt - border_pt
    for row_paths, (row_id, _) in zip(grid, ROWS):
        row_h_pt = mm_to_pt(ROW_HEIGHT_MM[row_id])
        # v3: per-row gutter (5a uses 3mm to fit 75mm panels in 183mm page).
        row_gutter_mm = ROW_GUTTER_MM.get(row_id, default_gutter_mm)
        gutter_pt = mm_to_pt(row_gutter_mm)
        n_cols = len(row_paths)
        cell_w = (inner_w - (n_cols - 1) * gutter_pt) / n_cols
        log.info(
            "row %s: %d cells, cell %.1f×%.1f pt (gutter=%.1fmm, top=%.1f)",
            row_id, n_cols, cell_w, row_h_pt, row_gutter_mm, cur_top,
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
    log.info("Wrote vector-composed PDF: %s", out)


def main() -> int:
    args = parse_args()
    grid = collect_panels(args.panel_dir)
    log.info("Collected panels: %d rows", len(grid))
    compose(
        grid, args.out,
        args.page_width_mm, args.page_height_mm, args.gutter_mm, args.border_mm,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
