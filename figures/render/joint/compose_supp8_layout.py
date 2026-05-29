"""compose_supp8_layout.py — assemble Supp 8 panels into a 3×3 grid PDF
(rows = compartments, cols = UMAP | confusion | heatmap).

Each row presents one compartment's annotation evidence:
- Col A: joint embedding labeled by L1.5 cascade label (8a)
- Col B: per-compartment FLEX-label × joint-label confusion matrix (8b) —
  the cross-platform alignment evidence
- Col C: per-compartment canonical+artifact marker heatmap (8c) —
  the artifact-state validation evidence

pypdf preserves vector content of each panel via merge_transformed_page
with a scale + translate transformation (the target container does not
have pikepdf installed). Pattern:
publication/figures/render/flex/compose_supp5_layout.py:117-174 (compose
function), adapted to pypdf's Transformation API and a 3×3 grid with
explicit per-cell placement.

Outputs:
    supp8_layout.pdf
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import List, Tuple

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("compose_supp8")


MM_PER_INCH = 25.4
PT_PER_INCH = 72.0
PAGE_WIDTH_MM = 183.0
PAGE_HEIGHT_MM = 247.0

# Column widths (left-to-right): UMAP, confusion, heatmap.
# Sized to match each panel type's intrinsic aspect ratio at row_height = 55,
# minimizing letterbox whitespace after aspect-preserve scaling in compose.
# Heatmaps render at ~1.0 aspect (anisotropic cells in render_supp8c) so
# the heatmap cell is roughly square; UMAPs are ~1.2 wide; confusion is
# ~1.0 (Epi/Imm) or taller (Str).
# Page math: 60 + 55 + 60 = 175; +2*2 gutters +2*2 borders = 183 = page width.
# Tighter gutters (2mm vs 5mm) buy ~10% larger cells while staying on-page.
# UMAP column expanded to 60 (was 55) by trimming confusion (60→55).
COL_WIDTHS_MM = [60.0, 55.0, 60.0]
COL_GUTTER_MM = 2.0
BORDER_MM = 2.0
# 3 × 60 + 2 × 2 gutter + 2 × 2 border = 188 mm (within 247).
ROW_HEIGHT_MM = 60.0
ROW_GUTTER_MM = 2.0


COMPARTMENTS = [
    ("epithelial", "Epithelial"),
    ("immune", "Immune"),
    ("stromal", "Stromal"),
]


def panel_filename(col: int, compartment_dir: str) -> str:
    """col 0 = UMAP, 1 = confusion, 2 = heatmap."""
    return [
        f"supp8a_joint_umap_{compartment_dir}.pdf",
        f"supp8b_confusion_{compartment_dir}.pdf",
        f"supp8c_heatmap_{compartment_dir}.pdf",
    ][col]


def mm_to_pt(mm: float) -> float:
    return mm / MM_PER_INCH * PT_PER_INCH


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--panel-dir", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    return p.parse_args()


def collect_panels(panel_dir: Path) -> List[List[Path]]:
    """Returns 3×3 grid of panel paths (rows = compartments, cols = panels)."""
    grid: List[List[Path]] = []
    for compartment_dir, _comp_full in COMPARTMENTS:
        row: List[Path] = []
        for col in range(3):
            p = panel_dir / panel_filename(col, compartment_dir)
            if not p.exists():
                sys.exit(f"FATAL: missing panel {p}")
            row.append(p)
        grid.append(row)
    return grid


def compose(grid: List[List[Path]], out: Path) -> None:
    """Compose the 3×3 grid using pypdf (pikepdf is not available in the
    target container; pypdf is). pypdf preserves vector content via
    merge_transformed_page with a scale + translate transformation.
    """
    from pypdf import PdfReader, PdfWriter, Transformation  # type: ignore

    writer = PdfWriter()
    page_w_pt = mm_to_pt(PAGE_WIDTH_MM)
    page_h_pt = mm_to_pt(PAGE_HEIGHT_MM)
    out_page = writer.add_blank_page(width=page_w_pt, height=page_h_pt)

    for row_idx, row_panels in enumerate(grid):
        y_top_mm = BORDER_MM + row_idx * (ROW_HEIGHT_MM + ROW_GUTTER_MM)
        x_mm = BORDER_MM
        for col_idx, panel_path in enumerate(row_panels):
            w_mm = COL_WIDTHS_MM[col_idx]
            h_mm = ROW_HEIGHT_MM

            reader = PdfReader(str(panel_path))
            src_page = reader.pages[0]
            mb = src_page.mediabox
            src_w_pt = float(mb.right) - float(mb.left)
            src_h_pt = float(mb.top) - float(mb.bottom)

            cell_w_pt = mm_to_pt(w_mm)
            cell_h_pt = mm_to_pt(h_mm)
            scale = min(cell_w_pt / src_w_pt, cell_h_pt / src_h_pt)
            placed_w_pt = src_w_pt * scale
            placed_h_pt = src_h_pt * scale

            # PDF coordinate origin is bottom-left.
            x_pt = mm_to_pt(x_mm) + (cell_w_pt - placed_w_pt) / 2.0
            y_top_pt = page_h_pt - mm_to_pt(y_top_mm)
            y_bottom_pt = y_top_pt - cell_h_pt + (cell_h_pt - placed_h_pt) / 2.0
            # Translation accounts for the source mediabox origin offset
            # (mediabox.left/bottom may be non-zero).
            tx = x_pt - float(mb.left) * scale
            ty = y_bottom_pt - float(mb.bottom) * scale

            op = Transformation().scale(scale).translate(tx, ty)
            out_page.merge_transformed_page(src_page, op)
            log.info(
                "  placed %s → %.1f×%.1fmm at x=%.1fmm, y_top=%.1fmm",
                panel_path.name,
                placed_w_pt / PT_PER_INCH * MM_PER_INCH,
                placed_h_pt / PT_PER_INCH * MM_PER_INCH,
                x_mm, y_top_mm,
            )
            x_mm += w_mm + COL_GUTTER_MM

    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "wb") as fh:
        writer.write(fh)
    log.info("Wrote vector-composed PDF: %s", out)


def main() -> int:
    args = parse_args()
    grid = collect_panels(args.panel_dir)
    log.info("Collected %d×%d panel grid", len(grid), len(grid[0]))
    compose(grid, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
