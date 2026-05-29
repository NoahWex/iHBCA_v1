#!/usr/bin/env python3
"""Compose the M6-only 1x2 hero panel (signature + intensity).

Sibling of build_3d_hero_M4_M6_composite.py — same aesthetics, just one
motif and one row instead of two. Reads the single-panel signature_M6
+ intensity_M6 PDFs produced by render_from_saved_views.py against
saved_views_3d_hero_m6.yaml, composes them side-by-side at the locked
cell dimensions (4.50 x 5.98 in), no bbox_inches='tight' to preserve
FOV alignment. PyMuPDF vector composite preserves rasterized polygon
layers + vector text.
"""
import argparse
import sys
from pathlib import Path

import fitz  # PyMuPDF


CELL_W_IN = 4.50
CELL_H_IN = 5.98
H_GAP_IN = 0.04
M_TOP = 0.05
M_BOT = 0.05
M_LEFT = 0.05
M_RIGHT = 0.05


def in_to_pt(inches: float) -> float:
    return inches * 72.0


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-dir", required=True, type=Path,
                    help="Dir containing signature_M6.pdf + intensity_M6.pdf")
    ap.add_argument("--out", required=True, type=Path,
                    help="Output composite PDF path")
    return ap.parse_args()


def main() -> int:
    args = parse_args()

    # Discover the two source PDFs by suffix
    sig = next(args.in_dir.glob("*__signature_M6.pdf"), None)
    intn = next(args.in_dir.glob("*__intensity_M6.pdf"), None)
    if sig is None or intn is None:
        print(f"ERROR: missing source PDFs in {args.in_dir}", file=sys.stderr)
        print(f"  signature_M6.pdf: {sig}", file=sys.stderr)
        print(f"  intensity_M6.pdf: {intn}", file=sys.stderr)
        return 1
    print(f"signature: {sig.name}")
    print(f"intensity: {intn.name}")

    page_w_in = M_LEFT + CELL_W_IN + H_GAP_IN + CELL_W_IN + M_RIGHT
    page_h_in = M_TOP + CELL_H_IN + M_BOT
    page_w_pt = in_to_pt(page_w_in)
    page_h_pt = in_to_pt(page_h_in)
    cell_w_pt = in_to_pt(CELL_W_IN)
    cell_h_pt = in_to_pt(CELL_H_IN)

    doc = fitz.open()
    page = doc.new_page(width=page_w_pt, height=page_h_pt)

    # Black backdrop matches the v6 gallery convention
    page.draw_rect(fitz.Rect(0, 0, page_w_pt, page_h_pt),
                   color=(0, 0, 0), fill=(0, 0, 0), overlay=False)

    # Cell positions (PDF coords: origin top-left)
    cell_top = in_to_pt(M_TOP)
    cell_bot = cell_top + cell_h_pt
    sig_left = in_to_pt(M_LEFT)
    sig_right = sig_left + cell_w_pt
    int_left = sig_right + in_to_pt(H_GAP_IN)
    int_right = int_left + cell_w_pt

    sig_rect = fitz.Rect(sig_left, cell_top, sig_right, cell_bot)
    int_rect = fitz.Rect(int_left, cell_top, int_right, cell_bot)

    sig_doc = fitz.open(sig)
    int_doc = fitz.open(intn)
    page.show_pdf_page(sig_rect, sig_doc, 0, keep_proportion=True)
    page.show_pdf_page(int_rect, int_doc, 0, keep_proportion=True)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    doc.save(str(args.out))
    print(f"wrote {args.out}")
    print(f"page: {page_w_in:.2f} x {page_h_in:.2f} in")
    return 0


if __name__ == "__main__":
    sys.exit(main())
