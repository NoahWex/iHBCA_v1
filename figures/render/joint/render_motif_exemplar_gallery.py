#!/usr/bin/env python3
"""Motif exemplar gallery v6 — PyMuPDF (fitz) composite, vector preserved.

Layout: 3 cols × 3 motif-groups; each group = title + signature row + intensity row.
Letter-height page (11in), width auto-derived from cell aspect.
"""
import os
import fitz

BASE = ("/Users/NoahWechter/Library/Application Support/CRSP Desktop/Volumes.noindex/"
        "CRSP Lab - dalawson.localized/nwechter/iHBCA_publication/"
        "coordination/staging/fig3_promotion_preview_20260520/figures/output/fig3/"
        "motif_exemplars_v1")

EXEMPLARS = [
    ("M0", "Fb_Activated", "Pat2_P3_P_LOQ",    f"{BASE}/m0/Pat2_P3_P_LOQ_xenium_1__view_10987_landscape_zoomout25"),
    ("M1", "EC",           "Pat2_P3_M_UOQ",    f"{BASE}/m1/Pat2_P3_M_UOQ_xenium_1__view_9207_landscape_zoomout25"),
    ("M2", "Fb",           "UCI220228_P2_LIQ", f"{BASE}/m2/UCI220228_P2_LIQ_xenium_1__view_9860_landscape_zoomout25"),
    ("M3", "PV",           "UCI604_P3_A_LOQ",  f"{BASE}/m3/UCI604_P3_A_LOQ_xenium_1__M3_unit000_landscape_zoomout25"),
    ("M4", "T+DC",         "Pat1_P1",          f"{BASE}/m4/Pat1_P1_xenium_1__view_7870_zoomout25"),
    ("M5", "LEC",          "Pat2_P1_xenium_2", f"{BASE}/m5/Pat2_P1_xenium_2__view_12457_M5_landscape_zoomout25"),
    ("M6", "B-cell pocket","UCI604_P3_A_UOQ",  f"{BASE}/m6/UCI604_P3_A_UOQ_xenium_1__view_12892_M6_landscape_zoomout25"),
    ("M7", "Mac",          "Pat2_P3_M_LOQ",    f"{BASE}/m7/Pat2_P3_M_LOQ_xenium_1__view_13242_M7_landscape_zoomout25"),
    ("M8", "Plas",         "UCI220228_P2_UIQ", f"{BASE}/m8/UCI220228_P2_UIQ_xenium_1__view_13414_M8_landscape_zoomout25"),
]

for m, _, _, stem in EXEMPLARS:
    for kind in ("signature", "intensity"):
        p = f"{stem}__{kind}_{m}.pdf"
        if not os.path.exists(p):
            raise SystemExit(f"MISSING: {p}")
print("All 18 source PDFs present.")

IN = 72.0
PAGE_H = 11.0 * IN
M_TOP, M_BOT, M_LEFT, M_RIGHT = 0.15*IN, 0.15*IN, 0.15*IN, 0.15*IN
N_GROUPS, N_COLS = 3, 3
H_GAP      = 0.04 * IN
V_GAP_PAIR = 0.06 * IN
V_GAP_GRP  = 0.28 * IN
TITLE_H    = 0.32 * IN
SRC_ASPECT = 0.753

avail_h = PAGE_H - M_TOP - M_BOT
group_h = (avail_h - (N_GROUPS - 1) * V_GAP_GRP) / N_GROUPS
CELL_H  = (group_h - TITLE_H - V_GAP_PAIR) / 2
CELL_W  = CELL_H / SRC_ASPECT
PAGE_W  = M_LEFT + M_RIGHT + N_COLS * CELL_W + (N_COLS - 1) * H_GAP

print(f"Page {PAGE_W/IN:.2f}×{PAGE_H/IN:.2f}in; cell {CELL_W/IN:.2f}×{CELL_H/IN:.2f}in")

out = fitz.open()
page = out.new_page(width=PAGE_W, height=PAGE_H)
page.draw_rect(fitz.Rect(0, 0, PAGE_W, PAGE_H), color=(0,0,0), fill=(0,0,0))

for g in range(N_GROUPS):
    group_top_y = M_TOP + g * (group_h + V_GAP_GRP)
    title_top_y = group_top_y
    sig_top_y   = group_top_y + TITLE_H
    int_top_y   = sig_top_y + CELL_H + V_GAP_PAIR
    for col in range(N_COLS):
        idx = g * N_COLS + col
        m, label, sample, stem = EXEMPLARS[idx]
        cell_left = M_LEFT + col * (CELL_W + H_GAP)
        line1 = f"{m}  {label}"
        line2 = sample
        l1_w = fitz.get_text_length(line1, fontname="hebo", fontsize=8)
        l2_w = fitz.get_text_length(line2, fontname="helv", fontsize=6)
        page.insert_text(
            fitz.Point(cell_left + (CELL_W - l1_w)/2, title_top_y + 9),
            line1, fontname="hebo", fontsize=8, color=(1,1,1),
        )
        page.insert_text(
            fitz.Point(cell_left + (CELL_W - l2_w)/2, title_top_y + 18),
            line2, fontname="helv", fontsize=6, color=(1,1,1),
        )
        for kind, top_y in (("signature", sig_top_y), ("intensity", int_top_y)):
            src_path = f"{stem}__{kind}_{m}.pdf"
            rect = fitz.Rect(cell_left, top_y, cell_left + CELL_W, top_y + CELL_H)
            src = fitz.open(src_path)
            page.show_pdf_page(rect, src, 0, keep_proportion=True)
            src.close()

OUT = "/tmp/s_motif_exemplar_gallery_v6.pdf"
out.save(OUT)
out.close()
print(f"Wrote {OUT}")
