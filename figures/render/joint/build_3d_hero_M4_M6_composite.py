"""3d main panel — Pat2_P3_A_UOQ view_88055 M4+M6 signature + intensity.

2 rows (M4, M6) × 2 cols (signature carriers, intensity score). Composes
pre-rendered per-FOV checklist panels into a single hero figure.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.image import imread

SHARE = ("/Users/NoahWechter/Library/Application Support/CRSP Desktop/"
         "Volumes.noindex/CRSP Lab - dalawson.localized/nwechter")
PREVIEW = (f"{SHARE}/iHBCA_publication/coordination/staging/"
           "fig3_promotion_preview_20260520")
SRC = (f"{PREVIEW}/figures/output/fig3/supplemental/distribution/"
       "per_fov_checklist")
OUT = f"{PREVIEW}/figures/output/fig3/main/3d_hero_M4_M6.pdf"

PREFIX = "Pat2_P3_A_UOQ_xenium_1__view_88055_L1-5_types__"

PANELS = {
    ("M4", "signature"): f"{SRC}/{PREFIX}signature_M4.png",
    ("M4", "intensity"): f"{SRC}/{PREFIX}intensity_M4.png",
    ("M6", "signature"): f"{SRC}/{PREFIX}signature_M6.png",
    ("M6", "intensity"): f"{SRC}/{PREFIX}intensity_M6.png",
}

ROWS = ["M4", "M6"]
COLS = ["signature", "intensity"]

fig, axes = plt.subplots(len(ROWS), len(COLS),
                          figsize=(len(COLS) * 3.0, len(ROWS) * 3.0),
                          squeeze=False)
for i, motif in enumerate(ROWS):
    for j, layer in enumerate(COLS):
        ax = axes[i, j]
        path = PANELS[(motif, layer)]
        if not os.path.exists(path):
            ax.text(0.5, 0.5, f"missing\n{os.path.basename(path)}",
                    transform=ax.transAxes, fontsize=6, ha="center",
                    va="center", color="red")
        else:
            img = imread(path)
            ax.imshow(img)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        if i == 0:
            ax.set_title(layer, fontsize=8, fontweight="bold", pad=4)
        if j == 0:
            ax.set_ylabel(motif, fontsize=8, fontweight="bold",
                           rotation=0, ha="right", va="center", labelpad=8)

fig.suptitle("Pat2_P3_A_UOQ view 88055 — M4 + M6",
              fontsize=9, y=0.995)
plt.tight_layout(rect=(0.02, 0.01, 0.99, 0.97), h_pad=0.2, w_pad=0.2)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT, format="pdf", bbox_inches="tight")
fig.savefig(OUT.replace(".pdf", ".png"), dpi=600, bbox_inches="tight")
plt.close(fig)
print(f"wrote {OUT}")
