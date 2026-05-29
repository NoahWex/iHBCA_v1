"""
Standalone renderer for S7 panel b: LASP-override per-marker row-z-score
heatmap of all res-5.0 LASP-region sub-clusters (blocked by majority_L2,
secondary-ordered by combined LASP-marker signal within block).

Extracts the heatmap geom from `render_s7_composite.py:render_panel_b_into`
into a first-class standalone script. Thin wrapper — no logic duplication.

CLI:
  --substrate   Path to s7b_lasp_violin_long.csv[.gz] (panel-b long-form CSV;
                gzipped variant is canonical at publication/figures/data/fig1/lasp_override/).
  --out-dir     Destination directory. Writes:
                  s_ihbca_lasp_subcluster_heatmap.pdf
                  s_ihbca_lasp_subcluster_heatmap.png
                  s_ihbca_lasp_subcluster_heatmap.panel.yaml

Local run (seconds):
  python render_s7_panel_b_heatmap.py \\
    --substrate publication/figures/data/fig1/lasp_override/s7b_lasp_violin_long.csv.gz \\
    --out-dir   ./output
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import yaml

# Ensure co-located modules are importable regardless of CWD.
sys.path.insert(0, str(Path(__file__).resolve().parent))

import render_s7_composite as composite  # noqa: E402 — render_panel_b_into + dims


MM_PER_IN = 25.4


PANEL_B_HEIGHT_MM = composite.PANEL_B_HEIGHT_MM   # 70.0 — matches composite envelope
COMPOSITE_WIDTH_MM = composite.COMPOSITE_WIDTH_MM # 183.0 — matches composite width

# Standalone figure dimensions: same width × panel-b height; small extra padding
# below to fit the inset colorbar that sits beneath the heatmap.
FIG_WIDTH_MM = COMPOSITE_WIDTH_MM
FIG_HEIGHT_MM = PANEL_B_HEIGHT_MM + 10.0  # +10mm for colorbar inset + bottom xlabel


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--substrate", required=True, type=Path,
        help="Path to s7b_lasp_violin_long.csv (cell_id, leiden_5.0, majority_L2, marker_gene, log1p_expression)",
    )
    p.add_argument(
        "--out-dir", required=True, type=Path,
        help="Destination directory for PDF + PNG + panel.yaml.",
    )
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading substrate: {args.substrate}", flush=True)
    df_long = pd.read_csv(args.substrate)
    print(f"  rows: {len(df_long):,}; "
          f"markers: {df_long['marker_gene'].nunique()}; "
          f"sub-clusters: {df_long['leiden_5.0'].nunique()}",
          flush=True)

    # Figure dims in inches; aligns to the 70mm composite envelope so the
    # standalone reads identically to the panel-b slot in the composite.
    fig_w_in = FIG_WIDTH_MM / MM_PER_IN
    fig_h_in = FIG_HEIGHT_MM / MM_PER_IN

    fig, ax = plt.subplots(figsize=(fig_w_in, fig_h_in))
    # Margins picked to leave room for block labels above, colorbar inset below,
    # and the cluster-count xlabel below the heatmap.
    fig.subplots_adjust(left=0.06, right=0.985, top=0.84, bottom=0.30)

    meta = composite.render_panel_b_into(fig, ax, df_long)
    print(f"  rendered: {meta['n_rows']} markers × {meta['n_cols']} sub-clusters; "
          f"row z range [{meta['row_zscore_range'][0]:.2f}, {meta['row_zscore_range'][1]:.2f}]",
          flush=True)

    panel_id = "s_ihbca_lasp_subcluster_heatmap"
    pdf_path = args.out_dir / f"{panel_id}.pdf"
    png_path = args.out_dir / f"{panel_id}.png"
    yaml_path = args.out_dir / f"{panel_id}.panel.yaml"

    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, bbox_inches="tight", dpi=600)
    plt.close(fig)
    print(f"Wrote {pdf_path} ({pdf_path.stat().st_size:,} bytes)")
    print(f"Wrote {png_path} ({png_path.stat().st_size:,} bytes)")

    panel_meta = {
        "panel_id": panel_id,
        "source": "render_s7_panel_b_heatmap.py (extracted from render_s7_composite.py:render_panel_b_into)",
        "substrate": str(args.substrate),
        "n_rows_markers": meta["n_rows"],
        "n_cols_subclusters": meta["n_cols"],
        "block_counts": meta["block_counts"],
        "row_zscore_range": meta["row_zscore_range"],
        "layout": meta["layout"],
        "dims_mm": {"width": FIG_WIDTH_MM, "height": FIG_HEIGHT_MM},
    }
    with open(yaml_path, "w") as fh:
        yaml.safe_dump(panel_meta, fh, sort_keys=False)
    print(f"Wrote {yaml_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
