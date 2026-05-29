#!/usr/bin/env python3
"""
build_fig3_flex_l2s_markers.py — assemble Fig 3 supplemental table from the
FLEX Spatial HBCA L2S (spatial behavioral group) marker catalog.

Source: Spatial_HBCA/project/05_Analyses/outputs/behavioral_group_characterization_20260318/
        results/markers_allvsall/allvsall_markers.csv

Method: presto::wilcoxauc on Seurat object with Idents = l2_label (40 L2S
groups across the FLEX 258K-cell dataset). One-vs-rest within the full
vocabulary. Source script: behavioral_group_characterization_20260318/
scripts/m1d_findmarkers.Rmd (chunk: allvsall_markers).

Output schema (column rename for reviewer audience, parallel to T3 L1 markers):
  gene_symbol, cell_type, label_level, avg_log_expression, log_fold_change,
  wilcox_statistic, auROC, p_val, FDR, pct_expressed_in, pct_expressed_out
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


COLUMN_MAP = {
    "feature":   "gene_symbol",
    "group":     "cell_type",
    "avgExpr":   "avg_log_expression",
    "logFC":     "log_fold_change",   # presto: natural-log effect size on log1p-normalized data
    "statistic": "wilcox_statistic",
    "auc":       "auROC",
    "pval":      "p_val",
    "padj":      "FDR",
    "pct_in":    "pct_expressed_in",
    "pct_out":   "pct_expressed_out",
}


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--src", required=True, type=Path,
                   help="Spatial_HBCA allvsall_markers.csv")
    p.add_argument("--out-csv", required=True, type=Path)
    args = p.parse_args()

    print(f"[fig3-l2s] loading {args.src}")
    df = pd.read_csv(args.src)
    print(f"[fig3-l2s]   rows: {len(df):,}  unique groups: {df['group'].nunique()}")

    df = df.rename(columns=COLUMN_MAP)
    df.insert(2, "label_level", "L2S")

    # Convert percentages so 0-100 -> 0-1 for parity with T3 (L1 markers used
    # pct_in as fraction 0-1, presto reports percent 0-100). Halve memory by
    # using float32 for the percent columns.
    for col in ["pct_expressed_in", "pct_expressed_out"]:
        df[col] = df[col].astype(float) / 100.0

    # Final column order
    cols = ["gene_symbol", "cell_type", "label_level",
            "avg_log_expression", "log_fold_change",
            "wilcox_statistic", "auROC",
            "p_val", "FDR",
            "pct_expressed_in", "pct_expressed_out"]
    df = df[cols].copy()

    df = df.sort_values(["cell_type", "FDR", "p_val"]).reset_index(drop=True)

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out_csv, index=False)
    print(f"[fig3-l2s] wrote {args.out_csv} ({len(df):,} rows × {len(df.columns)} cols)")
    print(f"[fig3-l2s] groups present ({df['cell_type'].nunique()}):")
    for ct, n in df.groupby("cell_type").size().items():
        print(f"           {ct:40s} {n:>5}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
