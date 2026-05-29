#!/usr/bin/env python3
"""
build_t4_brca1_de.py — assemble Fig 2 supplemental table T4 from per-cell-type
limma-voom BRCA1-vs-non-carrier DE result CSVs.

Two cohort × formula slices are pooled (long-format), distinguished by the
`cohort_id` and `formula_id` columns already present in each source CSV:

  - C_tested / A2_sampletype_tested  (4 fibroblast L2 subtypes; matches Fig 2 c)
  - A_tested / A2_tested             (22 L1 cell types; broader BRCA1-tested cohort)

Filter: raw P.Value <= 0.1 (mirrors T2 convention for keeping file size finite).

Output schema:
  gene_symbol, ensembl_id, cell_type, label_level, cohort_id, formula_id,
  log2FC, AveExpr, t_statistic, p_val, FDR, B_statistic, n_AR, n_BR1

(`label_level` is L2 for the 4 fibro subtypes, L1 for the broader cohort.)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


# Per-cell file -> label_level
L2_FILES = {"str__Fibro_IGF1", "str__Fibro_major",
            "str__Fibro_prematrix", "str__Fibro_SFRP4"}


def load_dir(dirpath: Path, label_level: str) -> pd.DataFrame:
    rows = []
    for p in sorted(dirpath.glob("*.csv")):
        df = pd.read_csv(p)
        df["label_level"] = label_level
        rows.append(df)
        print(f"[t4]   {p.name}: {len(df):,} rows")
    return pd.concat(rows, ignore_index=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--de-root", required=True, type=Path,
                        help="iHBCA_V1 risk_l2_voom_20260519/outputs/de_results")
    parser.add_argument("--out-csv", required=True, type=Path)
    parser.add_argument("--p-filter", type=float, default=0.1,
                        help="Drop rows with raw P.Value > this (default 0.1)")
    args = parser.parse_args()

    c_tested = args.de_root / "C_tested" / "A2_sampletype_tested"
    a_tested = args.de_root / "A_tested" / "A2_tested"
    for d in (c_tested, a_tested):
        if not d.is_dir():
            sys.exit(f"missing dir: {d}")

    print(f"[t4] loading C_tested/A2_sampletype_tested (4 fibro L2):")
    c = load_dir(c_tested, "L2")
    print(f"[t4]   subtotal: {len(c):,} rows  unique cell_types: {c['L2'].nunique()}")

    print(f"[t4] loading A_tested/A2_tested (L1 cell types):")
    a = load_dir(a_tested, "L1")
    print(f"[t4]   subtotal: {len(a):,} rows  unique cell_types: {a['L2'].nunique()}")

    df = pd.concat([c, a], ignore_index=True)
    print(f"[t4] combined rows: {len(df):,}")

    before = len(df)
    df = df[df["P.Value"] <= args.p_filter].copy()
    print(f"[t4] after P.Value <= {args.p_filter} filter: {len(df):,} "
          f"(dropped {before - len(df):,})")

    # Drop duplicate aliases (Austin's pipeline carried FDR=adj.P.Val + PValue=P.Value as copies)
    df = df.drop(columns=["FDR", "PValue"], errors="ignore")
    # Rename to canonical schema
    df = df.rename(columns={
        "symbol":    "gene_symbol",
        "gene_id":   "ensembl_id",
        "logFC":     "log2FC",
        "t":         "t_statistic",
        "P.Value":   "p_val",
        "adj.P.Val": "FDR",
        "B":         "B_statistic",
        "L2":        "cell_type",
    })

    # Final column order
    cols = ["gene_symbol", "ensembl_id", "cell_type", "label_level",
            "cohort_id", "formula_id",
            "log2FC", "AveExpr", "t_statistic",
            "p_val", "FDR", "B_statistic",
            "n_AR", "n_BR1"]
    cols = [c for c in cols if c in df.columns]
    df = df[cols].copy()

    # Sort: cell type, then FDR
    df = df.sort_values(["label_level", "cell_type", "FDR"]).reset_index(drop=True)

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out_csv, index=False)
    print(f"[t4] wrote {args.out_csv} ({len(df):,} rows × {len(df.columns)} cols)")
    print(f"[t4] cell types present:")
    for (lvl, ct), n in df.groupby(["label_level", "cell_type"]).size().items():
        print(f"           {lvl:3s}  {ct:30s} {n:>6}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
