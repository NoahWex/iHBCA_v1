#!/usr/bin/env python3
"""Build the cell -> L1 alignment CSV used by the L1-resolution pseudobulk run.

Mirror of 00_build_cell_to_l2.py but for L1: L1_label = {compartment}__{lineage}
from labels_full.csv. Excludes is_artifact=True cells (same gate as L2 build).

Output columns expected by pseudobulk_aggregate_full_object.py:
    numeric_id, L1_label
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def load_compartment_meta(meta_csv: str, comp_label: str) -> pd.DataFrame:
    cols = ["global_numeric_id", "cell_id", "patientID", "dataset"]
    df = pd.read_csv(meta_csv, low_memory=False, usecols=cols)
    df["compartment"] = comp_label
    print(f"  {comp_label}: {len(df):,} cells", flush=True)
    return df


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--epi-meta", required=True)
    p.add_argument("--imm-meta", required=True)
    p.add_argument("--str-meta", required=True)
    p.add_argument("--labels-full", required=True)
    p.add_argument("--output-csv", required=True,
                   help="Path for cell_to_l1.csv (numeric_id, L1_label)")
    p.add_argument("--donor-summary-csv", required=True)
    args = p.parse_args()

    print("=== build_cell_to_l1 (labels_full.csv source) ===", flush=True)

    # ---- Step 1: load per-compartment meta ----
    print("[step 1] load per-compartment enriched metadata", flush=True)
    cells = pd.concat([
        load_compartment_meta(args.epi_meta, "epi"),
        load_compartment_meta(args.imm_meta, "imm"),
        load_compartment_meta(args.str_meta, "str"),
    ], ignore_index=True)
    print(f"  total cells: {len(cells):,}", flush=True)

    # ---- Step 2: load labels_full.csv (lineage = L1) ----
    print("[step 2] load labels_full.csv (canonical L1)", flush=True)
    labels = pd.read_csv(args.labels_full, low_memory=False,
                         usecols=["cell_id", "compartment", "lineage", "is_artifact"])
    print(f"  labels_full: {len(labels):,} cell entries", flush=True)
    labels["is_artifact_bool"] = labels["is_artifact"].astype(str).str.lower().isin(["true", "1"])
    labels = labels[~labels["is_artifact_bool"]].copy()
    print(f"  after artifact drop: {len(labels):,} cells", flush=True)

    # Compose L1_label = {compartment}__{lineage}
    labels["L1_label"] = (labels["compartment"].astype(str) + "__"
                          + labels["lineage"].astype(str)
                              .str.replace(" ", "_", regex=False)
                              .str.replace("-", "_", regex=False))
    n_l1 = labels["L1_label"].nunique()
    print(f"  unique L1_label values: {n_l1}", flush=True)
    top_l1 = labels["L1_label"].value_counts()
    print("  L1_label counts (all):", flush=True)
    for k, v in top_l1.items():
        print(f"    {k}: {v:,}", flush=True)

    # ---- Step 3: join labels -> per-compartment meta on cell_id ----
    print("[step 3] join labels_full -> per-compartment metadata on cell_id", flush=True)
    aligned = cells.merge(
        labels[["cell_id", "L1_label"]],
        on="cell_id", how="left"
    )
    n_with = aligned["L1_label"].notna().sum()
    print(f"  cells with L1_label after join: {n_with:,}/{len(aligned):,} "
          f"({100*n_with/len(aligned):.1f}%)", flush=True)

    # ---- Step 4: write cell_to_l1.csv ----
    print("[step 4] write cell_to_l1.csv", flush=True)
    out_align = aligned.dropna(subset=["L1_label"])[["global_numeric_id", "L1_label"]]
    out_align = out_align.rename(columns={"global_numeric_id": "numeric_id"})
    Path(args.output_csv).parent.mkdir(parents=True, exist_ok=True)
    out_align.to_csv(args.output_csv, index=False)
    print(f"  wrote: {args.output_csv} ({len(out_align):,} rows)", flush=True)

    # ---- Step 5: per-donor x L1 cell counts ----
    print("[step 5] per-donor x L1 cell counts", flush=True)
    labeled = aligned.dropna(subset=["L1_label"])
    summary = (labeled.groupby(["L1_label", "patientID"], dropna=True)
                      .size()
                      .reset_index(name="n_cells"))
    summary.to_csv(args.donor_summary_csv, index=False)
    print(f"  wrote: {args.donor_summary_csv} ({len(summary):,} rows)", flush=True)

    print("\nDone.", flush=True)


if __name__ == "__main__":
    sys.exit(main())
