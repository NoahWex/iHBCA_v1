#!/usr/bin/env python3
"""Build the cell -> L2 alignment CSV used by the pseudobulk aggregator.

Label source: labels_full.csv (canonical L2.0c, ~25 types, full atlas).
Join key: cell_id (between labels_full and per-compartment enriched meta).

Why labels_full.csv and not per-compartment `level2`:
    The `level2` column in the enriched per-compartment metadata is a finer
    upstream annotation (~119 unique values, covers only 46.7% of cells).
    The canonical L2.0c vocabulary used by risk_main_effects_20260507 and
    the Fig 2 promotion lives in labels_full.csv. That's the right source
    for publication-aligned DE.

Output columns expected by pseudobulk_aggregate_full_object.py:
    numeric_id, L2_label
where numeric_id holds the global_numeric_id value (the script's join key)
and L2_label is `{compartment}__{label}` (matching the existing
risk_main_effects naming convention: epi__BMYO_basal, str__Fibro_SFRP4).

Also: writes a sample_type_coarse heterogeneity report for the 51 donors
with >1 value across their cells (per CP1.5 investigation).
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def load_compartment_meta(meta_csv: str, comp_label: str) -> pd.DataFrame:
    cols = ["global_numeric_id", "cell_id", "patientID", "dataset",
            "sample_type_coarse", "sample_type"]
    df = pd.read_csv(meta_csv, low_memory=False, usecols=cols)
    df["compartment"] = comp_label
    print(f"  {comp_label}: {len(df):,} cells", flush=True)
    return df


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--epi-meta", required=True)
    p.add_argument("--imm-meta", required=True)
    p.add_argument("--str-meta", required=True)
    p.add_argument("--labels-full", required=True,
                   help="iHBCA_publication labels_full.csv (cell_id, compartment, label, is_artifact, ...)")
    p.add_argument("--output-csv", required=True,
                   help="Path for cell_to_l2.csv (numeric_id, L2_label)")
    p.add_argument("--donor-summary-csv", required=True)
    p.add_argument("--donor-meta-out", required=True)
    p.add_argument("--existing-donor-meta", required=True)
    p.add_argument("--stc-heterogeneity-csv", required=True,
                   help="Per-donor sample_type_coarse value-count table for donors with >1 STC value")
    args = p.parse_args()

    print("=== build_cell_to_l2 (labels_full.csv source) ===", flush=True)

    # ---- Step 1: load per-compartment meta (small subset of cols) ----
    print("[step 1] load per-compartment enriched metadata", flush=True)
    dfs = [
        load_compartment_meta(args.epi_meta, "epi"),
        load_compartment_meta(args.imm_meta, "imm"),
        load_compartment_meta(args.str_meta, "str"),
    ]
    cells = pd.concat(dfs, ignore_index=True)
    print(f"  total cells: {len(cells):,}", flush=True)

    # ---- Step 2: load labels_full.csv ----
    print("[step 2] load labels_full.csv (canonical L2.0c)", flush=True)
    labels = pd.read_csv(args.labels_full, low_memory=False,
                          usecols=["cell_id", "compartment", "label", "is_artifact"])
    print(f"  labels_full: {len(labels):,} cell entries", flush=True)
    print(f"  unique labels: {labels['label'].nunique()}", flush=True)
    n_artifact = labels["is_artifact"].astype(str).str.lower().isin(["true", "1"]).sum()
    print(f"  is_artifact=True cells: {n_artifact:,} (will be excluded)", flush=True)

    # Normalize is_artifact -> boolean for safe filtering
    labels["is_artifact_bool"] = labels["is_artifact"].astype(str).str.lower().isin(["true", "1"])
    labels = labels[~labels["is_artifact_bool"]].copy()
    print(f"  after artifact drop: {len(labels):,} cells", flush=True)

    # Compose L2_label = {compartment}__{label}, matching existing naming
    # convention (e.g., str__Fibro_SFRP4) used by risk_main_effects F3 markers.
    labels["L2_label"] = (labels["compartment"].astype(str) + "__"
                          + labels["label"].astype(str).str.replace(" ", "_", regex=False)
                                                       .str.replace("-", "_", regex=False))
    print(f"  unique L2_label values: {labels['L2_label'].nunique()}", flush=True)
    top_l2 = labels["L2_label"].value_counts().head(15)
    print("  top 15 L2_label values:", flush=True)
    for k, v in top_l2.items():
        print(f"    {k}: {v:,}", flush=True)

    # ---- Step 3: join labels_full -> per-compartment meta on cell_id ----
    print("[step 3] join labels_full -> per-compartment metadata on cell_id", flush=True)
    aligned = cells.merge(
        labels[["cell_id", "L2_label"]],
        on="cell_id", how="left"
    )
    n_with_label = aligned["L2_label"].notna().sum()
    print(f"  cells with L2_label after join: {n_with_label:,}/{len(aligned):,} "
          f"({100*n_with_label/len(aligned):.1f}%)", flush=True)

    # ---- Step 4: write cell_to_l2.csv ----
    print("[step 4] write cell_to_l2.csv", flush=True)
    out_align = aligned.dropna(subset=["L2_label"])[["global_numeric_id", "L2_label"]]
    out_align = out_align.rename(columns={"global_numeric_id": "numeric_id"})
    Path(args.output_csv).parent.mkdir(parents=True, exist_ok=True)
    out_align.to_csv(args.output_csv, index=False)
    print(f"  wrote: {args.output_csv} ({len(out_align):,} rows)", flush=True)

    # ---- Step 5: per-donor x L2 cell counts ----
    print("[step 5] per-donor x L2 cell counts", flush=True)
    labeled = aligned.dropna(subset=["L2_label"])
    summary = (labeled.groupby(["L2_label", "patientID"], dropna=True)
                      .size()
                      .reset_index(name="n_cells"))
    summary.to_csv(args.donor_summary_csv, index=False)
    print(f"  wrote: {args.donor_summary_csv} ({len(summary):,} rows)", flush=True)

    # ---- Step 6: sample_type_coarse heterogeneity report ----
    # Recode the literal-string "NaN" -> real NaN first
    print("[step 6] sample_type_coarse heterogeneity check", flush=True)
    aligned["sample_type_coarse"] = aligned["sample_type_coarse"].replace("NaN", pd.NA)

    stc_per_donor = (aligned.dropna(subset=["sample_type_coarse"])
                            .groupby("patientID")["sample_type_coarse"]
                            .agg(lambda s: dict(s.value_counts())))
    stc_per_donor = stc_per_donor[stc_per_donor.apply(len) > 1]
    print(f"  donors with >1 sample_type_coarse value: {len(stc_per_donor)}", flush=True)

    # Flatten to a tabular report
    hetero_rows = []
    for donor, val_counts in stc_per_donor.items():
        total = sum(val_counts.values())
        for v, n in sorted(val_counts.items(), key=lambda x: -x[1]):
            hetero_rows.append({
                "patientID": donor,
                "sample_type_coarse": v,
                "n_cells": n,
                "pct_of_donor": round(100 * n / total, 2),
            })
    hetero_df = pd.DataFrame(hetero_rows)
    hetero_df.to_csv(args.stc_heterogeneity_csv, index=False)
    print(f"  wrote: {args.stc_heterogeneity_csv} ({len(hetero_df)} rows)", flush=True)
    if len(hetero_df) > 0:
        # Print first 10 lines
        print("  preview:", flush=True)
        print(hetero_df.head(20).to_string(index=False), flush=True)

    # ---- Step 7: augment donor metadata ----
    # Drop sample_type_coarse from final donor metadata since B4_austin is
    # being removed from the formula ladder. Just pass through the existing
    # donor metadata unchanged (no augmentation needed). Keep this step so
    # downstream scripts find an enriched_donor_metadata.csv at the expected path.
    print("[step 7] write enriched_donor_metadata.csv (passthrough)", flush=True)
    existing = pd.read_csv(args.existing_donor_meta)
    existing.to_csv(args.donor_meta_out, index=False)
    print(f"  wrote: {args.donor_meta_out} ({len(existing)} donors, "
          f"{existing.columns.size} cols)", flush=True)

    print("\nDone.", flush=True)


if __name__ == "__main__":
    sys.exit(main())
