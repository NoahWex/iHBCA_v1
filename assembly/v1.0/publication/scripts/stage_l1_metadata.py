#!/usr/bin/env python3
"""
stage_l1_metadata.py - Create L1 harmonized donor staging CSV
================================================================
Reads harmonized_donor_metadata.csv, selects/renames columns for obs
injection, normalizes missing values, and writes L1_harmonized_donor.csv.

Run once to produce the staging file. Assembly scripts consume the output.

Usage:
  python stage_l1_metadata.py --repo-root /path/to/iHBCAv1_upload
"""

import argparse
import csv
import os
import sys

import pandas as pd


# Join keys (kept in CSV for assembly scripts, but not injected into obs)
L1_JOIN_KEYS = ["ihbca_donor_id", "study"]

# Obs-ready columns from harmonized CSV
L1_OBS_COLUMNS = [
    "age_continuous",
    "age_binary",
    "parity_count",
    "parity_binary",
    "age_at_first_birth",
    "brca_genotype",
    "cancer_history",
    "tissue_indication",
    "risk_status_binary",
    "risk_genotype_only",
    "menopausal_status_detailed",
    "menopausal_status_binary",
    "ethnicity_verbatim",
    "ethnicity_grouped",
    "bmi_continuous",
    "bmi_category",
    "sample_preservation",
    "sample_type",
    "facs_status",
    "dissociation_minutes",
    "metadata_notes",
]

# Numeric columns (keep NaN for missing)
NUMERIC_COLUMNS = {
    "age_continuous",
    "parity_count",
    "age_at_first_birth",
    "bmi_continuous",
    "dissociation_minutes",
}

# All columns to select from harmonized CSV
ALL_L1_COLUMNS = L1_JOIN_KEYS + L1_OBS_COLUMNS


def stage_l1(repo_root):
    """Read harmonized CSV, select and normalize columns, write L1 CSV."""
    harm_path = os.path.join(
        repo_root,
        "external_studies/harmonization/outputs/harmonized_metadata/"
        "harmonized_donor_metadata.csv",
    )
    out_dir = os.path.join(repo_root, "publication/config/metadata_stages")
    out_path = os.path.join(out_dir, "L1_harmonized_donor.csv")

    # --- Load ---
    if not os.path.exists(harm_path):
        print(f"ERROR: Harmonized metadata not found: {harm_path}")
        sys.exit(1)

    df = pd.read_csv(harm_path, dtype=str)
    print(f"Loaded {len(df)} donors from {harm_path}")
    print(f"  Columns: {list(df.columns)}")

    # --- Validate expected columns ---
    missing = [c for c in ALL_L1_COLUMNS if c not in df.columns]
    if missing:
        print(f"ERROR: Missing columns in harmonized CSV: {missing}")
        sys.exit(1)

    # --- Select ---
    l1 = df[ALL_L1_COLUMNS].copy()
    print(f"\nSelected {len(ALL_L1_COLUMNS)} columns ({len(L1_JOIN_KEYS)} join keys + {len(L1_OBS_COLUMNS)} obs)")

    # --- Convert numeric columns to float ---
    for col in NUMERIC_COLUMNS:
        l1[col] = pd.to_numeric(l1[col], errors="coerce")

    # --- Normalize string columns: empty/NaN -> "unknown" ---
    string_cols = [c for c in L1_OBS_COLUMNS if c not in NUMERIC_COLUMNS]
    for col in string_cols:
        l1[col] = l1[col].fillna("").str.strip()
        # metadata_notes: use empty string for missing (not "unknown")
        if col == "metadata_notes":
            continue
        l1[col] = l1[col].replace("", "unknown")

    # --- Validate ---
    n_donors = len(l1)
    n_studies = l1["study"].nunique()
    dup_within_study = l1.duplicated(subset=["study", "ihbca_donor_id"], keep=False)
    if dup_within_study.any():
        dups = l1[dup_within_study][["study", "ihbca_donor_id"]]
        print(f"WARNING: {dup_within_study.sum()} duplicate donor IDs within study:")
        print(dups.to_string(index=False))

    # --- Summary ---
    print(f"\n--- L1 Summary ---")
    print(f"  Donors: {n_donors}")
    print(f"  Studies: {n_studies}")
    print(f"  Columns: {len(l1.columns)} ({len(L1_JOIN_KEYS)} join keys + {len(L1_OBS_COLUMNS)} obs)")

    print(f"\n  Missing values per column:")
    for col in L1_OBS_COLUMNS:
        if col in NUMERIC_COLUMNS:
            n_missing = l1[col].isna().sum()
            label = "NaN"
        else:
            if col == "metadata_notes":
                n_missing = (l1[col] == "").sum()
                label = '""'
            else:
                n_missing = (l1[col] == "unknown").sum()
                label = '"unknown"'
        if n_missing > 0:
            print(f"    {col}: {n_missing}/{n_donors} {label}")

    print(f"\n  Per-study donor counts:")
    for study, group in l1.groupby("study"):
        print(f"    {study}: {len(group)} donors")

    # --- Write (quote all string fields to preserve leading whitespace) ---
    os.makedirs(out_dir, exist_ok=True)
    l1.to_csv(out_path, index=False, quoting=csv.QUOTE_NONNUMERIC)
    print(f"\nWritten: {out_path}")
    print(f"  Size: {os.path.getsize(out_path):,} bytes")


def main():
    parser = argparse.ArgumentParser(
        description="Stage L1 harmonized donor metadata for h5ad assembly"
    )
    parser.add_argument(
        "--repo-root",
        required=True,
        help="Path to iHBCAv1_upload repo root",
    )
    args = parser.parse_args()
    stage_l1(args.repo_root)


if __name__ == "__main__":
    main()
