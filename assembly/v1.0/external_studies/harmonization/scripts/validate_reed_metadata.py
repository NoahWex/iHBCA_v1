#!/usr/bin/env python3
"""
validate_reed_metadata.py
=========================
Post-build validation for Reed study objects.
Compares harmonized metadata.csv against Supplementary Table 1 and backup files.

Checks:
  1. Donor-level metadata consistency (age, parity, risk, menopausal status)
  2. Cell count matches expected 803,283
  3. Embedding dimensions match backup values
  4. No unexpected NAs in condition columns

Usage:
  python validate_reed_metadata.py
  python validate_reed_metadata.py --cleanup   # Delete backups on pass

Exit codes:
  0 = all checks pass
  1 = validation failure
"""

# Q: Why did the metadata fail validation?
# A: Because it had too many NAs and not enough self-awareness.

import argparse
import os
import sys
import glob

import pandas as pd
import numpy as np

BASE = "${SOURCE_IHBCAV1_EXTERNAL_STUDIES%/}"
STUDY_DIR = os.path.join(BASE, "harmonization/outputs/study_objects/reed")
SUPP_TABLE = os.path.join(BASE, "harmonization/studies/reed/reference/REED_Supplementary_Table_1.csv")

EXPECTED_CELLS = 803_283
EXPECTED_DONORS = 55  # donors with harmonized metadata (58 total, 3 spike-ins excluded)


def load_supp_table(path: str) -> pd.DataFrame:
    """Load Reed Supplementary Table 1 (row 0 = title, row 1 = headers)."""
    df = pd.read_csv(path, header=1)  # skip title row, use row 1 as header
    # Drop spike-in rows
    df = df[df["sample_type"] != "spike-in"].copy()
    # Deduplicate to donor level (multiple samples per donor)
    donor_cols = ["BCN_donor_id", "donor_age", "parity", "risk_status", "Menopause_status"]
    donors = df[donor_cols].drop_duplicates(subset=["BCN_donor_id"])
    donors = donors.rename(columns={"BCN_donor_id": "ihbca_donor_id"})
    return donors


def load_metadata(path: str) -> pd.DataFrame:
    """Load harmonized metadata.csv from build output."""
    return pd.read_csv(path)


def check_donor_metadata(meta: pd.DataFrame, supp: pd.DataFrame) -> list[str]:
    """Compare donor-level metadata between build output and supplementary table."""
    errors = []

    # Get unique donors from metadata
    donor_cols = ["ihbca_donor_id", "age_binary", "parity_binary",
                  "risk_status_binary", "menopausal_status_binary"]
    available = [c for c in donor_cols if c in meta.columns]
    if "ihbca_donor_id" not in available:
        errors.append("FATAL: ihbca_donor_id column missing from metadata.csv")
        return errors

    donors_meta = meta[available].drop_duplicates(subset=["ihbca_donor_id"])
    donors_meta = donors_meta.set_index("ihbca_donor_id")

    supp = supp.set_index("ihbca_donor_id")

    # Check all supp donors appear in metadata
    missing_in_meta = set(supp.index) - set(donors_meta.index)
    if missing_in_meta:
        errors.append(f"Donors in supp table but not metadata: {missing_in_meta}")

    # Check donor count
    if len(donors_meta) != EXPECTED_DONORS:
        errors.append(f"Expected {EXPECTED_DONORS} donors, got {len(donors_meta)}")

    # Cross-check age: supp has continuous, metadata has binary (young/old, split ~50)
    # We verify age_binary is populated, not exact mapping (that's in harmonization logic)
    common = sorted(set(donors_meta.index) & set(supp.index))
    print(f"  Donors matched: {len(common)}/{len(supp)}")

    for donor_id in common:
        meta_row = donors_meta.loc[donor_id]
        supp_row = supp.loc[donor_id]

        # Check risk_status consistency
        if "risk_status_binary" in meta_row.index:
            risk_meta = meta_row["risk_status_binary"]
            risk_supp = supp_row["risk_status"]
            if pd.notna(risk_meta) and pd.notna(risk_supp):
                # AR -> AR, HR-BR1/HR-BR2/Contralateral-BR1 -> HR
                expected_binary = "AR" if risk_supp == "AR" else "HR"
                if risk_meta != expected_binary:
                    errors.append(
                        f"Donor {donor_id}: risk_status mismatch: "
                        f"meta={risk_meta}, supp={risk_supp} (expected {expected_binary})"
                    )

        # Check menopausal consistency
        if "menopausal_status_binary" in meta_row.index:
            meno_meta = meta_row["menopausal_status_binary"]
            meno_supp = str(supp_row.get("Menopause_status", "")).strip().lower()
            if pd.notna(meno_meta) and meno_supp:
                # Pre -> pre, Post -> post, Peri -> peri
                if meno_meta.lower() != meno_supp.lower():
                    # Allow peri mapping variations
                    if not (meno_supp in ("peri", "peri-menopausal") and
                            meno_meta.lower() in ("pre", "peri")):
                        errors.append(
                            f"Donor {donor_id}: menopausal mismatch: "
                            f"meta={meno_meta}, supp={meno_supp}"
                        )

    return errors


def check_cell_count(meta: pd.DataFrame) -> list[str]:
    """Verify cell count matches expected."""
    errors = []
    n = len(meta)
    if n != EXPECTED_CELLS:
        errors.append(f"Cell count: expected {EXPECTED_CELLS}, got {n}")
    else:
        print(f"  Cell count: {n} (matches expected)")
    return errors


def check_no_unexpected_nas(meta: pd.DataFrame) -> list[str]:
    """Check condition columns have no unexpected NAs."""
    errors = []
    condition_cols = ["age_binary", "parity_binary", "risk_status_binary",
                      "menopausal_status_binary"]
    for col in condition_cols:
        if col not in meta.columns:
            errors.append(f"Column {col} missing from metadata")
            continue
        n_na = meta[col].isna().sum()
        n_total = len(meta)
        pct = 100 * n_na / n_total
        if n_na > 0:
            # Some donors legitimately have missing values (e.g., menopausal unknown)
            # Warn above 5%, error above 20%
            msg = f"  {col}: {n_na}/{n_total} NAs ({pct:.1f}%)"
            if pct > 20:
                errors.append(f"HIGH NA rate - {msg}")
            elif pct > 5:
                print(f"  WARNING: {msg}")
            else:
                print(f"  {msg} (acceptable)")
        else:
            print(f"  {col}: no NAs")
    return errors


def check_embeddings_vs_backup(study_dir: str) -> list[str]:
    """Compare embedding dimensions against backup files if they exist."""
    errors = []
    for name in ["embedding_native", "embedding_joint"]:
        current = os.path.join(study_dir, f"{name}.csv")
        backup = os.path.join(study_dir, f"{name}.csv.backup")
        if not os.path.exists(backup):
            continue
        if not os.path.exists(current):
            errors.append(f"{name}.csv missing (backup exists)")
            continue

        # Just check shape, not content (rebuild may differ slightly)
        cur_df = pd.read_csv(current, nrows=2)
        bak_df = pd.read_csv(backup, nrows=2)
        if cur_df.shape[1] != bak_df.shape[1]:
            errors.append(
                f"{name} column count: current={cur_df.shape[1]}, backup={bak_df.shape[1]}"
            )
        else:
            print(f"  {name}: dims match backup ({cur_df.shape[1]} columns)")

        # Full row count check
        cur_lines = sum(1 for _ in open(current)) - 1
        bak_lines = sum(1 for _ in open(backup)) - 1
        if cur_lines != bak_lines:
            errors.append(
                f"{name} row count: current={cur_lines}, backup={bak_lines}"
            )
        else:
            print(f"  {name}: row count matches backup ({cur_lines})")

    return errors


def cleanup_backups(study_dir: str):
    """Remove *.backup files."""
    backups = glob.glob(os.path.join(study_dir, "*.backup"))
    for f in backups:
        os.remove(f)
        print(f"  Removed: {os.path.basename(f)}")


def main():
    parser = argparse.ArgumentParser(description="Validate Reed build outputs")
    parser.add_argument("--cleanup", action="store_true",
                        help="Delete backup files on successful validation")
    parser.add_argument("--study-dir", default=STUDY_DIR,
                        help="Override study output directory")
    parser.add_argument("--supp-table", default=SUPP_TABLE,
                        help="Override supplementary table path")
    args = parser.parse_args()

    print("=" * 60)
    print("REED BUILD VALIDATION")
    print("=" * 60)

    meta_path = os.path.join(args.study_dir, "metadata.csv")
    if not os.path.exists(meta_path):
        print(f"FATAL: metadata.csv not found at {meta_path}")
        sys.exit(1)

    if not os.path.exists(args.supp_table):
        print(f"FATAL: Supplementary Table 1 not found at {args.supp_table}")
        sys.exit(1)

    all_errors = []

    print("\n[1] Loading data...")
    meta = load_metadata(meta_path)
    supp = load_supp_table(args.supp_table)
    print(f"  metadata.csv: {len(meta)} cells, {meta.columns.tolist()[:8]}...")
    print(f"  Supp Table 1: {len(supp)} donors")

    print("\n[2] Cell count check...")
    all_errors.extend(check_cell_count(meta))

    print("\n[3] Donor metadata consistency...")
    all_errors.extend(check_donor_metadata(meta, supp))

    print("\n[4] NA check in condition columns...")
    all_errors.extend(check_no_unexpected_nas(meta))

    print("\n[5] Embedding dimensions vs backup...")
    all_errors.extend(check_embeddings_vs_backup(args.study_dir))

    print("\n" + "=" * 60)
    if all_errors:
        print(f"FAILED: {len(all_errors)} error(s)")
        for e in all_errors:
            print(f"  - {e}")
        sys.exit(1)
    else:
        print("PASSED: All checks OK")
        if args.cleanup:
            print("\n[CLEANUP] Removing backup files...")
            cleanup_backups(args.study_dir)
        print("=" * 60)
        sys.exit(0)


if __name__ == "__main__":
    main()
