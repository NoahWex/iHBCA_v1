#!/usr/bin/env python3
"""Patch obs columns in restructured h5ads in-place using h5py.

Fixes:
  1. facs_status: for murrow/pal/reed cells, overwrite with FACS_status
     (L1 CSV is donor-level; FACS_status from source is sample-level and correct)
  2. n_genes, n_counts, percent_mito: convert from categorical-of-strings to float32

Usage:
    python3 patch_obs_fields.py file1.h5ad [file2.h5ad ...]
"""

import sys
import h5py
import numpy as np


STUDIES_TO_FIX_FACS = {"murrow", "pal", "reed"}
NUMERIC_COLS = ["n_genes", "n_counts", "percent_mito"]


def decode_categorical(grp):
    """Return numpy array of decoded string values from an h5py categorical group."""
    codes = grp["codes"][:]
    categories = grp["categories"][:]
    # categories may be bytes or str
    if categories.dtype.kind in ("S", "O"):
        categories = np.array([c.decode() if isinstance(c, bytes) else c for c in categories])
    result = np.empty(len(codes), dtype=object)
    valid = codes >= 0
    result[valid] = categories[codes[valid]]
    result[~valid] = None  # -1 == NaN category
    return result


def decode_categorical_as_float(grp):
    """Decode categorical-of-strings column to float32 array."""
    str_vals = decode_categorical(grp)
    result = np.empty(len(str_vals), dtype=np.float32)
    for i, v in enumerate(str_vals):
        if v is None or v == "nan" or v == "":
            result[i] = np.nan
        else:
            try:
                result[i] = float(v)
            except (ValueError, TypeError):
                result[i] = np.nan
    return result


def get_category_index(categories, value):
    """Return index of value in categories array, or -1 if not found."""
    for i, c in enumerate(categories):
        cat = c.decode() if isinstance(c, bytes) else c
        if cat == value:
            return i
    return -1


def patch(path: str) -> None:
    print(f"\n--- {path} ---")

    with h5py.File(path, "r") as f:
        # Quick checks
        if "obs" not in f:
            print("  SKIP: no /obs group")
            return
        obs = f["obs"]
        required = ["dataset", "facs_status", "FACS_status"] + NUMERIC_COLS
        missing = [c for c in required if c not in obs]
        if missing:
            print(f"  SKIP: missing obs columns: {missing}")
            return

    with h5py.File(path, "a") as f:
        obs = f["obs"]

        # ----------------------------------------------------------------
        # 1. facs_status fix: use FACS_status for murrow/pal/reed cells
        # ----------------------------------------------------------------
        print("  [1] Fixing facs_status for murrow/pal/reed...")

        dataset_vals = decode_categorical(obs["dataset"])
        facs_status_cats = obs["facs_status/categories"][:]
        facs_status_codes = obs["facs_status/codes"][:]
        FACS_status_cats = obs["FACS_status/categories"][:]
        FACS_status_codes = obs["FACS_status/codes"][:]

        # Build mask of cells to fix
        fix_mask = np.array([
            (v.decode() if isinstance(v, bytes) else v) in STUDIES_TO_FIX_FACS
            if v is not None else False
            for v in dataset_vals
        ])
        n_fix = fix_mask.sum()
        print(f"    Cells to fix: {n_fix:,}")

        # For each fixed cell, look up FACS_status value and find its index in facs_status categories
        new_codes = facs_status_codes.copy()
        changed = 0
        for i in np.where(fix_mask)[0]:
            src_code = FACS_status_codes[i]
            if src_code < 0:
                continue  # NaN → leave as-is
            src_cat = FACS_status_cats[src_code]
            src_val = src_cat.decode() if isinstance(src_cat, bytes) else src_cat
            tgt_code = get_category_index(facs_status_cats, src_val)
            if tgt_code == -1:
                print(f"    WARNING: '{src_val}' not in facs_status categories — skipping cell {i}")
                continue
            if new_codes[i] != tgt_code:
                new_codes[i] = tgt_code
                changed += 1

        obs["facs_status/codes"][:] = new_codes
        print(f"    Updated {changed:,} codes")

        # Verify
        new_facs = decode_categorical(obs["facs_status"])
        new_FACS = decode_categorical(obs["FACS_status"])
        still_wrong = 0
        for i in np.where(fix_mask)[0]:
            if new_facs[i] != new_FACS[i]:
                still_wrong += 1
        if still_wrong:
            print(f"    WARNING: {still_wrong:,} cells still disagree after fix")
        else:
            print(f"    Verification: all murrow/pal/reed cells now agree")

        # ----------------------------------------------------------------
        # 2. Convert n_genes, n_counts, percent_mito to float32
        # ----------------------------------------------------------------
        print("  [2] Converting numeric cols to float32...")

        for col in NUMERIC_COLS:
            if col not in obs:
                print(f"    {col}: not found, skipping")
                continue

            grp = obs[col]
            if not hasattr(grp, "keys") or "codes" not in grp:
                print(f"    {col}: not categorical, skipping (dtype={grp.dtype})")
                continue

            vals = decode_categorical_as_float(grp)
            n_nan = np.isnan(vals).sum()
            n_valid = (~np.isnan(vals)).sum()
            print(f"    {col}: {n_valid:,} valid, {n_nan:,} NaN → float32")

            del obs[col]
            ds = obs.create_dataset(col, data=vals, dtype=np.float32)
            ds.attrs["encoding-type"] = "array"
            ds.attrs["encoding-version"] = "0.2.0"

        # column-order stays unchanged (same column names, just different encoding)
        print(f"  Done.")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} file1.h5ad [file2.h5ad ...]")
        sys.exit(1)

    for path in sys.argv[1:]:
        patch(path)

    print("\nAll files patched.")
