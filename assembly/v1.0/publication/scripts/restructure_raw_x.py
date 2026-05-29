#!/usr/bin/env python3
"""Restructure integrated h5ad to match HCA recommended layout.

Before:
  X = raw counts (sparse CSR float32)
  layers/log_normalized = normalized expression (sparse CSR float32)
  raw = (does not exist)

After:
  X = normalized expression (from layers/log_normalized)
  raw/X = raw counts (from X)
  raw/var = copy of var (minus feature_is_filtered, prohibited by HCA schema)
  layers = (removed)

Uses h5py move operations — matrix data is not copied, only HDF5
internal references are renamed. The actual data movement is O(1).
"""

import argparse
import time

import h5py


def restructure(path: str, dry_run: bool = False) -> None:
    """Restructure h5ad in-place."""

    # --- Pre-flight checks ---
    with h5py.File(path, "r") as f:
        assert "X" in f, "Missing /X group"
        assert "layers" in f, "Missing /layers group"
        assert "log_normalized" in f["layers"], "Missing /layers/log_normalized"
        assert "raw" not in f, "/raw already exists — already restructured?"
        assert "var" in f, "Missing /var group"

        # Record structure for verification
        x_keys = sorted(f["X"].keys())
        ln_keys = sorted(f["layers/log_normalized"].keys())
        x_shape_attr = dict(f["X"].attrs)
        ln_shape_attr = dict(f["layers/log_normalized"].attrs)

        print(f"Pre-flight OK:")
        print(f"  /X keys: {x_keys}")
        print(f"  /X attrs: {x_shape_attr}")
        print(f"  /layers/log_normalized keys: {ln_keys}")
        print(f"  /layers/log_normalized attrs: {ln_shape_attr}")
        print(f"  /var keys: {sorted(f['var'].keys())}")

    if dry_run:
        print("\n[DRY RUN] Would restructure — exiting.")
        return

    # --- Restructure ---
    t0 = time.time()
    with h5py.File(path, "a") as f:
        # 1. Create /raw group
        print("\n1. Creating /raw group...")
        f.create_group("raw")

        # 2. Move /X → /raw/X (O(1) rename)
        print("2. Moving /X → /raw/X...")
        f.move("X", "raw/X")

        # 3. Copy /var → /raw/var (raw needs its own var)
        print("3. Copying /var → /raw/var...")
        f.copy("var", "raw/var")

        # 3b. Drop feature_is_filtered from raw/var — HCA schema prohibits it there
        if "feature_is_filtered" in f["raw/var"]:
            print("3b. Removing raw/var/feature_is_filtered (prohibited by HCA schema)...")
            del f["raw/var/feature_is_filtered"]
            # Also remove from column-order attribute or anndata will crash on load
            col_order = list(f["raw/var"].attrs.get("column-order", []))
            if "feature_is_filtered" in col_order:
                col_order.remove("feature_is_filtered")
                import numpy as np
                f["raw/var"].attrs["column-order"] = np.array(col_order, dtype=object)

        # 4. Move /layers/log_normalized → /X (O(1) rename)
        print("4. Moving /layers/log_normalized → /X...")
        f.move("layers/log_normalized", "X")

        # 5. Clean up empty layers group
        if len(f["layers"]) == 0:
            print("5. Removing empty /layers group...")
            del f["layers"]
        else:
            remaining = sorted(f["layers"].keys())
            print(f"5. /layers still has keys: {remaining} — keeping.")

        # 6. Set encoding-type on /raw for anndata compatibility
        print("6. Setting /raw encoding metadata...")
        f["raw"].attrs["encoding-type"] = "raw"
        f["raw"].attrs["encoding-version"] = "0.1.0"

    elapsed = time.time() - t0
    print(f"\nRestructure complete in {elapsed:.1f}s")

    # --- Verify ---
    print("\n--- Verification ---")
    with h5py.File(path, "r") as f:
        # Check new structure exists
        assert "X" in f, "FAIL: /X missing after restructure"
        assert "raw" in f, "FAIL: /raw missing after restructure"
        assert "X" in f["raw"], "FAIL: /raw/X missing"
        assert "var" in f["raw"], "FAIL: /raw/var missing"
        assert "layers" not in f or "log_normalized" not in f.get("layers", {}), \
            "FAIL: /layers/log_normalized still exists"

        # Verify the swap: new /X should have log_normalized's structure,
        # /raw/X should have original X's structure
        new_x_keys = sorted(f["X"].keys())
        raw_x_keys = sorted(f["raw/X"].keys())
        new_x_attrs = dict(f["X"].attrs)
        raw_x_attrs = dict(f["raw/X"].attrs)

        print(f"  /X keys: {new_x_keys} (was log_normalized: {ln_keys})")
        print(f"  /X attrs: {new_x_attrs}")
        print(f"  /raw/X keys: {raw_x_keys} (was X: {x_keys})")
        print(f"  /raw/X attrs: {raw_x_attrs}")
        raw_var_keys = sorted(f["raw/var"].keys())
        print(f"  /raw/var keys: {raw_var_keys}")
        assert "feature_is_filtered" not in f["raw/var"], \
            "FAIL: feature_is_filtered still present in raw/var"

        # Verify data shapes match expectations
        if "indptr" in f["X"]:
            n_rows_x = f["X/indptr"].shape[0] - 1
            print(f"  /X n_rows (from indptr): {n_rows_x}")
        if "indptr" in f["raw/X"]:
            n_rows_raw = f["raw/X/indptr"].shape[0] - 1
            print(f"  /raw/X n_rows (from indptr): {n_rows_raw}")

        # Check no /layers or empty
        if "layers" in f:
            print(f"  WARNING: /layers still present with keys: {sorted(f['layers'].keys())}")
        else:
            print(f"  /layers: removed (clean)")

    print("\nAll checks passed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("h5ad", help="Path to h5ad file")
    parser.add_argument("--dry-run", action="store_true",
                        help="Check structure without modifying")
    args = parser.parse_args()
    restructure(args.h5ad, dry_run=args.dry_run)
