#!/usr/bin/env python3
"""Remove feature_is_filtered from raw/var in already-restructured h5ads.

The HCA schema prohibits feature_is_filtered in raw.var. The initial
restructure_raw_x.py run (job 50031446) copied it along with the rest
of /var. This patch removes it in-place.

Usage:
    python3 patch_remove_feature_is_filtered_raw.py file1.h5ad [file2.h5ad ...]
"""

import sys
import h5py


def patch(path: str) -> None:
    print(f"\n--- {path} ---")

    with h5py.File(path, "r") as f:
        has_raw_var = "raw" in f and "var" in f["raw"]
        if not has_raw_var:
            print("  SKIP: no raw/var group found")
            return
        has_dataset = "feature_is_filtered" in f["raw/var"]
        col_order = list(f["raw/var"].attrs.get("column-order", []))
        has_in_colorder = "feature_is_filtered" in col_order

    if not has_dataset and not has_in_colorder:
        print("  OK: feature_is_filtered not in raw/var dataset or column-order (already clean)")
        return

    print(f"  dataset present: {has_dataset}, in column-order: {has_in_colorder}")

    with h5py.File(path, "a") as f:
        import numpy as np

        # Delete the dataset if present
        if "feature_is_filtered" in f["raw/var"]:
            del f["raw/var/feature_is_filtered"]
            print(f"  Deleted raw/var/feature_is_filtered dataset")

        # Fix column-order attribute: delete and recreate to reliably update
        # variable-length string attributes (in-place overwrite is unreliable in h5py)
        col_order = list(f["raw/var"].attrs.get("column-order", []))
        if "feature_is_filtered" in col_order:
            col_order.remove("feature_is_filtered")
            del f["raw/var"].attrs["column-order"]
            f["raw/var"].attrs.create(
                "column-order",
                data=np.array(col_order, dtype=h5py.string_dtype()),
            )
            print(f"  Rebuilt column-order attr (removed feature_is_filtered)")
        else:
            print(f"  column-order attr already clean")

        remaining = sorted(f["raw/var"].keys())
        new_col_order = list(f["raw/var"].attrs.get("column-order", []))
        print(f"  raw/var keys now: {remaining}")
        print(f"  raw/var column-order now: {new_col_order}")
        assert "feature_is_filtered" not in new_col_order, \
            "FAIL: feature_is_filtered still in column-order after patch"

    print("  Done.")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} file1.h5ad [file2.h5ad ...]")
        sys.exit(1)

    for path in sys.argv[1:]:
        patch(path)

    print("\nAll files patched.")
