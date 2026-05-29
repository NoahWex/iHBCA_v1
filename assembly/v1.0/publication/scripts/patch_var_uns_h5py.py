#!/usr/bin/env python3
"""
patch_var_uns_h5py.py - Fix CxG 5.3.2 validation blockers via h5py
====================================================================
1. Delete uns["layer_descriptions"] (deprecated in CxG 5.3.2)
2. Rename var["feature_biotype"] -> var["feature_biotype_gencode"]
   (CxG reserves "feature_biotype" as an auto-populated column)

Memory: ~100 MB (no data loading, just HDF5 group operations)

Usage:
  python patch_var_uns_h5py.py --h5ad path/to/file.h5ad [--dry-run]

"""

import argparse

import h5py
import numpy as np


def patch_file(h5ad_path, dry_run=False):
    """Patch one h5ad file in-place."""
    print(f"\n{'[DRY RUN] ' if dry_run else ''}Patching: {h5ad_path}")

    mode = "r" if dry_run else "a"
    with h5py.File(h5ad_path, mode) as f:
        changes = []

        # --- Fix 1: delete uns["layer_descriptions"] ---
        if "uns" in f and "layer_descriptions" in f["uns"]:
            if not dry_run:
                del f["uns"]["layer_descriptions"]
            changes.append("deleted uns/layer_descriptions")
        else:
            print("  uns/layer_descriptions: not present, skip")

        # --- Fix 2: rename var["feature_biotype"] -> var["feature_biotype_gencode"] ---
        if "var" in f:
            var = f["var"]

            # Check if it's a dataset or group (categorical = group)
            if "feature_biotype" in var:
                if "feature_biotype_gencode" in var:
                    print("  var/feature_biotype_gencode: already exists, skip rename")
                else:
                    if not dry_run:
                        var.move("feature_biotype", "feature_biotype_gencode")
                    changes.append("renamed var/feature_biotype -> feature_biotype_gencode")

                # Update column-order attr
                if "column-order" in var.attrs:
                    col_order = list(var.attrs["column-order"].astype(str))
                    if "feature_biotype" in col_order:
                        idx = col_order.index("feature_biotype")
                        col_order[idx] = "feature_biotype_gencode"
                        if not dry_run:
                            var.attrs["column-order"] = np.array(
                                col_order, dtype=h5py.string_dtype()
                            )
                        changes.append("updated var column-order attr")
            else:
                print("  var/feature_biotype: not present, skip")

        if changes:
            for c in changes:
                print(f"  {c}")
        else:
            print("  no changes needed")

    return len(changes) > 0


def main():
    parser = argparse.ArgumentParser(
        description="Fix CxG 5.3.2 validation blockers in h5ad via h5py"
    )
    parser.add_argument("--h5ad", nargs="+", required=True, help="Path(s) to h5ad file(s)")
    parser.add_argument("--dry-run", action="store_true", help="Report changes without applying")
    args = parser.parse_args()

    n_changed = 0
    for path in args.h5ad:
        if patch_file(path, args.dry_run):
            n_changed += 1

    print(f"\n{'[DRY RUN] ' if args.dry_run else ''}Done: {n_changed}/{len(args.h5ad)} files modified")


if __name__ == "__main__":
    main()
