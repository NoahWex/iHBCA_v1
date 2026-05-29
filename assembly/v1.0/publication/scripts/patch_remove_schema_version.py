#!/usr/bin/env python3
"""Remove CxG-reserved 'schema_version' key from uns in h5ad files.

The CxG validator reserves this key and sets it automatically during upload.
Having it pre-set causes validation failure.

Usage:
    python patch_remove_schema_version.py file1.h5ad [file2.h5ad ...]
"""

import argparse
import h5py


def patch_file(path):
    print(f"\nPatching: {path}")
    with h5py.File(path, "r+") as f:
        uns = f["uns"]
        if "schema_version" in uns:
            val = uns["schema_version"][()].decode() if hasattr(uns["schema_version"][()], "decode") else uns["schema_version"][()]
            del uns["schema_version"]
            print(f"  Removed uns['schema_version'] (was: {val})")
        else:
            print("  uns['schema_version'] not present — nothing to do")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("files", nargs="+", help="h5ad files to patch")
    args = p.parse_args()

    for path in args.files:
        patch_file(path)

    print("\nDone.")


if __name__ == "__main__":
    main()
