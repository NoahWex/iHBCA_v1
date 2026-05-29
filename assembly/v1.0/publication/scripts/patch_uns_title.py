#!/usr/bin/env python3
"""
patch_uns_title.py — Set uns["title"] in one or more h5ad files (in-place, h5py).

Usage:
  python3 patch_uns_title.py --h5ad <file1> [<file2> ...] --title "New Title"
"""
import argparse
import sys
import h5py
import numpy as np


def inspect_uns(path: str):
    print(f"\n--- uns in {path} ---")
    with h5py.File(path, "r") as f:
        if "uns" not in f:
            print("  No uns group found.")
            return
        uns = f["uns"]
        for k in uns.keys():
            try:
                v = uns[k][()]
                if isinstance(v, (bytes, np.bytes_)):
                    v = v.decode()
                elif isinstance(v, np.ndarray) and v.dtype.kind == "S":
                    v = v.astype(str)
                print(f"  {k}: {repr(v)[:160]}")
            except Exception:
                print(f"  {k}: (group/unreadable)")


def patch_title(path: str, title: str, dry_run: bool = False):
    with h5py.File(path, "r" if dry_run else "r+") as f:
        if "uns" not in f:
            print(f"  ERROR: no 'uns' group in {path}")
            return False

        uns = f["uns"]
        existing = None
        if "title" in uns:
            raw = uns["title"][()]
            if isinstance(raw, (bytes, np.bytes_)):
                raw = raw.decode()
            existing = raw

        print(f"  current uns['title']: {repr(existing)}")
        print(f"  -> setting to:        {repr(title)}")

        if dry_run:
            print("  (dry run — no changes written)")
            return True

        if "title" in uns:
            del uns["title"]
        uns.create_dataset(
            "title",
            data=np.bytes_(title.encode("utf-8")),
            dtype=h5py.string_dtype(),
        )

        # Verify
        written = uns["title"][()]
        if isinstance(written, (bytes, np.bytes_)):
            written = written.decode()
        assert written == title, f"Mismatch after write: {repr(written)}"
        print(f"  VERIFIED: {repr(written)}")
        return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", nargs="+", required=True)
    parser.add_argument("--title", required=True)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--inspect-only", action="store_true")
    args = parser.parse_args()

    ok = True
    for path in args.h5ad:
        print(f"\n=== {path} ===")
        if args.inspect_only:
            inspect_uns(path)
            continue
        inspect_uns(path)
        result = patch_title(path, args.title, dry_run=args.dry_run)
        if not result:
            ok = False

    if not ok:
        sys.exit(1)
    print("\nDone.")


if __name__ == "__main__":
    main()
