#!/usr/bin/env python3
# Step 20 — Build canonical manifest
#
# Walks integration_intermediate/{compartment}/scvi_n100/ (excluding patches/)
# and deseq/{compartment}/, records all files with relative paths anchored to
# outputs_root, and writes outputs/manifest.yaml. Fails loudly on missing
# expected files — no silent fallbacks.
#
# Usage:
#   python 20_build_manifest.py \
#       --outputs-root /path/to/outputs \
#       --out          /path/to/manifest.yaml \
#       [--date YYYY-MM-DD]

import argparse
import json
import os
import sys
from datetime import date

import yaml


COMPARTMENTS = ["Epithelial", "Immune", "Stromal"]

REQUIRED_FILES = [
    "latent.csv", "umap.csv", "obs.csv",
    "counts.mtx.gz", "genes.tsv", "cells.tsv", "manifest.json",
]


def fail(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def relpath(abs_path, outputs_root):
    return os.path.relpath(abs_path, outputs_root)


def validate_compartment(comp, outputs_root):
    abs_base = os.path.join(outputs_root, "integration_intermediate", comp, "scvi_n100")

    if not os.path.isdir(abs_base):
        fail(f"{comp}: integration directory not found: {abs_base}")

    # Confirm required files exist
    missing = [f for f in REQUIRED_FILES if not os.path.exists(os.path.join(abs_base, f))]
    if missing:
        fail(f"{comp}: missing required files:\n  " + "\n  ".join(missing))

    # Walk all files, exclude patches/
    integration_files = {}
    for fname in sorted(os.listdir(abs_base)):
        fpath = os.path.join(abs_base, fname)
        if os.path.isfile(fpath):
            integration_files[fname] = relpath(fpath, outputs_root)

    # Cell and gene counts
    with open(os.path.join(abs_base, "obs.csv")) as f:
        header = f.readline().strip().split(",")
        n_cells = sum(1 for _ in f)

    with open(os.path.join(abs_base, "genes.tsv")) as f:
        n_genes = sum(1 for _ in f)

    leiden_resolutions = sorted(
        float(c.replace("leiden_", "")) for c in header if c.startswith("leiden_")
    )

    # Provenance from manifest.json
    with open(os.path.join(abs_base, "manifest.json")) as f:
        provenance = json.load(f)

    print(f"  {comp}: {n_cells:,} cells x {n_genes:,} genes, "
          f"leiden resolutions: {leiden_resolutions}")

    entry = {
        "n_cells": n_cells,
        "n_genes": n_genes,
        "leiden_resolutions": leiden_resolutions,
        "provenance": provenance,
        "integration": integration_files,
        "deseq": None,
    }

    # DESeq outputs
    deseq_abs = os.path.join(outputs_root, "deseq", comp)
    if os.path.isdir(deseq_abs):
        csvs = sorted(f for f in os.listdir(deseq_abs) if f.endswith(".csv"))
        if csvs:
            entry["deseq"] = {
                f: relpath(os.path.join(deseq_abs, f), outputs_root) for f in csvs
            }
            print(f"    DESeq: {len(csvs)} resolution CSV(s)")
        else:
            print(f"    DESeq: directory exists but no CSVs yet")
    else:
        print(f"    DESeq: not yet run")

    return entry


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outputs-root", required=True,
                        help="Absolute path to outputs/ directory")
    parser.add_argument("--out", required=True,
                        help="Output path for manifest.yaml")
    parser.add_argument("--date", default=str(date.today()))
    args = parser.parse_args()

    if not os.path.isdir(args.outputs_root):
        fail(f"outputs-root does not exist: {args.outputs_root}")

    print(f"Outputs root: {args.outputs_root}")
    print(f"Manifest out: {args.out}\n")

    manifest = {
        "version": "1.0",
        "date": args.date,
        "outputs_root": args.outputs_root,
        "compartments": {},
    }

    total_cells = 0
    for comp in COMPARTMENTS:
        entry = validate_compartment(comp, args.outputs_root)
        manifest["compartments"][comp] = entry
        total_cells += entry["n_cells"]

    manifest["total_cells"] = total_cells
    print(f"\nTotal cells: {total_cells:,}")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as f:
        yaml.dump(manifest, f, default_flow_style=False, sort_keys=False)

    print(f"Manifest written: {args.out}")


if __name__ == "__main__":
    main()
