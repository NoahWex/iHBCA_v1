#!/usr/bin/env python3
# Step 20a — Assemble full-object integration intermediate
#
# Concatenates per-compartment counts matrices and obs metadata into a single
# integration_intermediate/full/scvi_n100/ directory, matching the structure
# of the per-compartment directories. Pulls full-object latent and UMAP from
# the step 17 sweep outputs where they already exist.
#
# Fails loudly on gene universe mismatches or missing inputs.
#
# Usage:
#   python 20a_concat_full_object.py \
#       --integration-root /path/to/integration_intermediate \
#       --step17-root      /path/to/17_ScviIntegration \
#       --out-dir          /path/to/integration_intermediate/full/scvi_n100

import argparse
import json
import os
import shutil
import sys
from datetime import date

import pandas as pd
import scipy.io
import scipy.sparse


COMPARTMENTS = ["Epithelial", "Immune", "Stromal"]


def fail(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--integration-root", required=True,
                        help="Path to integration_intermediate/")
    parser.add_argument("--step17-root", required=True,
                        help="Path to 17_ScviIntegration/ sidecar directory")
    parser.add_argument("--out-dir", required=True,
                        help="Output directory (integration_intermediate/full/scvi_n100)")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # -----------------------------------------------------------------------
    # Load and validate per-compartment data
    # -----------------------------------------------------------------------
    mats: list[scipy.sparse.csc_matrix] = []
    obs_frames: list[pd.DataFrame] = []

    # Load reference gene list from first compartment
    first_dir = os.path.join(args.integration_root, COMPARTMENTS[0], "scvi_n100")
    genes_ref: list[str] = open(os.path.join(first_dir, "genes.tsv")).read().splitlines()

    for comp in COMPARTMENTS:
        comp_dir = os.path.join(args.integration_root, comp, "scvi_n100")
        print(f"Loading {comp}...")

        genes = open(os.path.join(comp_dir, "genes.tsv")).read().splitlines()
        if genes != genes_ref:
            fail(f"{comp} gene universe differs from {COMPARTMENTS[0]}")

        cells = open(os.path.join(comp_dir, "cells.tsv")).read().splitlines()
        _raw = scipy.io.mmread(os.path.join(comp_dir, "counts.mtx.gz"))
        mat: scipy.sparse.csc_matrix = _raw.tocsc()  # type: ignore[union-attr,assignment]

        # Matrix stored as cells x genes on disk
        if mat.shape != (len(cells), len(genes_ref)):
            fail(f"{comp}: matrix shape {mat.shape} doesn't match "
                 f"({len(cells)} cells x {len(genes_ref)} genes)")

        obs = pd.read_csv(os.path.join(comp_dir, "obs.csv"))
        obs.insert(1, "compartment", comp)

        mats.append(mat)
        obs_frames.append(obs)
        n_c, n_g = mat.shape[0], mat.shape[1]  # type: ignore[index]
        print(f"  {comp}: {n_c:,} cells x {n_g:,} genes")

    # -----------------------------------------------------------------------
    # Concatenate
    # -----------------------------------------------------------------------
    print("\nConcatenating...")
    # Matrix is stored cells x genes on disk; concatenate along cell axis (rows)
    full_mat: scipy.sparse.csc_matrix = scipy.sparse.vstack(mats, format="csc")  # type: ignore[assignment]
    full_obs = pd.concat(obs_frames, ignore_index=True)
    full_cells = full_obs["cell_id"].tolist()

    n_cells_full, n_genes_full = full_mat.shape[0], full_mat.shape[1]  # type: ignore[index]
    print(f"Full object: {n_cells_full:,} cells x {n_genes_full:,} genes")

    if len(full_cells) != n_cells_full:
        fail(f"Cell count mismatch: obs has {len(full_cells)}, matrix has {n_cells_full}")

    # -----------------------------------------------------------------------
    # Write counts
    # -----------------------------------------------------------------------
    print("Writing counts.mtx.gz...")
    tmp_mtx = os.path.join(args.out_dir, "_counts_tmp.mtx")
    scipy.io.mmwrite(tmp_mtx, full_mat)
    os.system(f"gzip -f {tmp_mtx} && mv {tmp_mtx}.gz {os.path.join(args.out_dir, 'counts.mtx.gz')}")

    # -----------------------------------------------------------------------
    # Write genes, cells, obs
    # -----------------------------------------------------------------------
    with open(os.path.join(args.out_dir, "genes.tsv"), "w") as f:
        f.write("\n".join(genes_ref) + "\n")

    with open(os.path.join(args.out_dir, "cells.tsv"), "w") as f:
        f.write("\n".join(full_cells) + "\n")

    full_obs.to_csv(os.path.join(args.out_dir, "obs.csv"), index=False)

    # -----------------------------------------------------------------------
    # Pull latent + UMAP from step 17 sweep outputs
    # -----------------------------------------------------------------------
    latent_src = os.path.join(args.step17_root, "embeddings", "latent_50d.csv")
    umap_src   = os.path.join(args.step17_root, "winner", "full_object", "umap.csv")

    for src, dst_name in [(latent_src, "latent.csv"), (umap_src, "umap.csv")]:
        if not os.path.exists(src):
            fail(f"Missing step 17 file: {src}")
        dst = os.path.join(args.out_dir, dst_name)
        shutil.copy(src, dst)
        print(f"Copied {dst_name} from step 17")

    # -----------------------------------------------------------------------
    # Write manifest.json
    # -----------------------------------------------------------------------
    obs_columns = list(full_obs.columns)
    leiden_cols = [c for c in obs_columns if c.startswith("leiden_")]

    manifest = {
        "source": "concatenated from per-compartment integration_intermediate",
        "compartments": COMPARTMENTS,
        "date": str(date.today()),
        "config_label": "scvi_n100",
        "n_cells": n_cells_full,
        "n_genes": n_genes_full,
        "latent_key": "X_emb",
        "umap_key": "X_umap",
        "obs_columns": obs_columns,
        "leiden_resolutions": sorted(
            float(c.replace("leiden_", "")) for c in leiden_cols
        ),
    }

    with open(os.path.join(args.out_dir, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)

    print(f"\nDone. Full object written to: {args.out_dir}")


if __name__ == "__main__":
    main()
