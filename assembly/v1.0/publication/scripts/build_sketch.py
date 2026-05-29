#!/usr/bin/env python3
"""Build a geometric sketch of the integrated iHBCA object.

Produces a ~280K-cell downsampled h5ad that preserves biological diversity
via geometric sketching (Hie et al. 2019) on the scVI latent space,
applied per-donor (1000 cells/donor, keeping all cells for small donors).

Usage:
    python build_sketch.py \
        --input /path/to/all-breast-cells.h5ad \
        --output /path/to/all-breast-cells-sketch.h5ad \
        [--cells-per-donor 1000] \
        [--seed 42]
"""

import argparse
import os
import shutil
import time

import anndata as ad
import geosketch
import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="Path to enriched all-breast-cells.h5ad")
    p.add_argument("--output", required=True, help="Path for sketch output h5ad")
    p.add_argument("--cells-per-donor", type=int, default=1000,
                   help="Target cells per donor (default: 1000)")
    p.add_argument("--seed", type=int, default=42,
                   help="Random seed for reproducibility (default: 42)")
    return p.parse_args()


def build_sketch(adata, cells_per_donor, seed):
    """Select cells via per-donor geometric sketch on X_scvi_100.

    Returns sorted array of selected integer indices into adata.obs.
    """
    embedding_key = "X_scvi_100"
    if embedding_key not in adata.obsm:
        raise KeyError(f"Embedding '{embedding_key}' not found in obsm. "
                       f"Available: {list(adata.obsm.keys())}")

    embedding = adata.obsm[embedding_key]
    donor_ids = adata.obs["donor_id"]

    # Sort donors alphabetically for determinism
    unique_donors = sorted(donor_ids.unique())
    print(f"  Total donors: {len(unique_donors)}")
    print(f"  Target cells/donor: {cells_per_donor}")
    print(f"  Seed: {seed}")

    selected_indices = []
    donors_below_threshold = []
    kept_all_count = 0
    sketched_count = 0

    for donor in unique_donors:
        donor_mask = (donor_ids == donor).values
        donor_positions = np.where(donor_mask)[0]
        n_cells = len(donor_positions)

        if n_cells <= cells_per_donor:
            # Keep all cells for small donors
            selected_indices.extend(donor_positions.tolist())
            donors_below_threshold.append(donor)
            kept_all_count += 1
        else:
            # Geometric sketch on this donor's embedding
            donor_embedding = embedding[donor_positions]
            sketch_idx = geosketch.gs(donor_embedding, cells_per_donor,
                                      replace=False, seed=seed)
            # Map back to global indices
            global_idx = donor_positions[sketch_idx]
            selected_indices.extend(global_idx.tolist())
            sketched_count += 1

    print(f"  Donors kept in full: {kept_all_count}")
    print(f"  Donors sketched: {sketched_count}")
    print(f"  Total selected cells: {len(selected_indices)}")

    # Sort indices for deterministic row ordering
    selected_indices = np.sort(np.array(selected_indices))

    return selected_indices, donors_below_threshold


def main():
    args = parse_args()
    t0 = time.time()

    # ── Load parent ──
    print(f"Loading parent h5ad: {args.input}")
    print(f"  File size: {os.path.getsize(args.input) / 1e9:.1f} GB")
    adata = ad.read_h5ad(args.input)
    t_load = time.time() - t0
    print(f"  Loaded in {t_load:.0f}s: {adata.shape[0]} cells x {adata.shape[1]} genes")
    print(f"  obs columns: {len(adata.obs.columns)}")
    print(f"  obsm keys: {list(adata.obsm.keys())}")
    print(f"  layers: {list(adata.layers.keys())}")

    # ── Geometric sketch ──
    print("\nRunning geometric sketch...")
    selected_indices, donors_below = build_sketch(
        adata, args.cells_per_donor, args.seed
    )

    # ── Subset ──
    print("\nSubsetting...")
    sketch = adata[selected_indices].copy()
    print(f"  Sketch shape: {sketch.shape[0]} cells x {sketch.shape[1]} genes")

    # ── Set is_primary_data = False ──
    sketch.obs["is_primary_data"] = False

    # ── Add sketch metadata to uns ──
    sketch.uns["sketch_metadata"] = {
        "parent_object": "all-breast-cells.h5ad",
        "parent_n_obs": int(adata.shape[0]),
        "sampling_method": "geometric_sketch",
        "sampling_library": "geosketch",
        "sampling_library_version": geosketch.__version__,
        "embedding_used": "X_scvi_100",
        "cells_per_donor": args.cells_per_donor,
        "random_state": args.seed,
        "actual_n_obs": int(sketch.shape[0]),
        "n_donors": int(sketch.obs["donor_id"].nunique()),
        "donors_below_threshold": donors_below,
    }

    # ── Verify integrity ──
    print("\nIntegrity checks:")
    assert sketch.shape[1] == adata.shape[1], "Gene count mismatch"
    print(f"  var shape matches parent: {sketch.shape[1]} genes")
    assert set(sketch.obs.columns) == set(adata.obs.columns), "obs columns mismatch"
    print(f"  obs columns match parent: {len(sketch.obs.columns)}")
    assert set(sketch.obsm.keys()) == set(adata.obsm.keys()), "obsm keys mismatch"
    print(f"  obsm keys match parent: {list(sketch.obsm.keys())}")
    assert set(sketch.layers.keys()) == set(adata.layers.keys()), "layers mismatch"
    print(f"  layers match parent: {list(sketch.layers.keys())}")
    n_donors = sketch.obs["donor_id"].nunique()
    assert n_donors == adata.obs["donor_id"].nunique(), "Donor count mismatch"
    print(f"  All {n_donors} donors represented")
    assert (sketch.obs["is_primary_data"] == False).all(), "is_primary_data not False"
    print(f"  is_primary_data = False for all cells")

    # ── Write to /tmp then copy (CRSP stale file handle prevention) ──
    tmp_path = f"/tmp/sketch_{os.getpid()}.h5ad"
    print(f"\nWriting sketch to {tmp_path}...")
    sketch.write_h5ad(tmp_path)
    tmp_size = os.path.getsize(tmp_path)
    print(f"  Written: {tmp_size / 1e9:.1f} GB")

    print(f"Copying to final output: {args.output}")
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    shutil.copy(tmp_path, args.output)
    os.remove(tmp_path)

    # Verify copy
    time.sleep(3)  # NFS cache lag
    final_size = os.path.getsize(args.output)
    assert final_size == tmp_size, f"Copy size mismatch: {final_size} vs {tmp_size}"
    print(f"  Verified: {final_size / 1e9:.1f} GB")

    t_total = time.time() - t0
    print(f"\nDone in {t_total:.0f}s ({t_total/60:.1f}m)")
    print(f"\nSummary:")
    print(f"  Parent: {adata.shape[0]:,} cells")
    print(f"  Sketch: {sketch.shape[0]:,} cells ({sketch.shape[0]/adata.shape[0]*100:.1f}%)")
    print(f"  Donors below threshold ({args.cells_per_donor}): {len(donors_below)}")
    print(f"  Output: {args.output}")


if __name__ == "__main__":
    main()
