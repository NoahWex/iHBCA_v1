#!/usr/bin/env python3
"""
step_17_integration_sweep / postprocess_winner.py

Postprocess the selected winner integration: compute a tighter neighbor
graph, UMAP, and multi-resolution Leiden, and export flat CSV sidecars
that downstream annotation (25b_singler_l1.Rmd, 25c_l1_harmonize.Rmd)
consumes.

Sidecars produced in {--output-dir}:
  umap.csv                  cell_id, UMAP_1, UMAP_2
  leiden_multi.csv          cell_id + one column per resolution
  metadata_with_umap.csv    full obs + UMAP coords

Leiden resolutions are LOCKED via sweep_configs.yaml.leiden.postprocess_resolutions
to [0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0] (coordinator decision 2026-04-11).
The sweep-stage script uses a single resolution (leiden_scoring_resolution
= 1.0) for scIB scoring; the full multi-resolution set is computed only
here, on the winner, so we do not pay the cost 9 times.

"""

try:
    import numba
    _o = numba.njit
    def _nc(*a, **k):  # noqa: E306
        k.pop("cache", None)
        return _o(*a, **k)
    numba.njit = _nc
    _ov = numba.vectorize
    def _ncv(*a, **k):  # noqa: E306
        k.pop("cache", None)
        return _ov(*a, **k)
    numba.vectorize = _ncv
except ImportError:
    pass

import argparse
import os
import sys
import time

import pandas as pd
import scanpy as sc
import yaml


DEFAULT_CONFIG_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "configs", "sweep_configs.yaml",
)


def _die(msg, hint=None):
    sys.stderr.write(f"FATAL: {msg}\n")
    if hint:
        sys.stderr.write(f"  hint: {hint}\n")
    sys.exit(2)


def _load_config(path):
    if not os.path.exists(path):
        _die(f"sweep_configs.yaml not found: {path}")
    with open(path) as fh:
        return yaml.safe_load(fh)


def _pick_embed_key(adata, explicit=None):
    """Find the integration embedding key in obsm."""
    if explicit:
        if explicit not in adata.obsm:
            _die(f"--embed-key {explicit!r} not in obsm "
                 f"({list(adata.obsm.keys())})")
        return explicit
    for k in ["X_emb", "X_scVI", "X_pca"]:
        if k in adata.obsm:
            return k
    _die("No known embedding found in obsm",
         hint=f"keys={list(adata.obsm.keys())}; expected X_emb/X_scVI/X_pca")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True,
                        help="Winner integration H5AD "
                             "(from {sweep}/{target}/{config}/integrated.h5ad)")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for sidecars. Typically "
                             "CFG_CANONICAL_SWEEP_WINNER/{target}/.")
    parser.add_argument("--embed-key", default=None,
                        help="Embedding key in obsm. Auto-detect if omitted.")
    parser.add_argument("--config-path", default=DEFAULT_CONFIG_PATH,
                        help="Path to sweep_configs.yaml")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    print("=" * 70)
    print(f"step_17 postprocess_winner")
    print(f"  h5ad:       {args.h5ad}")
    print(f"  output_dir: {args.output_dir}")
    print(f"  config:     {args.config_path}")
    print("=" * 70)

    if not os.path.exists(args.h5ad):
        _die(f"winner H5AD not found: {args.h5ad}")

    cfg = _load_config(args.config_path)
    leiden_cfg = cfg.get("leiden", {})
    resolutions = leiden_cfg.get("postprocess_resolutions", [0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0])
    n_neighbors = cfg.get("defaults", {}).get("postprocess_n_neighbors", 15)
    print(f"  resolutions: {resolutions}")
    print(f"  n_neighbors: {n_neighbors}")

    if args.dry_run:
        print("=== DRY RUN MODE ===")
        print(f"Would read: {args.h5ad}")
        print(f"Would write: umap.csv, leiden_multi.csv, metadata_with_umap.csv "
              f"in {args.output_dir}")
        print("VALIDATION PASSED")
        return

    os.makedirs(args.output_dir, exist_ok=True)

    t0 = time.time()
    print(f"Loading: {args.h5ad}")
    adata = sc.read_h5ad(args.h5ad)
    print(f"  {adata.n_obs} cells x {adata.n_vars} genes "
          f"({time.time()-t0:.1f}s)")
    print(f"  obsm keys: {list(adata.obsm.keys())}")

    embed_key = _pick_embed_key(adata, explicit=args.embed_key)
    print(f"  Using embedding: {embed_key} "
          f"({adata.obsm[embed_key].shape[1]} dims)")

    print("Computing neighbors...")
    sc.pp.neighbors(adata, use_rep=embed_key, n_neighbors=n_neighbors)

    print("Computing UMAP...")
    sc.tl.umap(adata)
    umap_df = pd.DataFrame(
        adata.obsm["X_umap"],
        index=adata.obs_names,
        columns=["UMAP_1", "UMAP_2"],
    )
    umap_df.index.name = "cell_id"
    umap_path = os.path.join(args.output_dir, "umap.csv")
    umap_df.to_csv(umap_path)
    print(f"  Saved: {umap_path}")

    print("Computing multi-resolution Leiden...")
    leiden_df = pd.DataFrame(index=adata.obs_names)
    leiden_df.index.name = "cell_id"
    for res in resolutions:
        key = f"leiden_{res}"
        sc.tl.leiden(adata, resolution=res, key_added=key)
        n_clusters = adata.obs[key].nunique()
        leiden_df[key] = adata.obs[key].values
        print(f"  res={res}: {n_clusters} clusters")
    leiden_path = os.path.join(args.output_dir, "leiden_multi.csv")
    leiden_df.to_csv(leiden_path)
    print(f"  Saved: {leiden_path}")

    # Full metadata export (obs + UMAP)
    meta = adata.obs.copy()
    meta["UMAP_1"] = umap_df["UMAP_1"].values
    meta["UMAP_2"] = umap_df["UMAP_2"].values
    meta.index.name = "cell_id"
    meta_path = os.path.join(args.output_dir, "metadata_with_umap.csv")
    meta.to_csv(meta_path)
    print(f"  Saved: {meta_path}")

    print(f"\nDone ({time.time()-t0:.1f}s total)")


if __name__ == "__main__":
    main()
