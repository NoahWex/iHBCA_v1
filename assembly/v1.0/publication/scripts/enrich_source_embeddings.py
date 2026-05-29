#!/usr/bin/env python3
"""
enrich_source_embeddings.py - Port joint scVI embeddings to source h5ads
========================================================================
Reads integrated all-breast-cells.h5ad, extracts per-study X_scvi_100
embeddings, computes per-study UMAP, and writes enriched source h5ads.

Adds:
  obsm["X_ihbca_scvi_100"]       - 100D joint scVI from integrated object
  obsm["X_umap_ihbca_scvi_100"]  - per-study UMAP from joint embedding
  obs["in_ihbca_integrated"]     - boolean flag
  uns["obsm_descriptions"]       - provenance dict for all obsm keys

Removes:
  obsm["X_scVI_joint"]           - misleading, superseded by X_ihbca_scvi_100
  obsm["X_scVI_native"]          - misleading, original studies have different integrations

Usage:
  python enrich_source_embeddings.py \
    --study gray \
    --repo-root /path/to/iHBCAv1_upload

"""

import argparse
import shutil
import sys
import tempfile
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

# Match rate gates: script exits non-zero if below threshold
MATCH_GATES = {
    "gray": 0.995,
    "kumar": 0.995,
    "murrow": 0.995,
    "nee": 0.995,
    "twigger": 0.995,
    "reed": 0.995,
    "pal": 0.85,
}

# Mapping from integrated obs["dataset"] values to registry study names
# Author-share mode uses lowercase; legacy CxG uses capitalized
DATASET_TO_STUDY = {
    "Gray": "gray",
    "Kumar": "kumar",
    "Murrow": "murrow",
    "Nee": "nee",
    "Twigger": "twigger",
    "Reed": "reed",
    "Pal": "pal",
}

STUDY_TO_FILENAME = {
    "gray": "gray2022.h5ad",
    "kumar": "kumar2023.h5ad",
    "murrow": "murrow2022.h5ad",
    "nee": "nee2023.h5ad",
    "twigger": "twigger2022.h5ad",
    "reed": "reed2024.h5ad",
    "pal": "pal2021.h5ad",
}

OBSM_DESCRIPTIONS = {
    "X_umap": "Per-study UMAP from iHBCA harmonization pipeline (batch-corrected embedding basis)",
    "X_ihbca_scvi_100": (
        "100D joint scVI embedding from iHBCA integrated object "
        "(all-breast-cells.h5ad obsm['X_scvi_100'])"
    ),
    "X_umap_ihbca_scvi_100": (
        "Per-study UMAP computed from X_ihbca_scvi_100 "
        "(scanpy, n_neighbors=15, min_dist=0.5, random_state=42)"
    ),
}


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def load_integrated_embeddings(integrated_path):
    """Load X_scvi_100, obs_names, and obs['dataset'] from integrated h5ad.

    Uses backed='r' mode to avoid loading the full 2.12M-cell count matrix.

    Returns:
        embeddings: np.ndarray (N x 100, float32)
        obs_names: pd.Index of cell IDs
        dataset_series: pd.Series of study labels
    """
    print(f"Loading integrated embeddings: {integrated_path}")
    t0 = time.time()

    try:
        adata = ad.read_h5ad(integrated_path, backed="r")
        embeddings = adata.obsm["X_scvi_100"][:]  # force into memory
        obs_names = adata.obs_names.copy()
        dataset_series = adata.obs["dataset"].copy()
        adata.file.close()
    except Exception as e:
        print(f"  Backed mode failed ({e}), falling back to full load")
        adata = ad.read_h5ad(integrated_path)
        embeddings = adata.obsm["X_scvi_100"]
        obs_names = adata.obs_names.copy()
        dataset_series = adata.obs["dataset"].copy()
        del adata

    print(f"  Loaded {len(obs_names):,} cells x {embeddings.shape[1]}D "
          f"in {time.time() - t0:.1f}s")
    print(f"  Dataset values: {sorted(dataset_series.unique())}")

    return embeddings, obs_names, dataset_series


def extract_study_embeddings(embeddings, obs_names, dataset_series, study):
    """Extract embedding rows for a single study.

    Normalizes dataset labels via DATASET_TO_STUDY mapping with lowercase
    fallback (author-share mode already uses lowercase).

    Returns:
        study_embeddings: np.ndarray (N_study x 100, float32)
        study_cell_ids: pd.Index of cell IDs for this study
    """
    # Normalize dataset labels to lowercase study names
    normalized = dataset_series.map(DATASET_TO_STUDY).fillna(dataset_series.str.lower())
    mask = normalized == study
    n_cells = mask.sum()

    if n_cells == 0:
        raise ValueError(
            f"No cells found for study '{study}' in integrated object. "
            f"Dataset values: {sorted(dataset_series.unique())}"
        )

    study_cell_ids = obs_names[mask]
    study_embeddings = embeddings[np.asarray(mask)]

    print(f"  Extracted {n_cells:,} cells for study '{study}'")
    return study_embeddings, study_cell_ids


# ---------------------------------------------------------------------------
# UMAP computation
# ---------------------------------------------------------------------------


def compute_study_umap(embeddings, n_neighbors=15, min_dist=0.5, random_state=42):
    """Compute 2D UMAP from joint scVI embeddings for matched cells.

    Uses scanpy.pp.neighbors + scanpy.tl.umap on a temporary AnnData.

    Returns: np.ndarray (N_matched x 2, float32)
    """
    import scanpy as sc

    print(f"  Computing UMAP from {embeddings.shape[0]:,} cells "
          f"({embeddings.shape[1]}D -> 2D)")
    print(f"    n_neighbors={n_neighbors}, min_dist={min_dist}, "
          f"random_state={random_state}")
    t0 = time.time()

    tmp = ad.AnnData(
        X=np.zeros((embeddings.shape[0], 1), dtype=np.float32),
        obs=pd.DataFrame(index=range(embeddings.shape[0])),
    )
    tmp.obsm["X_scvi"] = embeddings.astype(np.float32)

    sc.pp.neighbors(tmp, use_rep="X_scvi", n_neighbors=n_neighbors)
    sc.tl.umap(tmp, min_dist=min_dist, random_state=random_state)

    umap_coords = tmp.obsm["X_umap"].astype(np.float32)
    print(f"    UMAP shape: {umap_coords.shape} in {time.time() - t0:.1f}s")

    return umap_coords


# ---------------------------------------------------------------------------
# Enrichment
# ---------------------------------------------------------------------------


def enrich_source_h5ad(source_path, study_embeddings, study_cell_ids,
                       study_umap, study, output_path=None, dry_run=False):
    """Enrich a source h5ad with joint embeddings and UMAP.

    - Matches cell IDs between source and integrated
    - Drops X_scVI_joint and X_scVI_native (misleading, superseded)
    - Stores X_ihbca_scvi_100 (100D, float32, NaN for unmatched)
    - Stores X_umap_ihbca_scvi_100 (2D, float32, NaN for unmatched)
    - Adds obs['in_ihbca_integrated'] boolean flag
    - Adds uns['obsm_descriptions'] provenance dict
    - Writes back to same path (or output_path if specified)

    Returns: dict with match stats
    """
    print(f"\nEnriching: {source_path}")
    t0 = time.time()

    adata = ad.read_h5ad(source_path)
    n_cells = adata.n_obs
    print(f"  Source cells: {n_cells:,}")
    print(f"  Existing obsm keys: {list(adata.obsm.keys())}")

    # --- Cell ID matching (vectorized) ---
    source_in_integrated = adata.obs_names.isin(study_cell_ids)
    matched_idx = np.where(source_in_integrated)[0]
    n_matched = len(matched_idx)
    n_unmatched = n_cells - n_matched
    match_pct = n_matched / max(n_cells, 1) * 100

    print(f"  Matched: {n_matched:,} / {n_cells:,} ({match_pct:.2f}%)")
    print(f"  Unmatched: {n_unmatched:,}")

    # Log unmatched cell IDs for Pal (or any study with >0 unmatched)
    if n_unmatched > 0:
        unmatched_ids = adata.obs_names[~source_in_integrated]
        if n_unmatched <= 50:
            print(f"  Unmatched IDs: {list(unmatched_ids)}")
        else:
            print(f"  Unmatched IDs (first 50): {list(unmatched_ids[:50])}")

    # --- Gate check ---
    gate = MATCH_GATES.get(study, 0.995)
    if match_pct / 100 < gate:
        print(f"  FATAL: Match rate {match_pct:.2f}% below gate "
              f"{gate * 100:.1f}% for {study}")
        sys.exit(1)

    stats = {
        "study": study,
        "source_cells": n_cells,
        "matched": n_matched,
        "unmatched": n_unmatched,
        "match_pct": match_pct,
    }

    if dry_run:
        print("  DRY RUN — not writing changes")
        return stats

    # --- Build index mapping ---
    int_id_to_row = pd.Series(
        range(len(study_cell_ids)), index=study_cell_ids
    )
    int_rows = int_id_to_row.loc[adata.obs_names[matched_idx]].values

    # --- Populate new obsm arrays ---
    emb_100d = np.full((n_cells, 100), np.nan, dtype=np.float32)
    umap_2d = np.full((n_cells, 2), np.nan, dtype=np.float32)
    emb_100d[matched_idx] = study_embeddings[int_rows]
    umap_2d[matched_idx] = study_umap[int_rows]

    # --- Drop misleading embeddings ---
    for key in ["X_scVI_joint", "X_scVI_native"]:
        if key in adata.obsm:
            del adata.obsm[key]
            print(f"  Dropped obsm['{key}']")

    # --- Write new data ---
    adata.obsm["X_ihbca_scvi_100"] = emb_100d
    adata.obsm["X_umap_ihbca_scvi_100"] = umap_2d
    adata.obs["in_ihbca_integrated"] = source_in_integrated
    adata.uns["obsm_descriptions"] = OBSM_DESCRIPTIONS

    # --- Verification output ---
    print(f"\n  obsm keys after enrichment:")
    for key in sorted(adata.obsm.keys()):
        arr = adata.obsm[key]
        n_nan = np.isnan(arr).any(axis=1).sum() if arr.dtype.kind == "f" else 0
        tag = "[NEW]" if key.startswith("X_ihbca") or key.startswith("X_umap_ihbca") else "[preserved]"
        print(f"    {key}: {arr.shape} {arr.dtype}  {tag}"
              f"{'  (' + str(n_nan) + ' NaN rows)' if n_nan > 0 else ''}")

    print(f"  obs['in_ihbca_integrated']: "
          f"{adata.obs['in_ihbca_integrated'].sum():,} True / "
          f"{(~adata.obs['in_ihbca_integrated']).sum():,} False")

    # --- Write back (CRSP safety: /tmp first, then copy) ---
    if output_path is None:
        output_path = source_path

    old_size = source_path.stat().st_size if source_path.exists() else 0

    tmp_dir = Path(tempfile.mkdtemp())
    tmp_path = tmp_dir / Path(output_path).name

    print(f"\n  Writing to /tmp...")
    adata.write_h5ad(tmp_path)
    tmp_size = tmp_path.stat().st_size

    print(f"  Copying to {output_path}...")
    shutil.copy(str(tmp_path), str(output_path))
    tmp_path.unlink()
    tmp_dir.rmdir()

    new_size = output_path.stat().st_size
    print(f"  File size: {old_size / 1e9:.2f} GB -> {new_size / 1e9:.2f} GB")
    if abs(new_size - tmp_size) > 1024:
        print(f"  WARNING: Copy size mismatch: tmp={tmp_size}, final={new_size}")

    print(f"  Enrichment complete in {time.time() - t0:.1f}s")
    return stats


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Port joint scVI embeddings from integrated to source h5ads"
    )
    parser.add_argument(
        "--study", required=True,
        choices=list(STUDY_TO_FILENAME.keys()),
        help="Study name"
    )
    parser.add_argument(
        "--repo-root", required=True,
        help="Path to iHBCAv1_upload repo root"
    )
    parser.add_argument(
        "--integrated",
        help="Override path to integrated h5ad (default: repo-root/.../all-breast-cells.h5ad)"
    )
    parser.add_argument(
        "--output",
        help="Override output h5ad path (default: overwrite source in place)"
    )
    parser.add_argument(
        "--n-neighbors", type=int, default=15,
        help="UMAP n_neighbors (default: 15)"
    )
    parser.add_argument(
        "--min-dist", type=float, default=0.5,
        help="UMAP min_dist (default: 0.5)"
    )
    parser.add_argument(
        "--random-state", type=int, default=42,
        help="UMAP random state (default: 42)"
    )
    parser.add_argument(
        "--source",
        help="Override source h5ad path (default: repo-root/.../source_datasets/<study>.h5ad)"
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Report matches without writing"
    )
    args = parser.parse_args()

    repo_root = Path(args.repo_root)
    study = args.study

    print("=" * 70)
    print(f"SOURCE EMBEDDING ENRICHMENT — {study}")
    print("=" * 70)

    # Resolve paths
    integrated_path = (
        Path(args.integrated) if args.integrated
        else repo_root / "publication/outputs/integrated_objects/all-breast-cells.h5ad"
    )
    source_path = (
        Path(args.source) if args.source
        else repo_root / "publication/outputs/source_datasets" / STUDY_TO_FILENAME[study]
    )
    output_path = Path(args.output) if args.output else None

    if not integrated_path.exists():
        print(f"FATAL: Integrated h5ad not found: {integrated_path}")
        sys.exit(1)
    if not source_path.exists():
        print(f"FATAL: Source h5ad not found: {source_path}")
        sys.exit(1)

    print(f"  Integrated: {integrated_path}")
    print(f"  Source:      {source_path}")
    print(f"  Output:      {output_path or '(overwrite source)'}")

    # 1. Load integrated embeddings
    print(f"\n[1/4] Loading integrated embeddings...")
    embeddings, obs_names, dataset_series = load_integrated_embeddings(integrated_path)

    # 2. Extract study-specific embeddings
    print(f"\n[2/4] Extracting embeddings for {study}...")
    study_emb, study_ids = extract_study_embeddings(
        embeddings, obs_names, dataset_series, study
    )

    # 3. Compute per-study UMAP
    print(f"\n[3/4] Computing UMAP...")
    study_umap = compute_study_umap(
        study_emb,
        n_neighbors=args.n_neighbors,
        min_dist=args.min_dist,
        random_state=args.random_state,
    )

    # 4. Enrich source h5ad
    print(f"\n[4/4] Enriching source h5ad...")
    stats = enrich_source_h5ad(
        source_path, study_emb, study_ids, study_umap, study,
        output_path=output_path,
        dry_run=args.dry_run,
    )

    # Summary
    print(f"\n{'=' * 70}")
    print(f"ENRICHMENT {'(DRY RUN) ' if args.dry_run else ''}COMPLETE: {study}")
    print(f"  Matched: {stats['matched']:,} / {stats['source_cells']:,} "
          f"({stats['match_pct']:.2f}%)")
    print(f"  Unmatched: {stats['unmatched']:,}")
    print(f"{'=' * 70}")

    sys.exit(0)


if __name__ == "__main__":
    main()
