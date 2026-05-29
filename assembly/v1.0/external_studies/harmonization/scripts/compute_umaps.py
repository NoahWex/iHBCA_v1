#!/usr/bin/env python3
"""
compute_umaps.py — Compute joint UMAP from embedding CSV.

Computes UMAP from joint scVI embedding (100D) for studies that don't have
native UMAP. Uses scanpy for consistency with existing integration pipeline.

Usage:
    python compute_umaps.py --study gray
    python compute_umaps.py --study gray --base-path /path/to/external_studies

Parameters match 03_Integration pattern (5_generate_embeddings_v2.py):
    n_neighbors = 30
    min_dist = 0.3
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

# UMAP computes best in silence, unlike grad students
sc.settings.verbosity = 1


def compute_umap(
    embedding_path: Path,
    output_path: Path,
    n_neighbors: int = 30,
    min_dist: float = 0.3
) -> None:
    """Compute UMAP from embedding CSV and save to output path."""
    print(f"Loading embedding: {embedding_path}")
    emb_df = pd.read_csv(embedding_path, index_col=0)
    print(f"  Shape: {emb_df.shape[0]:,} cells x {emb_df.shape[1]} dims")

    # Create minimal AnnData with embedding in obsm
    adata = ad.AnnData(X=np.zeros((len(emb_df), 1)))
    adata.obs_names = emb_df.index
    adata.obsm['X_latent'] = emb_df.values

    # Compute neighbors and UMAP
    print(f"Computing neighbors (n_neighbors={n_neighbors})...")
    sc.pp.neighbors(adata, use_rep='X_latent', n_neighbors=n_neighbors)

    print(f"Computing UMAP (min_dist={min_dist})...")
    sc.tl.umap(adata, min_dist=min_dist)

    # Save UMAP coordinates
    umap_df = pd.DataFrame(
        adata.obsm['X_umap'],
        index=emb_df.index,
        columns=['UMAP_1', 'UMAP_2']
    )
    umap_df.to_csv(output_path)
    print(f"Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Compute joint UMAP from embedding CSV"
    )
    parser.add_argument(
        "--study",
        required=True,
        help="Study name (gray, kumar, murrow, nee, twigger, pal_norm_*, reed)"
    )
    parser.add_argument(
        "--base-path",
        default="${SOURCE_IHBCAV1_EXTERNAL_STUDIES%/}",
        help="Base path to external_studies directory"
    )
    parser.add_argument(
        "--n-neighbors",
        type=int,
        default=30,
        help="Number of neighbors for UMAP [default: 30]"
    )
    parser.add_argument(
        "--min-dist",
        type=float,
        default=0.3,
        help="Minimum distance for UMAP [default: 0.3]"
    )
    args = parser.parse_args()

    # Resolve paths
    # Input: from published (where embedding_joint.csv already is)
    pub_dir = Path(args.base_path) / "outputs" / args.study / "published"
    embedding_path = pub_dir / "embedding_joint.csv"

    # Output: to study_objects (consistent with R script flow)
    study_objects_dir = Path(args.base_path) / "harmonization" / "outputs" / "study_objects" / args.study
    study_objects_dir.mkdir(parents=True, exist_ok=True)
    output_path = study_objects_dir / "umap_joint.csv"

    # Validate input exists
    if not embedding_path.exists():
        print(f"ERROR: No joint embedding found for {args.study}")
        print(f"  Expected: {embedding_path}")
        sys.exit(1)

    print(f"=== compute_umaps.py ===")
    print(f"  Study: {args.study}")
    print(f"  Input: {embedding_path}")
    print(f"  Output: {output_path}")
    print()

    compute_umap(
        embedding_path,
        output_path,
        n_neighbors=args.n_neighbors,
        min_dist=args.min_dist
    )

    print()
    print(f"=== Done: {args.study} joint UMAP ===")


if __name__ == "__main__":
    main()
