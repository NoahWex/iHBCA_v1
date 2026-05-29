"""01_cluster_compartment.py — extract / compute per-compartment Leiden clustering across resolutions.

Hybrid mode: for each requested resolution, extract from obs.csv if a
`leiden_<res>` column already exists, otherwise compute via scanpy on the
50D scVI latent (n_neighbors from --n-neighbors).

The scvi_n100 bundle's obs.csv typically has leiden_{0.3, 0.5, 0.8, 1.0, 5.0}
precomputed. To match V1's broader exploration range we also sweep {1.5, 3.0}
(fills the gap between anchor and fine_cluster). Track C uses the full set so
downstream can compare per-cell-type best-resolution (V1 pattern: cluster-mode
labels at low res + fine_cluster overrides at high res).

Output: `{out_dir}/{Compartment}_clusters.csv` — one row per cell, columns:
    cell_id, sample_id, patient_id, position, leiden_<res> (one per requested resolution)

Pattern: V1 `pseudobulk_aggregate_full_object.py` reads precomputed leiden
mappings; recompute path mirrors `scanpy.tl.leiden` standard usage.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib_config")

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("cluster_compartment")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--compartment", required=True, choices=["Epithelial", "Immune", "Stromal"])
    p.add_argument("--bundle-dir", required=True, type=Path,
                   help="Path to scvi_n100 bundle dir (contains obs.csv, latent.csv)")
    p.add_argument("--resolutions", nargs="+", required=True,
                   help="Resolution suffixes; extracted from obs.csv if leiden_<res> exists, "
                        "otherwise computed via scanpy.tl.leiden on latent.csv")
    p.add_argument("--n-neighbors", type=int, default=30,
                   help="n_neighbors for scanpy.pp.neighbors when recomputing (default 30)")
    p.add_argument("--n-cells", type=int, default=None,
                   help="Subsample to N cells for subset test (random_state=42)")
    p.add_argument("--out-dir", required=True, type=Path)
    return p.parse_args()


def compute_leiden(latent_df: pd.DataFrame, resolution: float, n_neighbors: int) -> pd.Series:
    """Run leiden clustering on the latent embedding. Returns Series indexed by cell_id.

    Uses sklearn.NearestNeighbors + igraph + leidenalg directly (no scanpy/pynndescent
    to avoid numba caching against the read-only system python at /opt/apps).
    Equivalent in result to scanpy.tl.leiden default usage:
      kNN graph (undirected, simplified) → RBConfigurationVertexPartition
      with resolution_parameter and seed=42.
    """
    from sklearn.neighbors import NearestNeighbors
    import igraph as ig
    import leidenalg

    log.info("    leiden via sklearn+igraph+leidenalg: resolution=%.3f n_neighbors=%d",
             resolution, n_neighbors)
    X = latent_df.values.astype(np.float32)
    n_cells = X.shape[0]

    log.info("    building kNN graph (k=%d)...", n_neighbors)
    nn = NearestNeighbors(n_neighbors=n_neighbors + 1, n_jobs=-1).fit(X)
    _, indices = nn.kneighbors(X)

    # Edges: (i, j) for j in neighbors of i (skip self at column 0)
    src = np.repeat(np.arange(n_cells), n_neighbors)
    dst = indices[:, 1:].ravel()
    edges = list(zip(src.tolist(), dst.tolist()))

    log.info("    constructing igraph (n=%d cells, %d directed edges → undirected)",
             n_cells, len(edges))
    g = ig.Graph(n=n_cells, edges=edges, directed=True)
    g = g.as_undirected(mode="collapse")
    g.simplify(combine_edges="first")

    log.info("    running leidenalg.find_partition (RBConfiguration, seed=42)...")
    partition = leidenalg.find_partition(
        g,
        leidenalg.RBConfigurationVertexPartition,
        resolution_parameter=resolution,
        seed=42,
    )

    membership = [str(c) for c in partition.membership]
    out = pd.Series(membership, index=latent_df.index, name="cluster")
    log.info("    %d clusters found (modularity=%.4f)",
             len(set(membership)), partition.modularity)
    return out


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    obs_csv = args.bundle_dir / "obs.csv"
    if not obs_csv.exists():
        sys.exit(f"ERROR: obs.csv not found at {obs_csv}")

    log.info("Loading obs: %s", obs_csv)
    obs = pd.read_csv(obs_csv)
    log.info("  %d cells, columns: %s", len(obs), list(obs.columns))

    required = {"cell_id", "sample_id", "patient_id"}
    missing = required - set(obs.columns)
    if missing:
        sys.exit(f"ERROR: obs.csv missing required columns: {missing}")

    if args.n_cells is not None and args.n_cells < len(obs):
        log.info("Subsampling to %d cells (seed=42)", args.n_cells)
        obs = obs.sample(n=args.n_cells, random_state=42).reset_index(drop=True)

    base_cols = ["cell_id", "sample_id", "patient_id"]
    if "position" in obs.columns:
        base_cols.append("position")
    out = obs[base_cols].copy()

    # Partition resolutions into (extract, recompute)
    leiden_cols = [f"leiden_{r}" for r in args.resolutions]
    extract_pairs = [(r, c) for r, c in zip(args.resolutions, leiden_cols) if c in obs.columns]
    recompute_pairs = [(r, c) for r, c in zip(args.resolutions, leiden_cols) if c not in obs.columns]
    log.info("Resolutions to extract from obs.csv: %s",
             [r for r, _ in extract_pairs])
    log.info("Resolutions to recompute via scanpy: %s",
             [r for r, _ in recompute_pairs])

    for res, col in extract_pairs:
        out[col] = obs[col].astype(str).values

    if recompute_pairs:
        latent_csv = args.bundle_dir / "latent.csv"
        if not latent_csv.exists():
            sys.exit(f"ERROR: latent.csv required for recompute but not found at {latent_csv}")
        log.info("Loading latent: %s", latent_csv)
        latent = pd.read_csv(latent_csv, index_col=0)
        log.info("  latent shape: %s", latent.shape)
        # Restrict latent to subsampled cells (if applicable)
        latent = latent.loc[obs["cell_id"].values]
        for res, col in recompute_pairs:
            log.info("Computing %s ...", col)
            cluster = compute_leiden(latent, float(res), args.n_neighbors)
            cluster = cluster.reindex(obs["cell_id"].values).reset_index(drop=True)
            out[col] = cluster.values

    log.info("Cluster counts per resolution:")
    for c in leiden_cols:
        col = out[c]
        log.info("  %s: %d clusters (top 5 sizes: %s)",
                 c, int(col.nunique()),
                 col.value_counts().head(5).to_dict())

    out_csv = args.out_dir / f"{args.compartment}_clusters.csv"
    out.to_csv(out_csv, index=False)
    log.info("Wrote %s (%d rows × %d cols)", out_csv, len(out), out.shape[1])


if __name__ == "__main__":
    main()
