"""Leiden clustering on the full FLEX + Xenium joint latent (joint_v6).

Extends 01_flex_leiden_joint.py: same neighbor-graph and Leiden settings, but
no platform filter — clusters all ~1.3M cells jointly. Output cluster IDs are
used by 04_joint_markers.R and 05_render_joint_heatmaps.R to produce the
three-platform marker-and-heatmap preview (Xenium nuclear, FLEX on Xenium
panel, FLEX full transcriptome).

Input:  joint_v6 latent (100-dim) + obs (cell_id, platform, ...)
Build:  PCA-free kNN graph on joint-latent cosine (RAPIDS, GPU)
Leiden: CPU igraph at resolutions {0.3, 0.5, 1.0, 1.5}
Output: joint_leiden_assignments.csv  (cell_id + leiden_{res} columns)
"""
import argparse
import os

import anndata as ad
import numpy as np
import pandas as pd
import rapids_singlecell as rsc
import scanpy as sc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--latent", required=True)
    ap.add_argument("--obs", required=True)
    ap.add_argument("--out-path", required=True,
                    help="Destination CSV (e.g. pipeline/outputs/embedding/joint_leiden_assignments.csv)")
    ap.add_argument("--k", type=int, default=30)
    ap.add_argument("--resolutions", nargs="+", type=float,
                    default=[0.3, 0.5, 1.0, 1.5])
    args = ap.parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(args.out_path)), exist_ok=True)

    lat = pd.read_csv(args.latent, index_col=0)
    obs = pd.read_csv(args.obs); obs.index = obs["cell_id"]
    obs = obs.reindex(lat.index)
    print(f"joint cells: {lat.shape}", flush=True)

    a = ad.AnnData(obs=obs.copy())
    a.obsm["X_concord"] = lat.astype(np.float32).values
    print("rapids neighbors (GPU)...", flush=True)
    rsc.pp.neighbors(a, n_neighbors=args.k, use_rep="X_concord", metric="cosine")

    out = pd.DataFrame(index=lat.index)
    for r in args.resolutions:
        key = f"leiden_{r}"
        print(f"leiden res={r} (CPU, igraph)", flush=True)
        sc.tl.leiden(a, resolution=r, key_added=key,
                     flavor="igraph", n_iterations=2, directed=False)
        out[key] = a.obs[key].values
        print(f"  {key}: {out[key].nunique()} clusters", flush=True)

    out.index.name = "cell_id"
    out.to_csv(args.out_path)
    print(f"Wrote {args.out_path}", flush=True)


if __name__ == "__main__":
    main()
