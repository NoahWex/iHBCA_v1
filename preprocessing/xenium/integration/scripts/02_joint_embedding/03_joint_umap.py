"""UMAP on the joint_nuc_n100 winner latent (FLEX + Xenium).

RAPIDS neighbors (GPU, fast on 1.4M cells) + scanpy CPU UMAP (umap-learn,
outlier-robust — RAPIDS cuml.UMAP produced ±2000-unit outliers from weakly
connected components in the kNN graph).

Pattern: concord eval `03_Integration/flex_xenium_concord_20260413/scripts/02_evaluate.py`
  sc.pp.neighbors + sc.tl.umap with cosine.

Outputs:
  joint_umap.csv            (cell_id, UMAP1, UMAP2) — full 1.38M cells
  flex_from_joint_umap.csv  FLEX subset
  xenium_from_joint_umap.csv Xenium subset
"""
import argparse, os
import numpy as np
import pandas as pd
import anndata as ad
import rapids_singlecell as rsc
import scanpy as sc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--latent", required=True)
    ap.add_argument("--obs", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("-k", type=int, default=50)
    ap.add_argument("--init", default="random", choices=["random", "spectral"])
    ap.add_argument("--min-dist", type=float, default=0.3)
    ap.add_argument("--spread", type=float, default=1.0)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    lat = pd.read_csv(args.latent, index_col=0)
    obs = pd.read_csv(args.obs); obs.index = obs["cell_id"]
    obs = obs.reindex(lat.index)
    print(f"latent: {lat.shape}, platforms: {obs['platform'].value_counts().to_dict()}",
          flush=True)

    a = ad.AnnData(X=np.zeros((lat.shape[0], 1), dtype=np.float32), obs=obs.copy())
    a.obsm["X_concord"] = lat.values.astype(np.float32)

    print("rapids neighbors (GPU)...", flush=True)
    rsc.pp.neighbors(a, n_neighbors=args.k, use_rep="X_concord", metric="cosine")
    print(f"scanpy UMAP (CPU, umap-learn) init={args.init} min_dist={args.min_dist} spread={args.spread}...", flush=True)
    sc.tl.umap(a, init_pos=args.init, min_dist=args.min_dist, spread=args.spread)

    umap = pd.DataFrame(a.obsm["X_umap"], index=lat.index,
                        columns=["UMAP1", "UMAP2"])
    umap.index.name = "cell_id"
    umap.to_csv(os.path.join(args.out_dir, "joint_umap.csv"))
    print(f"  joint_umap.csv: {umap.shape}", flush=True)

    flex_mask = obs["platform"] == "flex"
    umap[flex_mask].to_csv(os.path.join(args.out_dir, "flex_from_joint_umap.csv"))
    umap[~flex_mask].to_csv(os.path.join(args.out_dir, "xenium_from_joint_umap.csv"))
    print("Done.", flush=True)


if __name__ == "__main__":
    main()
