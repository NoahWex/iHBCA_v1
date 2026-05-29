"""Leiden clustering on FLEX side of Phase 2c winner joint latent.

Input: joint_nuc_n100 latent + obs (filter to platform=='flex')
Build: PCA-free kNN graph (k=30 on joint-latent cosine), Leiden at res {0.1, 0.3, 0.5, 1.0}
Output: flex_leiden_joint.csv (cell_id, leiden_0.1, leiden_0.3, leiden_0.5, leiden_1.0)
"""
import argparse, os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import rapids_singlecell as rsc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--latent", required=True)
    ap.add_argument("--obs", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--k", type=int, default=30)
    ap.add_argument("--resolutions", nargs="+", type=float,
                    default=[0.1, 0.3, 0.5, 1.0])
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    lat = pd.read_csv(args.latent, index_col=0)
    obs = pd.read_csv(args.obs); obs.index = obs["cell_id"]
    obs = obs.reindex(lat.index)
    flex = obs["platform"] == "flex"
    lat_f = lat[flex].astype(np.float32)
    obs_f = obs[flex]
    print(f"FLEX cells: {lat_f.shape}", flush=True)

    a = ad.AnnData(obs=obs_f.copy())
    a.obsm["X_concord"] = lat_f.values
    print("rapids neighbors (GPU)...", flush=True)
    rsc.pp.neighbors(a, n_neighbors=args.k, use_rep="X_concord", metric="cosine")

    out = pd.DataFrame(index=lat_f.index)
    for r in args.resolutions:
        key = f"leiden_{r}"
        print(f"leiden res={r} (CPU, igraph)", flush=True)
        sc.tl.leiden(a, resolution=r, key_added=key,
                     flavor="igraph", n_iterations=2, directed=False)
        out[key] = a.obs[key].values
        print(f"  {key}: {out[key].nunique()} clusters", flush=True)

    out.index.name = "cell_id"
    out.to_csv(os.path.join(args.out_dir, "flex_leiden_joint.csv"))
    print(f"Wrote {args.out_dir}/flex_leiden_joint.csv", flush=True)


if __name__ == "__main__":
    main()
