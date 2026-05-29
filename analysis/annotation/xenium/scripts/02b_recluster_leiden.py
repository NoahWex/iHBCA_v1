"""Regenerate leiden_assignments.csv using CPU scanpy neighbors.

Uses sc.pp.neighbors (cosine, k=30) + sc.tl.leiden(random_state=42) to match
the neighbor graph used in 05_2x2_marker_comparison.py.  RAPIDS neighbors in
02_compartment_analysis.py produce a different graph → different cluster numbers.
This script replaces leiden_assignments.csv so cascade YAMLs (written from 05_2x2
HTML reports) apply correctly.

Usage:
    python 02b_recluster_leiden.py --compartment Epithelial \
        --latent <joint_latent.csv> --obs <joint_obs.csv> \
        --out-dir <Epithelial_clean/>
"""
import argparse, os
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

DEFAULT_RESOLUTIONS = [0.1, 0.3, 0.5, 0.7, 1.0, 1.5]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--compartment", required=True)
    ap.add_argument("--latent",      required=True, help="joint_latent.csv")
    ap.add_argument("--obs",         required=True, help="joint_obs.csv")
    ap.add_argument("--out-dir",     required=True, help="Compartment_clean/ root")
    ap.add_argument("-k", type=int, default=30)
    ap.add_argument("--resolutions", nargs="+", type=float, default=DEFAULT_RESOLUTIONS,
                    help="Leiden resolutions to compute (default: all 6)")
    args = ap.parse_args()

    comp = args.compartment
    print(f"[{comp}] Loading latent + obs...", flush=True)
    lat = pd.read_csv(args.latent, index_col=0)
    obs = pd.read_csv(args.obs).set_index("cell_id").reindex(lat.index)
    print(f"  {len(lat):,} cells", flush=True)

    a = ad.AnnData(obs=obs.copy())
    a.obsm["X_concord"] = lat.values.astype(np.float32)

    print(f"  CPU neighbors k={args.k} cosine...", flush=True)
    sc.pp.neighbors(a, n_neighbors=args.k, use_rep="X_concord", metric="cosine")

    for res in args.resolutions:
        key = f"leiden_{res}"
        print(f"  Leiden res={res}...", flush=True)
        sc.tl.leiden(a, resolution=res, key_added=key, random_state=42)
        a.obs[key] = a.obs[key].astype(str)
        print(f"    {a.obs[key].nunique()} clusters", flush=True)

    leiden_cols = [f"leiden_{r}" for r in args.resolutions]
    leiden_df = a.obs[leiden_cols].copy()
    leiden_df.index.name = "cell_id"
    out_path = os.path.join(args.out_dir, "leiden_assignments.csv")
    leiden_df.to_csv(out_path)
    print(f"  Written {len(leiden_df):,} rows → {out_path}", flush=True)


if __name__ == "__main__":
    main()
