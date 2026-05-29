"""Joint-latent kNN vote: assign each Xenium cell an L0.5 label by cosine kNN
to FLEX cells in the Phase 2c winner joint latent (joint_nuc_n100).

Output schema (cell_id ordered by Xenium obs):
  cell_id, knn_label, knn_confidence, knn_runner_up, knn_runner_frac

confidence = fraction of k neighbors voting the winning label
runner_frac = fraction voting the second-place label
"""
import argparse, os
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from collections import Counter


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--latent", required=True)
    ap.add_argument("--obs", required=True)
    ap.add_argument("--l0p5-csv", required=True, help="FLEX L0.5 per cell")
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("-k", type=int, default=30)
    args = ap.parse_args()
    os.makedirs(os.path.dirname(args.out_csv), exist_ok=True)

    print("Loading latent + obs...", flush=True)
    lat = pd.read_csv(args.latent, index_col=0)
    obs = pd.read_csv(args.obs); obs.index = obs["cell_id"]
    obs = obs.reindex(lat.index)

    flex_mask = obs["platform"] == "flex"
    lat_f = lat[flex_mask].astype(np.float32).values
    lat_x = lat[~flex_mask].astype(np.float32).values
    flex_ids = lat.index[flex_mask]
    xen_ids = lat.index[~flex_mask]
    print(f"FLEX: {lat_f.shape}  XENIUM: {lat_x.shape}", flush=True)

    # normalize for cosine via Euclidean
    lat_f /= (np.linalg.norm(lat_f, axis=1, keepdims=True) + 1e-12)
    lat_x /= (np.linalg.norm(lat_x, axis=1, keepdims=True) + 1e-12)

    l0p5 = pd.read_csv(args.l0p5_csv).set_index("cell_id")["l0p5"]
    labels = l0p5.reindex(flex_ids)
    miss = labels.isna().sum()
    print(f"FLEX cells missing L0.5: {miss}", flush=True)
    # fill missing with sentinel to avoid NaN propagation in vote
    labels = labels.fillna("_UNLABELED_").values

    print(f"Fitting NN on FLEX (k={args.k})...", flush=True)
    nn = NearestNeighbors(n_neighbors=args.k, metric="euclidean", n_jobs=-1)
    nn.fit(lat_f)
    print("Querying Xenium...", flush=True)
    _, idx = nn.kneighbors(lat_x)
    neigh_labels = labels[idx]  # (n_xen, k)

    print("Voting...", flush=True)
    lab_out = np.empty(len(xen_ids), dtype=object)
    conf_out = np.zeros(len(xen_ids), dtype=np.float32)
    ru_out = np.empty(len(xen_ids), dtype=object)
    ru_frac = np.zeros(len(xen_ids), dtype=np.float32)
    for i in range(len(xen_ids)):
        c = Counter(neigh_labels[i])
        top = c.most_common(2)
        lab_out[i] = top[0][0]
        conf_out[i] = top[0][1] / args.k
        if len(top) > 1:
            ru_out[i] = top[1][0]
            ru_frac[i] = top[1][1] / args.k
        else:
            ru_out[i] = ""
            ru_frac[i] = 0.0

    df = pd.DataFrame({
        "cell_id": xen_ids,
        "knn_label": lab_out,
        "knn_confidence": conf_out,
        "knn_runner_up": ru_out,
        "knn_runner_frac": ru_frac,
    })
    df.to_csv(args.out_csv, index=False)
    print(f"Wrote {args.out_csv}", flush=True)
    print(df["knn_label"].value_counts().to_string(), flush=True)
    print(f"mean confidence: {df['knn_confidence'].mean():.3f}", flush=True)


if __name__ == "__main__":
    main()
