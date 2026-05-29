"""Phase 6 Axis 1: FLEX-kNN proximity in joint_nuc_n100 latent.

Computes per-cell mean distance to k nearest FLEX neighbors (cosine).
- FLEX self-kNN: distribution used to set the threshold
- Xenium to FLEX: per-Xenium-cell proximity score

"""
import argparse, os
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors


def l2norm(X):
    n = np.linalg.norm(X, axis=1, keepdims=True)
    return X / np.where(n > 0, n, 1.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--latent", required=True)
    ap.add_argument("--obs", required=True)
    ap.add_argument("--k", type=int, default=30)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print("Loading latent...", flush=True)
    L = pd.read_csv(args.latent).set_index("cell_id")
    obs = pd.read_csv(args.obs).set_index("cell_id").reindex(L.index)
    print(f"  latent: {L.shape}, platforms: {obs['platform'].value_counts().to_dict()}", flush=True)

    X = L.values.astype(np.float32)
    Xn = l2norm(X)
    flex_mask = (obs["platform"] == "flex").values
    xen_mask = (obs["platform"] == "xenium").values
    Xf = Xn[flex_mask]
    Xx = Xn[xen_mask]
    flex_ids = L.index[flex_mask]
    xen_ids = L.index[xen_mask]
    print(f"  flex: {Xf.shape[0]}, xenium: {Xx.shape[0]}", flush=True)

    print(f"Fitting NN on FLEX (k={args.k})...", flush=True)
    nn = NearestNeighbors(n_neighbors=args.k + 1, metric="euclidean", n_jobs=-1)
    nn.fit(Xf)

    print("FLEX self-kNN...", flush=True)
    dflex, _ = nn.kneighbors(Xf)
    dflex = dflex[:, 1:]  # drop self
    flex_mean_knn = dflex.mean(axis=1)
    flex_max_knn = dflex.max(axis=1)

    pd.DataFrame({
        "cell_id": flex_ids,
        "mean_knn_dist": flex_mean_knn,
        "max_knn_dist": flex_max_knn,
    }).to_csv(os.path.join(args.out_dir, "flex_self_knn.csv"), index=False)

    print("Xenium → FLEX kNN...", flush=True)
    dxen, _ = nn.kneighbors(Xx, n_neighbors=args.k)
    xen_mean_knn = dxen.mean(axis=1)
    xen_max_knn = dxen.max(axis=1)

    pd.DataFrame({
        "cell_id": xen_ids,
        "mean_knn_dist_to_flex": xen_mean_knn,
        "max_knn_dist_to_flex": xen_max_knn,
    }).to_csv(os.path.join(args.out_dir, "xenium_to_flex_knn.csv"), index=False)

    qs = [1, 5, 25, 50, 75, 95, 99]
    fs = np.percentile(flex_mean_knn, qs)
    xs = np.percentile(xen_mean_knn, qs)
    summary = pd.DataFrame({"q": qs, "flex_self_mean_knn": fs, "xen_to_flex_mean_knn": xs})
    summary.to_csv(os.path.join(args.out_dir, "proximity_summary.csv"), index=False)

    print("\nFLEX self-kNN distance quantiles (cosine, L2-normed Euclidean):", flush=True)
    print(summary.to_string(index=False), flush=True)

    thr_q95 = float(np.percentile(flex_mean_knn, 95))
    thr_q99 = float(np.percentile(flex_mean_knn, 99))
    pass_q95 = int((xen_mean_knn <= thr_q95).sum())
    pass_q99 = int((xen_mean_knn <= thr_q99).sum())
    n_xen = len(xen_mean_knn)
    print(f"\nXenium pass at FLEX-q95 ({thr_q95:.4f}): {pass_q95}/{n_xen} ({100*pass_q95/n_xen:.1f}%)", flush=True)
    print(f"Xenium pass at FLEX-q99 ({thr_q99:.4f}): {pass_q99}/{n_xen} ({100*pass_q99/n_xen:.1f}%)", flush=True)


if __name__ == "__main__":
    main()
