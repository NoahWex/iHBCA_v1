"""L0.5 x Leiden contingency on FLEX joint-latent clusters.

For each resolution: build contingency table, compute per-cluster purity
(max row fraction) and entropy (Shannon), plus per-L0.5 recall (best-cluster
capture fraction).

Outputs:
  contingency_res{r}.csv       rows=L0.5, cols=Leiden cluster, values=cell count
  cluster_summary_res{r}.csv   cluster, size, top_l0p5, purity, entropy
  l0p5_recall_res{r}.csv       l0p5, size, best_cluster, n_in_best, recall
"""
import argparse, os
import numpy as np
import pandas as pd


def entropy(p):
    p = p[p > 0]
    return float(-np.sum(p * np.log2(p)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--leiden-csv", required=True)
    ap.add_argument("--l0p5-csv", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    lei = pd.read_csv(args.leiden_csv, index_col="cell_id")
    l0p5 = pd.read_csv(args.l0p5_csv).set_index("cell_id")["l0p5"]
    common = lei.index.intersection(l0p5.index)
    lei = lei.loc[common]; l0p5 = l0p5.loc[common]
    print(f"common cells: {len(common)}", flush=True)

    res_cols = [c for c in lei.columns if c.startswith("leiden_")]
    for col in res_cols:
        r = col.split("_", 1)[1]
        ct = pd.crosstab(l0p5, lei[col])
        ct.to_csv(os.path.join(args.out_dir, f"contingency_res{r}.csv"))

        # per-cluster: purity + entropy
        col_tot = ct.sum(axis=0)
        col_frac = ct.div(col_tot, axis=1)
        cluster_rows = []
        for cl in ct.columns:
            p = col_frac[cl].values
            top = ct[cl].idxmax()
            cluster_rows.append({
                "cluster": cl, "size": int(col_tot[cl]),
                "top_l0p5": top,
                "purity": float(p.max()),
                "entropy": entropy(p),
            })
        pd.DataFrame(cluster_rows).sort_values(
            "size", ascending=False
        ).to_csv(os.path.join(args.out_dir, f"cluster_summary_res{r}.csv"),
                 index=False)

        # per-L0.5: best-cluster capture
        row_tot = ct.sum(axis=1)
        rec_rows = []
        for lab in ct.index:
            best = ct.loc[lab].idxmax()
            rec_rows.append({
                "l0p5": lab, "size": int(row_tot[lab]),
                "best_cluster": best,
                "n_in_best": int(ct.loc[lab, best]),
                "recall": float(ct.loc[lab, best] / row_tot[lab]),
            })
        pd.DataFrame(rec_rows).sort_values(
            "size", ascending=False
        ).to_csv(os.path.join(args.out_dir, f"l0p5_recall_res{r}.csv"),
                 index=False)
        print(f"res={r}: {len(ct.columns)} clusters, "
              f"median purity={np.median([r['purity'] for r in cluster_rows]):.3f}",
              flush=True)

    print("Done.", flush=True)


if __name__ == "__main__":
    main()
