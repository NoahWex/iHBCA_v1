"""Phase 6: 3-axis Xenium filter.

Axis 1: FLEX-kNN proximity — mean_knn_dist_to_flex <= threshold (FLEX self-kNN qN)
Axis 2: NMP rule — nmp_nuclear > nmp_cyto → EXCLUDE (so pass_nmp = nuclear <= cyto)
Axis 3: noise floor — pass_qc_whole from Phase 1

Output: three_axis_filter.csv with per-axis booleans + passes_all
"""
import argparse, os
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--proximity", required=True,
                    help="xenium_to_flex_knn.csv from 01_flex_proximity")
    ap.add_argument("--flex-self", required=True,
                    help="flex_self_knn.csv for threshold")
    ap.add_argument("--proximity-q", type=float, default=95.0,
                    help="Percentile of FLEX-self mean_knn_dist used as threshold")
    ap.add_argument("--nmp-csv", required=True)
    ap.add_argument("--qc-csv", required=True)
    ap.add_argument("--out-csv", required=True)
    args = ap.parse_args()
    os.makedirs(os.path.dirname(args.out_csv), exist_ok=True)

    prox = pd.read_csv(args.proximity)
    flex_self = pd.read_csv(args.flex_self)
    thr = float(np.percentile(flex_self["mean_knn_dist"], args.proximity_q))
    print(f"Axis 1 threshold (FLEX self q{args.proximity_q}): {thr:.4f}", flush=True)

    nmp = pd.read_csv(args.nmp_csv, usecols=["cell_id", "nmp_nuclear", "nmp_cyto"])
    qc = pd.read_csv(args.qc_csv, usecols=["cell_id", "pass_qc_whole"])

    df = prox.merge(nmp, on="cell_id", how="left").merge(qc, on="cell_id", how="left")
    print(f"Merged rows: {len(df)}", flush=True)

    df["passes_proximity"] = df["mean_knn_dist_to_flex"] <= thr
    df["passes_nmp"] = df["nmp_nuclear"] <= df["nmp_cyto"]
    df["passes_qc_whole"] = df["pass_qc_whole"].astype("boolean").fillna(False)
    df["passes_all"] = df["passes_proximity"] & df["passes_nmp"] & df["passes_qc_whole"]

    out = df[["cell_id", "mean_knn_dist_to_flex", "nmp_nuclear", "nmp_cyto",
              "pass_qc_whole", "passes_proximity", "passes_nmp",
              "passes_qc_whole", "passes_all"]]
    out.to_csv(args.out_csv, index=False)
    print(f"Wrote {args.out_csv} shape={out.shape}", flush=True)

    n = len(out)
    print(f"\n=== Per-axis pass rates (n={n}) ===", flush=True)
    for a in ["passes_proximity", "passes_nmp", "passes_qc_whole", "passes_all"]:
        p = int(out[a].sum())
        print(f"  {a}: {p}/{n} ({100*p/n:.1f}%)", flush=True)

    print("\n=== Axis intersection (Venn-ish) ===", flush=True)
    combos = out.groupby(["passes_proximity", "passes_nmp", "passes_qc_whole"]).size()
    print(combos.to_string(), flush=True)


if __name__ == "__main__":
    main()
