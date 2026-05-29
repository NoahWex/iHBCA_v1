"""Axis 1 v4: per-L0.5 threshold = FLEX self-kNN 95th percentile within type.

Writes three_axis_filter_v4.csv with `passes_all` = prox_v4 & nmp & qc_whole
(column name matches Phase 7a train script; v1 passes_all preserved as
passes_all_v1 for audit).
"""
import argparse, os
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--flex-knn", required=True)
    ap.add_argument("--flex-l0p5", required=True)
    ap.add_argument("--xen-knn", required=True)
    ap.add_argument("--xen-l0p5", required=True)
    ap.add_argument("--three-axis", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--p", type=float, default=0.95)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    flex = (pd.read_csv(args.flex_knn)
              .merge(pd.read_csv(args.flex_l0p5)[["cell_id", "l0p5"]],
                     on="cell_id", how="left"))
    thr_by_type = (flex.groupby("l0p5")["mean_knn_dist"]
                       .quantile(args.p).rename("thr_v4"))
    global_fallback = float(flex["mean_knn_dist"].quantile(args.p))
    print(f"Per-L0.5 FLEX q{int(args.p*100)} thresholds:", flush=True)
    for k, v in thr_by_type.items():
        print(f"  {k}: {v:.4f}", flush=True)
    print(f"Global fallback (unmapped types): {global_fallback:.4f}", flush=True)

    xen = (pd.read_csv(args.xen_knn)
             .merge(pd.read_csv(args.xen_l0p5)[["cell_id", "l0p5_final"]],
                    on="cell_id", how="left"))
    xen["thr_v4"] = xen["l0p5_final"].map(thr_by_type).fillna(global_fallback)
    xen["passes_proximity_v4"] = (xen["mean_knn_dist_to_flex"] <=
                                  xen["thr_v4"]).astype(int)

    prior = pd.read_csv(args.three_axis)[["cell_id", "passes_proximity",
                                          "passes_nmp", "passes_qc_whole",
                                          "passes_all"]].rename(
        columns={"passes_all": "passes_all_v1",
                 "passes_proximity": "passes_proximity_v1"})
    out = prior.merge(xen[["cell_id", "mean_knn_dist_to_flex", "l0p5_final",
                           "thr_v4", "passes_proximity_v4"]],
                      on="cell_id", how="left")
    out["passes_proximity"] = out["passes_proximity_v4"].fillna(0).astype(int)
    out["passes_all"] = (out["passes_proximity"] &
                         out["passes_nmp"].astype(int) &
                         out["passes_qc_whole"].astype(int)).astype(int)
    out.to_csv(os.path.join(args.out_dir, "three_axis_filter_v4.csv"),
               index=False)

    n = len(out)
    v1p = int(out["passes_proximity_v1"].sum())
    v4p = int(out["passes_proximity"].sum())
    v1a = int(out["passes_all_v1"].sum())
    v4a = int(out["passes_all"].sum())
    print(f"\nXenium cells: {n}", flush=True)
    print(f"  proximity: v1={v1p} ({v1p/n:.1%}) -> v4={v4p} ({v4p/n:.1%})",
          flush=True)
    print(f"  all:       v1={v1a} ({v1a/n:.1%}) -> v4={v4a} ({v4a/n:.1%})",
          flush=True)

    byt = (out.groupby("l0p5_final", dropna=False)
              .agg(n=("cell_id", "size"),
                   pass_v1=("passes_proximity_v1", "sum"),
                   pass_v4=("passes_proximity", "sum"),
                   thr_v4=("thr_v4", "first"))
              .reset_index())
    byt["rate_v1"] = byt["pass_v1"] / byt["n"]
    byt["rate_v4"] = byt["pass_v4"] / byt["n"]
    byt.to_csv(os.path.join(args.out_dir, "per_l0p5_pass_rates_v4.csv"),
               index=False)
    print("\nPer-L0.5 pass rates (v4 = FLEX q95):\n" +
          byt.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
