#!/usr/bin/env python3
"""Compute per-cell DA logFC for the C5 epi_content contrast on the joint Milo.

Reuses the joint nhood_membership.csv (pre-computed for the fig 4 supp work).
Mirrors make_per_cell_da_flex.R behavior with the same output schema.

Output (per_cell_logfc.csv): cell_id, mean_logfc, n_nhoods, n_sig_pos,
                             n_sig_neg, status ("sig+", "sig-", "ns")
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--membership-csv", required=True,
                   help="All_nhood_membership.csv (cell_id, nhood_idx, ...)")
    p.add_argument("--da-csv", required=True,
                   help="da_results.csv for this scope")
    p.add_argument("--out-csv", required=True)
    p.add_argument("--sig-fdr", type=float, default=0.05)
    p.add_argument("--min-nhoods-for-sig", type=int, default=1)
    args = p.parse_args()

    print(f"[per-cell-da] loading membership: {args.membership_csv}")
    mem = pd.read_csv(args.membership_csv)
    # expected columns: cell_id, nhood_idx (1-based or 0-based; check first row)
    assert {"cell_id", "nhood_idx"}.issubset(mem.columns), \
        f"membership missing cell_id/nhood_idx: {mem.columns.tolist()}"

    print(f"[per-cell-da] loading DA: {args.da_csv}")
    da = pd.read_csv(args.da_csv)
    assert {"logFC", "SpatialFDR"}.issubset(da.columns)
    n_nhoods = len(da)

    # DA is 1-indexed when written by Miloreaching the bridge step
    nhood_min = mem["nhood_idx"].min()
    nhood_max = mem["nhood_idx"].max()
    print(f"  nhood_idx range: [{nhood_min}, {nhood_max}], da nhoods={n_nhoods}")
    # normalize to 0-based for vector indexing
    if nhood_min == 1 and nhood_max == n_nhoods:
        mem["nhood_idx"] = mem["nhood_idx"] - 1
    elif nhood_min == 0 and nhood_max == n_nhoods - 1:
        pass
    else:
        sys.exit(f"nhood_idx range doesn't match DA rows: got [{nhood_min},{nhood_max}], expected 0..{n_nhoods-1} or 1..{n_nhoods}")

    logfc = np.asarray(da["logFC"], dtype=float)
    fdr   = np.asarray(da["SpatialFDR"], dtype=float)
    sig   = (fdr < args.sig_fdr) & np.isfinite(fdr)
    sig_pos = sig & (logfc > 0)
    sig_neg = sig & (logfc < 0)

    # Vectorized per-cell stats via groupby
    nh_idx = np.asarray(mem["nhood_idx"], dtype=int)
    mem["lfc"]     = logfc[nh_idx]
    mem["is_sigp"] = sig_pos[nh_idx]
    mem["is_sign"] = sig_neg[nh_idx]

    out = mem.groupby("cell_id").agg(
        mean_logfc=("lfc", "mean"),
        n_nhoods=("nhood_idx", "count"),
        n_sig_pos=("is_sigp", "sum"),
        n_sig_neg=("is_sign", "sum"),
    ).reset_index()

    out["mean_logfc"] = out["mean_logfc"].round(5)
    out["n_sig_pos"]  = out["n_sig_pos"].astype(int)
    out["n_sig_neg"]  = out["n_sig_neg"].astype(int)

    pos_thresh = (out["n_sig_pos"] >= args.min_nhoods_for_sig) & \
                 (out["n_sig_pos"] > out["n_sig_neg"])
    neg_thresh = (out["n_sig_neg"] >= args.min_nhoods_for_sig) & \
                 (out["n_sig_neg"] > out["n_sig_pos"])
    out["status"] = np.where(pos_thresh, "sig+",
                    np.where(neg_thresh, "sig-", "ns"))

    Path(args.out_csv).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_csv, index=False)
    print(f"[per-cell-da] wrote {args.out_csv} ({len(out)} cells)")
    print(f"  status: sig+ {(out.status=='sig+').sum()}  "
          f"sig- {(out.status=='sig-').sum()}  "
          f"ns {(out.status=='ns').sum()}")
    print(f"  mean_logfc range: "
          f"[{out['mean_logfc'].min():.3f}, {out['mean_logfc'].max():.3f}]")


if __name__ == "__main__":
    main()
