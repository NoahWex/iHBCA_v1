#!/usr/bin/env python3
# 09f_smooth_full_l1.py — Cluster-vote smoothing for full-object L1 labels.
#
# Per-compartment harmonize covers cells that landed in Epi/Stro/Imm Seurat
# objects (~238K of the 258K canonical 5-way filtered set). The remaining
# ~19K cells passed canonical QC but were dropped by compartment classification
# / per-compartment Seurat builds. They are present in the full-object
# integration (latent + leiden) but lack an L1 label.
#
# This script smooths the gap by majority-voting within full-object leiden
# clusters. For each unlabeled cell, assigns the dominant L1 label among
# labeled cells in its leiden cluster.
#
# Output: harmonized_l1_labels_full.csv — superset of per-compartment file
# with smoothed labels for the gap cells, marked via L1.0_source.

import argparse
import pandas as pd
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--obs", required=True,
                    help="full obs.csv (cell_id, leiden_*, sample_id)")
    ap.add_argument("--labels", required=True,
                    help="harmonized_l1_labels.csv (per-comp concatenated)")
    ap.add_argument("--label-col", default="L1.0")
    ap.add_argument("--vote-leiden", default="leiden_2.0",
                    help="leiden column on full obs to vote within")
    ap.add_argument("--min-cluster-labeled", type=int, default=10,
                    help="min labeled cells in cluster to accept vote")
    ap.add_argument("--out", required=True,
                    help="output CSV (harmonized_l1_labels_full.csv)")
    args = ap.parse_args()

    obs = pd.read_csv(args.obs)
    lab = pd.read_csv(args.labels)
    print(f"obs: {len(obs):,}  labels: {len(lab):,}")
    print(f"vote leiden: {args.vote_leiden}  label col: {args.label_col}")

    # Restrict labels to cells present in obs.
    lab = lab[lab["cell_id"].isin(obs["cell_id"])].copy()
    print(f"labels after intersect with obs: {len(lab):,}")

    # Join leiden onto labels for per-cluster vote.
    obs_idx = obs.set_index("cell_id")[args.vote_leiden]
    lab["_leiden"] = lab["cell_id"].map(obs_idx)

    # Per-cluster majority vote.
    vote = (lab.dropna(subset=["_leiden"])
                .groupby("_leiden")[args.label_col]
                .agg(lambda s: s.value_counts().idxmax())
                .to_dict())
    counts = (lab.dropna(subset=["_leiden"])
                  .groupby("_leiden")
                  .size()
                  .to_dict())
    purity = (lab.dropna(subset=["_leiden"])
                  .groupby("_leiden")[args.label_col]
                  .agg(lambda s: s.value_counts().iloc[0] / len(s))
                  .to_dict())

    # Apply: any cell in obs not already labeled gets cluster vote.
    labeled_ids = set(lab["cell_id"])
    smoothed_rows = []
    for _, row in obs.iterrows():
        cid = row["cell_id"]
        if cid in labeled_ids:
            continue
        leiden = row[args.vote_leiden]
        if pd.isna(leiden) or counts.get(leiden, 0) < args.min_cluster_labeled:
            smoothed_rows.append({"cell_id": cid, args.label_col: pd.NA,
                                  "vote_purity": pd.NA,
                                  "L1.0_source": "unsmoothed",
                                  "leiden_cluster": leiden})
            continue
        smoothed_rows.append({
            "cell_id": cid,
            args.label_col: vote[leiden],
            "vote_purity": round(purity[leiden], 6),
            "L1.0_source": "leiden_vote",
            "leiden_cluster": leiden,
        })
    smoothed = pd.DataFrame(smoothed_rows)

    # Carry forward existing labels (preserve their source). If the input
    # labels lack vote_purity / L1.0_source, fill them as "compartment_harmonize".
    existing = lab[["cell_id", args.label_col]].copy()
    if "vote_purity" in lab.columns:
        existing["vote_purity"] = lab["vote_purity"]
    else:
        existing["vote_purity"] = pd.NA
    if "L1.0_source" in lab.columns:
        existing["L1.0_source"] = lab["L1.0_source"]
    else:
        existing["L1.0_source"] = "compartment_harmonize"
    existing["leiden_cluster"] = lab["_leiden"] if "_leiden" in lab.columns else pd.NA

    out = pd.concat([existing, smoothed], ignore_index=True)
    print(f"\nFinal label distribution:")
    print(out[args.label_col].value_counts(dropna=False))
    print(f"\nSource breakdown:")
    print(out["L1.0_source"].value_counts(dropna=False))

    smoothed_count = (out["L1.0_source"] == "leiden_vote").sum()
    unsmoothed = (out["L1.0_source"] == "unsmoothed").sum()
    coverage = (out[args.label_col].notna().sum() / len(obs)) * 100
    print(f"\nSmoothed: {smoothed_count:,}  unsmoothed: {unsmoothed:,}")
    print(f"Coverage: {out[args.label_col].notna().sum():,} / "
          f"{len(obs):,} ({coverage:.2f}%)")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    print(f"\nWrote: {args.out}")


if __name__ == "__main__":
    main()
