#!/usr/bin/env python3
"""
Aggregate per-config scIB sidecar CSVs into a ranked summary.

Reads all per-config CSVs from a scib output directory (one row per config,
emitted by the 05_* sidecars and the legacy 05_scib_benchmark.py monolith)
and produces:

  1. scib_summary.csv  — long table: all configs × all metrics + composite
                          scores + per-compartment rank
  2. scib_winners.csv  — one row per compartment: best method by composite
                          score (the input contract for downstream consumers)

Composite scoring follows the scIB paper convention:
    Batch correction = mean(ASW_batch, graph_conn, PCR_batch, iLISI)
    Bio conservation = mean(NMI, ARI, ASW_label, isolated_labels_F1, cLISI)
    Overall          = 0.4 * batch + 0.6 * bio
Missing metrics are dropped from the mean (NaN-tolerant); the composite is
only emitted when at least one metric is available in each block.

Usage:
    python 06_scib_aggregate.py \\
        --scib-dir <scib_output_dir> \\
        --output-dir <scib_output_dir> \\
        [--batch-weight 0.4] [--bio-weight 0.6]
"""

import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd


# Metric groupings (match scIB convention; column names follow scib.metrics output)
BATCH_METRICS = [
    "ASW_batch", "ASW_label/batch",  # scib uses ASW_label/batch as the batch ASW
    "graph_conn",
    "PCR_batch",
    "iLISI",
]

BIO_METRICS = [
    "NMI_cluster/label",  # scib's NMI key
    "ARI_cluster/label",  # scib's ARI key
    "ASW_label",          # bio ASW
    "isolated_label_F1",
    "isolated_label_silhouette",
    "cLISI",
]

# Some scib versions use shorter keys; accept both.
ALIAS_MAP = {
    "NMI": "NMI_cluster/label",
    "ARI": "ARI_cluster/label",
    "isolated_labels_F1": "isolated_label_F1",
    "isolated_labels_asw": "isolated_label_silhouette",
}


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Apply alias map and ensure expected columns exist (NaN if missing)."""
    df = df.rename(columns=ALIAS_MAP)
    return df


def collect_metrics(df: pd.DataFrame, metric_list: list) -> pd.Series:
    """Return numeric values for the named metrics that actually exist in df."""
    available = [m for m in metric_list if m in df.columns]
    if not available:
        return pd.Series(dtype="float64")
    return df[available].apply(pd.to_numeric, errors="coerce").iloc[0]


def compute_composite(row_df: pd.DataFrame, batch_weight: float,
                      bio_weight: float) -> dict:
    """Compute batch / bio / overall composite for a single config row."""
    batch_vals = collect_metrics(row_df, BATCH_METRICS)
    bio_vals = collect_metrics(row_df, BIO_METRICS)

    batch_arr = np.asarray(batch_vals.to_numpy(), dtype=float)
    bio_arr = np.asarray(bio_vals.to_numpy(), dtype=float)
    batch_score = float(np.nanmean(batch_arr)) if batch_arr.size else np.nan
    bio_score = float(np.nanmean(bio_arr)) if bio_arr.size else np.nan

    if np.isnan(batch_score) and np.isnan(bio_score):
        overall = np.nan
    elif np.isnan(batch_score):
        overall = bio_score
    elif np.isnan(bio_score):
        overall = batch_score
    else:
        overall = batch_weight * batch_score + bio_weight * bio_score

    return {
        "batch_score": batch_score,
        "bio_score": bio_score,
        "overall_score": overall,
        "n_batch_metrics_used": int((~batch_vals.isna()).sum())
                                if len(batch_vals) else 0,
        "n_bio_metrics_used": int((~bio_vals.isna()).sum())
                              if len(bio_vals) else 0,
    }


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--scib-dir", required=True,
                        help="Directory containing per-config scIB CSVs")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--batch-weight", type=float, default=0.4)
    parser.add_argument("--bio-weight", type=float, default=0.6)
    parser.add_argument("--pattern", default="*.csv",
                        help="Glob pattern for input CSVs (default: *.csv)")
    args = parser.parse_args()

    sys.stdout.reconfigure(line_buffering=True)  # type: ignore[attr-defined]
    os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 72, flush=True)
    print("06_scib_aggregate", flush=True)
    print("=" * 72, flush=True)
    print(f"scib_dir:    {args.scib_dir}", flush=True)
    print(f"output_dir:  {args.output_dir}", flush=True)
    print(f"weights:     batch={args.batch_weight}, bio={args.bio_weight}",
          flush=True)
    print("", flush=True)

    csv_paths = sorted(glob.glob(os.path.join(args.scib_dir, args.pattern)))
    # Exclude our own outputs if --output-dir == --scib-dir
    csv_paths = [p for p in csv_paths
                 if not os.path.basename(p).startswith(("scib_summary",
                                                        "scib_winners"))]
    if not csv_paths:
        raise FileNotFoundError(
            f"No scIB CSVs found in {args.scib_dir} matching {args.pattern}"
        )
    print(f"Found {len(csv_paths)} per-config CSVs", flush=True)

    rows = []
    for path in csv_paths:
        try:
            df = pd.read_csv(path)
        except Exception as e:
            print(f"  WARN: failed to read {path}: {e}", flush=True)
            continue
        if df.empty:
            print(f"  WARN: empty CSV {path}", flush=True)
            continue
        df = normalize_columns(df)

        # Per-config error rows (when 05 hit an exception, it writes a row
        # with an "error" column instead of metrics)
        if "error" in df.columns and pd.notna(df.iloc[0].get("error", None)):
            method = df.iloc[0].get("method", os.path.basename(path).replace(".csv", ""))
            comp = df.iloc[0].get("compartment", "?")
            print(f"  ERROR config {method} ({comp}): "
                  f"{str(df.iloc[0]['error'])[:120]}", flush=True)
            rows.append({
                "method": method,
                "compartment": comp,
                "embedding_key": df.iloc[0].get("embedding_key", ""),
                "n_cells": df.iloc[0].get("n_cells", np.nan),
                "n_genes": df.iloc[0].get("n_genes", np.nan),
                "scib_version": df.iloc[0].get("scib_version", ""),
                "error": str(df.iloc[0]["error"])[:200],
                "batch_score": np.nan,
                "bio_score": np.nan,
                "overall_score": np.nan,
            })
            continue

        # Healthy row
        row = df.iloc[0].to_dict()
        composite = compute_composite(df, args.batch_weight, args.bio_weight)
        row.update(composite)
        rows.append(row)

    if not rows:
        raise RuntimeError("No usable rows aggregated — all CSVs failed.")

    summary = pd.DataFrame(rows)

    # Per-compartment rank by overall_score (higher = better; NaNs at bottom)
    if "compartment" in summary.columns and "overall_score" in summary.columns:
        summary["compartment_rank"] = (
            summary.groupby("compartment")["overall_score"]
                   .rank(ascending=False, method="min", na_option="bottom")
                   .astype("Int64")
        )
    else:
        summary["compartment_rank"] = pd.NA

    # Reorder columns: identifiers first, then composite, then raw metrics
    id_cols = ["compartment", "method", "embedding_key", "compartment_rank",
               "overall_score", "batch_score", "bio_score",
               "n_batch_metrics_used", "n_bio_metrics_used",
               "n_cells", "n_genes", "scib_version"]
    other_cols = [c for c in summary.columns if c not in id_cols]
    summary = summary[[c for c in id_cols if c in summary.columns] + other_cols]

    summary_path = os.path.join(args.output_dir, "scib_summary.csv")
    summary.to_csv(summary_path, index=False)
    print(f"\nWrote summary: {summary_path} ({len(summary)} rows)", flush=True)

    # --- Winners table: best per compartment ---
    if "compartment" in summary.columns and "overall_score" in summary.columns:
        # Build the scored subset via numpy boolean mask. pyright stubs make
        # pd.to_numeric/Series.notna chain ambiguous; explicit cast quiets it.
        score_arr = np.asarray(summary["overall_score"], dtype=float)
        scored: pd.DataFrame = summary.loc[~np.isnan(score_arr)].copy()
        scored = scored.sort_values(
            by=["compartment", "overall_score"],
            ascending=[True, False],
        )
        winners = scored.groupby("compartment", as_index=False).first()
        winners_path = os.path.join(args.output_dir, "scib_winners.csv")
        winners.to_csv(winners_path, index=False)
        print(f"Wrote winners: {winners_path} ({len(winners)} rows)",
              flush=True)
        print("\n=== Winners per compartment ===", flush=True)
        for _, w in winners.iterrows():
            print(f"  {w['compartment']:>11s}: {w['method']:>14s}  "
                  f"overall={w['overall_score']:.3f}  "
                  f"batch={w.get('batch_score', float('nan')):.3f}  "
                  f"bio={w.get('bio_score', float('nan')):.3f}",
                  flush=True)

    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
