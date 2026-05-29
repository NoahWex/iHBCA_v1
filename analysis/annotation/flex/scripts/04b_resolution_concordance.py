"""04b_resolution_concordance.py — F1 / mean-score analysis across resolutions per V1 label.

For each V1 L2 label and each resolution, ask: how cleanly does some cluster at
this resolution capture this label?

Metric protocol:
    1. Per cell, compute combined UCell score for the label
       (identity * w_identity + canonical * w_canonical) / (w_identity + w_canonical)
    2. Define cell-level positives at a fixed UCell threshold (default 0.5).
    3. For each (resolution, cluster) pair, compute precision/recall/F1
       treating cluster membership as the predicted positive.
    4. Best F1 per (label, resolution) = maximum across clusters at that resolution.
    5. Recommend resolution range = resolutions where best-F1 is within ε of max
       (default ε = 0.05).

Also reports the cluster-mean combined score per (label, resolution) — useful when
the label is too rare for a F1-stable estimate (e.g., LASP-KIT at low resolution
where it gets lumped into LASP-major).

Inputs:
    --ucell-per-cell-csv  outputs/ucell/{Compartment}_ucell_per_cell.csv
    --clusters-csv        outputs/clusters/{Compartment}_clusters.csv (multi-res)
    --resolutions         e.g. 0.3 0.5 0.8 1.0 5.0
    --yaml                V1 yaml (for label list, lineage, source_resolution)
    --compartment         Epithelial / Immune / Stromal
    --weight-identity, --weight-canonical
    --positive-threshold  UCell combined score cutoff for cell-level positive (default 0.5)
    --recommend-eps       within-eps tolerance for recommended resolution range (default 0.05)
    --out-dir

Outputs:
    {Compartment}_resolution_concordance.csv
        long: label, resolution, best_cluster, n_cluster_cells, n_label_positives,
              precision, recall, f1, cluster_mean_combined
    {Compartment}_resolution_recommendations.csv
        per-label summary: label, lineage, v1_source_resolution,
                           best_resolution, best_f1, recommended_resolutions
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")

import numpy as np
import pandas as pd
import yaml

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("concordance")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ucell-per-cell-csv", required=True, type=Path)
    p.add_argument("--clusters-csv", required=True, type=Path)
    p.add_argument("--resolutions", nargs="+", required=True)
    p.add_argument("--yaml", required=True, type=Path)
    p.add_argument("--compartment", required=True, choices=["Epithelial", "Immune", "Stromal"])
    p.add_argument("--weight-identity", type=float, default=2.0)
    p.add_argument("--weight-canonical", type=float, default=1.0)
    p.add_argument("--positive-threshold", type=float, default=0.5)
    p.add_argument("--recommend-eps", type=float, default=0.05)
    p.add_argument("--out-dir", required=True, type=Path)
    return p.parse_args()


def load_yaml_label_meta(yaml_path: Path) -> pd.DataFrame:
    with open(yaml_path) as fh:
        yml = yaml.safe_load(fh)
    rows = []
    for label, spec in (yml.get("labels", {}) or {}).items():
        if spec.get("is_artifact"):
            continue
        rows.append({
            "label": label,
            "lineage": spec.get("lineage", ""),
            "v1_source_resolution": str(spec.get("source_resolution", "")),
            "v1_mode": spec.get("mode", ""),
        })
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    log.info("Loading per-cell UCell: %s", args.ucell_per_cell_csv)
    per_cell = pd.read_csv(args.ucell_per_cell_csv)
    if "cell_id" not in per_cell.columns:
        sys.exit("ERROR: per_cell CSV missing cell_id")

    sig_cols = [c for c in per_cell.columns if "__" in c]
    labels = sorted(set(c.rsplit("__", 1)[0] for c in sig_cols))
    log.info("  %d cells, %d labels", len(per_cell), len(labels))

    w_total = args.weight_identity + args.weight_canonical
    combined = pd.DataFrame({"cell_id": per_cell["cell_id"]})
    for lb in labels:
        canon = per_cell.get(f"{lb}__canonical", pd.Series(0.0, index=per_cell.index))
        ident = per_cell.get(f"{lb}__identity", pd.Series(0.0, index=per_cell.index))
        combined[lb] = (ident * args.weight_identity + canon * args.weight_canonical) / w_total

    log.info("Loading clusters: %s", args.clusters_csv)
    clusters = pd.read_csv(args.clusters_csv)
    if "cell_id" not in clusters.columns:
        sys.exit("ERROR: clusters CSV missing cell_id")
    leiden_cols = [f"leiden_{r}" for r in args.resolutions]
    missing = [c for c in leiden_cols if c not in clusters.columns]
    if missing:
        sys.exit(f"ERROR: clusters CSV missing leiden columns: {missing}")

    df = clusters[["cell_id"] + leiden_cols].merge(combined, on="cell_id", how="inner")
    log.info("  merged: %d cells", len(df))

    label_meta = load_yaml_label_meta(args.yaml)

    long_rows = []
    summary_rows = []

    for lb in labels:
        if lb not in df.columns:
            continue
        scores = df[lb].values
        positives = scores >= args.positive_threshold
        n_pos = int(positives.sum())
        if n_pos == 0:
            log.warning("  label %s: no cells above threshold %.2f; skipping",
                        lb, args.positive_threshold)
            continue

        best_per_res: dict[str, dict] = {}
        for res, leiden_col in zip(args.resolutions, leiden_cols):
            cluster_vals = df[leiden_col].astype(str).values
            cluster_unique = np.unique(cluster_vals)
            best = None
            for cl in cluster_unique:
                in_cluster = cluster_vals == cl
                tp = int((in_cluster & positives).sum())
                pred_pos = int(in_cluster.sum())
                if pred_pos == 0:
                    continue
                precision = tp / pred_pos
                recall = tp / n_pos if n_pos > 0 else 0.0
                f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
                cluster_mean = float(np.mean(scores[in_cluster]))
                row = {
                    "label": lb,
                    "resolution": res,
                    "cluster": cl,
                    "n_cluster_cells": pred_pos,
                    "n_label_positives": n_pos,
                    "tp": tp,
                    "precision": round(precision, 4),
                    "recall": round(recall, 4),
                    "f1": round(f1, 4),
                    "cluster_mean_combined": round(cluster_mean, 4),
                }
                long_rows.append(row)
                if best is None or f1 > best["f1"]:
                    best = row
            if best is not None:
                best_per_res[res] = best

        if not best_per_res:
            continue
        best_overall = max(best_per_res.values(), key=lambda r: r["f1"])
        f1_max = best_overall["f1"]
        recommended = [
            r for r, info in best_per_res.items()
            if info["f1"] >= f1_max - args.recommend_eps
        ]
        meta_row = (label_meta[label_meta["label"] == lb]
                    .iloc[0] if (label_meta["label"] == lb).any() else None)
        summary_rows.append({
            "label": lb,
            "lineage": meta_row["lineage"] if meta_row is not None else "",
            "v1_source_resolution": meta_row["v1_source_resolution"] if meta_row is not None else "",
            "v1_mode": meta_row["v1_mode"] if meta_row is not None else "",
            "n_label_positives": n_pos,
            "best_resolution": best_overall["resolution"],
            "best_cluster": best_overall["cluster"],
            "best_f1": f1_max,
            "best_precision": best_overall["precision"],
            "best_recall": best_overall["recall"],
            "best_cluster_mean_combined": best_overall["cluster_mean_combined"],
            "recommended_resolutions": ";".join(recommended),
        })

    long_df = pd.DataFrame(long_rows)
    long_path = args.out_dir / f"{args.compartment}_resolution_concordance.csv"
    long_df.to_csv(long_path, index=False)
    log.info("Wrote %s (%d rows)", long_path, len(long_df))

    summary_df = pd.DataFrame(summary_rows).sort_values("best_f1", ascending=False)
    summary_path = args.out_dir / f"{args.compartment}_resolution_recommendations.csv"
    summary_df.to_csv(summary_path, index=False)
    log.info("Wrote %s (%d labels)", summary_path, len(summary_df))

    log.info("Top 5 labels by best F1:")
    for _, r in summary_df.head(5).iterrows():
        log.info("  %s: best_res=%s f1=%.3f (recommended=%s; v1_source_res=%s)",
                 r["label"], r["best_resolution"], r["best_f1"],
                 r["recommended_resolutions"], r["v1_source_resolution"])


if __name__ == "__main__":
    main()
