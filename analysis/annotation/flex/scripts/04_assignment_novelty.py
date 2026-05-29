"""04_assignment_novelty.py — best-match L2S assignment + novelty flagging per cluster.

Run per (compartment × resolution). Reads the multi-resolution UCell long CSV
and filters to the requested resolution; reads limma top-markers for the same
resolution; combines canonical + identity scores via weighted sum; assigns
best-match V1 L2 label as L2S and flags novelty heuristics.

L2S = Track C's FLEX-native L2 annotation set, sourced from V1's L2 vocabulary
(annotation_v2_{epi,imm,str}.yaml) but tracked separately because the labels
attach to FLEX clusters via marker scoring (vs V1's transfer).

Inputs:
    --ucell-csv      outputs/ucell/{Compartment}_ucell_per_cluster.csv (long, multi-res)
    --limma-top-csv  outputs/limma/limma_top_markers_{Compartment}_leiden_{res}.csv
    --yaml           V1 yaml (canonical_markers per label, for marker_mismatch heuristic)
    --clusters-csv   outputs/clusters/{Compartment}_clusters.csv (multi-res)
    --compartment    Epithelial / Immune / Stromal
    --resolution     resolution suffix (e.g. 1.0)
    --weight-identity, --weight-canonical    UCell weighting (default 2:1)
    --novelty-low-confidence-max-score   default 0.30
    --novelty-tight-margin-max           default 0.05
    --novelty-marker-overlap-max         default 2
    --out-dir

Output: {Compartment}_leiden_{res}_assignments.csv with columns:
    resolution, cluster, n_cells, compartment,
    L2S_assignment, L2S_runner_up,
    ucell_combined, ucell_canonical, ucell_identity,
    score_margin, is_low_confidence, is_tight_margin, is_marker_mismatch,
    marker_overlap_with_canonical, is_novel, suggested_novel_label,
    limma_top_markers
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")

import pandas as pd
import yaml

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("assignment")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ucell-csv", required=True, type=Path)
    p.add_argument("--limma-top-csv", required=True, type=Path)
    p.add_argument("--yaml", required=True, type=Path)
    p.add_argument("--clusters-csv", required=True, type=Path)
    p.add_argument("--compartment", required=True, choices=["Epithelial", "Immune", "Stromal"])
    p.add_argument("--resolution", required=True)
    p.add_argument("--weight-identity", type=float, default=2.0)
    p.add_argument("--weight-canonical", type=float, default=1.0)
    p.add_argument("--novelty-low-confidence-max-score", type=float, default=0.30)
    p.add_argument("--novelty-tight-margin-max", type=float, default=0.05)
    p.add_argument("--novelty-marker-overlap-max", type=int, default=2)
    p.add_argument("--out-dir", required=True, type=Path)
    return p.parse_args()


def load_yaml_canonical_markers(yaml_path: Path) -> dict[str, list[str]]:
    with open(yaml_path) as fh:
        yml = yaml.safe_load(fh)
    out = {}
    for label, spec in yml.get("labels", {}).items():
        if spec.get("is_artifact"):
            continue
        out[label] = list(spec.get("canonical_markers", []) or [])
    return out


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    log.info("Loading UCell per-cluster: %s", args.ucell_csv)
    ucell = pd.read_csv(args.ucell_csv)
    required = {"resolution", "cluster", "label", "score_type", "median_score", "n_cells"}
    if not required.issubset(ucell.columns):
        sys.exit(f"ERROR: ucell CSV missing columns {required - set(ucell.columns)}")
    ucell = ucell[ucell["resolution"].astype(str) == args.resolution].copy()
    if ucell.empty:
        sys.exit(f"ERROR: no UCell rows for resolution={args.resolution}")
    log.info("  %d rows for resolution=%s, %d clusters, %d labels",
             len(ucell), args.resolution,
             ucell["cluster"].nunique(), ucell["label"].nunique())

    pivot = ucell.pivot_table(
        index=["cluster", "label", "n_cells"],
        columns="score_type",
        values="median_score",
        aggfunc="first",
    ).reset_index()
    pivot.columns.name = None
    for col in ("canonical", "identity"):
        if col not in pivot.columns:
            pivot[col] = 0.0
    pivot[["canonical", "identity"]] = pivot[["canonical", "identity"]].fillna(0.0)
    w_total = args.weight_identity + args.weight_canonical
    pivot["combined"] = (
        pivot["identity"] * args.weight_identity
        + pivot["canonical"] * args.weight_canonical
    ) / w_total

    log.info("Loading limma top markers: %s", args.limma_top_csv)
    limma = pd.read_csv(args.limma_top_csv)
    if "cluster" not in limma.columns or "gene" not in limma.columns:
        sys.exit(f"ERROR: limma CSV missing cluster/gene columns; have {list(limma.columns)}")
    limma["cluster"] = limma["cluster"].astype(str)
    cluster_top10 = (limma.groupby("cluster", group_keys=False)
                          .head(10)
                          .groupby("cluster")["gene"].apply(list).to_dict())
    cluster_top_str = (limma.groupby("cluster")["gene"]
                       .apply(lambda s: ",".join(s.head(15).astype(str)))
                       .to_dict())

    log.info("Loading clusters: %s", args.clusters_csv)
    clusters_df = pd.read_csv(args.clusters_csv)
    leiden_col = f"leiden_{args.resolution}"
    if leiden_col not in clusters_df.columns:
        sys.exit(f"ERROR: clusters CSV missing column {leiden_col}")
    n_per_cluster = clusters_df[leiden_col].astype(str).value_counts().to_dict()

    log.info("Loading yaml canonical markers: %s", args.yaml)
    label_canonical = load_yaml_canonical_markers(args.yaml)

    pivot["cluster"] = pivot["cluster"].astype(str)
    rows = []
    for cl in sorted(pivot["cluster"].unique(), key=lambda s: (len(s), s)):
        sub = (pivot[pivot["cluster"] == cl]
               .sort_values(by="combined", ascending=False)
               .reset_index(drop=True))
        if sub.empty:
            log.warning("  cluster %s: no UCell scores; skipping", cl)
            continue
        best = sub.iloc[0]
        runner = sub.iloc[1] if len(sub) >= 2 else None
        margin = float(best["combined"] - (runner["combined"] if runner is not None else 0.0))

        is_low = bool(best["combined"] < args.novelty_low_confidence_max_score)
        is_tight = bool(margin < args.novelty_tight_margin_max)

        top10_genes = cluster_top10.get(cl, [])
        best_canonical = label_canonical.get(str(best["label"]), [])
        overlap = len(set(top10_genes) & set(best_canonical))
        is_mismatch = bool(overlap <= args.novelty_marker_overlap_max)

        is_novel = bool(is_low or is_tight or is_mismatch)
        tags = []
        if is_low: tags.append("low_confidence")
        if is_tight: tags.append("tight_margin")
        if is_mismatch: tags.append("marker_mismatch")
        suggested = ";".join(tags)

        rows.append({
            "resolution": args.resolution,
            "cluster": cl,
            "n_cells": int(n_per_cluster.get(cl, 0)),
            "compartment": args.compartment,
            "L2S_assignment": str(best["label"]),
            "L2S_runner_up": str(runner["label"]) if runner is not None else "",
            "ucell_combined": round(float(best["combined"]), 4),
            "ucell_canonical": round(float(best["canonical"]), 4),
            "ucell_identity": round(float(best["identity"]), 4),
            "score_margin": round(margin, 4),
            "is_low_confidence": is_low,
            "is_tight_margin": is_tight,
            "is_marker_mismatch": is_mismatch,
            "marker_overlap_with_canonical": overlap,
            "is_novel": is_novel,
            "suggested_novel_label": suggested,
            "limma_top_markers": cluster_top_str.get(cl, ""),
        })

    out_df = pd.DataFrame(rows)
    out_path = args.out_dir / f"{args.compartment}_leiden_{args.resolution}_assignments.csv"
    out_df.to_csv(out_path, index=False)
    log.info("Wrote %s (%d clusters)", out_path, len(out_df))

    log.info("Summary: %d/%d clusters flagged novel (low_conf=%d, tight_margin=%d, marker_mismatch=%d)",
             int(out_df["is_novel"].sum()), len(out_df),
             int(out_df["is_low_confidence"].sum()),
             int(out_df["is_tight_margin"].sum()),
             int(out_df["is_marker_mismatch"].sum()))


if __name__ == "__main__":
    main()
