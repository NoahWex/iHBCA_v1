#!/usr/bin/env python3
"""
Per-sweep SVD diagnostic on a scIB metric matrix.

Catches failure modes that single-number composite scores hide:

  1. PC1 dominance — if PC1 > 90% of variance, the metric set is uninformative
     (the sweep is effectively 1D; methods are too similar to discriminate).
  2. Single-metric collapse — if any individual metric loads > 0.9 on PC1,
     that metric IS the score; the benchmark has degenerated to one column.
  3. Method-family separability — project configs onto PC1/PC2 and report
     cluster centers per method family. If families overlap on PC1, the
     dominant axis does not separate the methods.
  4. Effective dimensionality — participation ratio = 1 / sum(var_i^2) on the
     normalized PC variance vector. PR ~ 1 indicates 1D collapse, PR ~
     n_metrics indicates uniform noise, intermediate values indicate
     structured multi-dim variance.

Inputs:
    --csv <scib_consolidated.csv>
        Long-format CSV with at least one row per (compartment, config) and
        one column per metric. Compartment column is optional — if absent,
        all configs are treated as one sweep.
    --compartment-col name (default: 'compartment'; pass empty string to skip)
    --config-col name (default: 'config')
    --metric-cols 'm1,m2,m3,...' (default: NMI,ARI,ASW_label,graph_conn,
                                  ASW_batch,kBET,isolated_labels)
    --output path (optional) — write loadings + projections to JSON

Outputs:
    Stdout: per-compartment SVD report (variance, loadings, projections,
    warnings). Optional JSON file with full numerical results.
"""

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd


PR_LOW_WARN = 1.5      # PR < 1.5 means very near 1D collapse
PC1_HIGH_WARN = 0.90   # PC1 > 90% means single dominant axis
LOAD_HIGH_WARN = 0.90  # any |loading| > 0.9 on PC1 means single metric IS the score


def svd_one_sweep(metric_df: pd.DataFrame, name: str = ""):
    """Run SVD on a (n_configs, n_metrics) frame. Returns dict with results.

    Standardization note: each metric is z-scored before SVD. This puts every
    metric on equal footing regardless of absolute dynamic range — a metric
    with spread 0.001 contributes the same unit variance as one with spread
    0.5 after z-scoring. This is the standard scIB-benchmark convention but
    means the loadings reflect within-metric ranking structure, not absolute
    quality differences. A metric with very small absolute spread can carry
    a substantial PC1 loading even when its real contribution to integration
    quality is negligible. Cross-check loadings against raw metric ranges
    when interpreting the dominant axis.
    """
    if len(metric_df) < 3:
        return {"name": name, "error": f"too few configs ({len(metric_df)}) for SVD"}

    # Standardize per metric (z-score)
    means = metric_df.mean()
    stds = metric_df.std()
    if (stds == 0).any():
        zero_metrics = stds[stds == 0].index.tolist()
        return {"name": name, "error": f"zero variance in {zero_metrics}"}
    X = (metric_df - means) / stds

    U, s, Vt = np.linalg.svd(X.values, full_matrices=False)
    var = (s ** 2) / (s ** 2).sum()
    pr = float(1.0 / (var ** 2).sum())

    n_pcs = min(len(s), 5)
    # loadings: dict[PC label -> dict[metric -> loading]]
    loadings = {f"PC{k+1}": {m: round(float(Vt[k, j]), 4)
                              for j, m in enumerate(metric_df.columns)}
                for k in range(n_pcs)}
    # projection: dict[PC label -> dict[config_name -> score]]
    proj_arr = U[:, :n_pcs] * s[:n_pcs]
    projection = {f"PC{k+1}": {c: round(float(proj_arr[i, k]), 4)
                                for i, c in enumerate(metric_df.index)}
                  for k in range(n_pcs)}

    return {
        "name": name,
        "n_configs": int(len(metric_df)),
        "n_metrics": int(metric_df.shape[1]),
        "n_pcs": int(n_pcs),
        "variance_explained": var.tolist(),
        "participation_ratio": pr,
        "loadings": loadings,
        "projection": projection,
        "configs": list(metric_df.index),
    }


def project_by_family(projection: pd.DataFrame, configs: pd.Series):
    """Group config projections by method family for separability check."""
    def fam(c):
        if c == "pca_n50":
            return "pca"
        if c == "harmony_n50":
            return "harmony"
        if "austin" in c:
            return "scvi_austin"
        # Order matters: scanvi check must come BEFORE scvi check, since
        # "scanvi_n50".startswith("scvi") is False (good) but the explicit
        # ordering keeps the family categorization unambiguous.
        if c.startswith("scanvi"):
            return "scanvi"
        if c.startswith("scvi"):
            return "scvi"
        return "other"

    fam_series = configs.apply(fam)
    rows = []
    for family in sorted(fam_series.unique()):
        mask = fam_series == family
        n = int(mask.sum())
        rows.append({
            "family": family,
            "n": n,
            "PC1_mean": float(projection.loc[mask, "PC1"].mean()),
            "PC1_min": float(projection.loc[mask, "PC1"].min()),
            "PC1_max": float(projection.loc[mask, "PC1"].max()),
            "PC2_mean": float(projection.loc[mask, "PC2"].mean()),
        })
    return pd.DataFrame(rows)


def warn_check(result: dict, metric_cols: list) -> list[str]:
    """Apply the failure-mode checks. Returns list of human-readable warnings."""
    warnings = []
    if "error" in result:
        return [f"FAILED: {result['error']}"]

    var = result["variance_explained"]
    if var[0] > PC1_HIGH_WARN:
        warnings.append(
            f"WARN [PC1 dominance]: PC1 = {var[0]:.3f} > {PC1_HIGH_WARN}. "
            f"The metric set has collapsed to ~1D. Methods in this sweep are too "
            f"similar to differentiate, OR the metric set is too redundant. "
            f"Add a divergent baseline (e.g. unintegrated PCA) or drop redundant metrics."
        )

    pr = result["participation_ratio"]
    if pr < PR_LOW_WARN:
        warnings.append(
            f"WARN [effective dim]: participation ratio = {pr:.2f} < {PR_LOW_WARN}. "
            f"The metric matrix is essentially low-dimensional. Composite scoring is "
            f"a roundabout way of reporting one number."
        )

    # loadings is dict[PC label -> dict[metric -> loading]]
    pc1_loadings = result["loadings"].get("PC1", {})
    over = [(m, abs(v)) for m, v in pc1_loadings.items() if abs(v) > LOAD_HIGH_WARN]
    if over:
        offending = [m for m, _ in over]
        warnings.append(
            f"WARN [single-metric PC1]: {offending} carries |loading| > {LOAD_HIGH_WARN} "
            f"on PC1. The composite score is effectively that single metric. "
            f"Either drop the redundant metrics or rebalance."
        )

    return warnings


def render_report(name, result, metric_cols, configs):
    """Pretty-print one sweep's SVD report."""
    print()
    print("=" * 80)
    print(f"SVD diagnostic: {name}")
    print("=" * 80)
    if "error" in result:
        print(f"  ERROR: {result['error']}")
        return

    print(f"  n_configs: {result['n_configs']}")
    print(f"  n_metrics: {result['n_metrics']}")
    print(f"  participation ratio (effective dim): {result['participation_ratio']:.2f}")

    var = result["variance_explained"]
    print(f"  variance explained: " + ", ".join(
        f"PC{i+1}={v:.3f}" for i, v in enumerate(var[:5])))

    # Loadings: dict[PC label -> dict[metric -> loading]]
    n_pcs_show = min(result.get("n_pcs", 3), 3)
    loadings_df = pd.DataFrame({
        f"PC{k+1}": [result["loadings"][f"PC{k+1}"].get(m, float("nan")) for m in metric_cols]
        for k in range(n_pcs_show)
    }, index=metric_cols)
    print()
    print("  Loadings (rows=metrics, cols=PCs):")
    print("  " + loadings_df.round(3).to_string().replace("\n", "\n  "))

    if configs is not None and "PC1" in result["projection"]:
        # projection is dict[PC label -> dict[config_name -> score]]
        config_list = list(configs.values) if hasattr(configs, "values") else list(configs)
        proj = pd.DataFrame({
            "PC1": [result["projection"]["PC1"].get(c, 0.0) for c in config_list],
            "PC2": [result["projection"]["PC2"].get(c, 0.0) for c in config_list]
                   if "PC2" in result["projection"]
                   else [0.0] * len(config_list),
        })
        family_df = project_by_family(proj, pd.Series(config_list))
        print()
        print("  Method-family separation on PC1/PC2:")
        print("  " + family_df.round(3).to_string(index=False).replace("\n", "\n  "))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--csv", required=True)
    p.add_argument("--compartment-col", default="compartment",
                   help="Pass empty string to disable per-compartment splitting")
    p.add_argument("--config-col", default="config")
    p.add_argument("--metric-cols",
                   default="NMI,ARI,ASW_label,graph_conn,ASW_batch,kBET,isolated_labels")
    p.add_argument("--output", default=None,
                   help="Optional JSON output for full numerical results")
    args = p.parse_args()

    df = pd.read_csv(args.csv)
    metric_cols = [m.strip() for m in args.metric_cols.split(",") if m.strip()]

    missing = [m for m in metric_cols if m not in df.columns]
    if missing:
        sys.exit(f"FAIL: metric columns missing from CSV: {missing}")
    if args.config_col not in df.columns:
        sys.exit(f"FAIL: config column '{args.config_col}' missing from CSV")

    print(f"=== scIB metric SVD diagnostic ===")
    print(f"  CSV: {args.csv}")
    print(f"  Metrics: {metric_cols}")

    all_warnings = []
    full_results = {}

    if args.compartment_col and args.compartment_col in df.columns:
        compartments = sorted(df[args.compartment_col].unique())
        print(f"  Compartments: {compartments}")
        for comp in compartments:
            sub = df[df[args.compartment_col] == comp].dropna(subset=metric_cols).copy()
            metric_df = sub[metric_cols]
            metric_df.index = sub[args.config_col].values
            result = svd_one_sweep(metric_df, name=comp)
            full_results[comp] = result
            render_report(comp, result, metric_cols, sub[args.config_col])
            warnings = warn_check(result, metric_cols)
            for w in warnings:
                all_warnings.append(f"[{comp}] {w}")
    else:
        print(f"  Compartment col disabled; treating as single sweep")
        sub = df.dropna(subset=metric_cols).copy()
        metric_df = sub[metric_cols]
        metric_df.index = sub[args.config_col].values
        result = svd_one_sweep(metric_df, name="all")
        full_results["all"] = result
        render_report("all", result, metric_cols, sub[args.config_col])
        all_warnings.extend(warn_check(result, metric_cols))

    print()
    print("=" * 80)
    print("Warnings summary")
    print("=" * 80)
    if all_warnings:
        for w in all_warnings:
            print(f"  {w}")
    else:
        print("  No warnings raised. Metric matrix has structured multi-dim variance.")

    if args.output:
        with open(args.output, "w") as f:
            json.dump(full_results, f, indent=2, default=float)
        print()
        print(f"Wrote full numerical results: {args.output}")


if __name__ == "__main__":
    main()
