#!/usr/bin/env python3
"""
step_17_integration_sweep / svd_diagnostic.py

SVD on the wide scIB metric matrix (configs x metrics) to catch failure
modes that a single composite score hides:

1. PC1 dominance: if PC1 > PC1_HIGH_WARN (default 0.90), the metric set
   has collapsed to ~1D. Methods are too similar to differentiate, or the
   metric set itself is too redundant to be informative.
2. Single-metric collapse: if any individual metric loads > LOAD_HIGH_WARN
   on PC1, that metric IS the score. The composite is a rename.
3. Effective dimensionality: participation ratio = 1 / sum(var_i^2) on the
   normalized PC variance vector. PR < PR_LOW_WARN means near-1D collapse.
4. Method-family separation: project configs onto PC1/PC2 and report
   means per method family. If scVI and PCA overlap on PC1, the dominant
   axis does not separate the methods.

Default behaviour with --fail-on-pc1 (locked decision in sweep_configs.yaml
svd_thresholds.fail_on_pc1=true): exit non-zero if PC1 exceeds the threshold,
forcing the coordinator to investigate before promoting the winner.

Inputs:
  --csv                 scib_comparison.csv produced by scib_aggregate.py
  --config-path         sweep_configs.yaml (for thresholds)
  --compartment-col     column name identifying compartment runs (optional)
  --config-col          column name identifying configs (default 'config')
  --metric-cols         comma-separated metric columns (default pulls all
                        scored + sidecar metrics from sweep_configs.yaml)
  --output              optional JSON file with full numerical results

"""

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import yaml


DEFAULT_CONFIG_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "configs", "sweep_configs.yaml",
)


def _die(msg, hint=None):
    sys.stderr.write(f"FATAL: {msg}\n")
    if hint:
        sys.stderr.write(f"  hint: {hint}\n")
    sys.exit(2)


def _load_config(path):
    if not os.path.exists(path):
        _die(f"sweep_configs.yaml not found: {path}")
    with open(path) as fh:
        return yaml.safe_load(fh)


def svd_one_sweep(metric_df: pd.DataFrame, name: str = ""):
    """Run SVD on a (n_configs, n_metrics) matrix. Returns result dict.

    Standardization: each metric is z-scored before SVD so metrics with
    wildly different absolute spreads are put on equal footing. This is
    the standard scIB-benchmark convention; cross-check loadings against
    raw metric ranges if the dominant axis is surprising.

    """
    if len(metric_df) < 3:
        return {"name": name, "error": f"too few configs ({len(metric_df)}) for SVD"}

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
    loadings = {
        f"PC{k+1}": {m: round(float(Vt[k, j]), 4)
                     for j, m in enumerate(metric_df.columns)}
        for k in range(n_pcs)
    }
    proj_arr = U[:, :n_pcs] * s[:n_pcs]
    projection = {
        f"PC{k+1}": {c: round(float(proj_arr[i, k]), 4)
                     for i, c in enumerate(metric_df.index)}
        for k in range(n_pcs)
    }

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
    """Group config projections by method family for separability check.

    """
    def fam(c):
        if c == "pca_n50":
            return "pca"
        if c == "harmony_n50":
            return "harmony"
        if "austin" in c:
            return "scvi_austin"
        # Order matters: scanvi before scvi (scanvi doesn't startswith scvi,
        # but explicit is cheap insurance against future config label drift).
        if c.startswith("scanvi"):
            return "scanvi"
        if c.startswith("scvi"):
            return "scvi"
        return "other"

    fam_series = configs.apply(fam)
    rows = []
    for family in sorted(fam_series.unique()):
        mask = fam_series == family
        rows.append({
            "family": family,
            "n": int(mask.sum()),
            "PC1_mean": float(projection.loc[mask, "PC1"].mean()),
            "PC1_min": float(projection.loc[mask, "PC1"].min()),
            "PC1_max": float(projection.loc[mask, "PC1"].max()),
            "PC2_mean": float(projection.loc[mask, "PC2"].mean()),
        })
    return pd.DataFrame(rows)


def warn_check(result, thresholds):
    """Apply failure-mode checks. Returns list of warning strings.

    """
    warnings = []
    if "error" in result:
        return [f"FAILED: {result['error']}"]

    pc1_warn = float(thresholds.get("pc1_high_warn", 0.90))
    pr_warn = float(thresholds.get("pr_low_warn", 1.5))
    load_warn = float(thresholds.get("load_high_warn", 0.90))

    var = result["variance_explained"]
    if var[0] > pc1_warn:
        warnings.append(
            f"WARN [PC1 dominance]: PC1 = {var[0]:.3f} > {pc1_warn}. "
            f"The metric set has collapsed to ~1D. Methods in this sweep are "
            f"too similar to differentiate, OR the metric set is too redundant. "
            f"Add a divergent baseline (e.g. unintegrated PCA) or drop redundant metrics."
        )

    pr = result["participation_ratio"]
    if pr < pr_warn:
        warnings.append(
            f"WARN [effective dim]: participation ratio = {pr:.2f} < {pr_warn}. "
            f"The metric matrix is essentially low-dimensional. Composite scoring "
            f"is a roundabout way of reporting one number."
        )

    pc1_loadings = result["loadings"].get("PC1", {})
    over = [(m, abs(v)) for m, v in pc1_loadings.items() if abs(v) > load_warn]
    if over:
        offending = [m for m, _ in over]
        warnings.append(
            f"WARN [single-metric PC1]: {offending} carries |loading| > {load_warn} "
            f"on PC1. The composite score is effectively that single metric."
        )

    return warnings


def render_report(name, result, metric_cols, configs):
    """Pretty-print one sweep's SVD report.

    """
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

    n_pcs_show = min(result.get("n_pcs", 3), 3)
    loadings_df = pd.DataFrame({
        f"PC{k+1}": [result["loadings"][f"PC{k+1}"].get(m, float("nan"))
                     for m in metric_cols]
        for k in range(n_pcs_show)
    }, index=metric_cols)
    print()
    print("  Loadings (rows=metrics, cols=PCs):")
    print("  " + loadings_df.round(3).to_string().replace("\n", "\n  "))

    if configs is not None and "PC1" in result["projection"]:
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


def _default_metric_cols(cfg):
    """Resolve metric column list from sweep_configs.yaml."""
    scored = cfg.get("metrics", {}).get("scored", [])
    sidecar = cfg.get("metrics", {}).get("sidecar", [])
    return list(scored) + list(sidecar)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--csv", required=True,
                   help="scib_comparison.csv from scib_aggregate.py")
    p.add_argument("--compartment-col", default="",
                   help="Column name for compartment-resolved sweeps. "
                        "Empty string = treat all rows as a single sweep.")
    p.add_argument("--config-col", default="config")
    p.add_argument("--metric-cols", default=None,
                   help="Comma-separated metric columns. Default = all scored "
                        "+ sidecar metrics from sweep_configs.yaml.")
    p.add_argument("--config-path", default=DEFAULT_CONFIG_PATH)
    p.add_argument("--output", default=None,
                   help="Optional JSON output for full numerical results")
    p.add_argument("--fail-on-pc1", action="store_true",
                   help="Exit non-zero if PC1 > pc1_high_warn threshold. "
                        "Matches sweep_configs.yaml:svd_thresholds.fail_on_pc1.")
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()

    print("=" * 70)
    print("=== step_17 scIB metric SVD diagnostic ===")
    print(f"  CSV:    {args.csv}")
    print(f"  Config: {args.config_path}")
    print("=" * 70)

    if args.dry_run:
        print("=== DRY RUN MODE ===")
        if not os.path.exists(args.csv):
            _die(f"--csv not found: {args.csv}")
        print("VALIDATION PASSED")
        return

    if not os.path.exists(args.csv):
        _die(f"--csv not found: {args.csv}")

    cfg = _load_config(args.config_path)
    thresholds = cfg.get("svd_thresholds", {})
    fail_on_pc1 = args.fail_on_pc1 or bool(thresholds.get("fail_on_pc1", False))

    df = pd.read_csv(args.csv)
    if args.metric_cols:
        metric_cols = [m.strip() for m in args.metric_cols.split(",") if m.strip()]
    else:
        metric_cols = _default_metric_cols(cfg)

    # Only use metrics that are actually present in the CSV (sidecar PCR/kBET
    # may be missing if the SLURM array didn't request them).
    available = [m for m in metric_cols if m in df.columns]
    dropped = [m for m in metric_cols if m not in df.columns]
    if dropped:
        print(f"  NOTE: dropping metrics not in CSV: {dropped}")
    metric_cols = available

    # Drop metrics where ALL values are NaN (e.g. bio metrics when no labels
    # CSV was provided). Such columns carry zero information for SVD and would
    # cause dropna() to eliminate every row.
    all_nan = [m for m in metric_cols if df[m].isna().all()]
    if all_nan:
        print(f"  NOTE: dropping all-NaN metrics (no labels available): {all_nan}")
        metric_cols = [m for m in metric_cols if m not in all_nan]

    # Drop zero-variance metrics (e.g. graph_conn = 1.0 for all configs).
    # Z-scoring a constant column produces NaN, which causes SVD failure.
    zero_var = [m for m in metric_cols if m in df.columns and df[m].std() == 0]
    if zero_var:
        print(f"  NOTE: dropping zero-variance metrics (constant across configs): {zero_var}")
        metric_cols = [m for m in metric_cols if m not in zero_var]

    if len(metric_cols) < 2:
        _die(f"need >=2 metric columns for SVD, got {metric_cols}")
    if args.config_col not in df.columns:
        _die(f"config column {args.config_col!r} missing from CSV")

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
            for w in warn_check(result, thresholds):
                all_warnings.append(f"[{comp}] {w}")
    else:
        print(f"  Single-sweep mode (no compartment column)")
        sub = df.dropna(subset=metric_cols).copy()
        metric_df = sub[metric_cols]
        metric_df.index = sub[args.config_col].values
        result = svd_one_sweep(metric_df, name="all")
        full_results["all"] = result
        render_report("all", result, metric_cols, sub[args.config_col])
        all_warnings.extend(warn_check(result, thresholds))

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
        print(f"Wrote full results: {args.output}")

    # Fail-loud gate
    if fail_on_pc1:
        pc1_warn = float(thresholds.get("pc1_high_warn", 0.90))
        bad = []
        for name, result in full_results.items():
            if "error" in result:
                continue
            if result["variance_explained"][0] > pc1_warn:
                bad.append(f"{name} (PC1={result['variance_explained'][0]:.3f})")
        if bad:
            _die(f"PC1 dominance detected (threshold={pc1_warn}): {bad}",
                 hint="Inspect the warnings above. Consider: (1) adding divergent "
                      "baselines to the sweep; (2) dropping redundant metrics; "
                      "(3) overriding winner selection manually with "
                      "scib_aggregate.py --manual-override.")


if __name__ == "__main__":
    main()
