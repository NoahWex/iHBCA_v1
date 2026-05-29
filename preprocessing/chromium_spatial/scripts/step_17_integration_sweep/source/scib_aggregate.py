#!/usr/bin/env python3
"""
step_17_integration_sweep / scib_aggregate.py

Aggregate per-metric JSONs across configs into a wide comparison CSV,
apply the PINNED scoring formula from sweep_configs.yaml, and report
the winner (config with the highest `overall` score).

Reads:
  {scib_dir}/{config}/{metric}.json  for each config and each metric

Writes:
  {scib_dir}/scib_comparison.csv     wide CSV, one row per config, sorted
                                     by `overall` descending

Scoring formula (locked; defined in sweep_configs.yaml:scoring):
  bio_mean   = mean(ASW_label, NMI, ARI)
  batch_mean = mean(ASW_batch, graph_conn)
  overall    = 0.6 * bio_mean + 0.4 * batch_mean

Sidecar columns (NOT in scoring): isolated_labels, PCR, kBET.

"""

import argparse
import glob
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


def load_config_metrics(config_dir):
    """Read all metric JSONs for one config dir.

    """
    json_files = sorted(glob.glob(os.path.join(config_dir, "*.json")))
    metrics = {}
    for jf in json_files:
        with open(jf) as f:
            data = json.load(f)
        metrics[data["metric"]] = data["value"]
        metrics[f"{data['metric']}_n_cells"] = data["n_cells"]
    return metrics


def aggregate(scib_dir, cfg, expected_configs=None, allow_nan_bio=False):
    """Walk scib_dir/{config}/ and build a scored dataframe.


    allow_nan_bio: when True, NaN bio metrics are tolerated (e.g. no labels
      CSV supplied). bio_mean is computed as nanmean of available values; if
      all bio metrics are NaN, overall = batch_mean alone (batch-only scoring).
      NaN batch metrics still cause a hard failure — they indicate a real
      computation problem, not just missing labels.
    """
    bio_keys = cfg["scoring"]["bio_metrics"]
    batch_keys = cfg["scoring"]["batch_metrics"]
    sidecar_keys = cfg["metrics"].get("sidecar", [])
    bio_w = float(cfg["scoring"]["bio_weight"])
    batch_w = float(cfg["scoring"]["batch_weight"])

    config_dirs = sorted(d for d in glob.glob(os.path.join(scib_dir, "*"))
                         if os.path.isdir(d))
    if not config_dirs:
        _die(f"no config subdirs in {scib_dir}")

    rows = []
    for cd in config_dirs:
        cfg_label = os.path.basename(cd)
        if expected_configs and cfg_label not in expected_configs:
            continue
        m = load_config_metrics(cd)

        # Batch metrics must always be finite — NaN here is a real failure.
        for k in batch_keys:
            if k not in m:
                _die(f"{cfg_label}/{k}.json missing")
            if m[k] is None or not np.isfinite(m[k]):
                _die(f"{cfg_label}/{k}.json value is null/non-finite ({m[k]})")

        # Bio metrics: fail loud unless --allow-nan-bio
        for k in bio_keys:
            if k not in m:
                _die(f"{cfg_label}/{k}.json missing")
            if m[k] is not None and not np.isfinite(m[k]):
                _die(f"{cfg_label}/{k}.json value is non-finite ({m[k]})")
            if m[k] is None and not allow_nan_bio:
                _die(f"{cfg_label}/{k}.json value is null/non-finite ({m[k]})",
                     hint="If no labels CSV was provided, bio metrics will be "
                          "null. Re-run with --allow-nan-bio to score on "
                          "batch metrics only.")

        batch_vals = [m[k] for k in batch_keys]
        batch_mean = float(np.mean(batch_vals))

        bio_vals_finite = [m[k] for k in bio_keys if m.get(k) is not None]
        if bio_vals_finite:
            bio_mean = float(np.mean(bio_vals_finite))
            overall = bio_w * bio_mean + batch_w * batch_mean
        else:
            # No bio labels available — rank by batch only.
            bio_mean = float("nan")
            overall = batch_mean

        row = {"config": cfg_label}
        for k in bio_keys + batch_keys + sidecar_keys:
            row[k] = m.get(k)
        row["bio_mean"] = round(bio_mean, 6) if np.isfinite(bio_mean) else None
        row["batch_mean"] = round(batch_mean, 6)
        row["overall"] = round(overall, 6)
        row["n_cells"] = m.get(f"{bio_keys[0]}_n_cells")
        rows.append(row)

    if not rows:
        _die(f"no configs aggregated from {scib_dir}")

    df = pd.DataFrame(rows).sort_values("overall", ascending=False).reset_index(drop=True)
    return df


def _write_winner_record(df, scib_dir, manual_override=None):
    """Record winner + runner-up + margin to a JSON audit file."""
    winner_row = df.iloc[0]
    runner_row = df.iloc[1] if len(df) > 1 else None
    margin = None
    if runner_row is not None:
        margin = float(winner_row["overall"]) - float(runner_row["overall"])

    def _maybe_float(v):
        return float(v) if v is not None and not (isinstance(v, float) and np.isnan(v)) else None

    winner = {
        "selected_config": manual_override or str(winner_row["config"]),
        "selection_method": "manual_override" if manual_override else "top_overall",
        "winner_overall": float(winner_row["overall"]),
        "winner_bio_mean": _maybe_float(winner_row["bio_mean"]),
        "winner_batch_mean": float(winner_row["batch_mean"]),
        "runner_up_config": str(runner_row["config"]) if runner_row is not None else None,
        "runner_up_overall": float(runner_row["overall"]) if runner_row is not None else None,
        "margin": margin,
        "n_configs_evaluated": int(len(df)),
    }

    # If manual override, look up the overridden row
    if manual_override:
        sel = df[df["config"] == manual_override]
        if sel.empty:
            _die(f"--manual-override {manual_override!r} not in aggregated configs")
        sel_row = sel.iloc[0]
        winner["selected_overall"] = float(sel_row["overall"])
        winner["selected_rank"] = int(sel.index[0] + 1)

    out_path = os.path.join(scib_dir, "winner_selection.json")
    with open(out_path, "w") as f:
        json.dump(winner, f, indent=2)
    print(f"Wrote: {out_path}")
    return winner


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--scib-dir", required=True,
                   help="Directory containing per-config subdirs of metric JSONs. "
                        "Typically CFG_CANONICAL_SWEEP_SCIB or a target subdir.")
    p.add_argument("--output-name", default="scib_comparison.csv",
                   help="CSV filename written into --scib-dir")
    p.add_argument("--configs", nargs="*", default=None,
                   help="Optional whitelist (default: all subdirs)")
    p.add_argument("--manual-override", default=None,
                   help="Override automatic winner selection with this config "
                        "label. Writes winner_selection.json with "
                        "selection_method=manual_override and preserves the "
                        "automatic rank in selected_rank.")
    p.add_argument("--allow-nan-bio", action="store_true",
                   help="Allow null bio metrics (e.g. no labels CSV supplied). "
                        "Falls back to batch-only scoring when all bio metrics "
                        "are NaN. NaN batch metrics still cause a hard failure.")
    p.add_argument("--config-path", default=DEFAULT_CONFIG_PATH,
                   help="Path to sweep_configs.yaml")
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()

    print(f"=== step_17 scIB aggregation ===")
    print(f"  scib_dir:    {args.scib_dir}")
    print(f"  config:      {args.config_path}")

    if args.dry_run:
        print("=== DRY RUN MODE ===")
        if not os.path.isdir(args.scib_dir):
            _die(f"--scib-dir not a directory: {args.scib_dir}")
        cfg = _load_config(args.config_path)
        print(f"  bio metrics:   {cfg['scoring']['bio_metrics']}")
        print(f"  batch metrics: {cfg['scoring']['batch_metrics']}")
        print(f"  weights:       bio={cfg['scoring']['bio_weight']} "
              f"batch={cfg['scoring']['batch_weight']}")
        print("VALIDATION PASSED")
        return

    cfg = _load_config(args.config_path)
    if args.allow_nan_bio:
        print("  NOTE: --allow-nan-bio set — bio metrics may be NaN; "
              "falling back to batch-only scoring where needed.")
    df = aggregate(args.scib_dir, cfg, expected_configs=args.configs,
                   allow_nan_bio=args.allow_nan_bio)

    out_path = os.path.join(args.scib_dir, args.output_name)
    df.to_csv(out_path, index=False)

    print(f"  Configs:  {len(df)}")
    print(f"  Cells:    {df['n_cells'].iloc[0]}")
    print()
    print(df.to_string(index=False))
    print()
    print(f"Wrote: {out_path}")

    winner = _write_winner_record(df, args.scib_dir,
                                  manual_override=args.manual_override)
    print()
    print("=== Winner selection ===")
    print(f"  selected: {winner['selected_config']} "
          f"(method={winner['selection_method']})")
    if winner.get("runner_up_config"):
        print(f"  runner-up: {winner['runner_up_config']} "
              f"(margin={winner['margin']:.6f})")


if __name__ == "__main__":
    main()
