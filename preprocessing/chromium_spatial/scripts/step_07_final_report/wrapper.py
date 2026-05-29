#!/usr/bin/env python3
# ============================================================================
# Step 07: Final Report - Wrapper
# ============================================================================
# Purpose
#   Aggregate per-sample summaries from step 01-06 sidecars, flag
#   statistical outliers across the cohort, and render an HTML QC
#   report used as the coordinator's review gate before unlocking
#   downstream integration.
#
# Methodology
#   - Metric aggregation walks each step's sample_manifests/{sample}.yaml
#     and pulls the summary block, producing a wide per-sample DataFrame.
#   - Outlier flagging uses cohort-wide MAD on priority metrics
#     (cell_removal_rate, vf_retention_rate). mad_threshold = 2.0 is the
#     standard FLEX convention: samples beyond +-2 MAD on a priority
#     metric are flagged for manual review, not automatically excluded.
#   - The HTML report is rendered via a Jinja template that embeds the
#     generated PNGs; it is the primary human-readable artifact from L1.
#
# Inputs
#   - CFG_PREPROCESSING_ROOT (scans each step's sample_manifests dir)
#   - Module configs with step_07_final_report block
#
# Outputs
#   - ${CFG_PREPROCESSING_STEP_07}/reports/preprocessing_qc_report.html
#   - ${CFG_PREPROCESSING_STEP_07}/data/flagged_samples.csv
#   - ${CFG_PREPROCESSING_STEP_07}/data/sample_metrics.csv
#   - ${CFG_PREPROCESSING_STEP_07}/plots/*.png
# ============================================================================

from __future__ import annotations

import argparse
import os
import sys
from datetime import datetime
from pathlib import Path

import yaml

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR / "source"))


def _log(action, value=None, level="INFO"):
    ts = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    if value is None:
        print(f"[{ts}] {level} step_07 {action}", flush=True)
    else:
        print(f"[{ts}] {level} step_07 {action}={value}", flush=True)


def _fail_loud(msg, hint=None):
    _log("FAIL", msg, level="ERROR")
    if hint:
        print(f"  hint: {hint}", flush=True)
    sys.exit(2)


def _parse_args():
    p = argparse.ArgumentParser(description="Step 07: Final Report")
    p.add_argument("--preprocessing-root",
                   default=os.environ.get("CFG_PREPROCESSING_ROOT", ""),
                   help="Root of all preprocessing step outputs.")
    p.add_argument("--central-cell-status",
                   default=os.environ.get("CFG_PREPROCESSING_CENTRAL_CELL_STATUS", ""))
    p.add_argument("--raw-data-manifest",
                   default=os.environ.get("CFG_RAW_DATA_SAMPLE_MANIFEST", ""))
    p.add_argument("--output-root",
                   default=os.environ.get("CFG_PREPROCESSING_STEP_07", ""))
    p.add_argument("--module-configs",
                   default=os.environ.get(
                       "CFG_MODULE_CONFIGS",
                       os.path.join(os.environ.get("CFG_PROJECT_ROOT", ""),
                                    "config", "module_configs.yaml")))
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def _load_yaml(path):
    with open(path, "r") as fh:
        return yaml.safe_load(fh)


def _aggregate_sample_manifests(preprocessing_root: Path, expected_samples: list[str]):
    """Walk per-step sample_manifests/ dirs and merge into one record per sample.

    Fails loudly if any expected sample is missing from any L1 step's manifests
    (steps 01/02/03/05/06). Step 04 is excluded because it is superseded.
    """
    steps = [
        ("01_BaselineVFs", "step_01_baseline_vfs"),
        ("02_CellFiltering", "step_02_cell_filtering"),
        ("03_DoubletDetection", "step_03_doublet"),
        ("05_PostQCVFs", "step_05_post_qc_vfs"),
        ("06_FilteredPreview", "step_06_filtered_preview"),
    ]
    aggregated = {sid: {"sample_id": sid} for sid in expected_samples}

    for step_dir, step_key in steps:
        mfst_dir = preprocessing_root / step_dir / "sample_manifests"
        if not mfst_dir.is_dir():
            _fail_loud(
                f"{step_dir}/sample_manifests missing",
                "run the step's aggregator to materialize per-sample yaml manifests.",
            )
        present = set()
        for yml in sorted(mfst_dir.glob("*.yaml")):
            doc = _load_yaml(yml)
            sid = doc.get("sample_info", {}).get("sample_id")
            if sid is None:
                continue
            present.add(sid)
            aggregated.setdefault(sid, {"sample_id": sid})
            # Pull summary + outputs into flattened keys
            step_block = doc.get(step_key, {})
            aggregated[sid][f"{step_key}_status"] = step_block.get("status", "unknown")
            summary = step_block.get("summary", {}) or {}
            for k, v in summary.items():
                aggregated[sid][f"{step_key}_{k}"] = v

        missing = set(expected_samples) - present
        if missing:
            _fail_loud(
                f"{len(missing)} sample(s) missing from {step_dir} manifests",
                f"first missing: {sorted(missing)[0]}",
            )

    return aggregated


def main():
    args = _parse_args()
    _log("start")
    _log("dry_run", args.dry_run)

    for label, val in [
        ("preprocessing_root", args.preprocessing_root),
        ("central_cell_status", args.central_cell_status),
        ("raw_data_manifest", args.raw_data_manifest),
        ("output_root", args.output_root),
        ("module_configs", args.module_configs),
    ]:
        if not val:
            _fail_loud(f"{label} is empty")
        _log(label, val)

    for label, p in [
        ("central_cell_status", args.central_cell_status),
        ("raw_data_manifest", args.raw_data_manifest),
        ("module_configs", args.module_configs),
    ]:
        if not os.path.exists(p):
            _fail_loud(f"{label} not found: {p}")

    preprocessing_root = Path(args.preprocessing_root)
    if not preprocessing_root.is_dir():
        _fail_loud(f"preprocessing_root not a directory: {preprocessing_root}")

    import pandas as pd

    module_config = _load_yaml(args.module_configs)
    step_07_cfg = module_config["modules"].get("step_07_final_report")
    if step_07_cfg is None:
        _fail_loud("module_configs missing modules.step_07_final_report")
    algo_params = step_07_cfg["algorithm_params"]
    output_params = step_07_cfg.get("output_params", {"report_name": "preprocessing_qc_report.html",
                                                      "plot_dpi": 300})

    raw_manifest = pd.read_csv(args.raw_data_manifest, sep="\t",
                               dtype=str, keep_default_na=False)
    expected_samples = sorted(raw_manifest["sample_id"].unique())
    _log("expected_samples", len(expected_samples))

    output_base = Path(args.output_root)
    for sub in ("reports", "plots", "data", "logs"):
        (output_base / sub).mkdir(parents=True, exist_ok=True)

    if args.dry_run:
        print("=== DRY RUN MODE ===")
        print("Step: step_07_final_report")
        print("Container: scgpt_gpu")
        print("Inputs validated:")
        print(f"  {args.preprocessing_root}: OK")
        print(f"  {args.central_cell_status}: OK")
        print(f"  {args.module_configs}: OK")
        print(f"  expected samples (from raw-data manifest): {len(expected_samples)}")
        print("Outputs planned:")
        print(f"  {output_base / 'reports' / output_params['report_name']}")
        print(f"  {output_base / 'data' / 'flagged_samples.csv'}")
        print(f"  {output_base / 'plots' / '*.png'}")
        print("Resources: 1 CPU, 16 GB, 30m (single job)")
        print("VALIDATION PASSED")
        return 0

    # Real execution
    import matplotlib
    matplotlib.use("Agg")

    from aggregate_metrics import compute_derived_metrics
    from flag_outliers import apply_flagging_pipeline, create_flagged_records
    from generate_plots import generate_section3_plots, set_plotting_style
    from render_report import render_html_report

    aggregated = _aggregate_sample_manifests(preprocessing_root, expected_samples)
    df_metrics = pd.DataFrame(aggregated.values())
    df_metrics = compute_derived_metrics(df_metrics)
    _log("n_samples_aggregated", len(df_metrics))

    metrics_csv_path = output_base / "data" / "sample_metrics.csv"
    df_metrics.to_csv(metrics_csv_path, index=False)
    _log("wrote", metrics_csv_path)

    df_metrics, _ = apply_flagging_pipeline(
        df_metrics,
        metrics=algo_params["priority_metrics"],
        threshold=algo_params["mad_threshold"],
    )
    df_flagged = create_flagged_records(df_metrics, algo_params["priority_metrics"])
    flagged_csv_path = output_base / "data" / "flagged_samples.csv"
    df_flagged.to_csv(flagged_csv_path, index=False)
    _log("wrote", flagged_csv_path)
    _log("n_flagged", len(df_flagged))

    set_plotting_style()
    metrics_config = {
        "color_removal": "#3498db",
        "color_vf": "#2ecc71",
        "color_step02": "#f39c12",
        "color_step03": "#e74c3c",
        "title_removal_dist": "Cell Removal Rate Distribution",
        "title_removal_bar": "Cell Removal by Patient",
        "title_vf_dist": "VF Retention Rate Distribution",
        "title_vf_scatter": "Variable Feature Retention: Baseline vs Post-QC",
    }
    section3_plots = generate_section3_plots(df_metrics, metrics_config)
    plots_dir = output_base / "plots"
    for plot_name, fig in section3_plots.items():
        fig.savefig(plots_dir / f"{plot_name}.png",
                    dpi=output_params.get("plot_dpi", 300),
                    bbox_inches="tight")

    report_config = {
        "report_title": "Preprocessing QC Report",
        "mad_threshold": algo_params["mad_threshold"],
        "priority_metrics": algo_params["priority_metrics"],
        "generation_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    }
    html_content = render_html_report(
        df_metrics=df_metrics,
        df_flagged=df_flagged,
        section3_figs=section3_plots,
        appendix_plot_paths={},  # fresh canonical tree; appendix plots
                                  # can be re-wired against step_06 output
                                  # dir in a follow-up rev
        config=report_config,
    )
    report_path = output_base / "reports" / output_params.get(
        "report_name", "preprocessing_qc_report.html"
    )
    with open(report_path, "w", encoding="utf-8") as fh:
        fh.write(html_content)
    _log("wrote", report_path)
    _log("execute_end")
    return 0


if __name__ == "__main__":
    sys.exit(main())
