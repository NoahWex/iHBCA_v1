#!/usr/bin/env python3
# ============================================================================
# Step 06: Filtered Preview (per-patient before/after QC visualization)
# ============================================================================
# Purpose
#   Rebuild a patient-level AnnData from the post-QC cell set using the
#   filtered VFs from step 05, compute UMAPs, and render the before/after
#   QC evidence panels used in the step 07 HTML report and the
#   supplemental preprocessing figures. This step is visualization-only;
#   it does NOT change the cell set.
#
# Methodology
#   - Patient-level merge of all passing samples (needed for consistent
#     cluster context across positions).
#   - Uses Step 02's union VFs for the "unfiltered" reference panel and
#     Step 05's filtered VFs for the "filtered" rebuilt panel. The two
#     side-by-side UMAPs let a reviewer see which cell populations were
#     removed by MAD (step 02) vs doublet detection (step 03).
#   - Plot dimensions come from the module config (umap_plot_size,
#     plot_dpi) so that aesthetics changes do not require editing the
#     script.
#
# Inputs
#   - CFG_PREPROCESSING_CENTRAL_CELL_STATUS
#   - CFG_PREPROCESSING_STEP_02 (merged H5AD, union VFs, cluster labels)
#   - CFG_PREPROCESSING_STEP_03 (doublet scores and pass flags)
#   - CFG_PREPROCESSING_STEP_05 (filtered VF list per sample)
#   - Raw H5 from CFG_RAW_DATA_SAMPLE_MANIFEST
#
# Outputs
#   - ${CFG_PREPROCESSING_STEP_06}/patient_summaries/{patient}_*.png (8 per patient)
#   - ${CFG_PREPROCESSING_STEP_06}/sample_detail_faceted/{patient}/{sample}_*.png
#   - ${CFG_PREPROCESSING_STEP_06}/sample_manifests/{sample}.yaml
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
        print(f"[{ts}] {level} step_06 {action}", flush=True)
    else:
        print(f"[{ts}] {level} step_06 {action}={value}", flush=True)


def _fail_loud(msg, hint=None):
    _log("FAIL", msg, level="ERROR")
    if hint:
        print(f"  hint: {hint}", flush=True)
    sys.exit(2)


def _parse_args():
    p = argparse.ArgumentParser(description="Step 06: Filtered Preview (per-patient)")
    p.add_argument("--patient-id", default=os.environ.get("PATIENT_ID", ""))
    p.add_argument("--central-cell-status",
                   default=os.environ.get("CFG_PREPROCESSING_CENTRAL_CELL_STATUS", ""))
    p.add_argument("--raw-data-manifest",
                   default=os.environ.get("CFG_RAW_DATA_SAMPLE_MANIFEST", ""))
    p.add_argument("--step-02-root",
                   default=os.environ.get("CFG_PREPROCESSING_STEP_02", ""))
    p.add_argument("--step-03-root",
                   default=os.environ.get("CFG_PREPROCESSING_STEP_03", ""))
    p.add_argument("--step-04-root",
                   default=os.environ.get("CFG_PREPROCESSING_STEP_04", ""),
                   help="Optional; used to pull HBCA label columns if present.")
    p.add_argument("--step-05-root",
                   default=os.environ.get("CFG_PREPROCESSING_STEP_05", ""))
    p.add_argument("--output-root",
                   default=os.environ.get("CFG_PREPROCESSING_STEP_06", ""))
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


def _write_yaml_atomic(data, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + f".tmp.{os.getpid()}")
    with open(tmp, "w") as fh:
        yaml.dump(data, fh, default_flow_style=False, sort_keys=False)
    tmp.replace(path)


def main():
    args = _parse_args()
    _log("start")
    _log("dry_run", args.dry_run)

    if not args.patient_id:
        _fail_loud("patient-id empty",
                   "set --patient-id or PATIENT_ID; SLURM wrapper maps array idx -> patient.")

    for label, val in [
        ("central_cell_status", args.central_cell_status),
        ("raw_data_manifest", args.raw_data_manifest),
        ("step_02_root", args.step_02_root),
        ("step_03_root", args.step_03_root),
        ("step_05_root", args.step_05_root),
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

    for root in (args.step_02_root, args.step_03_root, args.step_05_root):
        if not os.path.isdir(root):
            _fail_loud(f"sidecar root not a directory: {root}")

    import pandas as pd

    module_config = _load_yaml(args.module_configs)
    step_06_cfg = module_config["modules"].get("step_06_filtered_preview")
    if step_06_cfg is None:
        _fail_loud("module_configs missing modules.step_06_filtered_preview")
    algo_params = step_06_cfg["algorithm_params"]
    processing_params = step_06_cfg.get("processing_params", {})
    _log("umap_plot_size", algo_params.get("umap_plot_size"))
    _log("plot_dpi", algo_params.get("plot_dpi"))

    raw_manifest = pd.read_csv(args.raw_data_manifest, sep="\t",
                               dtype=str, keep_default_na=False)
    patient_raw = raw_manifest[raw_manifest["patient_id"] == args.patient_id]
    if len(patient_raw) == 0:
        _fail_loud(f"no samples for patient {args.patient_id} in raw-data manifest")

    central_manifest = pd.read_csv(args.central_cell_status)
    for col in ("step_01_umi_pass", "step_02_mad_pass", "step_03_doublet_pass"):
        if col not in central_manifest.columns:
            _fail_loud(f"central_cell_status missing {col}")
        central_manifest[col] = central_manifest[col].astype(bool)

    patient_cells = central_manifest[central_manifest["patient_id"] == args.patient_id].copy()
    if len(patient_cells) == 0:
        _fail_loud(f"no cells for patient {args.patient_id}")

    sample_ids = sorted(patient_cells["sample_id"].unique())
    _log("n_samples", len(sample_ids))

    # Sidecar presence validation (fail loud)
    step_02_merged = Path(args.step_02_root) / "merged_objects" / f"{args.patient_id}_merged_post_module2.h5ad"
    missing = []
    for sid in sample_ids:
        if not (Path(args.step_03_root) / "metadata" / f"{sid}.csv").exists():
            missing.append(f"step_03 metadata for {sid}")
        if not (Path(args.step_05_root) / "features" / "filtered_vfs" / f"{sid}_vfs_filtered.txt").exists():
            missing.append(f"step_05 filtered_vfs for {sid}")
    if missing:
        _fail_loud(
            f"{len(missing)} sidecars missing for patient {args.patient_id}",
            f"first missing: {missing[0]}",
        )

    output_base = Path(args.output_root)
    dirs = {
        "reprocessed_objects": output_base / "reprocessed_objects",
        "patient_summaries": output_base / "patient_summaries",
        "sample_detail_faceted": output_base / "sample_detail_faceted" / args.patient_id,
        "sample_manifests": output_base / "sample_manifests",
        "logs": output_base / "logs",
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)

    if args.dry_run:
        print("=== DRY RUN MODE ===")
        print("Step: step_06_filtered_preview")
        print("Container: scgpt_gpu")
        print("Inputs validated:")
        print(f"  {args.central_cell_status}: OK")
        print(f"  {step_02_merged}: " + ("OK" if step_02_merged.exists() else "PRESENT-AT-EXECUTION"))
        print(f"  step_03/step_05 sidecars: OK ({len(sample_ids)} samples)")
        print("Outputs planned:")
        print(f"  {dirs['patient_summaries']}/{args.patient_id}_*.png (8 panels)")
        print(f"  {dirs['sample_detail_faceted']}/{{sample}}_*.png (x{len(sample_ids)})")
        print("Resources: 8 CPUs, 64 GB, 1h (SLURM array 0-3)")
        print("VALIDATION PASSED")
        return 0

    # -----------------------------------------------------------------
    # Real execution
    # -----------------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import scanpy as sc

    from object_rebuilder import rebuild_filtered_object
    from mad_calculator import compute_mad_standardized_metrics, compute_qc_metrics_if_missing
    from removal_tracker import assign_removal_status, compute_final_pass_flag
    from plot_patient_summary import generate_patient_summary_plots
    from plot_sample_detail import generate_sample_detail_plot
    from plotting_utils import set_plotting_style

    log_path = dirs["logs"] / f"{args.patient_id}.log"
    log_fh = open(log_path, "w")

    try:
        set_plotting_style()

        # Load merged step-02 object (the source of the "unfiltered" reference)
        if not step_02_merged.exists():
            _fail_loud(f"step-02 merged H5AD missing: {step_02_merged}",
                       "step 02 must write merged_objects/ before step 06.")
        adata_step02 = sc.read_h5ad(step_02_merged)
        _log("adata_step02_cells", adata_step02.n_obs)

        required_cols = ["cluster_coarse", "sample_id", "patient_id"]
        missing_obs = [c for c in required_cols if c not in adata_step02.obs.columns]
        if missing_obs:
            _fail_loud(f"step-02 H5AD missing obs columns: {missing_obs}")
        if "X_umap" not in adata_step02.obsm:
            _fail_loud("step-02 H5AD missing UMAP coords")

        # Join QC filter columns from central manifest (canonical source)
        patient_cells_indexed = patient_cells.set_index("cell_id")
        adata_step02.obs = adata_step02.obs.join(
            patient_cells_indexed[["step_02_mad_pass", "step_03_doublet_pass"]],
            how="left",
        )
        adata_step02.obs["cell_filter_pass_final"] = adata_step02.obs["step_02_mad_pass"]
        adata_step02.obs["doublet_filter_pass"] = adata_step02.obs["step_03_doublet_pass"]

        # Doublet scores from step-03 sidecars
        doublet_dfs = []
        for sid in sample_ids:
            csv_path = Path(args.step_03_root) / "metadata" / f"{sid}.csv"
            doublet_dfs.append(pd.read_csv(csv_path))
        doublet_merged = pd.concat(doublet_dfs, ignore_index=True)
        adata_step02.obs = adata_step02.obs.join(
            doublet_merged.set_index("cell_barcode")[["scDblFinder_score"]], how="left"
        )

        # Optional step-04 labels
        if args.step_04_root and os.path.isdir(args.step_04_root):
            celltype_dfs = []
            for sid in sample_ids:
                csv = Path(args.step_04_root) / "cell_metadata" / f"{sid}_labels.csv"
                if csv.exists():
                    celltype_dfs.append(pd.read_csv(csv, index_col=0))
            if celltype_dfs:
                merged = pd.concat(celltype_dfs, axis=0)
                cols = [c for c in ["predicted.celltype.SC", "predicted.celltype.SN",
                                    "HBCATransferredLabels.Kumar_2023"] if c in merged.columns]
                if cols:
                    adata_step02.obs = adata_step02.obs.join(merged[cols], how="left")

        adata_step02 = compute_final_pass_flag(adata_step02)
        adata_step02 = assign_removal_status(adata_step02)
        n_final_pass = int(adata_step02.obs["final_pass"].sum())
        _log("n_final_pass", n_final_pass)

        # Rebuild filtered per-sample AnnDatas from raw H5 + filtered VFs
        h5_lookup = {r.sample_id: r.h5_path for r in patient_raw.itertuples(index=False)}
        adata_list = []
        for sid in sample_ids:
            h5_path = h5_lookup[sid]
            if not os.path.exists(h5_path):
                _fail_loud(f"raw H5 missing for {sid}: {h5_path}")
            adata = sc.read_10x_h5(h5_path)
            adata.var_names_make_unique()
            passing_prefixed = adata_step02.obs[
                (adata_step02.obs["sample_id"] == sid)
                & (adata_step02.obs["final_pass"] == True)  # noqa: E712
            ].index
            passing_bare = [bc.replace(f"{sid}_", "") for bc in passing_prefixed]
            adata = adata[adata.obs.index.isin(passing_bare)].copy()
            adata.obs.index = sid + "_" + adata.obs.index
            adata.obs["sample_id"] = sid
            adata.obs["patient_id"] = args.patient_id
            adata_list.append(adata)

        vf_sets = []
        for sid in sample_ids:
            vf_path = Path(args.step_05_root) / "features" / "filtered_vfs" / f"{sid}_vfs_filtered.txt"
            with open(vf_path, "r") as fh:
                vf_sets.append({line.strip() for line in fh if line.strip()})
        union_vfs = sorted(set.union(*vf_sets))

        passing_obs = adata_step02.obs[adata_step02.obs["final_pass"] == True]  # noqa: E712
        cluster_labels = passing_obs["cluster_coarse"].copy()

        adata_rebuilt = rebuild_filtered_object(
            adata_list=adata_list,
            filtered_vfs=union_vfs,
            cluster_labels=cluster_labels,
            processing_params=processing_params,
        )
        adata_rebuilt = compute_qc_metrics_if_missing(adata_rebuilt)

        # Build the full (pre-filter) plotting dataset from step-02 object
        adata_full = adata_step02.copy()
        common = adata_full.var_names.intersection(union_vfs)
        adata_full = adata_full[:, common].copy()
        adata_full = compute_qc_metrics_if_missing(adata_full)
        adata_full = compute_mad_standardized_metrics(
            adata_full, min_group_size=algo_params.get("min_cluster_sample_size", 100)
        )

        plot_params = {
            "plot_size": algo_params.get("umap_plot_size", [12, 8]),
            "dpi": algo_params.get("plot_dpi", 300),
        }

        patient_figs = generate_patient_summary_plots(
            adata_full=adata_full,
            adata_filtered=adata_rebuilt,
            patient_id=args.patient_id,
            plot_params=plot_params,
        )
        patient_plot_paths = {}
        for plot_name, fig in patient_figs.items():
            p = dirs["patient_summaries"] / f"{args.patient_id}_{plot_name}.png"
            fig.savefig(p, dpi=plot_params["dpi"], bbox_inches="tight")
            plt.close(fig)
            patient_plot_paths[plot_name] = str(p)

        sample_detail_paths = {}
        for sid in sample_ids:
            sample_figs = generate_sample_detail_plot(
                adata_full=adata_full,
                adata_filtered=adata_rebuilt,
                patient_id=args.patient_id,
                sample_id=sid,
                plot_params=plot_params,
            )
            for plot_name, fig in sample_figs.items():
                p = dirs["sample_detail_faceted"] / f"{sid}_{plot_name}.png"
                fig.savefig(p, dpi=plot_params["dpi"], bbox_inches="tight")
                plt.close(fig)
                sample_detail_paths[f"{sid}_{plot_name}"] = str(p)

        # Sample manifests
        for sid in sample_ids:
            position_id = patient_raw[patient_raw["sample_id"] == sid]["position_id"].iloc[0]
            _write_yaml_atomic(
                {
                    "sample_info": {
                        "patient_id": args.patient_id,
                        "position_id": position_id,
                        "sample_id": sid,
                        "position_status": "step_06_complete",
                    },
                    "step_06_filtered_preview": {
                        "timestamp": datetime.now().isoformat(),
                        "status": "pass",
                        "outputs": {
                            "patient_summary_plots": patient_plot_paths,
                            "sample_detail_plots": {
                                k: v for k, v in sample_detail_paths.items() if k.startswith(sid + "_")
                            },
                            "log": str(log_path),
                        },
                    },
                },
                dirs["sample_manifests"] / f"{sid}.yaml",
            )

        _log("execute_end")
        return 0
    finally:
        log_fh.close()


if __name__ == "__main__":
    sys.exit(main())
