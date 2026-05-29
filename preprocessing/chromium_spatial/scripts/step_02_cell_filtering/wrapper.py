#!/usr/bin/env python3
# ============================================================================
# Step 02: Cell Filtering (per-patient MAD outlier removal) - Wrapper
# ============================================================================
# Purpose
#   Remove per-cluster outlier cells from each patient using MAD on log1p
#   gene/UMI counts and upper-bound MAD on mito percentage. Clustering is
#   performed at the patient level (merging all positions for that patient)
#   so MAD thresholds are computed within biologically consistent groups.
#   Outputs are disaggregated back to per-sample metadata CSVs so the
#   position-centric downstream layout is preserved.
#
# Methodology
#   - MAD intervals on log1p counts:
#       genes: [-3, 3]  -> flag very-low (debris) and very-high (doublet-ish)
#       umi:   [-3]     -> lower bound only (high UMI cells are legitimate)
#       mito:  [3]      -> upper bound only (stressed / dying cells)
#   - Mito safety floor (mito_min_threshold = 3%): never remove a cell
#     with < 3% mito regardless of MAD, because the MAD calculation is
#     unstable at very low mito percentages and would otherwise flag
#     healthy cells as high-mito outliers.
#   - Clustering resolution 0.1 (coarse), because the MAD grouping only
#     needs gross cell-type separation, not fine sub-states. A small
#     cluster-sample group (< 100 cells) is left unfiltered — MAD is
#     unreliable on tiny groups, and the cells are dominated by doublet
#     filtering in step 03.
#   - The cell-set-changing gate prints pre-MAD and post-MAD counts and
#     fails loudly if removal is outside [5%, 30%]. This catches both
#     "MAD did nothing" and "MAD ate the patient" failure modes.
#
# Single-sample rerun caveat
#   Because clustering is patient-level, rerunning one sample after step
#   01 requires re-executing the entire patient array task. Do NOT try
#   to splice a new sample into an existing patient output — the cluster
#   labels will not match. Always rerun the full patient.
#
# Inputs
#   - CFG_PREPROCESSING_CENTRAL_CELL_STATUS (central_cell_status.csv,
#     produced by the step 01 aggregator).
#   - CFG_PREPROCESSING_STEP_01/features/{sample}_vfs.txt and
#     cell_metadata/{sample}_metadata.csv for each sample in the patient.
#   - Raw CellRanger H5 paths resolved via CFG_RAW_DATA_SAMPLE_MANIFEST.
#   - Algorithm params from CFG_MODULE_CONFIGS
#     (modules.step_02_cell_filtering.algorithm_params).
#
# Outputs
#   - ${CFG_PREPROCESSING_STEP_02}/cell_metadata/{sample}_metadata.csv
#   - ${CFG_PREPROCESSING_STEP_02}/union_vfs/{patient}_union_vfs.txt
#   - ${CFG_PREPROCESSING_STEP_02}/manifest_updates/{patient}_step_02_update.csv
#   - ${CFG_PREPROCESSING_STEP_02}/visualizations/{patient}_*.pdf
#   - ${CFG_PREPROCESSING_STEP_02}/sample_manifests/{sample}.yaml
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


def _log(action: str, value=None, level: str = "INFO") -> None:
    ts = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    if value is None:
        print(f"[{ts}] {level} step_02 {action}", flush=True)
    else:
        print(f"[{ts}] {level} step_02 {action}={value}", flush=True)


def _fail_loud(msg: str, hint: str | None = None) -> "NoReturn":
    _log("FAIL", msg, level="ERROR")
    if hint:
        print(f"  hint: {hint}", flush=True)
    sys.exit(2)


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Step 02: Cell Filtering (patient-level)")
    p.add_argument(
        "--patient-id",
        default=os.environ.get("PATIENT_ID", ""),
        help="Patient ID to process (e.g. Pat1). Required. SLURM wrapper maps array index -> patient.",
    )
    p.add_argument(
        "--central-cell-status",
        default=os.environ.get("CFG_PREPROCESSING_CENTRAL_CELL_STATUS", ""),
        help="Path to central_cell_status.csv (defaults to CFG_PREPROCESSING_CENTRAL_CELL_STATUS).",
    )
    p.add_argument(
        "--raw-data-manifest",
        default=os.environ.get("CFG_RAW_DATA_SAMPLE_MANIFEST", ""),
        help="TSV with patient_id, position_id, sample_id, h5_path.",
    )
    p.add_argument(
        "--step-01-root",
        default=os.environ.get("CFG_PREPROCESSING_STEP_01", ""),
        help="Step 01 output root (reads per-sample VFs and metadata).",
    )
    p.add_argument(
        "--output-root",
        default=os.environ.get("CFG_PREPROCESSING_STEP_02", ""),
        help="Step 02 output root.",
    )
    p.add_argument(
        "--module-configs",
        default=os.environ.get(
            "CFG_MODULE_CONFIGS",
            os.path.join(os.environ.get("CFG_PROJECT_ROOT", ""), "config", "module_configs.yaml"),
        ),
        help="Path to module_configs.yaml.",
    )
    p.add_argument(
        "--min-removal-pct", type=float, default=5.0,
        help="Fail-loud lower bound on cell removal percentage (safety gate).",
    )
    p.add_argument(
        "--max-removal-pct", type=float, default=30.0,
        help="Fail-loud upper bound on cell removal percentage (safety gate).",
    )
    p.add_argument("--dry-run", action="store_true", help="Validate inputs and exit 0.")
    return p.parse_args()


def _load_yaml(path: str):
    with open(path, "r") as fh:
        return yaml.safe_load(fh)


def _write_yaml_atomic(data, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + f".tmp.{os.getpid()}")
    with open(tmp, "w") as fh:
        yaml.dump(data, fh, default_flow_style=False, sort_keys=False)
    tmp.replace(path)


def main() -> int:
    args = _parse_args()
    _log("start")
    _log("dry_run", args.dry_run)

    # ------------------------------------------------------------------
    # Validation preamble
    # ------------------------------------------------------------------
    if not args.patient_id:
        _fail_loud(
            "patient-id is empty",
            "set --patient-id or PATIENT_ID (SLURM wrapper maps task index -> patient).",
        )
    _log("patient_id", args.patient_id)

    for label, val in [
        ("central_cell_status", args.central_cell_status),
        ("raw_data_manifest", args.raw_data_manifest),
        ("step_01_root", args.step_01_root),
        ("output_root", args.output_root),
        ("module_configs", args.module_configs),
    ]:
        if not val:
            _fail_loud(f"{label} path is empty", f"set the corresponding CFG_* env or CLI flag.")
        _log(label, val)

    for label, path in [
        ("central_cell_status", args.central_cell_status),
        ("raw_data_manifest", args.raw_data_manifest),
        ("module_configs", args.module_configs),
    ]:
        if not os.path.exists(path):
            _fail_loud(f"{label} not found: {path}")
        size = os.path.getsize(path)
        if size == 0:
            _fail_loud(f"{label} is empty (0 bytes): {path}")
        _log(f"{label}_bytes", size)

    if not os.path.isdir(args.step_01_root):
        _fail_loud(f"step_01_root is not a directory: {args.step_01_root}")

    output_base = Path(args.output_root)
    output_base.mkdir(parents=True, exist_ok=True)
    if not os.access(output_base, os.W_OK):
        _fail_loud(f"output root not writable: {output_base}")

    # ------------------------------------------------------------------
    # Load inputs that the dry-run preamble needs to know about
    # ------------------------------------------------------------------
    import pandas as pd  # import here so the dry-run preamble is fast

    module_config = _load_yaml(args.module_configs)
    step_02_cfg = module_config["modules"].get("step_02_cell_filtering")
    if step_02_cfg is None:
        _fail_loud("module_configs.yaml missing modules.step_02_cell_filtering")
    algo_params = step_02_cfg["algorithm_params"]
    _log("mad_thresholds_genes", algo_params["mad_thresholds"]["genes"])
    _log("mad_thresholds_umi", algo_params["mad_thresholds"]["umi"])
    _log("mad_thresholds_mito", algo_params["mad_thresholds"]["mito"])
    _log("mito_min_threshold", algo_params["mito_min_threshold"])

    raw_manifest = pd.read_csv(
        args.raw_data_manifest, sep="\t", dtype=str, keep_default_na=False
    )
    required_cols = {"patient_id", "position_id", "sample_id", "h5_path"}
    missing = required_cols - set(raw_manifest.columns)
    if missing:
        _fail_loud(
            f"raw_data_manifest missing columns: {sorted(missing)}",
            "schema must be patient_id, position_id, sample_id, h5_path",
        )
    patient_raw = raw_manifest[raw_manifest["patient_id"] == args.patient_id].copy()
    if len(patient_raw) == 0:
        _fail_loud(f"no samples for patient {args.patient_id} in raw_data_manifest")
    _log("patient_samples", len(patient_raw))

    central_manifest = pd.read_csv(args.central_cell_status)
    if "step_01_umi_pass" not in central_manifest.columns:
        _fail_loud(
            "central_cell_status missing step_01_umi_pass column",
            "run the step 01 aggregator before step 02.",
        )
    central_manifest["step_01_umi_pass"] = central_manifest["step_01_umi_pass"].astype(bool)

    patient_mask = (
        (central_manifest["patient_id"] == args.patient_id)
        & (central_manifest["step_01_umi_pass"])
    )
    patient_cells_step01 = central_manifest[patient_mask].copy()
    n_cells_pre_mad = len(patient_cells_step01)
    if n_cells_pre_mad == 0:
        _fail_loud(
            f"no step_01-passing cells for patient {args.patient_id}",
            "confirm step 01 completed and the aggregator ran.",
        )
    _log("n_cells_pre_mad", n_cells_pre_mad)

    # Pre-validate that every sample has step 01 sidecars
    sample_ids = patient_cells_step01["sample_id"].unique().tolist()
    step_01_root = Path(args.step_01_root)
    missing_sidecars = []
    for sid in sample_ids:
        vf_path = step_01_root / "features" / f"{sid}_vfs.txt"
        meta_path = step_01_root / "cell_metadata" / f"{sid}_metadata.csv"
        if not vf_path.exists():
            missing_sidecars.append(str(vf_path))
        if not meta_path.exists():
            missing_sidecars.append(str(meta_path))
    if missing_sidecars:
        _fail_loud(
            f"{len(missing_sidecars)} step-01 sidecars missing for patient {args.patient_id}",
            f"first missing: {missing_sidecars[0]}",
        )

    # Prepare output subdirs up front (so dry-run prints realistic plan)
    dirs = {
        "cell_metadata": output_base / "cell_metadata",
        "union_vfs": output_base / "union_vfs",
        "manifest_updates": output_base / "manifest_updates",
        "visualizations": output_base / "visualizations",
        "sample_manifests": output_base / "sample_manifests",
        "logs": output_base / "logs",
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)

    union_vf_path = dirs["union_vfs"] / f"{args.patient_id}_union_vfs.txt"
    update_csv_path = dirs["manifest_updates"] / f"{args.patient_id}_step_02_update.csv"

    if args.dry_run:
        print("=== DRY RUN MODE ===")
        print("Step: step_02_cell_filtering")
        print("Container: scgpt_gpu")
        print("Inputs validated:")
        print(f"  {args.central_cell_status}: OK")
        print(f"  {args.raw_data_manifest}: OK")
        print(f"  {args.module_configs}: OK")
        print(f"  step_01_root sidecars: OK ({len(sample_ids)} samples)")
        print("Outputs planned:")
        print(f"  {update_csv_path}")
        print(f"  {union_vf_path}")
        print(f"  {dirs['cell_metadata']}/{{sample_id}}_metadata.csv (x{len(sample_ids)})")
        print("Resources: 4 CPUs, 32 GB, 6h (SLURM array 0-3)")
        print("VALIDATION PASSED")
        return 0

    # ------------------------------------------------------------------
    # Heavy imports
    # ------------------------------------------------------------------
    import numpy as np
    import scanpy as sc
    from cell_filter import run_cell_filtering
    from visualizations import generate_all_visualizations

    # ------------------------------------------------------------------
    # Load per-sample step-01 sidecars + raw H5
    # ------------------------------------------------------------------
    log_path = dirs["logs"] / f"{args.patient_id}.log"
    log_fh = open(log_path, "w")

    def tee(msg: str) -> None:
        print(msg, flush=True)
        log_fh.write(msg + "\n")

    try:
        passing_samples = []
        vf_sets = []
        adata_list = []

        position_lookup = {
            row.sample_id: row.position_id for row in patient_raw.itertuples(index=False)
        }
        h5_lookup = {
            row.sample_id: row.h5_path for row in patient_raw.itertuples(index=False)
        }

        for sid in sample_ids:
            if sid not in position_lookup:
                _fail_loud(
                    f"sample {sid} in central_cell_status but not in raw_data_manifest",
                    "raw_data_manifest must be a superset of central_cell_status samples.",
                )
            position_id = position_lookup[sid]
            h5_path = h5_lookup[sid]
            if not os.path.exists(h5_path):
                _fail_loud(f"raw H5 not found for {sid}: {h5_path}")

            vf_path = step_01_root / "features" / f"{sid}_vfs.txt"
            with open(vf_path, "r") as fh:
                vfs = {line.strip() for line in fh if line.strip()}
            vf_sets.append(vfs)

            sample_cells = patient_cells_step01[patient_cells_step01["sample_id"] == sid].copy()
            unprefixed = sample_cells["cell_id"].str.replace(f"{sid}_", "", regex=False)

            adata = sc.read_10x_h5(h5_path)
            adata.var_names_make_unique()
            adata = adata[adata.obs.index.isin(set(unprefixed))].copy()
            adata.obs.index = sid + "_" + adata.obs.index
            adata.obs["sample_id"] = sid
            adata.obs["patient_id"] = args.patient_id

            n_umi_lookup = sample_cells.set_index("cell_id")["step_01_n_umi"]
            adata.obs["n_umi"] = n_umi_lookup.loc[adata.obs.index].values

            adata_list.append(adata)
            passing_samples.append(
                {
                    "sample_id": sid,
                    "position_id": position_id,
                    "n_cells_input": adata.n_obs,
                }
            )
            tee(f"  {sid}: loaded {adata.n_obs} cells from {h5_path}")

        # Subset each AnnData to the union VF list (same pattern as source)
        union_vfs = sorted(set.union(*vf_sets))
        for i, adata in enumerate(adata_list):
            common = adata.var_names.intersection(union_vfs)
            adata_list[i] = adata[:, common].copy()

        # --------------------------------------------------------------
        # Algorithm
        # --------------------------------------------------------------
        tee("Running run_cell_filtering...")
        algo_config = {"module_2": algo_params}
        adata_patient = run_cell_filtering(
            vf_sets=vf_sets,
            adata_list=adata_list,
            patient_id=args.patient_id,
            config=algo_config,
        )
        tee(f"  patient object: {adata_patient.n_obs} cells x {adata_patient.n_vars} genes")

        # --------------------------------------------------------------
        # Cell-set-changing gate
        # --------------------------------------------------------------
        n_passing = int(adata_patient.obs["cell_filter_pass_final"].sum())
        n_removed = n_cells_pre_mad - n_passing
        removal_pct = 100.0 * n_removed / max(n_cells_pre_mad, 1)
        _log("n_cells_pre_mad", n_cells_pre_mad)
        _log("n_cells_post_mad", n_passing)
        _log("removed_pct", f"{removal_pct:.2f}")
        tee(f"  pre-MAD cells: {n_cells_pre_mad}")
        tee(f"  post-MAD cells: {n_passing}")
        tee(f"  removal %: {removal_pct:.2f}")

        if removal_pct < args.min_removal_pct:
            _fail_loud(
                f"MAD removal {removal_pct:.2f}% below lower safety bound {args.min_removal_pct}%",
                "this usually means MAD thresholds did not fire -- check cluster sizes.",
            )
        if removal_pct > args.max_removal_pct:
            _fail_loud(
                f"MAD removal {removal_pct:.2f}% above upper safety bound {args.max_removal_pct}%",
                "check mito_min_threshold and cluster_resolution_coarse -- MAD may be eating healthy cells.",
            )

        # --------------------------------------------------------------
        # Create exclusion flags (comma-separated per cell)
        # --------------------------------------------------------------
        def flag_row(row):
            flags = []
            if row.get("debris_prefilter_removed", False):
                flags.append("mad_outlier_umi_low")
            if row.get("mad_outlier_genes_low", False) and "mad_outlier_genes_low" not in flags:
                flags.append("mad_outlier_genes_low")
            if row.get("mad_outlier_genes_high", False):
                flags.append("mad_outlier_genes_high")
            if row.get("mad_outlier_umi_low", False) and "mad_outlier_umi_low" not in flags:
                flags.append("mad_outlier_umi_low")
            if row.get("mad_outlier_mito_high", False):
                flags.append("mad_outlier_mito_high")
            return ",".join(flags)

        adata_patient.obs["step_02_exclusion_flags"] = adata_patient.obs.apply(flag_row, axis=1)

        # --------------------------------------------------------------
        # Visualizations (non-critical -- fail loud if the library does,
        # but do not let a plotting hiccup destroy the filtered metadata)
        # --------------------------------------------------------------
        try:
            generate_all_visualizations(
                adata_patient, dirs["visualizations"], args.patient_id
            )
        except Exception as e:
            _fail_loud(
                f"visualization generation failed: {e}",
                "fix the plotting module before rerunning -- do not silently skip.",
            )

        # --------------------------------------------------------------
        # Write manifest update CSV (atomic)
        # --------------------------------------------------------------
        update_df = pd.DataFrame(
            {
                "cell_id": adata_patient.obs.index,
                "sample_id": adata_patient.obs["sample_id"],
                "patient_id": adata_patient.obs["patient_id"],
                "step_02_mad_pass": adata_patient.obs["cell_filter_pass_final"],
                "step_02_exclusion_flags": adata_patient.obs["step_02_exclusion_flags"],
            }
        )
        tmp_update = update_csv_path.with_suffix(update_csv_path.suffix + f".tmp.{os.getpid()}")
        update_df.to_csv(tmp_update, index=False)
        tmp_update.replace(update_csv_path)
        _log("wrote", update_csv_path)
        _log("bytes", os.path.getsize(update_csv_path))

        # Union VFs
        tmp_vfs = union_vf_path.with_suffix(union_vf_path.suffix + f".tmp.{os.getpid()}")
        tmp_vfs.write_text("\n".join(union_vfs) + "\n")
        tmp_vfs.replace(union_vf_path)
        _log("wrote", union_vf_path)

        # --------------------------------------------------------------
        # Per-sample metadata CSVs + sample manifests
        # --------------------------------------------------------------
        for entry in passing_samples:
            sid = entry["sample_id"]
            mask = adata_patient.obs["sample_id"] == sid
            sample_meta = adata_patient.obs[mask].copy()
            metadata_csv_path = dirs["cell_metadata"] / f"{sid}_metadata.csv"
            sample_meta.to_csv(metadata_csv_path, index=True)

            sample_manifest = {
                "sample_info": {
                    "patient_id": args.patient_id,
                    "position_id": entry["position_id"],
                    "sample_id": sid,
                    "position_status": "step_02_complete",
                },
                "step_02_cell_filtering": {
                    "timestamp": datetime.now().isoformat(),
                    "status": "pass",
                    "patient_id": args.patient_id,
                    "outputs": {
                        "cell_metadata": str(metadata_csv_path),
                        "union_vfs": str(union_vf_path),
                        "log": str(log_path),
                    },
                    "summary": {
                        "n_cells_input": int(len(sample_meta)),
                        "n_cells_passing": int(sample_meta["cell_filter_pass_final"].sum()),
                    },
                },
            }
            _write_yaml_atomic(sample_manifest, dirs["sample_manifests"] / f"{sid}.yaml")

        _log("execute_end")
        return 0
    finally:
        log_fh.close()


if __name__ == "__main__":
    sys.exit(main())
