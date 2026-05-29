#!/usr/bin/env python3
"""
step_17_integration_sweep / wrapper.py

Unified orchestrator for the multi-config integration sweep. The
`--target` flag selects which atlas slice to integrate; one wrapper
serves both the full-object and per-compartment entry points.

Targets:
  full : FIVE-WAY-filtered full atlas (all cells from
         CFG_PREPROCESSING_CENTRAL_CELL_STATUS where step_15_pass == True)
  epi  : Epithelial compartment only (from CFG_PREPROCESSING_CELL_COMPARTMENTS)
  str  : Stromal compartment only
  imm  : Immune compartment only
  all  : runs all of {full, epi, str, imm} in sequence

Usage (invoked by run/run_step_17_sweep.sh SLURM array):
  python wrapper.py \
      --target full \
      --config-idx 0 \
      --project-root "$CFG_PROJECT_ROOT"

Layout under CFG_PREPROCESSING_STEP_17:
  sweep/full_object/{config}/integrated.h5ad
  sweep/compartments/{compartment}/{config}/integrated.h5ad
  sweep/scib_results/                          (per-config metric JSONs)
  winner/                                      (postprocess sidecars)
  diagnostics/svd/                              (SVD report + JSON)
"""

import argparse
import os
import subprocess
import sys

import yaml


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SOURCE_DIR = os.path.join(SCRIPT_DIR, "source")
DEFAULT_CONFIG_PATH = os.path.join(SOURCE_DIR, "configs", "sweep_configs.yaml")

TARGETS = ["full", "epi", "str", "imm"]
COMPARTMENT_NAME = {
    "epi": "Epithelial",
    "str": "Stromal",
    "imm": "Immune",
}


def _die(msg, hint=None):
    sys.stderr.write(f"FATAL: {msg}\n")
    if hint:
        sys.stderr.write(f"  hint: {hint}\n")
    sys.exit(2)


def _load_yaml(path):
    if not os.path.exists(path):
        _die(f"sweep_configs.yaml not found: {path}")
    with open(path) as fh:
        return yaml.safe_load(fh)


def _parse_paths_yaml(paths_yaml_path):
    """Parse config/paths.yaml with ${ref} expansion.

    Uses the same expansion semantics as run/_common.sh _expand_refs: flat
    one-level substitution against already-loaded keys. We re-implement it
    in Python so the wrapper can be invoked directly (for testing) without
    going through the bash entry point.
    """
    if not os.path.exists(paths_yaml_path):
        _die(f"paths.yaml not found: {paths_yaml_path}")
    raw = {}
    with open(paths_yaml_path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line or line.lstrip().startswith("#"):
                continue
            if ":" not in line:
                continue
            # Only process top-level (non-indented) keys — matches _yaml_get
            if line[0].isspace():
                continue
            key, _, val = line.partition(":")
            raw[key.strip()] = val.strip()

    # Expand ${ref} iteratively until fixed point or failure
    expanded = {}
    for k, v in raw.items():
        expanded[k] = v
    changed = True
    max_iter = 10
    while changed and max_iter > 0:
        changed = False
        for k, v in list(expanded.items()):
            if "${" in v:
                new_v = v
                # Replace each ${ref} where ref is a known key
                for ref_key, ref_val in expanded.items():
                    token = "${" + ref_key + "}"
                    if token in new_v and "${" not in ref_val:
                        new_v = new_v.replace(token, ref_val)
                if new_v != v:
                    expanded[k] = new_v
                    changed = True
        max_iter -= 1
    return expanded


def _resolve_paths(project_root):
    """Load paths.yaml and return the keys the wrapper needs.

    Returns a dict with: step_17 root, step 16 VF dir, step 14 compartments,
    step 15 contamination cells, central cell status, raw data manifest,
    canonical sweep subdirs.
    """
    paths_yaml = os.path.join(project_root, "config", "paths.yaml")
    p = _parse_paths_yaml(paths_yaml)

    needed = [
        "preprocessing_step_17",
        "preprocessing_step16_vf_dir",
        "preprocessing_cell_compartments",
        "preprocessing_central_cell_status",
        "raw_data_sample_manifest",
        "canonical_sweep_full_object",
        "canonical_sweep_compartments",
        "canonical_sweep_scib",
        "canonical_sweep_winner",
        "canonical_sweep_svd",
    ]
    missing = [k for k in needed if k not in p]
    if missing:
        _die(f"paths.yaml missing keys: {missing}")

    return {k: p[k] for k in needed}


def _select_config(cfg, config_idx):
    configs = cfg["configs"]
    if config_idx < 0 or config_idx >= len(configs):
        _die(f"--config-idx {config_idx} out of range [0, {len(configs)-1}]")
    return configs[config_idx]


def _build_cell_list_for_target(target, paths, project_root, temp_dir):
    """Return a cell_list path appropriate to the target.

    For full: the central_cell_status.csv gated to step_15_pass == True.
    For a compartment: derive a per-compartment cell list by intersecting
    step 14 cell_compartments.csv with step 15 contamination flags.

    We do NOT precompute compartment cell lists in a separate step — Phase B
    generates them on the fly in {temp_dir}, which is cleaned up by _common.sh
    when the SLURM job exits. This keeps the canonical tree from accumulating
    derived cell lists that duplicate upstream info.

    Returns the path to the cell list on disk.

    NOTE: this routine runs ONCE per SLURM array task. It is safe to read
    the upstream CSVs and write small filtered CSVs for that task.
    """
    import pandas as pd

    central = paths["preprocessing_central_cell_status"]
    if not os.path.exists(central):
        _die(f"central_cell_status.csv not found: {central}",
             hint="step 16 has not written the step_15_pass column yet. "
                  "See L3_spec.md step_16 central manifest update.")

    os.makedirs(temp_dir, exist_ok=True)

    # Canonical 5-way filter: intersection of all per-step pass columns.
    # step_15_pass alone admits cells that earlier steps rejected (~321K);
    # the intersection yields the true canonical set (258,317 cells).
    pass_cols = ["step_01_umi_pass", "step_02_mad_pass", "step_03_doublet_pass",
                 "step_10b_manifold_pass", "step_15_pass"]

    if target == "full":
        cell_list = os.path.join(temp_dir, "cell_list_full.csv")
        df = pd.read_csv(central)
        missing_cols = [c for c in pass_cols if c not in df.columns]
        if missing_cols:
            _die(f"central_cell_status.csv missing pass columns: {missing_cols}",
                 hint="step 16 wrapper merges per-step pass flags into "
                      "central_cell_status. See L3_spec.md step_16.")
        cell_col = "cell_id" if "cell_id" in df.columns else df.columns[0]
        mask = df[pass_cols].astype(str).apply(
            lambda c: c.str.upper() == "TRUE").all(axis=1)
        keep = df.loc[mask, [cell_col]].rename(columns={cell_col: "cell_id"})
        keep.to_csv(cell_list, index=False)
        print(f"  full target: {len(keep)} cells -> {cell_list}")
        return cell_list

    # Compartment target: intersect cell_compartments with step_15_pass
    compartment_name = COMPARTMENT_NAME[target]
    comp_path = paths["preprocessing_cell_compartments"]
    if not os.path.exists(comp_path):
        _die(f"cell_compartments.csv not found: {comp_path}",
             hint="step 14 has not written compartment classifications yet.")

    cell_list = os.path.join(temp_dir, f"cell_list_{target}.csv")
    comp_df = pd.read_csv(comp_path)
    cell_col = "cell_id" if "cell_id" in comp_df.columns else comp_df.columns[0]
    # Expected column from step 14: `compartment` or `L0`
    label_col = None
    for candidate in ("compartment", "L0", "compartment_label"):
        if candidate in comp_df.columns:
            label_col = candidate
            break
    if label_col is None:
        _die(f"no compartment label column in {comp_path}",
             hint=f"columns: {list(comp_df.columns)}")
    comp_df = comp_df[comp_df[label_col] == compartment_name][[cell_col]]
    comp_df = comp_df.rename(columns={cell_col: "cell_id"})

    # Gate to canonical 5-way pass intersection (not just step_15).
    central_df = pd.read_csv(central)
    central_cell_col = "cell_id" if "cell_id" in central_df.columns else central_df.columns[0]
    central_df = central_df.rename(columns={central_cell_col: "cell_id"})
    missing_cols = [c for c in pass_cols if c not in central_df.columns]
    if missing_cols:
        _die(f"central_cell_status.csv missing pass columns: {missing_cols}")
    mask = central_df[pass_cols].astype(str).apply(
        lambda c: c.str.upper() == "TRUE").all(axis=1)
    passing = central_df.loc[mask, ["cell_id"]]

    merged = comp_df.merge(passing, on="cell_id", how="inner")
    if len(merged) == 0:
        _die(f"compartment {compartment_name} has 0 step_15-passing cells",
             hint="check upstream step 14 and step 15 outputs")
    merged.to_csv(cell_list, index=False)
    print(f"  {target} ({compartment_name}): {len(merged)} cells -> {cell_list}")
    return cell_list


def _target_output_dir(target, config_label, paths):
    """Resolve {sweep_root}/{full_object|compartments/<name>}/{config}/."""
    if target == "full":
        return os.path.join(paths["canonical_sweep_full_object"], config_label)
    compartment_name = COMPARTMENT_NAME[target]
    return os.path.join(paths["canonical_sweep_compartments"],
                        compartment_name, config_label)


def _merged_counts_path(target, paths):
    """Return expected path for the pre-assembled merged_counts.h5ad.

    Built by source/build_merged_counts.py via run/run_step_17_build_merged.sh.
    Must exist before the sweep array is submitted.
    """
    if target == "full":
        return os.path.join(paths["canonical_sweep_full_object"], "merged_counts.h5ad")
    compartment_name = COMPARTMENT_NAME[target]
    return os.path.join(paths["canonical_sweep_compartments"],
                        compartment_name, "merged_counts.h5ad")


def _invoke_integration_sweep(config_spec, target, cell_list, paths, output_dir,
                              test=False, dry_run=False):
    """Run source/integration_sweep.py with the resolved args for this config.

    Uses pre-assembled merged_counts.h5ad (built by build_merged_counts.py) so
    each task reads one H5AD instead of 62 per-sample CellRanger H5 files.
    The sweep array must NOT be submitted until merged_counts.h5ad exists for
    every target — run run_step_17_build_merged.sh first.
    """
    script = os.path.join(SOURCE_DIR, "integration_sweep.py")
    hp = config_spec.get("hyperparams", {}) or {}

    merged_path = _merged_counts_path(target, paths)
    if not os.path.exists(merged_path) and not dry_run:
        _die(f"merged_counts.h5ad not found for target '{target}': {merged_path}",
             hint="Run run_step_17_build_merged.sh for this target first.")

    cmd = [
        "python3", script,
        "--method", config_spec["method"],
        "--config-label", config_spec["label"],
        "--n-latent", str(config_spec["n_latent"]),
        "--vf-dir", paths["preprocessing_step16_vf_dir"],
        "--vf-subdir", "filtered_vfs",
        "--cell-list", cell_list,
        "--merged-input", merged_path,
        "--sample-list", paths["raw_data_sample_manifest"],  # required by argparse even with --merged-input
        "--batch-key", "patient_id",
        "--output-dir", output_dir,
        "--output-embedding",
    ]

    # scVI hyperparameters (hyperparams dict is empty for harmony/pca/concord)
    if "n_layers" in hp:
        cmd += ["--n-layers", str(hp["n_layers"])]
    if "dropout_rate" in hp:
        cmd += ["--dropout-rate", str(hp["dropout_rate"])]
    if "batch_size" in hp:
        cmd += ["--batch-size", str(hp["batch_size"])]
    if "early_stopping_patience" in hp:
        cmd += ["--early-stopping-patience", str(hp["early_stopping_patience"])]
    if "max_epochs" in hp:
        cmd += ["--max-epochs", str(hp["max_epochs"])]
    # Harmony hyperparameters
    if "theta" in hp:
        cmd += ["--theta", str(hp["theta"])]

    if test:
        cmd.append("--test")
    if dry_run:
        cmd.append("--dry-run")

    print(">>", " ".join(cmd))
    result = subprocess.run(cmd)
    if result.returncode != 0:
        _die(f"integration_sweep.py exited {result.returncode}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", required=True,
                        choices=TARGETS + ["all"],
                        help="Integration target: full atlas or a compartment")
    parser.add_argument("--config-idx", type=int, default=None,
                        help="Index into sweep_configs.yaml:configs. When set, "
                             "runs ONE config (for SLURM array mode). Omit to "
                             "run all configs sequentially (for interactive use).")
    parser.add_argument("--config-label", default=None,
                        help="Alternative to --config-idx: select by label.")
    parser.add_argument("--project-root", required=True,
                        help="Pipeline root (CFG_PROJECT_ROOT)")
    parser.add_argument("--config-path", default=DEFAULT_CONFIG_PATH)
    parser.add_argument("--temp-dir", default=None,
                        help="Directory for ephemeral cell-list CSVs. Defaults "
                             "to $TMPDIR or /tmp/step_17_wrapper.")
    parser.add_argument("--test", action="store_true",
                        help="Pass --test to integration_sweep.py (5000 cells, 10 epochs)")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    cfg = _load_yaml(args.config_path)
    paths = _resolve_paths(args.project_root)
    temp_dir = args.temp_dir or os.environ.get(
        "JOB_TEMP_DIR", "/tmp/step_17_wrapper")

    # Resolve which configs to run
    all_configs = cfg["configs"]
    if args.config_idx is not None and args.config_label is not None:
        _die("pass only one of --config-idx or --config-label")
    if args.config_idx is not None:
        configs_to_run = [_select_config(cfg, args.config_idx)]
    elif args.config_label is not None:
        matches = [c for c in all_configs if c["label"] == args.config_label]
        if not matches:
            _die(f"--config-label {args.config_label!r} not in sweep_configs.yaml")
        configs_to_run = matches
    else:
        configs_to_run = all_configs

    # Resolve which targets to run
    targets_to_run = TARGETS if args.target == "all" else [args.target]

    print("=" * 70)
    print("step_17 wrapper: integration sweep orchestrator")
    print(f"  targets:   {targets_to_run}")
    print(f"  configs:   {[c['label'] for c in configs_to_run]}")
    print(f"  project:   {args.project_root}")
    print(f"  temp_dir:  {temp_dir}")
    print("=" * 70)

    for target in targets_to_run:
        cell_list = _build_cell_list_for_target(
            target, paths, args.project_root, temp_dir)
        for config_spec in configs_to_run:
            output_dir = _target_output_dir(target, config_spec["label"], paths)
            os.makedirs(output_dir, exist_ok=True)
            print()
            print(f"--- target={target}  config={config_spec['label']} ---")
            print(f"    output: {output_dir}")
            _invoke_integration_sweep(
                config_spec, target, cell_list, paths, output_dir,
                test=args.test, dry_run=args.dry_run)

    print()
    print("=" * 70)
    print("wrapper.py done.")


if __name__ == "__main__":
    main()
