#!/usr/bin/env python3
"""
step_17_integration_sweep / build_merged_counts.py

Pre-assemble a per-target merged raw count matrix so that each sweep task
reads a single H5AD rather than re-loading 62 CellRanger H5 files.

Run ONCE per target via run/run_step_17_build_merged.sh BEFORE submitting
the sweep array (run/run_step_17_sweep.sh).

Outputs:
  {canonical_sweep_compartments}/{Compartment}/merged_counts.h5ad  (epi/str/imm)
  {canonical_sweep_full_object}/merged_counts.h5ad                 (full)

The saved H5AD contains raw counts (X) for the full gene space across all
passing cells in the target. VF subsetting is left to integration_sweep.py
per-config (fast — it is a column selection, not a re-load).

provenance: NEW. Wraps integration_sweep.load_per_sample_data with the same
cell-selection logic as wrapper._build_cell_list_for_target.
"""

import argparse
import os
import sys

import pandas as pd

# load_per_sample_data lives in the same source/ directory.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from integration_sweep import load_per_sample_data  # noqa: E402

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


def _parse_paths_yaml(paths_yaml_path):
    """Parse config/paths.yaml with ${ref} expansion.

    provenance: wrapper._parse_paths_yaml (identical logic, co-located copy
    to avoid cross-directory import of wrapper.py).
    Iterates until fixed point (up to 10 passes) to resolve chains like
    canonical_sweep_compartments → preprocessing_step_17 → preprocessing_root.
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
            if line[0].isspace():
                continue
            key, _, val = line.partition(":")
            raw[key.strip()] = val.strip()

    expanded = dict(raw)
    changed = True
    max_iter = 10
    while changed and max_iter > 0:
        changed = False
        for k, v in list(expanded.items()):
            if "${" in v:
                new_v = v
                for ref_key, ref_val in expanded.items():
                    token = "${" + ref_key + "}"
                    if token in new_v and "${" not in ref_val:
                        new_v = new_v.replace(token, ref_val)
                if new_v != v:
                    expanded[k] = new_v
                    changed = True
        max_iter -= 1
    return expanded


def _build_cell_list(target, paths, temp_dir):
    """Derive per-target passing cell list.

    provenance: wrapper._build_cell_list_for_target (same logic).
    """
    central = paths["preprocessing_central_cell_status"]
    if not os.path.exists(central):
        _die(f"central_cell_status.csv not found: {central}",
             hint="Run step 16 first (step_15_pass column must be present).")

    os.makedirs(temp_dir, exist_ok=True)
    cell_list_path = os.path.join(temp_dir, f"cell_list_{target}.csv")

    central_df = pd.read_csv(central)
    cell_col = "cell_id" if "cell_id" in central_df.columns else central_df.columns[0]
    central_df = central_df.rename(columns={cell_col: "cell_id"})
    # Canonical 5-way filter: cell must pass every QC gate, not just the
    # final step. step_15_pass alone admits ~321K cells (includes earlier-
    # filtered cells re-marked True at step 15); intersection of all 5
    # pass columns yields the true canonical set (258,317 cells).
    pass_cols = ["step_01_umi_pass", "step_02_mad_pass", "step_03_doublet_pass",
                 "step_10b_manifold_pass", "step_15_pass"]
    missing_cols = [c for c in pass_cols if c not in central_df.columns]
    if missing_cols:
        _die(f"central_cell_status.csv missing pass columns: {missing_cols}")
    mask = central_df[pass_cols].astype(str).apply(
        lambda col: col.str.upper() == "TRUE").all(axis=1)
    passing = central_df.loc[mask, ["cell_id"]]

    if target == "full":
        passing.to_csv(cell_list_path, index=False)
        print(f"  full: {len(passing)} cells")
        return cell_list_path

    compartment_name = COMPARTMENT_NAME[target]
    comp_path = paths["preprocessing_cell_compartments"]
    if not os.path.exists(comp_path):
        _die(f"cell_compartments.csv not found: {comp_path}",
             hint="Run step 14 first.")

    comp_df = pd.read_csv(comp_path)
    comp_cell_col = "cell_id" if "cell_id" in comp_df.columns else comp_df.columns[0]
    label_col = next(
        (c for c in ("compartment", "L0", "compartment_label") if c in comp_df.columns),
        None,
    )
    if label_col is None:
        _die(f"no compartment label column in {comp_path}",
             hint=f"columns: {list(comp_df.columns)}")

    comp_df = comp_df[comp_df[label_col] == compartment_name][[comp_cell_col]].rename(
        columns={comp_cell_col: "cell_id"})
    merged = comp_df.merge(passing, on="cell_id", how="inner")
    if len(merged) == 0:
        _die(f"compartment {compartment_name} has 0 step_15-passing cells",
             hint="Check step 14 and step 15 outputs.")
    merged.to_csv(cell_list_path, index=False)
    print(f"  {target} ({compartment_name}): {len(merged)} cells")
    return cell_list_path


def main():
    parser = argparse.ArgumentParser(
        description="Pre-assemble merged raw counts for one integration sweep target")
    parser.add_argument("--target", required=True, choices=["full", "epi", "str", "imm"],
                        help="Integration target (full atlas or one compartment)")
    parser.add_argument("--project-root", required=True,
                        help="Pipeline root (CFG_PROJECT_ROOT)")
    parser.add_argument("--temp-dir", default=None,
                        help="Scratch dir for cell-list CSVs. Defaults to "
                             "$JOB_TEMP_DIR or /tmp/step_17_build_merged")
    parser.add_argument("--dry-run", action="store_true",
                        help="Validate inputs, print plan, exit 0 before H5 loading")
    args = parser.parse_args()

    paths_yaml = os.path.join(args.project_root, "config", "paths.yaml")
    paths = _parse_paths_yaml(paths_yaml)

    needed = [
        "preprocessing_central_cell_status",
        "preprocessing_cell_compartments",
        "raw_data_sample_manifest",
        "preprocessing_step16_vf_dir",
        "canonical_sweep_full_object",
        "canonical_sweep_compartments",
    ]
    missing = [k for k in needed if k not in paths]
    if missing:
        _die(f"paths.yaml missing keys: {missing}")

    temp_dir = args.temp_dir or os.environ.get(
        "JOB_TEMP_DIR", "/tmp/step_17_build_merged")

    # Resolve output path
    if args.target == "full":
        out_dir = paths["canonical_sweep_full_object"]
    else:
        compartment_name = COMPARTMENT_NAME[args.target]
        out_dir = os.path.join(paths["canonical_sweep_compartments"], compartment_name)
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "merged_counts.h5ad")

    print("=" * 70)
    print("step_17 build_merged_counts")
    print(f"  target:  {args.target}")
    print(f"  output:  {out_path}")
    print(f"  manifest:{paths['raw_data_sample_manifest']}")
    print("=" * 70)

    if args.dry_run:
        if not os.path.exists(paths["preprocessing_central_cell_status"]):
            _die(f"central_cell_status.csv not found: "
                 f"{paths['preprocessing_central_cell_status']}")
        if not os.path.exists(paths["raw_data_sample_manifest"]):
            _die(f"sample manifest not found: {paths['raw_data_sample_manifest']}")
        print("[DRY RUN] inputs validated — skipping H5 loads")
        return

    if os.path.exists(out_path):
        print(f"Overwriting existing: {out_path}")

    cell_list_path = _build_cell_list(args.target, paths, temp_dir)

    print("Loading per-sample H5 files...")
    merged, _ = load_per_sample_data(
        sample_list_path=paths["raw_data_sample_manifest"],
        cell_list_path=cell_list_path,
        vf_dir=paths["preprocessing_step16_vf_dir"],
        vf_subdir="filtered_vfs",
    )
    # merged contains raw counts for all genes (full gene space).
    # VF subsetting is intentionally deferred to integration_sweep.py so
    # each config can use the same merged H5AD without a gene-space dependency.

    print(f"Saving: {out_path}")
    merged.write_h5ad(out_path)
    print(f"Done: {merged.n_obs} cells x {merged.n_vars} genes -> {out_path}")


if __name__ == "__main__":
    main()
