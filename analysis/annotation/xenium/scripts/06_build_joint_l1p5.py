"""Build joint_l1p5.csv: cohort-wide cell-level table with L1.5 labels + metadata.

Joins per-compartment annotation outputs with cohort-level cell metadata and
library covariates. The output is the canonical L1.5 cohort table consumed by
the joint DA pipeline (DifferentialAbundance/...) and figure render scripts.

Output schema (19 columns):
    cell_id, platform, patient_id, library_id, p_level, quadrant, is_uoq,
    age, menopause, brca,
    compartment, l1p5_label, l1p5_short, l1p5_lineage, is_artifact,
    leiden_0.3, leiden_0.5, leiden_1.0, leiden_1.5

Inputs:
    --joint-obs           Cohort cell list with platform + compartment assignment
                            (pipeline/outputs/embedding/joint_obs.csv)
    --compartment-dir     Parent dir holding {epithelial,stromal,immune}/{cell_annotations,
                            leiden_assignments}.csv (pipeline/outputs/annotations/compartment)
    --library-covariates  Per-library metadata (library_id → p_level, quadrant, is_uoq,
                            age, menopause, brca). Currently lives under DifferentialAbundance/...
                            — flagged for promotion to a canonical pipeline output location.
    --out                 Output CSV path

Usage:
    python3 06_build_joint_l1p5.py \\
        --joint-obs           {paths.yaml: joint_obs} \\
        --compartment-dir     {paths.yaml: compartment_dir} \\
        --library-covariates  {DA/.../shared/library_covariates.csv} \\
        --out                 {paths.yaml: pipeline/outputs/annotations/joint_l1p5.csv}

Library_id derivation (from cell_id pattern, no external manifest needed):
    FLEX cell_id    "Pat1_P1_AAACCAATCCC...-1"          → library_id "Pat1_P1"
    Xenium cell_id  "Pat1_P1_xenium_1__aaalflil-1"      → library_id "Pat1_P1_xenium_1"

    Rule: if cell_id contains "__" (Xenium delimiter), library_id is the part before it;
    otherwise (FLEX), strip trailing "-N" suffix and drop the last underscore-segment
    (which is the 16-char barcode).
"""
import argparse
import os
import re
import sys

import pandas as pd

COMPARTMENTS = ["Epithelial", "Stromal", "Immune"]

LEIDEN_RESOLUTIONS = ["leiden_0.3", "leiden_0.5", "leiden_1.0", "leiden_1.5"]

OUTPUT_COLUMNS = [
    "cell_id", "platform", "patient_id", "library_id",
    "p_level", "quadrant", "is_uoq", "age", "menopause", "brca",
    "compartment", "l1p5_label", "l1p5_short", "l1p5_lineage", "is_artifact",
] + LEIDEN_RESOLUTIONS

LIBRARY_COVARIATE_COLS = [
    "library_id", "p_level", "quadrant", "is_uoq", "age", "menopause", "brca",
]


def derive_library_id(cell_id: str) -> str:
    # Xenium cells use "__" between library_id and barcode; FLEX uses single "_"
    # with a "-N" suffix on the 16-char barcode. Detect by presence of "__".
    if "__" in cell_id:
        return cell_id.split("__", 1)[0]
    stem = re.sub(r"-\d+$", "", cell_id)
    parts = stem.split("_")
    if len(parts) < 2:
        return cell_id
    return "_".join(parts[:-1])


def load_compartment(comp: str, comp_dir: str) -> pd.DataFrame:
    comp_lc = comp.lower()
    anno_path = os.path.join(comp_dir, comp_lc, "cell_annotations.csv")
    leiden_path = os.path.join(comp_dir, comp_lc, "leiden_assignments.csv")
    if not os.path.exists(anno_path):
        sys.exit(f"ERROR: missing {anno_path}")
    if not os.path.exists(leiden_path):
        sys.exit(f"ERROR: missing {leiden_path}")
    anno = pd.read_csv(anno_path, usecols=["cell_id", "label", "label_short", "lineage", "is_artifact"])
    anno = anno.rename(columns={"label": "l1p5_label", "label_short": "l1p5_short", "lineage": "l1p5_lineage"})
    leiden = pd.read_csv(leiden_path)
    have = [c for c in LEIDEN_RESOLUTIONS if c in leiden.columns]
    if not have:
        print(f"  {comp}: WARNING — no target leiden resolutions present in {leiden_path}", flush=True)
    leiden = leiden[["cell_id"] + have]
    df = anno.merge(leiden, on="cell_id", how="left")
    for col in LEIDEN_RESOLUTIONS:
        if col not in df.columns:
            df[col] = pd.NA
    df["compartment"] = comp
    print(f"  {comp}: {len(df):,} cells, leiden cols available: {have}", flush=True)
    return df


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--joint-obs", required=True)
    ap.add_argument("--compartment-dir", required=True)
    ap.add_argument("--library-covariates", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    print(f"Loading joint_obs: {args.joint_obs}", flush=True)
    obs = pd.read_csv(args.joint_obs)
    n_obs = len(obs)
    print(f"  {n_obs:,} cells", flush=True)
    needed_obs = {"cell_id", "patient_id", "platform", "compartment"}
    missing = needed_obs - set(obs.columns)
    if missing:
        sys.exit(f"ERROR: joint_obs missing columns: {missing}")

    print("Deriving library_id from cell_id pattern", flush=True)
    obs["library_id"] = [derive_library_id(c) for c in obs["cell_id"]]
    # Drop obs.compartment — unreliable upstream (Xenium cells get "unknown" since
    # compartment assignment happens later in phase 09). Use the per-compartment
    # subfolder source as the canonical compartment for each cell.
    if "compartment" in obs.columns:
        obs = obs.drop(columns=["compartment"])

    print(f"Loading library_covariates: {args.library_covariates}", flush=True)
    covs = pd.read_csv(args.library_covariates, usecols=LIBRARY_COVARIATE_COLS)
    print(f"  {len(covs):,} libraries", flush=True)
    n_libs_matched = obs["library_id"].isin(covs["library_id"]).sum()
    print(f"  obs library coverage: {n_libs_matched:,}/{n_obs:,} cells "
          f"({obs['library_id'].nunique()} unique library_ids)", flush=True)
    if n_libs_matched < n_obs:
        unmatched = sorted(set(obs.loc[~obs['library_id'].isin(covs['library_id']), 'library_id']))
        print(f"  WARNING: {len(unmatched)} library_ids in obs not in covariates: {unmatched[:5]}{'...' if len(unmatched)>5 else ''}", flush=True)

    obs = obs.merge(covs, on="library_id", how="left", validate="m:1")

    print("Loading per-compartment annotations + leiden", flush=True)
    label_parts = [load_compartment(c, args.compartment_dir) for c in COMPARTMENTS]
    labels = pd.concat(label_parts, ignore_index=True)
    print(f"  {len(labels):,} labeled cells across compartments", flush=True)

    n_labels_in_obs = obs["cell_id"].isin(labels["cell_id"]).sum()
    print(f"  obs label coverage: {n_labels_in_obs:,}/{n_obs:,} cells", flush=True)

    print("Joining cohort metadata + labels", flush=True)
    out = obs.merge(labels, on="cell_id", how="left", validate="1:1")

    missing_cols = [c for c in OUTPUT_COLUMNS if c not in out.columns]
    if missing_cols:
        sys.exit(f"ERROR: output missing columns: {missing_cols}")
    out = out[OUTPUT_COLUMNS]

    n_unlabeled = out["l1p5_label"].isna().sum()
    if n_unlabeled:
        print(f"  WARNING: {n_unlabeled:,} cells with no l1p5_label", flush=True)

    print(f"Writing {args.out}", flush=True)
    out.to_csv(args.out, index=False)
    print(f"Done: {len(out):,} rows × {len(out.columns)} cols", flush=True)


if __name__ == "__main__":
    main()
