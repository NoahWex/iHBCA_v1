"""Assemble cell_annotations.csv — the master per-cell table for the Xenium pipeline.

Joint scope: all FLEX + Xenium cells from the joint_v6 embedding.
FLEX cells have NMP, filter, and L0.5 columns set to NA (Xenium-only outputs).

Inputs (from pipeline/outputs/):
  embedding/joint_obs.csv      — cell_id, patient_id, platform, compartment
  embedding/joint_umap.csv     — cell_id, UMAP1, UMAP2
  annotations/l0p5_xenium.csv  — Xenium label arbitration (kNN primary)
  annotations/three_axis_filter.csv — per-axis pass/fail + NMP scores

Output:
  annotations/cell_annotations.csv
"""
import argparse, os
import pandas as pd
import numpy as np
import yaml

PATIENT_META = {
    "Pat1":      {"age": 72, "menopause": "Post",  "brca": "negative"},
    "Pat2":      {"age": 45, "menopause": "Pre",   "brca": "negative"},
    "UCI604":    {"age": 37, "menopause": "Pre",   "brca": "negative"},
    "UCI220228": {"age": 39, "menopause": "Pre",   "brca": "BRCA1"},
}

KNOWN_QUADRANTS = {"UOQ", "UIQ", "LOQ", "LIQ"}


def parse_position(cell_id, patient_id):
    """Extract position string from Xenium cell_id."""
    rest = cell_id[len(patient_id) + 1:]
    parts = rest.split("_xenium_")
    return parts[0] if len(parts) >= 2 else "unknown"


def extract_position_cols(obs):
    positions = [
        parse_position(cid, pid) if "_xenium_" in cid else "FLEX"
        for cid, pid in zip(obs.index, obs["patient_id"].values)
    ]
    obs["position"] = positions
    obs["p_level"] = obs["position"].str.split("_").str[0]
    obs["quadrant"] = obs["position"].str.split("_").str[-1].where(
        obs["position"].str.split("_").str[-1].isin(KNOWN_QUADRANTS), other=None
    )
    obs["is_uoq"] = obs["quadrant"] == "UOQ"
    return obs


def load_paths(project_root, config_path):
    with open(config_path) as f:
        return {k: os.path.join(project_root, v)
                for k, v in yaml.safe_load(f).items()
                if isinstance(v, str) and not v.startswith("/")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True,
                    help="Sub-track root holding config/paths.yaml.")
    ap.add_argument("--test", action="store_true",
                    help="Subset to first 5000 rows for validation")
    args = ap.parse_args()

    cfg = load_paths(args.project_root,
                     os.path.join(args.project_root, "config/paths.yaml"))

    print("Loading joint_obs...", flush=True)
    obs = pd.read_csv(cfg["joint_obs"]).set_index("cell_id")
    if args.test:
        obs = obs.iloc[:5000]
    print(f"  {len(obs):,} cells (platform: {obs['platform'].value_counts().to_dict()})",
          flush=True)

    print("Loading UMAP...", flush=True)
    umap = pd.read_csv(cfg["joint_umap"]).set_index("cell_id").reindex(obs.index)

    print("Loading L0.5 labels...", flush=True)
    l0 = pd.read_csv(cfg["l0p5_xenium"],
                     usecols=["cell_id", "l0p5_final", "vote_count",
                               "method", "knn_confidence"]).set_index("cell_id")

    print("Loading three-axis filter...", flush=True)
    flt = pd.read_csv(cfg["three_axis_filter"],
                      usecols=["cell_id", "passes_all", "passes_proximity",
                                "passes_nmp", "passes_qc_whole",
                                "nmp_nuclear", "nmp_cyto"]).set_index("cell_id")

    # Join — FLEX cells get NA for Xenium-only columns
    obs["l0p5"]         = l0["l0p5_final"].reindex(obs.index)
    obs["vote_count"]   = l0["vote_count"].reindex(obs.index)
    obs["label_method"] = l0["method"].reindex(obs.index)
    obs["knn_confidence"] = l0["knn_confidence"].reindex(obs.index)

    obs["passes_all"]       = flt["passes_all"].reindex(obs.index)
    obs["passes_proximity"] = flt["passes_proximity"].reindex(obs.index)
    obs["passes_nmp"]       = flt["passes_nmp"].reindex(obs.index)
    obs["passes_qc_whole"]  = flt["passes_qc_whole"].reindex(obs.index)
    obs["nmp_nuclear"]      = flt["nmp_nuclear"].reindex(obs.index)
    obs["nmp_cyto"]         = flt["nmp_cyto"].reindex(obs.index)
    obs["nmp_exclude_flag"] = (obs["nmp_nuclear"] > obs["nmp_cyto"])

    # Position metadata parsed from cell_id
    obs = extract_position_cols(obs)

    # Patient-level clinical metadata
    obs["age"]       = obs["patient_id"].map({p: m["age"]       for p, m in PATIENT_META.items()})
    obs["menopause"] = obs["patient_id"].map({p: m["menopause"] for p, m in PATIENT_META.items()})
    obs["brca"]      = obs["patient_id"].map({p: m["brca"]      for p, m in PATIENT_META.items()})

    # UMAP coordinates
    obs["UMAP1"] = umap["UMAP1"]
    obs["UMAP2"] = umap["UMAP2"]

    col_order = [
        "patient_id", "platform", "compartment",
        "l0p5", "vote_count", "label_method", "knn_confidence",
        "passes_all", "passes_proximity", "passes_nmp", "passes_qc_whole",
        "nmp_nuclear", "nmp_cyto", "nmp_exclude_flag",
        "position", "p_level", "quadrant", "is_uoq",
        "age", "menopause", "brca",
        "UMAP1", "UMAP2",
    ]
    obs = obs[[c for c in col_order if c in obs.columns]]

    out_path = cfg["cell_annotations"]
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    obs.to_csv(out_path)
    print(f"Wrote {len(obs):,} rows → {out_path}", flush=True)

    xen = obs[obs["platform"] == "xenium"]
    passes = xen["passes_all"].sum()
    print(f"  Xenium passes_all: {int(passes):,} / {len(xen):,} "
          f"({100*passes/len(xen):.1f}%)", flush=True)


if __name__ == "__main__":
    main()
