"""Write pipeline/outputs/manifest.json — provenance record for this pipeline run.

Counts rows in key output files and records pipeline version, model identifiers,
excluded samples, and output paths.
"""
import argparse, json, os
from datetime import date
import pandas as pd
import yaml


def count_rows(path):
    """Count data rows in a CSV (excludes header)."""
    try:
        return sum(1 for _ in open(path)) - 1
    except FileNotFoundError:
        return None


def load_paths(project_root, config_path):
    with open(config_path) as f:
        raw = yaml.safe_load(f)
    return {k: os.path.join(project_root, v) if isinstance(v, str) and not v.startswith("/") else v
            for k, v in raw.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    args = ap.parse_args()

    cfg = load_paths(args.project_root,
                     os.path.join(args.project_root, "config/paths.yaml"))

    ann_path = cfg["cell_annotations"]
    n_total = count_rows(ann_path)

    n_passes = None
    if n_total is not None:
        ann = pd.read_csv(ann_path, usecols=["platform", "passes_all"])
        xen = ann[ann["platform"] == "xenium"]
        n_passes = int(xen["passes_all"].sum())

    manifest = {
        "pipeline_version": "v1",
        "date": str(date.today()),
        "xenium_manifest": "config/xenium_manifest_clean.tsv",
        "n_samples": 65,
        "excluded_samples": [
            "UCI604_Tumor_xenium_1",
            "UCI604_Ipsilateral_xenium_1",
        ],
        "n_cells_joint": n_total,
        "n_xenium_passes_all": n_passes,
        "joint_model": "joint_nuc_n100_clean",
        "canonical_model": "joint_v6",
        "latent_dims": 100,
        "label_transfer": "majority-of-3 arbitration (kNN-joint primary)",
        "filter_version": "three_axis_filter_clean (proximity Q95, NMP rule, QC whole)",
        "outputs": {
            "joint_latent": "pipeline/outputs/embedding/joint_latent.csv",
            "joint_umap": "pipeline/outputs/embedding/joint_umap.csv",
            "joint_obs": "pipeline/outputs/embedding/joint_obs.csv",
            "cell_annotations": "pipeline/outputs/annotations/cell_annotations.csv",
            "l0p5_xenium": "pipeline/outputs/annotations/l0p5_xenium.csv",
            "three_axis_filter": "pipeline/outputs/annotations/three_axis_filter.csv",
        },
    }

    out_path = os.path.join(args.project_root, "pipeline/outputs/manifest.json")
    with open(out_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"Wrote manifest.json", flush=True)
    print(f"  n_cells_joint: {n_total:,}" if n_total else "  n_cells_joint: unknown",
          flush=True)
    print(f"  n_xenium_passes_all: {n_passes:,}" if n_passes else "", flush=True)


if __name__ == "__main__":
    main()
