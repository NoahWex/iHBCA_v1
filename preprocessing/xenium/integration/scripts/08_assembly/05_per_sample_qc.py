"""Per-sample QC summary: cell counts and pass rates from three_axis_filter.

Reads cell_annotations.csv (Xenium rows only) and xenium_manifest_clean.tsv.
Aggregates per xenium_id: total cells, passes_all count, pass rate.

Outputs:
  previews/per_sample_qc.csv   — tabular summary
  previews/per_sample_qc.png/.pdf — bar chart (n_cells + pass rate)
"""
import argparse, os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import yaml


def load_paths(project_root, config_path):
    with open(config_path) as f:
        return {k: os.path.join(project_root, v)
                for k, v in yaml.safe_load(f).items()
                if isinstance(v, str) and not v.startswith("/")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    args = ap.parse_args()

    cfg = load_paths(args.project_root,
                     os.path.join(args.project_root, "config/paths.yaml"))
    os.makedirs(cfg["previews_dir"], exist_ok=True)

    print("Loading cell_annotations (Xenium rows)...", flush=True)
    ann = pd.read_csv(cfg["cell_annotations"], index_col="cell_id",
                      usecols=["cell_id", "platform", "patient_id",
                                "position", "passes_all"])
    xen = ann[ann["platform"] == "xenium"].copy()
    xen["xenium_id"] = xen["patient_id"] + "_" + xen["position"] + "_xenium"
    print(f"  {len(xen):,} Xenium cells", flush=True)

    qc = (xen.groupby("xenium_id")
             .agg(n_cells=("passes_all", "count"),
                  n_pass=("passes_all", "sum"))
             .reset_index())
    qc["pass_rate"] = qc["n_pass"] / qc["n_cells"]
    qc = qc.sort_values("n_cells", ascending=False)

    csv_path = os.path.join(cfg["previews_dir"], "per_sample_qc.csv")
    qc.to_csv(csv_path, index=False)
    print(f"  Wrote {len(qc)} samples → per_sample_qc.csv", flush=True)

    # Bar charts
    fig, axes = plt.subplots(2, 1, figsize=(max(12, len(qc) * 0.25), 8))
    labels = qc["xenium_id"].str.replace("_xenium", "", regex=False)

    axes[0].bar(range(len(qc)), qc["n_cells"], color="#1f77b4", width=0.8)
    axes[0].set_ylabel("n cells")
    axes[0].set_title("Cells per sample", fontsize=8)
    axes[0].set_xticks([])

    axes[1].bar(range(len(qc)), qc["pass_rate"] * 100, color="#2ca02c", width=0.8)
    axes[1].axhline(50, color="red", linewidth=0.5, linestyle="--")
    axes[1].set_ylim(0, 100)
    axes[1].set_ylabel("% passes_all")
    axes[1].set_title("Pass rate per sample", fontsize=8)
    axes[1].set_xticks(range(len(qc)))
    axes[1].set_xticklabels(labels, rotation=90, fontsize=4)

    fig.tight_layout()
    base = os.path.join(cfg["previews_dir"], "per_sample_qc")
    fig.savefig(base + ".pdf", dpi=300, bbox_inches="tight")
    fig.savefig(base + ".png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  Wrote per_sample_qc.png/.pdf", flush=True)


if __name__ == "__main__":
    main()
