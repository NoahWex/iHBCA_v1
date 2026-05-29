#!/usr/bin/env python3
"""Plot X_umap_ihbca_scvi_100 for enriched source h5ads."""
import argparse
from pathlib import Path
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

STUDY_TO_FILENAME = {
    "gray": "gray2022.h5ad", "kumar": "kumar2023.h5ad",
    "murrow": "murrow2022.h5ad", "nee": "nee2023.h5ad",
    "twigger": "twigger2022.h5ad", "reed": "reed2024.h5ad",
    "pal": "pal2021.h5ad",
}
CELL_TYPE_COLS = ["cell_type", "cell_type_label", "level1.5", "level1", "cell_type_ontology_term_id"]

def find_cell_type_col(obs):
    for col in CELL_TYPE_COLS:
        if col in obs.columns:
            return col
    return None

def plot_study_umap(adata, study, output_dir):
    umap_key = "X_umap_ihbca_scvi_100"
    if umap_key not in adata.obsm:
        print(f"  WARNING: {umap_key} not found in {study}")
        return
    coords = adata.obsm[umap_key]
    valid = ~np.isnan(coords).any(axis=1)
    n_valid, n_total = valid.sum(), len(valid)
    coords_valid = coords[valid]
    ct_col = find_cell_type_col(adata.obs)
    if ct_col:
        labels = adata.obs[ct_col].values[valid]
        label_name = ct_col
    else:
        labels = np.full(n_valid, "unknown")
        label_name = "N/A"
    unique_labels = sorted(set(str(l) for l in labels))
    n_colors = len(unique_labels)
    cmap = plt.cm.get_cmap("tab20", max(n_colors, 1)) if n_colors <= 20 else plt.cm.get_cmap("gist_ncar", n_colors)
    label_to_color = {l: cmap(i) for i, l in enumerate(unique_labels)}
    fig, axes = plt.subplots(1, 2, figsize=(20, 8))
    ax = axes[0]
    colors = [label_to_color[str(l)] for l in labels]
    ax.scatter(coords_valid[:, 0], coords_valid[:, 1], c=colors, s=0.5, alpha=0.5, rasterized=True)
    ax.set_title(f"{study} - X_umap_ihbca_scvi_100\n({n_valid:,}/{n_total:,} cells, color={label_name})")
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2"); ax.set_aspect("equal")
    ax2 = axes[1]
    if "X_umap" in adata.obsm:
        orig = adata.obsm["X_umap"]
        orig_valid = ~np.isnan(orig).any(axis=1)
        orig_coords = orig[orig_valid]
        if ct_col:
            orig_labels = adata.obs[ct_col].values[orig_valid]
            orig_colors = [label_to_color.get(str(l), (0.5,0.5,0.5,1)) for l in orig_labels]
        else:
            orig_colors = "gray"
        ax2.scatter(orig_coords[:, 0], orig_coords[:, 1], c=orig_colors, s=0.5, alpha=0.5, rasterized=True)
        ax2.set_title(f"{study} - X_umap (original)\n({orig_valid.sum():,} cells)")
    else:
        ax2.text(0.5, 0.5, "No X_umap", ha="center", va="center", transform=ax2.transAxes)
    ax2.set_xlabel("UMAP 1"); ax2.set_ylabel("UMAP 2"); ax2.set_aspect("equal")
    if n_colors <= 30:
        from matplotlib.lines import Line2D
        handles = [Line2D([0],[0],marker="o",color="w",markerfacecolor=label_to_color[l],markersize=6,label=l) for l in unique_labels]
        fig.legend(handles=handles, loc="center right", bbox_to_anchor=(1.15,0.5), fontsize=7, title=label_name, title_fontsize=8)
    plt.tight_layout()
    out_path = output_dir / f"umap_enriched_{study}.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", required=True)
    parser.add_argument("--study", help="Single study (default: all)")
    args = parser.parse_args()
    repo_root = Path(args.repo_root)
    source_dir = repo_root / "publication/outputs/source_datasets"
    output_dir = repo_root / "publication/outputs/figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    studies = [args.study] if args.study else list(STUDY_TO_FILENAME.keys())
    for study in studies:
        h5ad_path = source_dir / STUDY_TO_FILENAME[study]
        if not h5ad_path.exists():
            print(f"Skipping {study}: not found"); continue
        print(f"\nPlotting {study}...")
        adata = ad.read_h5ad(h5ad_path)
        if "X_umap_ihbca_scvi_100" not in adata.obsm:
            print(f"  Not enriched yet, skipping"); continue
        plot_study_umap(adata, study, output_dir)
        del adata
    print(f"\nAll plots saved to: {output_dir}")

if __name__ == "__main__":
    main()
