#!/usr/bin/env python3
"""Build level1.5_annotation -> CL term mapping table from integrated h5ad.

Reads all-breast-cells.h5ad in backed mode, extracts every unique
(level1.5_annotation, cell_type_ontology_term_id) pair with cell counts
and per-label distribution percentages.

Output: publication/mappings/level15_to_cl_mapping_raw.csv
  Columns: level15_annotation, cl_term, n_cells, pct_of_label, total_label_cells

Review the raw output to produce the final reviewed mapping:
  publication/mappings/level15_to_cl_mapping.csv
  Columns: level15_annotation, assigned_cl_term, decision_rationale

Usage:
    python build_level15_mapping.py --repo-root /path/to/iHBCAv1_upload

"""

import argparse
import sys
from pathlib import Path

import anndata as ad


def build_mapping(repo_root):
    """Extract level1.5 -> CL term distribution from integrated h5ad."""
    h5ad_path = repo_root / "publication/outputs/integrated_objects/all-breast-cells.h5ad"
    print(f"Loading {h5ad_path} (backed mode)...")

    adata = ad.read_h5ad(h5ad_path, backed="r")
    print(f"  Shape: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # Extract the two columns we need
    obs_cols = list(adata.obs.columns)
    print(f"  obs columns ({len(obs_cols)}): {obs_cols[:20]}...")

    if "level1.5_annotation" not in adata.obs.columns:
        # Try alternate names
        for alt in ["level1.5", "level15_annotation", "Level1.5_annotation"]:
            if alt in adata.obs.columns:
                print(f"  Using alternate column: {alt}")
                l15_col = alt
                break
        else:
            print("ERROR: No level1.5_annotation column found!")
            print(f"  Available: {obs_cols}")
            sys.exit(1)
    else:
        l15_col = "level1.5_annotation"

    if "cell_type_ontology_term_id" not in adata.obs.columns:
        print("ERROR: cell_type_ontology_term_id not in obs!")
        sys.exit(1)

    # Read just the two columns into memory
    print("  Reading obs columns into memory...")
    df = adata.obs[[l15_col, "cell_type_ontology_term_id"]].copy()
    df.columns = ["level15_annotation", "cl_term"]
    print(f"  Loaded {len(df):,} cells")

    # Group by (level15_annotation, cl_term) and count
    print("  Computing per-label CL term distributions...")
    grouped = (
        df.groupby(["level15_annotation", "cl_term"])
        .size()
        .reset_index(name="n_cells")
    )

    # Add per-label totals and percentages
    label_totals = df.groupby("level15_annotation").size().reset_index(name="total_label_cells")
    grouped = grouped.merge(label_totals, on="level15_annotation")
    grouped["pct_of_label"] = (grouped["n_cells"] / grouped["total_label_cells"] * 100).round(2)

    # Sort by label, then by descending cell count within label
    grouped = grouped.sort_values(
        ["level15_annotation", "n_cells"],
        ascending=[True, False]
    )

    return grouped


def print_summary(df):
    """Print human-readable summary of the mapping."""
    labels = df["level15_annotation"].unique()
    print(f"\n{'='*70}")
    print(f"LEVEL1.5 -> CL TERM MAPPING SUMMARY")
    print(f"{'='*70}")
    print(f"  Total unique labels: {len(labels)}")
    print(f"  Total cells: {df['n_cells'].sum():,}")
    print()

    # Per-label summary
    for label in sorted(labels):
        label_df = df[df["level15_annotation"] == label]
        total = label_df["total_label_cells"].iloc[0]
        n_terms = len(label_df)
        top_term = label_df.iloc[0]

        # Classification
        if n_terms == 1 and top_term["cl_term"] == "unknown":
            category = "ALL_UNKNOWN"
        elif top_term["pct_of_label"] >= 95:
            category = "CLEAR"
        elif top_term["pct_of_label"] >= 80:
            category = "DOMINANT"
        else:
            category = "SPLIT"

        print(f"  {label} ({total:,} cells) [{category}]")
        for _, row in label_df.iterrows():
            print(f"    {row['cl_term']}: {row['n_cells']:,} ({row['pct_of_label']:.1f}%)")
        print()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root", required=True, type=Path,
        help="Path to iHBCAv1_upload root"
    )
    args = parser.parse_args()

    repo_root = args.repo_root.resolve()
    if not (repo_root / "publication").exists():
        print(f"ERROR: {repo_root} doesn't look like iHBCAv1_upload root")
        sys.exit(1)

    # Build raw mapping
    mapping_df = build_mapping(repo_root)

    # Print summary
    print_summary(mapping_df)

    # Write raw output
    out_path = repo_root / "publication/mappings/level15_to_cl_mapping_raw.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    mapping_df.to_csv(out_path, index=False)
    print(f"Written: {out_path} ({len(mapping_df)} rows)")

    # Count labels needing review
    labels = mapping_df["level15_annotation"].unique()
    unknown_only = 0
    clear = 0
    needs_review = 0
    for label in labels:
        ldf = mapping_df[mapping_df["level15_annotation"] == label]
        top = ldf.iloc[0]
        if len(ldf) == 1 and top["cl_term"] == "unknown":
            unknown_only += 1
        elif top["pct_of_label"] >= 95:
            clear += 1
        else:
            needs_review += 1

    print(f"\nReview summary:")
    print(f"  CLEAR (>= 95% dominant): {clear} labels")
    print(f"  ALL_UNKNOWN (unmappable): {unknown_only} labels")
    print(f"  NEEDS REVIEW (split): {needs_review} labels")
    print(f"\nNext: review raw CSV and produce level15_to_cl_mapping.csv")


if __name__ == "__main__":
    main()
