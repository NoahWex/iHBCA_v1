#!/usr/bin/env python3
"""Build CL term crosswalk CSVs for non-CxG source datasets.

Phases 1+2 of cl_term_backfill plan:
  Phase 1: Investigate cell ID alignment between integrated and source h5ads
  Phase 2: Build crosswalk CSVs mapping source_cell_id -> cell_type_ontology_term_id

Usage:
    python build_cl_crosswalk.py --repo-root /path/to/iHBCAv1_upload [--investigate-only]
"""

import argparse
import sys
from pathlib import Path

import anndata as ad
import pandas as pd


# ---------------------------------------------------------------------------
# Studies to process
# ---------------------------------------------------------------------------

TARGET_STUDIES = ["murrow", "nee", "pal"]

# Source h5ad filenames (must match source_dataset_registry.yaml)
SOURCE_H5AD = {
    "murrow": "murrow2022.h5ad",
    "nee": "nee2023.h5ad",
    "pal": "pal2021.h5ad",
}

# Cell ID mapping file locations (relative to repo_root)
# murrow and nee have single mapping files; pal has per-sub-study
MAPPING_PATHS = {
    "murrow": ["external_studies/outputs/murrow/published/cell_id_mapping.csv"],
    "nee": ["external_studies/outputs/nee/published/cell_id_mapping.csv"],
    "pal": [
        "external_studies/outputs/pal_norm_epi/published/cell_id_mapping.csv",
        "external_studies/outputs/pal_norm_total/published/cell_id_mapping.csv",
        "external_studies/outputs/pal_norm_b1/published/cell_id_mapping.csv",
    ],
}


# ---------------------------------------------------------------------------
# Phase 1: Investigate cell ID alignment
# ---------------------------------------------------------------------------


def load_level15_mapping(repo_root):
    """Load reviewed level1.5 -> CL term mapping if available.

    Returns: dict {level15_annotation: assigned_cl_term} or empty dict.
    """
    path = repo_root / "publication/mappings/level15_to_cl_mapping.csv"
    if not path.exists():
        return {}
    df = pd.read_csv(path)
    if "level15_annotation" not in df.columns or "assigned_cl_term" not in df.columns:
        print(f"  WARNING: level15_to_cl_mapping.csv missing required columns")
        return {}
    mapping = dict(zip(df["level15_annotation"], df["assigned_cl_term"]))
    print(f"  Loaded level1.5 -> CL mapping: {len(mapping)} labels")
    return mapping


def load_integrated_cl_terms(repo_root):
    """Load cell_type_ontology_term_id from integrated h5ad (backed mode).

    If a reviewed level15_to_cl_mapping.csv exists, resolves "unknown" CL terms
    via level1.5_annotation fallback before returning.

    Returns: dict {obs_name: cl_term} for all 2.12M cells.
    """
    path = repo_root / "publication/outputs/integrated_objects/all-breast-cells.h5ad"
    print(f"\n{'='*70}")
    print(f"Loading integrated h5ad (backed mode): {path}")
    print(f"{'='*70}")

    adata = ad.read_h5ad(path, backed="r")
    print(f"  Shape: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")
    print(f"  obs columns: {list(adata.obs.columns)[:20]}...")

    # Extract cell_type_ontology_term_id
    if "cell_type_ontology_term_id" not in adata.obs.columns:
        print("  ERROR: cell_type_ontology_term_id not in obs!")
        sys.exit(1)

    cl_series = adata.obs["cell_type_ontology_term_id"]
    obs_names = list(adata.obs_names)

    print(f"  First 5 obs_names: {obs_names[:5]}")
    print(f"  CL term distribution (top 20):")
    vc = cl_series.value_counts()
    for term, count in vc.head(20).items():
        print(f"    {term}: {count:,}")
    print(f"  Total unique CL terms: {vc.shape[0]}")

    # Check for study-identifying column
    for col in ["dataset_id", "collection_id", "study", "Study", "dataset"]:
        if col in adata.obs.columns:
            print(f"\n  Study column found: '{col}'")
            print(f"  Values: {adata.obs[col].value_counts().to_dict()}")
            break
    else:
        print("\n  No obvious study column found in obs.")

    # Build lookup dict — convert to str to allow new CL terms from level1.5 mapping
    cl_values = cl_series.astype(str).values.copy()
    n_unknown_before = (cl_values == "unknown").sum()

    # Apply level1.5 -> CL fallback for "unknown" cells
    l15_map = load_level15_mapping(repo_root)
    if l15_map:
        # Find level1.5_annotation column
        l15_col = None
        for candidate in ["level1.5_annotation", "level1.5", "level15_annotation"]:
            if candidate in adata.obs.columns:
                l15_col = candidate
                break

        if l15_col:
            l15_series = adata.obs[l15_col].values
            unknown_mask = cl_values == "unknown"
            n_resolved = 0
            for label, resolved_term in l15_map.items():
                if resolved_term == "unknown":
                    continue
                match = unknown_mask & (l15_series == label)
                n_match = match.sum()
                if n_match > 0:
                    cl_values[match] = resolved_term
                    n_resolved += n_match
            print(f"  Level1.5 fallback: resolved {n_resolved:,} / "
                  f"{n_unknown_before:,} unknown cells")
        else:
            print("  WARNING: No level1.5_annotation column — cannot apply fallback")

    cl_lookup = dict(zip(obs_names, cl_values))
    integrated_set = set(obs_names)
    print(f"  Integrated obs_names set: {len(integrated_set):,} cells")

    n_unknown_after = sum(1 for v in cl_values if v == "unknown")
    print(f"  Unknown CL terms: {n_unknown_before:,} -> {n_unknown_after:,}")

    return cl_lookup, integrated_set


def load_source_obs_names(repo_root, study):
    """Load obs_names from a source h5ad (backed mode)."""
    path = repo_root / "publication/outputs/source_datasets" / SOURCE_H5AD[study]
    print(f"\n  Loading source h5ad: {path.name}")
    adata = ad.read_h5ad(path, backed="r")
    obs_names = list(adata.obs_names)
    print(f"    Shape: {adata.shape[0]:,} cells")
    print(f"    First 5 obs_names: {obs_names[:5]}")
    return obs_names


def load_mapping_files(repo_root, study):
    """Load cell_id_mapping.csv(s) for a study.

    Returns: DataFrame with columns [ihbca_cell_id, component_cell_id].
    For pal, concatenates all 3 sub-study mapping files.
    """
    paths = MAPPING_PATHS[study]
    dfs = []
    for relpath in paths:
        p = repo_root / relpath
        if not p.exists():
            print(f"    WARNING: mapping file not found: {p}")
            continue
        df = pd.read_csv(p)
        # Find the ihbca_cell_id and component_cell_id columns
        comp_cols = [c for c in df.columns if "component_cell_id" in c]
        if not comp_cols:
            print(f"    WARNING: no component_cell_id column in {p.name}")
            continue
        comp_col = comp_cols[0]
        sub_df = df[["ihbca_cell_id", comp_col]].copy()
        sub_df.columns = ["ihbca_cell_id", "component_cell_id"]
        dfs.append(sub_df)
        print(f"    Loaded {p.name}: {len(sub_df):,} rows")

    if not dfs:
        return pd.DataFrame(columns=["ihbca_cell_id", "component_cell_id"])
    return pd.concat(dfs, ignore_index=True)


def investigate_study(repo_root, study, integrated_set):
    """Phase 1: Test cell ID overlap for one study.

    Returns: dict with investigation results and recommended strategy.
    """
    print(f"\n{'='*70}")
    print(f"INVESTIGATING: {study}")
    print(f"{'='*70}")

    source_obs = load_source_obs_names(repo_root, study)
    source_set = set(source_obs)

    # Test 1: Direct overlap (source obs_names vs integrated obs_names)
    direct_overlap = source_set & integrated_set
    print(f"\n  Direct overlap (source vs integrated): {len(direct_overlap):,}/{len(source_set):,}")

    if len(direct_overlap) == len(source_set):
        print(f"  -> FULL MATCH: direct join is sufficient")
        return {
            "study": study,
            "source_cells": len(source_set),
            "strategy": "direct",
            "overlap": len(direct_overlap),
            "coverage_pct": 100.0,
        }
    elif len(direct_overlap) > 0:
        pct = 100 * len(direct_overlap) / len(source_set)
        print(f"  -> PARTIAL MATCH ({pct:.1f}%): direct join covers most cells")
        return {
            "study": study,
            "source_cells": len(source_set),
            "strategy": "direct",
            "overlap": len(direct_overlap),
            "coverage_pct": pct,
        }

    # Test 2: Mapping-mediated overlap
    print(f"\n  No direct overlap. Trying mapping-mediated join...")
    mapping = load_mapping_files(repo_root, study)
    if mapping.empty:
        print(f"  -> NO MAPPING FILE: cannot resolve cell IDs")
        return {
            "study": study,
            "source_cells": len(source_set),
            "strategy": "none",
            "overlap": 0,
            "coverage_pct": 0.0,
        }

    # Test: ihbca_cell_id overlap with integrated
    ihbca_set = set(mapping["ihbca_cell_id"])
    ihbca_to_int = ihbca_set & integrated_set
    print(f"    ihbca_cell_id vs integrated: {len(ihbca_to_int):,}/{len(ihbca_set):,}")

    # Test: component_cell_id overlap with source
    comp_set = set(mapping["component_cell_id"])
    comp_to_src = comp_set & source_set
    print(f"    component_cell_id vs source: {len(comp_to_src):,}/{len(source_set):,}")

    # Test: ihbca_cell_id overlap with source (in case source uses ihbca format)
    ihbca_to_src = ihbca_set & source_set
    print(f"    ihbca_cell_id vs source: {len(ihbca_to_src):,}/{len(source_set):,}")

    # Determine best strategy
    if len(ihbca_to_src) > len(comp_to_src) and len(ihbca_to_src) > 0:
        # Source h5ad uses ihbca format — ihbca_cell_id matches both source and integrated
        strategy = "direct"
        overlap = len(ihbca_to_src & ihbca_to_int)
        pct = 100 * overlap / len(source_set) if source_set else 0
        print(f"  -> Source uses ihbca format. Direct join: {overlap:,} cells ({pct:.1f}%)")
    elif len(comp_to_src) > 0 and len(ihbca_to_int) > 0:
        # Source h5ad uses component format — need mapping crosswalk
        # Effective overlap = cells where component maps to source AND ihbca maps to integrated
        comp_in_src = mapping["component_cell_id"].isin(list(source_set))
        ihbca_in_int = mapping["ihbca_cell_id"].isin(list(integrated_set))
        merged = mapping[comp_in_src & ihbca_in_int]
        strategy = "mapping"
        overlap = len(merged)
        pct = 100 * overlap / len(source_set) if source_set else 0
        print(f"  -> Mapping-mediated join: {overlap:,} cells ({pct:.1f}%)")
    else:
        strategy = "none"
        overlap = 0
        pct = 0.0
        print(f"  -> NO VIABLE JOIN STRATEGY")

    return {
        "study": study,
        "source_cells": len(source_set),
        "strategy": strategy,
        "overlap": overlap,
        "coverage_pct": pct,
    }


# ---------------------------------------------------------------------------
# Phase 2: Build crosswalk CSVs
# ---------------------------------------------------------------------------


def build_crosswalk_direct(cl_lookup, source_obs):
    """Build crosswalk for studies where source obs_names == integrated obs_names."""
    rows = []
    for cid in source_obs:
        cl_term = cl_lookup.get(cid)
        if cl_term is not None:
            rows.append({"source_cell_id": cid, "cell_type_ontology_term_id": cl_term})

    df = pd.DataFrame(rows)
    return df


def build_crosswalk_mapping(repo_root, study, cl_lookup, source_obs):
    """Build crosswalk using cell_id_mapping.csv to bridge source -> integrated."""
    mapping = load_mapping_files(repo_root, study)

    # Build component -> ihbca lookup
    comp_to_ihbca = dict(
        zip(mapping["component_cell_id"], mapping["ihbca_cell_id"])
    )

    rows = []
    for cid in source_obs:
        ihbca_id = comp_to_ihbca.get(cid)
        if ihbca_id is not None:
            cl_term = cl_lookup.get(ihbca_id)
            if cl_term is not None:
                rows.append(
                    {"source_cell_id": cid, "cell_type_ontology_term_id": cl_term}
                )

    df = pd.DataFrame(rows)
    return df


def build_crosswalk(repo_root, study, cl_lookup, result):
    """Build and save crosswalk CSV for one study."""
    print(f"\n{'='*70}")
    print(f"BUILDING CROSSWALK: {study}")
    print(f"  Strategy: {result['strategy']}")
    print(f"{'='*70}")

    if result["strategy"] == "none":
        print(f"  SKIPPED: no viable join strategy")
        return

    source_obs = load_source_obs_names(repo_root, study)

    if result["strategy"] == "direct":
        df = build_crosswalk_direct(cl_lookup, source_obs)
    elif result["strategy"] == "mapping":
        df = build_crosswalk_mapping(repo_root, study, cl_lookup, source_obs)
    else:
        print(f"  ERROR: unknown strategy '{result['strategy']}'")
        return

    if df.empty:
        print(f"  WARNING: crosswalk is empty!")
        return

    # Stats
    n_source = len(source_obs)
    n_mapped = len(df)
    n_unknown = n_source - n_mapped
    print(f"  Mapped: {n_mapped:,}/{n_source:,} ({100*n_mapped/n_source:.1f}%)")
    print(f"  Remaining unknown: {n_unknown:,}")

    # CL term distribution
    print(f"  CL term distribution (top 15):")
    vc = df["cell_type_ontology_term_id"].value_counts()
    for term, count in vc.head(15).items():
        print(f"    {term}: {count:,}")
    print(f"  Total unique CL terms: {vc.shape[0]}")

    # Save
    out_path = repo_root / f"publication/mappings/cl_term_crosswalk_{study}.csv"
    df.to_csv(out_path, index=False)
    print(f"  Saved: {out_path} ({len(df):,} rows)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root", required=True, type=Path, help="Path to iHBCAv1_upload root"
    )
    parser.add_argument(
        "--investigate-only",
        action="store_true",
        help="Phase 1 only: investigate without building crosswalks",
    )
    parser.add_argument(
        "--studies",
        nargs="*",
        default=TARGET_STUDIES,
        help="Studies to process (default: murrow nee pal)",
    )
    args = parser.parse_args()

    repo_root = args.repo_root.resolve()
    if not (repo_root / "publication").exists():
        print(f"ERROR: {repo_root} doesn't look like iHBCAv1_upload root")
        sys.exit(1)

    # Phase 1: Load integrated and investigate each study
    cl_lookup, integrated_set = load_integrated_cl_terms(repo_root)

    results = {}
    for study in args.studies:
        results[study] = investigate_study(repo_root, study, integrated_set)

    # Summary
    print(f"\n{'='*70}")
    print(f"INVESTIGATION SUMMARY")
    print(f"{'='*70}")
    for study, r in results.items():
        print(
            f"  {study}: strategy={r['strategy']}, "
            f"overlap={r['overlap']:,}/{r['source_cells']:,} "
            f"({r['coverage_pct']:.1f}%)"
        )

    if args.investigate_only:
        print("\n--investigate-only: stopping before crosswalk generation.")
        return

    # Phase 2: Build crosswalks
    for study in args.studies:
        build_crosswalk(repo_root, study, cl_lookup, results[study])

    print(f"\n{'='*70}")
    print(f"DONE")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
