#!/usr/bin/env python3
"""
patch_obs_h5py.py - Patch obs columns directly in h5ad via h5py
================================================================
Adds/overwrites L1 harmonized donor metadata columns in an h5ad file
without loading the full object into memory. Operates directly on the
HDF5 structure that anndata uses.

Memory: ~2-5 GB (only obs dataframe, not X/layers/obsm)

Usage:
  python patch_obs_h5py.py \
    --h5ad path/to/all-breast-cells.h5ad \
    --repo-root /path/to/repo \
    --reorder  # also reorder obs columns

"""

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))

from ihbca.constants import L1_NUMERIC, L1_OBS_COLUMNS
from ihbca.loaders import load_donor_translations, load_l1_metadata


def read_obs_from_h5ad(h5ad_path):
    """Read obs dataframe from h5ad without loading X/layers/obsm."""
    with h5py.File(h5ad_path, "r") as f:
        obs_group = f["obs"]

        # Get index dataset name from attrs (anndata stores it as _index attr)
        index_col_name = obs_group.attrs.get("_index", "_index")
        if isinstance(index_col_name, bytes):
            index_col_name = index_col_name.decode()

        # Read the index dataset
        if index_col_name in obs_group:
            index = obs_group[index_col_name][()].astype(str)
        elif "_index" in obs_group:
            # fallback: literal _index dataset
            index = obs_group["_index"][()].astype(str)
        elif "__categories" in obs_group.attrs:
            # older format
            index = obs_group["index"][()].astype(str)
        else:
            index = np.arange(obs_group.attrs.get("_n_obs", 0))

        # Get column order from encoding
        if "column-order" in obs_group.attrs:
            col_order = list(obs_group.attrs["column-order"].astype(str))
        elif "__categories" in obs_group.attrs:
            col_order = [k for k in obs_group.keys() if k not in ("_index", "index", "__categories")]
        else:
            col_order = list(obs_group.attrs.get("column-order", []))

        # Exclude the index column from data columns
        index_exclude = {index_col_name, "_index"}
        col_order = [c for c in col_order if c not in index_exclude]

        # Read each column
        data = {}
        for col in col_order:
            if col not in obs_group:
                continue
            ds = obs_group[col]

            # Check if categorical
            encoding_type = ds.attrs.get("encoding-type", "")
            if encoding_type == "categorical":
                codes = ds["codes"][()]
                categories = ds["categories"][()].astype(str)
                ordered = ds.attrs.get("ordered", False)
                values = pd.Categorical.from_codes(
                    codes, categories=categories, ordered=ordered
                )
                data[col] = values
            elif encoding_type == "nullable-integer":
                data[col] = pd.array(ds[()], dtype=pd.Int64Dtype())
            else:
                arr = ds[()]
                if arr.dtype.kind in ("S", "U", "O"):
                    arr = arr.astype(str)
                data[col] = arr

        obs = pd.DataFrame(data, index=index)
        obs.index.name = obs_group.attrs.get("_index", None)

    return obs


def write_obs_to_h5ad(h5ad_path, obs):
    """Write obs dataframe back to h5ad, replacing existing obs group."""
    with h5py.File(h5ad_path, "a") as f:
        # Preserve obs attrs
        old_attrs = dict(f["obs"].attrs)

        # Delete existing obs group
        del f["obs"]

        # Create new obs group
        obs_group = f.create_group("obs")

        # Write index under its original dataset name
        index_col_name = old_attrs.get("_index", "_index")
        if isinstance(index_col_name, bytes):
            index_col_name = index_col_name.decode()
        index_data = np.array(obs.index, dtype=h5py.string_dtype())
        obs_group.create_dataset(index_col_name, data=index_data)

        # Update attrs
        col_order = np.array(list(obs.columns), dtype=h5py.string_dtype())
        obs_group.attrs["column-order"] = col_order
        obs_group.attrs["_index"] = index_col_name
        obs_group.attrs["encoding-type"] = "dataframe"
        obs_group.attrs["encoding-version"] = old_attrs.get("encoding-version", "0.2.0")

        # Write each column
        for col in obs.columns:
            series = obs[col]

            if isinstance(series.dtype, pd.CategoricalDtype):
                # Categorical encoding
                grp = obs_group.create_group(col)
                cats = np.array(series.cat.categories, dtype=h5py.string_dtype())
                codes = series.cat.codes.values
                grp.create_dataset("categories", data=cats)
                grp.create_dataset("codes", data=codes)
                grp.attrs["encoding-type"] = "categorical"
                grp.attrs["encoding-version"] = "0.2.0"
                grp.attrs["ordered"] = series.cat.ordered
            elif isinstance(series.dtype, (pd.Float64Dtype, pd.Int64Dtype)):
                # Nullable numeric — convert to regular numpy with NaN
                vals = series.to_numpy(dtype="float64", na_value=np.nan)
                obs_group.create_dataset(col, data=vals)
            elif series.dtype == bool or series.dtype == np.bool_:
                obs_group.create_dataset(col, data=series.values)
                obs_group[col].attrs["encoding-type"] = "array"
                obs_group[col].attrs["encoding-version"] = "0.2.0"
            elif series.dtype == object or series.dtype.kind in ("U", "S"):
                str_data = np.array(series.fillna(""), dtype=h5py.string_dtype())
                obs_group.create_dataset(col, data=str_data)
            else:
                obs_group.create_dataset(col, data=series.values)


def build_reverse_map(l1, xlat):
    """Build h5ad donor_id -> L1 row mapping."""
    reverse_map = {}
    for _, row in l1.iterrows():
        study_cap = row["study"].capitalize()
        h5ad_key = f"{study_cap}_{row['ihbca_donor_id']}"
        reverse_map[h5ad_key] = row

    for l1_key, cxg_id in xlat.items():
        parts = l1_key.split("_", 1)
        study_cap = parts[0]
        ihbca_id = parts[1]
        match = l1[(l1["study"] == study_cap.lower()) & (l1["ihbca_donor_id"] == ihbca_id)]
        if len(match) > 0:
            reverse_map[cxg_id] = match.iloc[0]
            if not cxg_id.startswith(f"{study_cap}_"):
                reverse_map[f"{study_cap}_{cxg_id}"] = match.iloc[0]

    return reverse_map


OBS_COLUMN_ORDER = [
    "organism_ontology_term_id", "assay_ontology_term_id",
    "tissue_ontology_term_id", "tissue_type", "disease_ontology_term_id",
    "donor_id", "sex_ontology_term_id", "development_stage_ontology_term_id",
    "self_reported_ethnicity_ontology_term_id", "is_primary_data",
    "suspension_type", "cell_type_ontology_term_id",
    "cell_type_label", "assay_label", "disease_label", "sex_label",
    "tissue_label", "self_reported_ethnicity_label", "development_stage_label",
    "ihbca_donor_id", "dataset",
    "age_continuous", "age_binary", "parity_count", "parity_binary",
    "age_at_first_birth", "brca_genotype", "cancer_history",
    "tissue_indication", "risk_status_binary", "risk_genotype_only",
    "menopausal_status_detailed", "menopausal_status_binary",
    "ethnicity_verbatim", "ethnicity_grouped",
    "bmi_continuous", "bmi_category", "sample_preservation",
    "sample_type", "facs_status", "dissociation_minutes", "metadata_notes",
    "sample_id", "manner_of_death", "sample_source", "sampled_site_condition",
    "sample_collection_method", "sample_preservation_method", "cell_enrichment",
    "institute", "library_id", "library_sequencing_run",
    "library_preparation_batch", "alignment_software", "sequenced_fragment",
    "reference_genome", "gene_annotation_version",
    "level0_annotation", "level1_annotation", "level1.5_annotation",
    "level0", "level1",
    "cellTypist_annotation_reed", "cellTypist_annotation_kumar",
    "n_genes", "percent_mito", "n_counts",
    "prob_spikein", "prob_spikein_dblt", "pred_spikein",
]


def patch_obs(h5ad_path, repo_root, reorder=False):
    """Patch obs in h5ad with all 21 L1 columns + ihbca_donor_id."""
    repo_root = Path(repo_root)

    print(f"Patching obs in: {h5ad_path}")
    print(f"Reading obs from h5ad...")
    obs = read_obs_from_h5ad(h5ad_path)
    print(f"  Shape: {obs.shape}")
    print(f"  Columns: {len(obs.columns)}")

    # Load L1
    l1 = load_l1_metadata(repo_root)
    if l1 is None:
        print("ERROR: L1 metadata not found")
        sys.exit(1)

    xlat = load_donor_translations(repo_root)
    reverse_map = build_reverse_map(l1, xlat)

    donor_ids = obs["donor_id"] if "donor_id" in obs.columns else obs.index
    unique_donors = donor_ids.unique() if hasattr(donor_ids, "unique") else set(donor_ids)
    n_matched = sum(1 for did in unique_donors if did in reverse_map)
    print(f"  Reverse map: {n_matched}/{len(unique_donors)} unique donor_ids matched")

    # --- ihbca_donor_id ---
    print("  Adding ihbca_donor_id...")
    ihbca_ids = []
    for did in donor_ids:
        did_str = str(did)
        if did_str in reverse_map:
            ihbca_ids.append(reverse_map[did_str]["ihbca_donor_id"])
        else:
            ihbca_ids.append(did_str.split("_", 1)[1] if "_" in did_str else did_str)
    obs["ihbca_donor_id"] = ihbca_ids

    # --- Broadcast all 21 L1 columns ---
    print(f"  Broadcasting {len(L1_OBS_COLUMNS)} L1 columns...")
    for col in L1_OBS_COLUMNS:
        is_numeric = col in L1_NUMERIC
        is_notes = col == "metadata_notes"

        values = []
        for did in donor_ids:
            did_str = str(did)
            if did_str in reverse_map:
                val = reverse_map[did_str].get(col, "")
                if pd.isna(val) or (isinstance(val, str) and val.strip() == ""):
                    if is_numeric:
                        values.append(np.nan)
                    elif is_notes:
                        values.append("")
                    else:
                        values.append("unknown")
                else:
                    if is_numeric:
                        try:
                            values.append(float(val))
                        except (ValueError, TypeError):
                            values.append(np.nan)
                    else:
                        values.append(str(val))
            else:
                values.append(np.nan if is_numeric else ("" if is_notes else "unknown"))

        if is_numeric:
            obs[col] = pd.array(values, dtype=pd.Float64Dtype())
            n_present = sum(1 for v in values if not pd.isna(v))
            print(f"    {col}: {n_present}/{len(obs)} non-NaN")
        else:
            obs[col] = values
            if is_notes:
                n_filled = sum(1 for v in values if v != "")
                print(f"    {col}: {n_filled}/{len(obs)} with notes")
            else:
                n_unknown = sum(1 for v in values if v == "unknown")
                print(f"    {col}: {len(obs) - n_unknown}/{len(obs)} mapped, {n_unknown} unknown")

    # --- Reorder columns ---
    if reorder:
        ordered = [c for c in OBS_COLUMN_ORDER if c in obs.columns]
        remaining = sorted(c for c in obs.columns if c not in set(ordered))
        obs = obs[ordered + remaining]
        print(f"  Reordered: {len(ordered)} priority + {len(remaining)} remaining")

    # --- Convert string columns to categorical ---
    n_cat = 0
    for col in obs.columns:
        if obs[col].dtype == object:
            if col == "is_primary_data":
                continue
            if col in L1_NUMERIC:
                continue
            n_unique = obs[col].nunique()
            if n_unique > 10000:
                continue
            obs[col] = pd.Categorical(obs[col])
            n_cat += 1
    print(f"  Categoricals: {n_cat} string columns converted")

    # --- Write back ---
    print(f"  Writing obs back to h5ad ({len(obs.columns)} columns)...")
    write_obs_to_h5ad(h5ad_path, obs)
    print(f"  Done.")

    print(f"\n  Final obs: {obs.shape}")
    print(f"  First 25 columns: {list(obs.columns[:25])}")


def main():
    parser = argparse.ArgumentParser(
        description="Patch obs columns in h5ad via h5py (no full load)"
    )
    parser.add_argument("--h5ad", required=True, help="Path to h5ad file")
    parser.add_argument("--repo-root", required=True, help="Path to repo root")
    parser.add_argument("--reorder", action="store_true", help="Reorder obs columns")
    args = parser.parse_args()
    patch_obs(args.h5ad, args.repo_root, args.reorder)


if __name__ == "__main__":
    main()
