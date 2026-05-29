#!/usr/bin/env python3
"""
assemble_h5ad.py - Build CxG-compliant h5ad from extracted intermediates
=========================================================================
Reads count matrix intermediates (from extract_counts_for_h5ad.R) plus
published metadata and UMAP coordinates, maps gene symbols to Ensembl IDs,
populates CxG-required ontology fields, and writes a single h5ad file.

For Pal (3 sub-studies -> 1 h5ad), handles concatenation and deduplication.

Usage:
  python assemble_h5ad.py --study gray --repo-root /path/to/iHBCAv1_upload
  python assemble_h5ad.py --study pal --repo-root /path/to/iHBCAv1_upload

"""

# Another day, another file format conversion. At least this one has a schema.

import argparse
import gzip
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from provenance import build_provenance_manifest

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

from ihbca.constants import (
    DATASET_META_FIELDS,
    FACS_TO_ENRICHMENT,
    L1_NUMERIC,
    L1_OBS_COLUMNS,
    PRESERVATION_NORMALIZE,
    RESERVED_OBS_RENAMES,
    TIER1_FIELDS_SOURCE as TIER1_FIELDS,
    UNS_TO_OBS_FIELDS,
    get_cxg_studies,
    get_pal_sub_studies,
)
from ihbca.ethnicity import downgrade_hancestro_terms, normalize_multi_ethnicity
from ihbca.hca_fields import populate_dataset_metadata, populate_derived_obs_fields
from ihbca.validation import (
    validate_dataset_metadata,
    validate_ensembl_coverage,
    validate_tier1_fields,
    validate_umap,
)
from ihbca.loaders import (
    load_cl_crosswalk,
    load_cxg_approved_genes,
    load_dataset_metadata,
    load_efo_assay_mapping,
    load_gene_mapping,
    load_hancestro_mapping,
    load_l1_metadata,
    load_pipeline_config,
    load_registry,
    get_study_inputs,
    get_study_output,
)

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def load_count_matrix(intermediates_dir):
    """Load sparse count matrix from 10x-format intermediates.

    R writes genes x cells (Market Matrix). We transpose to cells x genes.
    Returns: (csr_matrix, gene_names, cell_ids)
    """
    mtx_path = intermediates_dir / "counts.mtx.gz"
    if not mtx_path.exists():
        raise FileNotFoundError(f"Count matrix not found: {mtx_path}")

    mat = scipy.io.mmread(str(mtx_path))  # genes x cells, COO
    mat = sp.csr_matrix(mat.T)  # cells x genes, CSR

    feat_path = intermediates_dir / "features.tsv.gz"
    with gzip.open(str(feat_path), "rt") as f:
        gene_names = [line.strip() for line in f if line.strip()]

    bc_path = intermediates_dir / "barcodes.tsv.gz"
    with gzip.open(str(bc_path), "rt") as f:
        cell_ids = [line.strip() for line in f if line.strip()]

    if mat.shape != (len(cell_ids), len(gene_names)):
        raise ValueError(
            f"Matrix shape {mat.shape} != ({len(cell_ids)}, {len(gene_names)})"
        )

    return mat, gene_names, cell_ids


def enrich_donor_metadata(obs, study, repo_root, config=None):
    """Broadcast L1 harmonized donor metadata to per-cell obs.

    Reads L1_harmonized_donor.csv, filters to the current study,
    and joins on donor_id (= ihbca_donor_id for source datasets).
    """
    if repo_root is None:
        return obs

    l1 = load_l1_metadata(repo_root, config=config)
    if l1 is None:
        return obs

    # Filter to current study
    study_l1 = l1[l1["study"] == study]
    if study_l1.empty:
        print(f"  L1 donor metadata: no entries for study '{study}'")
        return obs

    print(f"  L1 donor metadata: {len(study_l1)} donors for {study}")

    # Build donor_id -> column value mappings
    n_matched = 0
    obs_donors = set(obs["donor_id"].unique()) if "donor_id" in obs.columns else set()

    for col in L1_OBS_COLUMNS:
        if col not in study_l1.columns:
            continue
        donor_map = dict(zip(study_l1["ihbca_donor_id"], study_l1[col]))
        obs[col] = obs["donor_id"].map(donor_map)
        if col in L1_NUMERIC:
            obs[col] = pd.to_numeric(obs[col], errors="coerce")

    # Count matched donors
    l1_donors = set(study_l1["ihbca_donor_id"])
    matched = obs_donors & l1_donors
    n_matched = len(matched)
    unmatched = obs_donors - l1_donors
    if unmatched:
        print(f"  WARNING: {len(unmatched)} h5ad donors not in L1: "
              f"{sorted(unmatched)[:10]}")

    n_cells = len(obs)
    matched_cells = obs["donor_id"].isin(l1_donors).sum() if "donor_id" in obs.columns else 0
    print(f"  Matched: {n_matched}/{len(obs_donors)} donors "
          f"({matched_cells:,}/{n_cells:,} cells)")
    print(f"  Added {len(L1_OBS_COLUMNS)} L1 columns to obs")

    return obs


def load_metadata(published_dir, cell_ids, fallback_dir=None):
    """Load cell metadata from published/metadata.csv, aligned to cell_ids.

    Falls back to fallback_dir/metadata.csv if published one is missing or empty
    (e.g. kumar's published metadata.csv is 0 bytes).

    Returns DataFrame indexed by cell_id with rows matching cell_ids order.
    """
    meta_path = published_dir / "metadata.csv"

    # Check if published metadata exists and is non-empty
    if not meta_path.exists() or meta_path.stat().st_size == 0:
        if fallback_dir is not None:
            fallback_path = fallback_dir / "metadata.csv"
            if fallback_path.exists() and fallback_path.stat().st_size > 0:
                print(f"  Published metadata missing/empty, using: {fallback_path}")
                meta_path = fallback_path
            else:
                raise FileNotFoundError(
                    f"Metadata not found in published ({published_dir}) "
                    f"or intermediates ({fallback_dir})"
                )
        else:
            raise FileNotFoundError(f"Metadata not found or empty: {meta_path}")

    meta = pd.read_csv(meta_path, low_memory=False)

    # Identify the cell ID column
    if "cell_id" in meta.columns:
        meta = meta.set_index("cell_id")
    else:
        # Some CSVs may use the first unnamed column as row names
        first_col = meta.columns[0]
        if first_col.startswith("Unnamed"):
            meta = meta.set_index(first_col)
            meta.index.name = "cell_id"
        else:
            raise ValueError(
                f"Cannot identify cell ID column. Columns: {list(meta.columns[:10])}"
            )

    # Align to barcodes order
    missing = set(cell_ids) - set(meta.index)
    if missing:
        pct = len(missing) / len(cell_ids) * 100
        print(f"  WARNING: {len(missing)} / {len(cell_ids)} ({pct:.1f}%) "
              f"cell IDs from barcodes not found in metadata")

    meta = meta.reindex(cell_ids)
    return meta


def load_umap(published_dir, umap_type, cell_ids):
    """Load UMAP coordinates from published directory.

    Tries umap_{type}.csv first. If not found:
      - native -> falls back to joint
      - joint -> returns None (Pal sub-studies have no UMAP files)

    Args:
        published_dir: Path to published/ directory
        umap_type: 'native' or 'joint'
        cell_ids: cell IDs for row alignment

    Returns: numpy array (n_cells, 2) or None if no UMAP available
    """
    umap_path = published_dir / f"umap_{umap_type}.csv"

    if not umap_path.exists():
        if umap_type == "native":
            print(f"  umap_native.csv not found, falling back to umap_joint.csv")
            return load_umap(published_dir, "joint", cell_ids)
        print(f"  No UMAP file found ({umap_path})")
        return None

    umap_df = pd.read_csv(umap_path, index_col=0)

    # Verify this is actually 2D UMAP (not a high-dim embedding)
    if "UMAP_1" not in umap_df.columns or "UMAP_2" not in umap_df.columns:
        print(f"  WARNING: {umap_path.name} lacks UMAP_1/UMAP_2 columns "
              f"(has: {list(umap_df.columns[:5])}...). Skipping.")
        return None

    # Align to cell_ids
    umap_df = umap_df.reindex(cell_ids)
    n_missing = umap_df.isna().any(axis=1).sum()
    if n_missing > 0:
        pct = n_missing / len(cell_ids) * 100
        print(f"  WARNING: {n_missing} ({pct:.1f}%) cells missing UMAP coordinates")

    # Fill missing with 0 (better than NaN for downstream tools)
    coords = umap_df[["UMAP_1", "UMAP_2"]].fillna(0.0).values
    return coords.astype(np.float32)


def load_embedding(published_dir, embedding_type, cell_ids):
    """Load high-dimensional embedding from published directory.

    Args:
        published_dir: Path to published/ directory
        embedding_type: 'joint' or 'native'
        cell_ids: cell IDs for row alignment

    Returns: numpy array (n_cells, n_dims) or None if CSV not found
    """
    csv_path = published_dir / f"embedding_{embedding_type}.csv"
    if not csv_path.exists():
        print(f"  embedding_{embedding_type}.csv not found")
        return None

    df = pd.read_csv(csv_path, index_col=0)
    df = df.reindex(cell_ids)

    n_missing = df.isna().any(axis=1).sum()
    if n_missing > 0:
        pct = n_missing / len(cell_ids) * 100
        print(f"  WARNING: {n_missing} ({pct:.1f}%) cells missing "
              f"{embedding_type} embedding")

    arr = df.fillna(0.0).values.astype(np.float32)
    print(f"  embedding_{embedding_type}: {arr.shape}")
    return arr


# ---------------------------------------------------------------------------
# Gene ID mapping
# ---------------------------------------------------------------------------


def map_gene_ids(gene_names, gene_mapping):
    """Map gene symbols to Ensembl IDs.

    Unmapped genes keep their original symbol as the ID.
    Returns: (ensembl_ids list, stats dict)
    """
    ensembl_ids = []
    unmapped = []

    for symbol in gene_names:
        eid = gene_mapping.get(symbol)
        if eid:
            ensembl_ids.append(eid)
        else:
            ensembl_ids.append(symbol)
            unmapped.append(symbol)

    stats = {
        "total": len(gene_names),
        "mapped": len(gene_names) - len(unmapped),
        "unmapped": len(unmapped),
        "pct_mapped": (len(gene_names) - len(unmapped)) / max(len(gene_names), 1) * 100,
        "unmapped_examples": unmapped[:20],
    }

    return ensembl_ids, stats


# ---------------------------------------------------------------------------
# CxG field population
# ---------------------------------------------------------------------------


def populate_cxg_fields(obs, study, efo_mapping, hancestro_mapping, repo_root=None,
                        target="hca", config=None):
    """Add CxG-required ontology fields to obs DataFrame.

    CxG studies (gray, twigger, reed): passthrough existing columns, add
    organism_ontology_term_id, fix types.

    Non-CxG studies: populate all 11 Tier 1 fields from mapping tables.
    CL terms looked up from pre-computed crosswalk if available.

    target: 'hca' (default) keeps HCA-native HANCESTRO terms (:0601/:0602).
            'cxg' applies downgrade to :0004 ancestry terms for CxG 5.3.2.
    """
    cxg_studies = get_cxg_studies(config) if config else {"gray", "kumar", "twigger", "reed"}
    is_cxg = study in cxg_studies

    # --- All studies ---
    obs["organism_ontology_term_id"] = "NCBITaxon:9606"
    obs["tissue_type"] = "tissue"

    # Canonical donor_id = ihbca_donor_id
    if "ihbca_donor_id" in obs.columns:
        obs["donor_id"] = obs["ihbca_donor_id"]

    if is_cxg:
        print(f"  CxG study — passthrough existing ontology columns")

        # Fix is_primary_data type (string -> bool)
        if "is_primary_data" in obs.columns:
            obs["is_primary_data"] = obs["is_primary_data"].map(
                {"True": True, "False": False, True: True, False: False}
            ).fillna(True).astype(bool)

        # A4: Drop CxG passthrough column reserved by HCA
        if "observation_joinid" in obs.columns:
            obs = obs.drop(columns=["observation_joinid"])
            print(f"  Dropped reserved column: observation_joinid")

    else:
        print(f"  Non-CxG study — populating ontology fields from mappings")

        # Constants
        obs["tissue_ontology_term_id"] = "UBERON:0000310"
        obs["sex_ontology_term_id"] = "PATO:0000383"
        obs["disease_ontology_term_id"] = "PATO:0000461"
        obs["is_primary_data"] = True
        obs["suspension_type"] = "cell"
        obs["development_stage_ontology_term_id"] = "HsapDv:0000258"
        # CL term backfill from integrated object crosswalk
        cl_crosswalk = load_cl_crosswalk(repo_root, study, config=config) if repo_root else {}
        if cl_crosswalk:
            obs["cell_type_ontology_term_id"] = obs.index.map(
                lambda cid: cl_crosswalk.get(cid, "unknown")
            )
            n_mapped = (obs["cell_type_ontology_term_id"] != "unknown").sum()
            print(f"  CL terms: {n_mapped:,}/{len(obs):,} mapped from crosswalk")
        else:
            obs["cell_type_ontology_term_id"] = "unknown"
            if repo_root:
                print(f"  CL terms: no crosswalk for {study}")

        # Assay: lookup by study name
        efo_study = "pal" if study == "pal" else study
        assay_term = efo_mapping.get(efo_study)
        if assay_term:
            obs["assay_ontology_term_id"] = assay_term
        else:
            print(f"  WARNING: No EFO assay term for '{efo_study}'")
            obs["assay_ontology_term_id"] = "unknown"

        # Ethnicity: lookup by ethnicity_verbatim
        if "ethnicity_verbatim" in obs.columns:
            obs["self_reported_ethnicity_ontology_term_id"] = (
                obs["ethnicity_verbatim"]
                .fillna("")
                .map(lambda v: hancestro_mapping.get(v, "unknown") if v else "unknown")
            )
        else:
            obs["self_reported_ethnicity_ontology_term_id"] = "unknown"

    # --- Post-processing for ALL studies ---

    # Normalize multi-value ethnicity: sorted, ascending order
    # HCA uses ' || ' delimiter, CxG uses ','
    if "self_reported_ethnicity_ontology_term_id" in obs.columns:
        sep = "," if target == "cxg" else " || "
        obs["self_reported_ethnicity_ontology_term_id"] = normalize_multi_ethnicity(
            obs["self_reported_ethnicity_ontology_term_id"], sep=sep
        )

    # CxG 5.3.2 compatibility: downgrade HCA-native HANCESTRO terms to :0004 branch
    # Only applied when --target cxg. HCA builds keep :0601/:0602 terms as-is.
    if target == "cxg" and "self_reported_ethnicity_ontology_term_id" in obs.columns:
        obs["self_reported_ethnicity_ontology_term_id"] = downgrade_hancestro_terms(
            obs["self_reported_ethnicity_ontology_term_id"]
        )
        # Re-normalize after downgrade — term replacement can break sort order
        obs["self_reported_ethnicity_ontology_term_id"] = normalize_multi_ethnicity(
            obs["self_reported_ethnicity_ontology_term_id"], sep=","
        )

    # A10: Drop deprecated ethnicity columns
    for col in ["ethnicity", "ethnicity_ontology_term_id"]:
        if col in obs.columns:
            obs = obs.drop(columns=[col])
            print(f"  Dropped deprecated column: {col}")

    # A3: Rename reserved obs columns (HCA reserves these label names)
    renames = {k: v for k, v in RESERVED_OBS_RENAMES.items() if k in obs.columns}
    if renames:
        obs = obs.rename(columns=renames)
        print(f"  Renamed {len(renames)} reserved obs columns: "
              f"{list(renames.keys())}")

    return obs


# ---------------------------------------------------------------------------
# HCA obs field population
# ---------------------------------------------------------------------------

# UNS_TO_OBS_FIELDS and FACS_TO_ENRICHMENT imported from ihbca.constants


def load_sra_run_table(repo_root, study):
    """Load per-study SRA run table if available.

    Returns: DataFrame with at least 'sample_id' column, or None.
    """
    path = repo_root / f"publication/mappings/sra_run_table_{study}.csv"
    if not path.exists():
        return None
    return pd.read_csv(path, dtype=str)


def populate_hca_obs_fields(adata, study, repo_root):
    """Populate HCA-required obs fields not covered by CxG Tier 1.

    B1: Broadcast 4 fields from uns to obs
    B2: Derive 8 constant/lookup fields
    B3: Library metadata from SRA run tables (3 fields)
    """
    obs = adata.obs
    print(f"  Populating HCA obs fields for {study}...")

    # --- B1: Broadcast uns -> obs ---
    for field in UNS_TO_OBS_FIELDS:
        if field in adata.uns:
            obs[field] = adata.uns[field]

    # --- B2: Derived constant/lookup fields ---
    obs = populate_derived_obs_fields(obs, FACS_TO_ENRICHMENT, PRESERVATION_NORMALIZE)

    # institute: per-study from dataset_metadata.yaml
    dmeta = load_dataset_metadata(repo_root)
    study_meta = dmeta.get("datasets", {}).get(study, {})
    obs["institute"] = study_meta.get("institute", "unknown")

    # --- B3: Library metadata from SRA run tables ---
    sra = load_sra_run_table(repo_root, study)
    if sra is not None and "donor_id" in obs.columns:
        # Map library metadata via donor_id -> SRA sample
        if "donor_id" in sra.columns:
            sra_map = sra.set_index("donor_id")
            for field in ["library_id", "library_sequencing_run", "library_preparation_batch"]:
                if field in sra_map.columns:
                    obs[field] = obs["donor_id"].map(
                        sra_map[field].to_dict()
                    ).fillna("unknown")
                else:
                    obs[field] = "unknown"
        else:
            for field in ["library_id", "library_sequencing_run", "library_preparation_batch"]:
                obs[field] = "unknown"
    else:
        for field in ["library_id", "library_sequencing_run", "library_preparation_batch"]:
            obs[field] = "unknown"

    n_fields = len(UNS_TO_OBS_FIELDS) + 8 + 3  # B1 + B2 + B3
    print(f"  Populated {n_fields} HCA obs fields")
    return obs


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def validate_h5ad(adata, study, registry):
    """Light validation checks on assembled h5ad. Returns list of issues."""
    issues = []

    # Cell count vs registry
    datasets = registry.get("datasets", {})
    study_info = datasets.get(study, {})
    expected = study_info.get("expected_cells")
    if expected is not None:
        diff = abs(adata.n_obs - expected)
        pct_diff = diff / max(expected, 1) * 100
        if pct_diff > 5:
            issues.append(
                f"Cell count: {adata.n_obs:,} vs expected {expected:,} "
                f"({pct_diff:.1f}% difference)"
            )
        else:
            print(f"  Cell count: {adata.n_obs:,} (expected {expected:,}, "
                  f"diff {pct_diff:.1f}%)")

    # Tier 1 fields
    issues.extend(validate_tier1_fields(adata, TIER1_FIELDS))

    # Ensembl coverage
    issues.extend(validate_ensembl_coverage(adata))

    # UMAP
    issues.extend(validate_umap(adata, check_shape=True))

    # Dataset metadata
    issues.extend(validate_dataset_metadata(adata, DATASET_META_FIELDS))

    return issues


# ---------------------------------------------------------------------------
# Assembly: single study
# ---------------------------------------------------------------------------


def assemble_single_study(study, repo_root, gene_mapping, efo_mapping,
                          hancestro_mapping, registry,
                          intermediates_dir=None, published_dir=None,
                          output_path=None, cxg_approved=None, target="hca",
                          config=None):
    """Assemble h5ad for a single (non-Pal) study."""

    if intermediates_dir is None:
        if config:
            intermediates_dir = repo_root / config["studies"][study]["inputs"]["intermediates"]
        else:
            intermediates_dir = (
                repo_root / "publication/outputs/source_datasets"
                / "intermediates" / study
            )
    if published_dir is None:
        if config:
            published_dir = repo_root / config["studies"][study]["inputs"]["published"]
        else:
            published_dir = (
                repo_root / "external_studies/outputs" / study / "published"
            )
    if output_path is None:
        if config:
            output_path = get_study_output(repo_root, study, config)
        else:
            filename = registry["datasets"][study]["output_filename"]
            output_path = (
                repo_root / "publication/outputs/source_datasets" / filename
            )

    print(f"\n{'=' * 70}")
    print(f"ASSEMBLING: {study}")
    print(f"  Intermediates: {intermediates_dir}")
    print(f"  Published:     {published_dir}")
    print(f"  Output:        {output_path}")
    print(f"{'=' * 70}")

    # 1. Load counts
    print(f"\n[1/6] Loading count matrix...")
    mat, gene_names, cell_ids = load_count_matrix(intermediates_dir)
    print(f"  Shape: {mat.shape[0]:,} cells x {mat.shape[1]:,} genes")
    print(f"  Non-zero: {mat.nnz:,}")

    # 2. Load metadata (falls back to intermediates for kumar etc.)
    print(f"\n[2/6] Loading metadata...")
    obs = load_metadata(published_dir, cell_ids, fallback_dir=intermediates_dir)
    print(f"  Columns: {obs.shape[1]}")

    # 3. Load UMAP (optional — Pal sub-studies have none)
    print(f"\n[3/6] Loading UMAP coordinates...")
    umap_coords = load_umap(published_dir, "native", cell_ids)
    if umap_coords is not None:
        print(f"  Shape: {umap_coords.shape}")
    else:
        print(f"  No UMAP available for {study}")

    # 3b. Load high-dimensional embeddings (optional)
    emb_joint = load_embedding(published_dir, "joint", cell_ids)
    emb_native = load_embedding(published_dir, "native", cell_ids)

    # 4. Map gene IDs
    print(f"\n[4/6] Mapping gene symbols -> Ensembl IDs...")
    ensembl_ids, mapping_stats = map_gene_ids(gene_names, gene_mapping)
    print(f"  Mapped: {mapping_stats['mapped']:,}/{mapping_stats['total']:,} "
          f"({mapping_stats['pct_mapped']:.1f}%)")
    if mapping_stats["unmapped_examples"]:
        print(f"  Unmapped examples: {mapping_stats['unmapped_examples'][:10]}")

    # 5. Populate CxG fields
    print(f"\n[5/6] Populating CxG ontology fields...")
    obs = populate_cxg_fields(obs, study, efo_mapping, hancestro_mapping, repo_root,
                              target=target, config=config)

    # 5b. Enrich donor metadata from L1 harmonized CSV
    print(f"\n[5b] Enriching donor metadata from L1 harmonized CSV...")
    obs = enrich_donor_metadata(obs, study, repo_root, config=config)

    # 6. Build AnnData
    print(f"\n[6/6] Building AnnData object...")

    var = pd.DataFrame({"gene_symbol": gene_names}, index=ensembl_ids)
    var.index.name = None

    # Filter genes against CxG approved set (GENCODE v44 / Ensembl 110)
    if cxg_approved is not None:
        approved_mask = var.index.isin(cxg_approved)
        if not approved_mask.all():
            n_drop = (~approved_mask).sum()
            dropped = var.index[~approved_mask].tolist()
            print(f"  Dropping {n_drop} genes not in CxG approved set")
            print(f"  Examples: {dropped[:20]}")
            keep = np.asarray(approved_mask)
            mat = mat[:, keep]
            var = var[keep]

    # Merge duplicate Ensembl IDs by summing counts
    mat, var = _merge_duplicate_genes(mat, var)

    adata = ad.AnnData(X=mat, obs=obs, var=var)
    if umap_coords is not None:
        adata.obsm["X_umap"] = umap_coords
    if emb_joint is not None:
        adata.obsm["X_scVI_joint"] = emb_joint
    if emb_native is not None:
        adata.obsm["X_scVI_native"] = emb_native

    # CxG requires raw.X to be float32
    if sp.issparse(adata.X) and adata.X.dtype != np.float32:
        print(f"  Casting X from {adata.X.dtype} to float32")
        adata.X = adata.X.astype(np.float32)
    elif not sp.issparse(adata.X) and hasattr(adata.X, 'dtype') and adata.X.dtype != np.float32:
        print(f"  Casting X from {adata.X.dtype} to float32")
        adata.X = adata.X.astype(np.float32)

    adata.raw = adata.copy()

    # A1: feature_is_filtered in var only (AFTER raw copy — CxG prohibits it in raw.var)
    adata.var["feature_is_filtered"] = False

    study_info = registry["datasets"].get(study, {})
    adata.uns["title"] = study_info.get("citation", f"{study} source dataset")
    populate_dataset_metadata(adata, study, repo_root)

    # 6b. Populate HCA obs fields (B1 uns->obs, B2 derived, B3 library)
    print(f"\n[6b] Populating HCA obs fields...")
    adata.obs = populate_hca_obs_fields(adata, study, repo_root)

    # Validate
    print(f"\n--- Validation ---")
    issues = validate_h5ad(adata, study, registry)
    if issues:
        for issue in issues:
            print(f"  WARNING: {issue}")
    else:
        print(f"  All checks passed")

    # Provenance
    print(f"\n--- Provenance ---")
    prov_inputs = {
        "counts": intermediates_dir / "counts.mtx.gz",
        "features": intermediates_dir / "features.tsv.gz",
        "barcodes": intermediates_dir / "barcodes.tsv.gz",
        "metadata": published_dir / "metadata.csv",
        "umap": published_dir / "umap_coordinates.csv",
    }
    emb_joint_path = published_dir / "embedding_joint.csv"
    if emb_joint_path.exists():
        prov_inputs["embedding_joint"] = emb_joint_path
    emb_native_path = published_dir / "embedding_native.csv"
    if emb_native_path.exists():
        prov_inputs["embedding_native"] = emb_native_path
    prov_config = {
        "registry": repo_root / "publication/config/source_dataset_registry.yaml",
        "dataset_metadata": repo_root / "publication/config/dataset_metadata.yaml",
        "gene_mapping": repo_root / "publication/mappings/gene_symbol_to_ensembl_full.tsv",
        "hancestro": repo_root / "publication/mappings/hancestro_ethnicity_mapping.tsv",
        "efo": repo_root / "publication/mappings/efo_assay_mapping.tsv",
    }
    sidecar = output_path.parent / f"{output_path.stem}.manifest.yaml"
    build_provenance_manifest(
        adata, output_path, study, repo_root,
        inputs=prov_inputs, config_files=prov_config,
        build_script="assemble_h5ad.py",
        extra={"mapping_stats": mapping_stats},
        sidecar_path=sidecar,
    )

    # Write
    print(f"\nWriting h5ad...")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)
    size_gb = output_path.stat().st_size / (1024 ** 3)
    print(f"  Written: {output_path} ({size_gb:.2f} GB)")

    return adata, issues, mapping_stats


# ---------------------------------------------------------------------------
# Assembly: Pal consolidated
# ---------------------------------------------------------------------------


def assemble_pal(repo_root, gene_mapping, efo_mapping, hancestro_mapping,
                 registry, intermediates_dir=None, published_dir=None,
                 output_path=None, cxg_approved=None, target="hca",
                 config=None):
    """Assemble consolidated pal2021.h5ad from 3 sub-studies."""

    pal_sub_studies = get_pal_sub_studies(config) if config else ["pal_norm_epi", "pal_norm_total", "pal_norm_b1"]

    if output_path is None:
        if config:
            output_path = get_study_output(repo_root, "pal", config)
        else:
            output_path = (
                repo_root / "publication/outputs/source_datasets/pal2021.h5ad"
            )

    print(f"\n{'=' * 70}")
    print(f"ASSEMBLING: pal (consolidated from {len(pal_sub_studies)} sub-studies)")
    print(f"  Output: {output_path}")
    print(f"{'=' * 70}")

    all_mats = []
    all_obs = []
    all_umaps = []
    all_emb_joints = []
    all_emb_natives = []
    all_gene_names = []

    # Get per-sub-study paths from config or derive from defaults
    if config:
        sub_inputs = get_study_inputs(repo_root, "pal", config)
    else:
        sub_inputs = None

    for i, sub in enumerate(pal_sub_studies):
        print(f"\n--- Sub-study: {sub} ---")

        if sub_inputs is not None:
            sub_intermediates = sub_inputs[i]["intermediates"]
            sub_published = sub_inputs[i]["published"]
        elif intermediates_dir is not None:
            sub_intermediates = intermediates_dir / sub
            sub_published = published_dir / sub if published_dir is not None else (
                repo_root / "external_studies/outputs" / sub / "published"
            )
        else:
            sub_intermediates = (
                repo_root / "publication/outputs/source_datasets"
                / "intermediates" / sub
            )
            sub_published = (
                repo_root / "external_studies/outputs" / sub / "published"
            )

        # Load counts
        print(f"  Loading counts...")
        mat, gene_names, cell_ids = load_count_matrix(sub_intermediates)
        print(f"  Shape: {mat.shape[0]:,} cells x {mat.shape[1]:,} genes")

        # Load metadata
        print(f"  Loading metadata...")
        obs = load_metadata(sub_published, cell_ids, fallback_dir=sub_intermediates)
        obs["sub_study"] = sub
        print(f"  Cells: {obs.shape[0]:,}")

        # Load UMAP
        print(f"  Loading UMAP...")
        umap = load_umap(sub_published, "native", cell_ids)
        if umap is not None:
            print(f"  UMAP shape: {umap.shape}")
        else:
            print(f"  No UMAP available for {sub}")

        # Load embeddings
        emb_joint = load_embedding(sub_published, "joint", cell_ids)
        emb_native = load_embedding(sub_published, "native", cell_ids)

        all_mats.append(mat)
        all_obs.append(obs)
        all_umaps.append(umap)
        all_emb_joints.append(emb_joint)
        all_emb_natives.append(emb_native)
        all_gene_names.append(gene_names)

    # --- Gene intersection across sub-studies ---
    # Sub-studies may have different gene sets (filtered independently).
    # Use intersection to avoid fabricating zero counts for truly missing genes.
    gene_sets = [set(g) for g in all_gene_names]
    common_genes = sorted(gene_sets[0].intersection(*gene_sets[1:]))
    print(f"\n--- Gene intersection ---")
    for i, sub in enumerate(pal_sub_studies):
        n_orig = len(all_gene_names[i])
        n_lost = n_orig - len(common_genes)
        print(f"  {sub}: {n_orig:,} genes, {n_lost:,} excluded")
    print(f"  Common genes: {len(common_genes):,}")
    reference_genes = common_genes

    # Subset each matrix to common genes (maintaining consistent column order)
    for i in range(len(all_mats)):
        gene_to_idx = {g: j for j, g in enumerate(all_gene_names[i])}
        col_indices = [gene_to_idx[g] for g in common_genes]
        all_mats[i] = all_mats[i][:, col_indices]

    # --- Consolidation ---
    print(f"\n--- Consolidation ---")
    combined_mat = sp.vstack(all_mats, format="csr")
    combined_obs = pd.concat(all_obs, axis=0)

    # UMAP: only if all sub-studies had it
    has_umap = all(u is not None for u in all_umaps)
    combined_umap = np.vstack(all_umaps) if has_umap else None
    # Embeddings: only if all sub-studies had them
    has_emb_joint = all(e is not None for e in all_emb_joints)
    combined_emb_joint = np.vstack(all_emb_joints) if has_emb_joint else None
    has_emb_native = all(e is not None for e in all_emb_natives)
    combined_emb_native = np.vstack(all_emb_natives) if has_emb_native else None
    print(f"  Combined: {combined_mat.shape[0]:,} cells x {combined_mat.shape[1]:,} genes")
    if not has_umap:
        print(f"  UMAP: not available (Pal sub-studies lack UMAP coordinates)")

    # Check for duplicate cell IDs
    dup_mask = combined_obs.index.duplicated(keep=False)
    if dup_mask.any():
        n_dup_total = dup_mask.sum()
        n_unique_dup = combined_obs.index[dup_mask].nunique()
        print(f"  WARNING: {n_unique_dup} cell IDs appear in multiple sub-studies "
              f"({n_dup_total} total rows)")
        print(f"  Keeping first occurrence only")
        keep_mask = ~combined_obs.index.duplicated(keep="first")
        # keep_mask may be numpy array (from Index.duplicated), not Series
        keep_arr = np.asarray(keep_mask)
        combined_mat = combined_mat[keep_arr]
        combined_obs = combined_obs[keep_arr]
        if combined_umap is not None:
            combined_umap = combined_umap[keep_arr]
        if combined_emb_joint is not None:
            combined_emb_joint = combined_emb_joint[keep_arr]
        if combined_emb_native is not None:
            combined_emb_native = combined_emb_native[keep_arr]
        print(f"  After dedup: {combined_mat.shape[0]:,} cells")

    # Map gene IDs
    print(f"\n  Mapping gene symbols -> Ensembl IDs...")
    ensembl_ids, mapping_stats = map_gene_ids(reference_genes, gene_mapping)
    print(f"  Mapped: {mapping_stats['mapped']:,}/{mapping_stats['total']:,} "
          f"({mapping_stats['pct_mapped']:.1f}%)")

    # Populate CxG fields
    print(f"\n  Populating CxG ontology fields...")
    combined_obs = populate_cxg_fields(
        combined_obs, "pal", efo_mapping, hancestro_mapping, repo_root,
        target=target, config=config
    )

    # Enrich donor metadata from L1 harmonized CSV
    print(f"\n  Enriching donor metadata from L1 harmonized CSV...")
    combined_obs = enrich_donor_metadata(combined_obs, "pal", repo_root, config=config)

    # Build AnnData
    print(f"\n  Building AnnData object...")
    var = pd.DataFrame({"gene_symbol": reference_genes}, index=ensembl_ids)
    var.index.name = None

    # Filter genes against CxG approved set (GENCODE v44 / Ensembl 110)
    if cxg_approved is not None:
        approved_mask = var.index.isin(cxg_approved)
        if not approved_mask.all():
            n_drop = (~approved_mask).sum()
            dropped = var.index[~approved_mask].tolist()
            print(f"  Dropping {n_drop} genes not in CxG approved set")
            print(f"  Examples: {dropped[:20]}")
            keep = np.asarray(approved_mask)
            combined_mat = combined_mat[:, keep]
            var = var[keep]

    # Merge duplicate Ensembl IDs by summing counts
    combined_mat, var = _merge_duplicate_genes(combined_mat, var)

    adata = ad.AnnData(X=combined_mat, obs=combined_obs, var=var)
    if combined_umap is not None:
        adata.obsm["X_umap"] = combined_umap
    if combined_emb_joint is not None:
        adata.obsm["X_scVI_joint"] = combined_emb_joint
    if combined_emb_native is not None:
        adata.obsm["X_scVI_native"] = combined_emb_native

    # CxG requires raw.X to be float32
    if sp.issparse(adata.X) and adata.X.dtype != np.float32:
        print(f"  Casting X from {adata.X.dtype} to float32")
        adata.X = adata.X.astype(np.float32)
    elif not sp.issparse(adata.X) and hasattr(adata.X, 'dtype') and adata.X.dtype != np.float32:
        print(f"  Casting X from {adata.X.dtype} to float32")
        adata.X = adata.X.astype(np.float32)

    adata.raw = adata.copy()

    # A1: feature_is_filtered in var only (AFTER raw copy — CxG prohibits it in raw.var)
    adata.var["feature_is_filtered"] = False

    study_info = registry["datasets"].get("pal", {})
    adata.uns["title"] = study_info.get("citation", "Pal et al. 2021 source dataset")
    populate_dataset_metadata(adata, "pal", repo_root)

    # Populate HCA obs fields (B1 uns->obs, B2 derived, B3 library)
    print(f"\n  Populating HCA obs fields...")
    adata.obs = populate_hca_obs_fields(adata, "pal", repo_root)

    # Validate
    print(f"\n--- Validation ---")
    issues = validate_h5ad(adata, "pal", registry)
    if issues:
        for issue in issues:
            print(f"  WARNING: {issue}")
    else:
        print(f"  All checks passed")

    # Provenance
    print(f"\n--- Provenance ---")
    prov_inputs = {}
    for idx, sub in enumerate(pal_sub_studies):
        if sub_inputs is not None:
            sub_int = sub_inputs[idx]["intermediates"]
            sub_pub = sub_inputs[idx]["published"]
        elif intermediates_dir is not None:
            sub_int = intermediates_dir / sub
            sub_pub = published_dir / sub if published_dir is not None else (
                repo_root / "external_studies/outputs" / sub / "published"
            )
        else:
            sub_int = (
                repo_root / "publication/outputs/source_datasets"
                / "intermediates" / sub
            )
            sub_pub = repo_root / "external_studies/outputs" / sub / "published"
        prov_inputs[f"{sub}_counts"] = sub_int / "counts.mtx.gz"
        prov_inputs[f"{sub}_features"] = sub_int / "features.tsv.gz"
        prov_inputs[f"{sub}_barcodes"] = sub_int / "barcodes.tsv.gz"
        prov_inputs[f"{sub}_metadata"] = sub_pub / "metadata.csv"
    prov_config = {
        "registry": repo_root / "publication/config/source_dataset_registry.yaml",
        "dataset_metadata": repo_root / "publication/config/dataset_metadata.yaml",
        "gene_mapping": repo_root / "publication/mappings/gene_symbol_to_ensembl_full.tsv",
        "hancestro": repo_root / "publication/mappings/hancestro_ethnicity_mapping.tsv",
        "efo": repo_root / "publication/mappings/efo_assay_mapping.tsv",
    }
    sidecar = output_path.parent / f"{output_path.stem}.manifest.yaml"
    build_provenance_manifest(
        adata, output_path, "pal", repo_root,
        inputs=prov_inputs, config_files=prov_config,
        build_script="assemble_h5ad.py",
        extra={"mapping_stats": mapping_stats},
        sidecar_path=sidecar,
    )

    # Write
    print(f"\nWriting h5ad...")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)
    size_gb = output_path.stat().st_size / (1024 ** 3)
    print(f"  Written: {output_path} ({size_gb:.2f} GB)")

    return adata, issues, mapping_stats


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _merge_duplicate_genes(mat, var):
    """Merge duplicate Ensembl IDs by summing count vectors.

    Uses sparse matrix multiplication for memory efficiency:
    builds a grouping matrix G (n_orig x n_unique) and computes mat @ G.
    """
    if not var.index.duplicated().any():
        return mat, var

    n_dup = var.index.duplicated().sum()
    print(f"  Merging {n_dup} duplicate Ensembl IDs by summing counts")

    unique_ids = var.index.unique()
    id_to_group = {eid: i for i, eid in enumerate(unique_ids)}
    group_indices = np.array([id_to_group[eid] for eid in var.index])

    n_orig = len(var.index)
    n_unique = len(unique_ids)

    # Grouping matrix: (n_orig x n_unique), entry [j, k] = 1 if gene j maps to group k
    group_mat = sp.csc_matrix(
        (np.ones(n_orig, dtype=np.float32), (np.arange(n_orig), group_indices)),
        shape=(n_orig, n_unique)
    )
    new_mat = sp.csr_matrix(mat @ group_mat)

    new_var = pd.DataFrame(index=unique_ids)
    new_var.index.name = None

    return new_mat, new_var


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Assemble CxG-compliant h5ad from extracted intermediates"
    )
    parser.add_argument(
        "--study", required=True,
        help="Study name: gray, kumar, murrow, nee, twigger, reed, or pal"
    )
    parser.add_argument(
        "--repo-root", required=True,
        help="Path to iHBCAv1_upload repo root"
    )
    parser.add_argument(
        "--intermediates-dir",
        help="Override auto-detected intermediates directory (Pal: expects sub-study subdirs)"
    )
    parser.add_argument(
        "--published-dir",
        help="Override auto-detected published directory (Pal: expects sub-study subdirs)"
    )
    parser.add_argument(
        "--output",
        help="Override auto-detected output h5ad path"
    )
    parser.add_argument(
        "--target", choices=["hca", "cxg"], default="hca",
        help="Target schema: hca (default) uses :0601/:0602 HANCESTRO terms; "
             "cxg downgrades to :0004 ancestry terms for CxG 5.3.2 compatibility"
    )
    args = parser.parse_args()

    repo_root = Path(args.repo_root)
    study = args.study

    # Load pipeline config
    config = load_pipeline_config(repo_root)

    # Load shared resources
    print("=" * 70)
    print(f"CxG H5AD ASSEMBLY — {study}")
    print("=" * 70)
    print("\nLoading shared resources...")
    registry = load_registry(repo_root, config=config)
    gene_mapping = load_gene_mapping(repo_root, config=config)
    hancestro_mapping = load_hancestro_mapping(repo_root, config=config)
    efo_mapping = load_efo_assay_mapping(repo_root, config=config)
    print(f"  Registry:  {len(registry['datasets'])} datasets")
    print(f"  Gene map:  {len(gene_mapping):,} symbols")
    print(f"  HANCESTRO: {len(hancestro_mapping)} ethnicity entries")
    print(f"  EFO assay: {len(efo_mapping)} study entries")

    # Load CxG approved gene set (GENCODE v44 / Ensembl 110)
    cxg_approved = load_cxg_approved_genes(repo_root, config=config)
    if cxg_approved is not None:
        print(f"  CxG genes: {len(cxg_approved):,} approved Ensembl IDs")

    # Optional path overrides
    output_path = Path(args.output) if args.output else None
    intermediates_dir = Path(args.intermediates_dir) if args.intermediates_dir else None
    published_dir = Path(args.published_dir) if args.published_dir else None

    target = args.target
    print(f"  Target:    {target}")

    # Assemble
    if study == "pal":
        adata, issues, stats = assemble_pal(
            repo_root, gene_mapping, efo_mapping, hancestro_mapping, registry,
            intermediates_dir=intermediates_dir,
            published_dir=published_dir,
            output_path=output_path,
            cxg_approved=cxg_approved,
            target=target,
            config=config,
        )
    else:
        if study not in registry["datasets"]:
            print(f"FATAL: Study '{study}' not found in registry. "
                  f"Known: {list(registry['datasets'].keys())}")
            sys.exit(1)

        adata, issues, stats = assemble_single_study(
            study, repo_root, gene_mapping, efo_mapping, hancestro_mapping,
            registry,
            intermediates_dir=intermediates_dir,
            published_dir=published_dir,
            output_path=output_path,
            cxg_approved=cxg_approved,
            target=target,
            config=config,
        )

    # Final summary
    print(f"\n{'=' * 70}")
    print(f"ASSEMBLY COMPLETE: {study}")
    print(f"  Cells:        {adata.n_obs:,}")
    print(f"  Genes:        {adata.n_vars:,}")
    print(f"  Ensembl map:  {stats['pct_mapped']:.1f}%")
    print(f"  Tier 1 issues: {len([i for i in issues if 'Missing Tier 1' in i])}")
    print(f"  Total issues:  {len(issues)}")
    print(f"{'=' * 70}")

    # Non-zero exit if critical issues (missing Tier 1 fields)
    if any("Missing Tier 1" in i for i in issues):
        print("\nFATAL: Missing required Tier 1 fields — h5ad is non-compliant")
        sys.exit(1)

    sys.exit(0)


if __name__ == "__main__":
    main()
