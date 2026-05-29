#!/usr/bin/env python3
"""
assemble_integrated.py - Assemble iHBCA integrated h5ad
=======================================================
Two build modes:

  AUTHOR SHARE MODE (preferred):
    Constructs AnnData from author share primitives — every field has clean
    provenance. Uses NPZ counts, gene_data.csv, annotations CSV, and
    embedding CSVs directly.

  LEGACY MODE:
    Reads the published integration_iHBCA.h5ad (stale CxG snapshot), overlays
    author share annotations and embeddings. Retained for backward compatibility.

Both modes enrich ethnicity via HANCESTRO mapping, populate CxG/HCA Tier 1 fields,
broadcast donor metadata, and write all-breast-cells + per-lineage h5ads.

Usage (author share):
  python assemble_integrated.py \
    --gene-data /path/to/gene_data.csv \
    --counts-npz /path/to/counts_matrix.npz \
    --umap /path/to/X_scVI100_UMAP.csv \
    --annotations /path/to/ihbca_level1.5_annotations.csv \
    --embeddings /path/to/X_scVI100.csv \
    --repo-root /path/to/iHBCAv1_upload \
    --output-dir /path/to/outputs/integrated_objects

Usage (legacy):
  python assemble_integrated.py \
    --h5ad /path/to/integration_iHBCA.h5ad \
    --annotations /path/to/ihbca_level1.5_annotations.csv \
    --embeddings /path/to/X_scVI100.csv \
    --repo-root /path/to/iHBCAv1_upload \
    --output-dir /path/to/outputs/integrated_objects

"""

# Two million cells walk into a bar. The bartender says "we're going to need a bigger RAM."

import argparse
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from provenance import build_provenance_manifest

import gc

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse
import yaml

from ihbca.constants import (
    ANNOTATION_COLUMNS_TO_MERGE,
    FACS_TO_ENRICHMENT as FACS_TO_ENRICHMENT_INT,
    L1_NUMERIC,
    L1_OBS_COLUMNS,
    LINEAGE_SPLITS,
    PRESERVATION_NORMALIZE,
    RESERVED_OBS_RENAMES,
    RESERVED_VAR_RENAMES,
    SUPERSEDED_COLUMNS,
    TIER1_FIELDS,
    UNS_TO_OBS_FIELDS_INT,
    get_dataset_to_study,
)
from ihbca.ethnicity import downgrade_hancestro_terms, normalize_multi_ethnicity
from ihbca.validation import (
    validate_ensembl_coverage,
    validate_schema_version,
    validate_tier1_fields,
    validate_umap,
)
from ihbca.loaders import (
    load_dataset_metadata,
    load_donor_translations,
    load_hancestro_mapping,
    load_l1_metadata,
    load_pipeline_config,
)
from ihbca.hca_fields import (
    collect_study_accessions,
    collect_study_pis,
    get_study_from_donor_id,
    populate_derived_obs_fields,
)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def load_annotations(annotations_path):
    """Load Austin's level1.5 annotations CSV.

    Returns DataFrame indexed by cellID.
    """
    print(f"Loading annotations: {annotations_path}")
    t0 = time.time()
    ann = pd.read_csv(annotations_path, dtype=str)
    # Use cellID as index for alignment with h5ad obs_names
    ann = ann.set_index("cellID")
    # Drop duplicate cellID.1 column
    if "cellID.1" in ann.columns:
        ann = ann.drop(columns=["cellID.1"])
    print(f"  {len(ann):,} cells x {len(ann.columns)} columns in {time.time() - t0:.1f}s")
    return ann


def remap_mh0023_cell_ids(ann):
    """Remap MH0023_mix_BR1_ prefixes to MH0023_mix_ to match CxG h5ad naming.

    The CxG h5ad has two MH0023 cell ID prefixes:
      - MH0023_mix_BARCODE-1     (5,019 cells, donor_id Pal_MH0023_p23)
      - MH0023_epi_BARCODE-1     (2,831 cells, donor_id Pal_MH0023_p123)

    The annotations CSV has three MH0023 prefixes:
      - MH0023_mix_BARCODE-1     (6,498 cells — catch-all, only 58 match h5ad)
      - MH0023_mix_BR1_BARCODE-1 (5,019 cells — maps 1:1 to h5ad MH0023_mix_)
      - MH0023_epi_BARCODE-1     (2,831 cells — already matches h5ad MH0023_epi_)

    Only MH0023_mix_BR1_ needs remapping. MH0023_epi_ already matches the h5ad
    and must NOT be remapped (doing so would destroy existing correct matches).

    If remapping creates duplicate cell IDs (barcode overlap between
    MH0023_mix_ catch-all and remapped MH0023_mix_BR1_), deduplicate keeping
    the last occurrence (BR1 sub-batch data preferred).

    Returns:
        (remapped_ann, keep_mask): remapped DataFrame and boolean mask for
        synchronized deduplication of row-aligned data (embeddings CSV).
    """
    print("Remapping MH0023 sub-batch cell IDs...")

    original_len = len(ann)
    new_index = ann.index.str.replace(
        r'^MH0023_mix_BR1_', 'MH0023_mix_', regex=True
    )

    n_remapped = (new_index != ann.index).sum()
    print(f"  Remapped {n_remapped:,} cell IDs (MH0023_mix_BR1_ -> MH0023_mix_)")

    ann.index = new_index

    # Check for duplicate cell IDs after remap
    dup_mask = ann.index.duplicated(keep='last')
    n_dups = dup_mask.sum()

    if n_dups > 0:
        print(f"  Deduplicating: {n_dups:,} duplicate cell IDs (keeping BR1 sub-batch)")
        keep_mask = ~dup_mask
        ann = ann[keep_mask]
        print(f"  Annotations: {original_len:,} -> {len(ann):,} rows")
    else:
        keep_mask = np.ones(original_len, dtype=bool)
        print("  No duplicate cell IDs after remap")

    return ann, keep_mask


def load_embeddings(embeddings_path, ann_index, keep_mask=None):
    """Load scVI 100D embeddings CSV, indexed by annotations cellID.

    The CSV has integer row index (not cell IDs) and row order matches the
    annotations CSV. We attach annotations cell IDs to enable h5ad alignment.

    Args:
        embeddings_path: Path to X_scVI100.csv
        ann_index: pandas Index of annotations cellIDs (same row order as CSV,
            after any remapping/deduplication)
        keep_mask: Optional boolean array to drop rows synchronized with
            annotations deduplication (e.g., from remap_mh0023_cell_ids).
            Applied after loading, before index assignment.

    Returns:
        DataFrame with cellID index and scVI_0..scVI_99 columns, float32.
    """
    print(f"Loading embeddings: {embeddings_path}")
    t0 = time.time()

    emb = pd.read_csv(embeddings_path)

    # Drop the Unnamed: 0 integer index column
    if "Unnamed: 0" in emb.columns:
        emb = emb.drop(columns=["Unnamed: 0"])

    # Apply deduplication mask if provided (synchronized with annotations remap)
    if keep_mask is not None:
        if len(emb) != len(keep_mask):
            raise ValueError(
                f"Embedding rows ({len(emb)}) != keep_mask length ({len(keep_mask)}). "
                "Row alignment with annotations broken."
            )
        n_dropped = (~keep_mask).sum()
        if n_dropped > 0:
            emb = emb[keep_mask].reset_index(drop=True)
            print(f"  Dropped {n_dropped:,} rows (synchronized with annotations dedup)")

    # Verify shape matches annotations
    if len(emb) != len(ann_index):
        raise ValueError(
            f"Embedding rows ({len(emb)}) != annotation rows ({len(ann_index)}). "
            "Row order alignment assumption broken."
        )

    # Assign annotations cell IDs as index
    emb.index = ann_index

    # Convert to float32 for memory efficiency
    emb = emb.astype(np.float32)

    print(f"  {len(emb):,} cells x {emb.shape[1]} dims in {time.time() - t0:.1f}s")
    return emb


def load_umap(umap_path, ann_index, keep_mask=None):
    """Load UMAP CSV (2D), indexed by annotations cellID.

    Same row-alignment logic as load_embeddings — CSV rows match annotations CSV.
    """
    print(f"Loading UMAP: {umap_path}")
    t0 = time.time()
    umap = pd.read_csv(umap_path)

    if "Unnamed: 0" in umap.columns:
        umap = umap.drop(columns=["Unnamed: 0"])

    if keep_mask is not None:
        if len(umap) != len(keep_mask):
            raise ValueError(
                f"UMAP rows ({len(umap)}) != keep_mask length ({len(keep_mask)})"
            )
        n_dropped = (~keep_mask).sum()
        if n_dropped > 0:
            umap = umap[keep_mask].reset_index(drop=True)
            print(f"  Dropped {n_dropped:,} rows (synchronized with annotations dedup)")

    if len(umap) != len(ann_index):
        raise ValueError(
            f"UMAP rows ({len(umap)}) != annotation rows ({len(ann_index)})"
        )

    umap.index = ann_index
    umap = umap.astype(np.float32)
    print(f"  {len(umap):,} cells x {umap.shape[1]} dims in {time.time() - t0:.1f}s")
    return umap


def construct_from_author_share(
    annotations_path, gene_data_path, embeddings_path, umap_path,
    counts_npz_path=None,
):
    """Construct integrated AnnData from author share primitives.

    Replaces the CxG h5ad base with direct author share files:
    - Cell metadata from annotations CSV
    - Gene metadata from gene_data.csv (scVI output)
    - Embeddings from X_scVI100.csv and X_scVI100_UMAP.csv
    - Counts from NPZ (optional — not needed for DA pipeline)

    Returns (adata, ann) where ann is the full annotations DataFrame.
    """
    print("Constructing AnnData from author share primitives...")
    t_total = time.time()

    # 1. Load annotations (cell metadata — 67 columns, 2.12M cells)
    ann = load_annotations(annotations_path)
    ann, keep_mask = remap_mh0023_cell_ids(ann)

    # Construct donor_id from annotations columns (matches CxG h5ad format)
    # CxG format: "{Study}_{patientID}" e.g. "Gray_HBCA_Donor_1"
    # Annotations CSV has lowercase dataset (gray, kumar, ...) — capitalize to match CxG.
    if "donor_id" not in ann.columns and "dataset" in ann.columns and "patientID" in ann.columns:
        ann["donor_id"] = ann["dataset"].str.capitalize() + "_" + ann["patientID"].astype(str)
        print(f"  Constructed donor_id: {ann['donor_id'].nunique()} unique donors")

    # 2. Load gene metadata
    print(f"Loading gene data: {gene_data_path}")
    gene_data = pd.read_csv(gene_data_path, index_col=0)
    print(f"  {len(gene_data):,} genes x {len(gene_data.columns)} columns")

    # 3. Load counts NPZ (optional)
    X = None
    if counts_npz_path and os.path.exists(counts_npz_path):
        print(f"Loading counts NPZ: {counts_npz_path}")
        t0 = time.time()
        X = scipy.sparse.load_npz(counts_npz_path)
        print(f"  Shape: {X.shape}, format: {X.format}, "
              f"dtype: {X.dtype} in {time.time() - t0:.1f}s")

        # Apply keep_mask for MH0023 dedup if needed
        if keep_mask is not None and (~keep_mask).sum() > 0:
            X = X[keep_mask]
            print(f"  After dedup: {X.shape}")

        # Verify dimensions
        if X.shape[0] != len(ann):
            raise ValueError(
                f"NPZ rows ({X.shape[0]}) != annotation rows ({len(ann)}). "
                "Cell alignment broken."
            )
        if X.shape[1] != len(gene_data):
            raise ValueError(
                f"NPZ cols ({X.shape[1]}) != gene_data rows ({len(gene_data)}). "
                "Gene alignment broken."
            )

        # Ensure CSR format
        if not scipy.sparse.issparse(X) or X.format != "csr":
            X = scipy.sparse.csr_matrix(X)
    else:
        if counts_npz_path:
            print(f"  WARNING: Counts NPZ not found: {counts_npz_path}")
        print("  Constructing without counts (X will be empty sparse)")
        X = scipy.sparse.csr_matrix((len(ann), len(gene_data)), dtype=np.float32)

    # 4. Build var DataFrame
    var = gene_data.copy()
    # Rename reserved var columns (same logic as legacy mode)
    var_renames = {k: v for k, v in RESERVED_VAR_RENAMES.items()
                   if k in var.columns}
    if var_renames:
        var = var.rename(columns=var_renames)
        print(f"  Renamed {len(var_renames)} reserved var columns")
    var["feature_is_filtered"] = False

    # 5. Construct AnnData
    adata = ad.AnnData(X=X, obs=ann, var=var)
    print(f"  AnnData: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # 6. Load embeddings and UMAP (before count-vector dedup so keep_mask aligns)
    emb = load_embeddings(embeddings_path, ann.index, keep_mask=keep_mask)
    umap_df = load_umap(umap_path, ann.index, keep_mask=keep_mask)

    adata.obsm["X_scvi_100"] = emb.reindex(adata.obs_names).values
    adata.obsm["X_umap"] = umap_df.reindex(adata.obs_names).values

    # 6b. Count-vector deduplication
    # Austin's NPZ contains Pal Normal controls in both NormB1Total and NormTotal
    # sub-studies (different cell IDs, identical count vectors). Austin filters
    # these in 10e_pal_data_preparation.R:112 before merging, but the raw NPZ
    # retains both copies. Deduplicate by sparse row hash to catch any identical
    # count vectors regardless of source.
    if scipy.sparse.issparse(adata.X):
        X_csr = adata.X.tocsr() if adata.X.format != "csr" else adata.X
        print("  Deduplicating count vectors...")
        t0 = time.time()
        seen_hashes = {}
        dup_indices = []
        for i in range(X_csr.shape[0]):
            start, end = X_csr.indptr[i], X_csr.indptr[i + 1]
            row_key = hash((tuple(X_csr.indices[start:end]),
                            tuple(X_csr.data[start:end])))
            if row_key in seen_hashes:
                dup_indices.append(i)
            else:
                seen_hashes[row_key] = i
        if dup_indices:
            keep = np.ones(adata.n_obs, dtype=bool)
            keep[dup_indices] = False
            adata = adata[keep].copy()
            ann = ann.loc[adata.obs_names]
            print(f"  Removed {len(dup_indices):,} duplicate count vectors "
                  f"({adata.n_obs:,} cells remain) in {time.time() - t0:.1f}s")
        else:
            print(f"  No duplicate count vectors found in {time.time() - t0:.1f}s")
    print(f"  obsm keys: {list(adata.obsm.keys())}")

    elapsed = time.time() - t_total
    print(f"  Author share construction complete in {elapsed:.1f}s")

    return adata, ann


# ---------------------------------------------------------------------------
# Tier 1 Constants (for author share mode)
# ---------------------------------------------------------------------------


def _map_age_to_hsapdv(age):
    """Map numeric age to the most specific non-deprecated HsapDv term.

    The old hierarchy (0083-0093) is entirely deprecated. Current terms
    (verified via EBI OLS, ontology version 2025-01-23):
      15-40: HsapDv:0000266  (young adult stage)
      40-60: HsapDv:0000267  (middle aged stage)
      60+:   HsapDv:0000227  (late adult stage)
      <15:   HsapDv:0000258  (adult stage — broadest non-deprecated adult term)
    Fallback: HsapDv:0000258 (adult stage, 15+)
    """
    try:
        age = float(age)
    except (ValueError, TypeError):
        return "unknown"
    if age < 40:
        return "HsapDv:0000266"
    elif age < 60:
        return "HsapDv:0000267"
    else:
        return "HsapDv:0000227"


def set_tier1_constants(adata):
    """Set constant Tier 1 fields that don't vary per cell or per study.

    In legacy (CxG h5ad) mode, these are inherited from the published object.
    In author share mode, they must be explicitly set since the annotations CSV
    doesn't include CxG schema fields.

    Per-study fields (assay, disease) are handled by populate_hca_obs_fields_integrated().
    """
    print("Setting Tier 1 constants (author share mode)...")

    constants = {
        "organism_ontology_term_id": "NCBITaxon:9606",     # Homo sapiens
        "tissue_ontology_term_id": "UBERON:0000310",       # breast
        "sex_ontology_term_id": "PATO:0000383",            # female
        "self_reported_ethnicity_ontology_term_id": "unknown",  # enriched later
        "suspension_type": "cell",
        "is_primary_data": False,                           # integrated = not primary
        "cell_type_ontology_term_id": "unknown",            # CL fallback fills later
    }

    for field, value in constants.items():
        adata.obs[field] = value
        print(f"  {field} = {value}")

    # development_stage_ontology_term_id: per-donor age-based mapping
    # HsapDv:0000087 is deprecated — use age-specific terms instead
    age_col = None
    for candidate in ["patient_age", "age", "age_continuous"]:
        if candidate in adata.obs.columns:
            age_col = candidate
            break

    if age_col:
        adata.obs["development_stage_ontology_term_id"] = (
            adata.obs[age_col].map(_map_age_to_hsapdv)
        )
        term_counts = adata.obs["development_stage_ontology_term_id"].value_counts()
        print(f"  development_stage_ontology_term_id mapped from {age_col}:")
        for term, count in term_counts.items():
            print(f"    {term}: {count:,} cells")
    else:
        # Fallback: broadest non-deprecated adult term (covers 15+)
        adata.obs["development_stage_ontology_term_id"] = "HsapDv:0000258"
        print("  development_stage_ontology_term_id = HsapDv:0000258 (fallback, no age column)")

    # donor_id: use from annotations if present, keep as-is
    if "donor_id" not in adata.obs.columns:
        print("  WARNING: donor_id not found in annotations")

    return adata


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------


def merge_annotations(adata, ann, columns_to_merge):
    """Merge annotation columns into h5ad obs by cell ID alignment.

    Left-joins on adata.obs_names: h5ad cells without annotations get NaN.
    """
    available = [c for c in columns_to_merge if c in ann.columns]
    missing = [c for c in columns_to_merge if c not in ann.columns]
    if missing:
        print(f"  WARNING: {len(missing)} requested columns not in annotations: {missing[:5]}...")

    # Already in obs — skip to avoid collision
    skip = [c for c in available if c in adata.obs.columns]
    if skip:
        print(f"  Skipping {len(skip)} columns already in obs: {skip}")
        available = [c for c in available if c not in skip]

    print(f"  Merging {len(available)} annotation columns into obs...")
    t0 = time.time()

    # Reindex annotations to match h5ad obs_names (left join)
    ann_subset = ann[available].reindex(adata.obs_names)

    # Count alignment
    matched = ann_subset.notna().any(axis=1).sum()
    unmatched = len(adata) - matched
    print(f"  Matched: {matched:,} / {len(adata):,} cells ({100 * matched / len(adata):.2f}%)")
    if unmatched > 0:
        print(f"  Unmatched (NaN fill): {unmatched:,} cells")

    # Merge into obs
    for col in available:
        adata.obs[col] = ann_subset[col].values

    print(f"  Done in {time.time() - t0:.1f}s")
    return adata


def add_embeddings(adata, emb):
    """Add scVI 100D embeddings to obsm['X_scvi_100'].

    Aligns embeddings to h5ad obs_names. Cells without embeddings get NaN.
    """
    print("Adding scVI 100D embeddings to obsm['X_scvi_100']...")
    t0 = time.time()

    # Reindex to h5ad obs_names
    emb_aligned = emb.reindex(adata.obs_names)
    n_matched = emb_aligned.notna().any(axis=1).sum()
    n_nan = len(adata) - n_matched

    print(f"  Matched: {n_matched:,} / {len(adata):,}")
    if n_nan > 0:
        print(f"  NaN-filled: {n_nan:,} cells (no embedding available)")
        # Fill NaN with 0 for obsm (numpy array can't have NaN in float32 obsm
        # without issues in downstream tools)
        emb_aligned = emb_aligned.fillna(0.0)

    adata.obsm["X_scvi_100"] = emb_aligned.values.astype(np.float32)
    print(f"  obsm['X_scvi_100'] shape: {adata.obsm['X_scvi_100'].shape}")
    print(f"  Done in {time.time() - t0:.1f}s")
    return adata


def enrich_ethnicity(adata, ann, hancestro_mapping):
    """Update self_reported_ethnicity_ontology_term_id from HANCESTRO mapping.

    Maps per-donor: annotations CSV 'ethnicity' column -> HANCESTRO term_id.
    Only overwrites cells where the existing value is unknown or missing —
    preserves any pre-existing ethnicity from the CxG input h5ad.
    """
    print("Enriching ethnicity via HANCESTRO mapping...")

    if "ethnicity" not in ann.columns:
        print("  WARNING: No 'ethnicity' column in annotations. Skipping.")
        return adata

    # Build per-cell ethnicity from annotations (aligned to h5ad)
    eth_series = ann["ethnicity"].reindex(adata.obs_names)

    # Map verbatim -> HANCESTRO term
    mapped = eth_series.fillna("").map(
        lambda v: hancestro_mapping.get(v, "unknown") if v else "unknown"
    )

    # Convert from Categorical to string if needed (CxG h5ads use Categorical
    # dtype, and pandas refuses .loc[] assignment of values not in the category set)
    for col in ["self_reported_ethnicity_ontology_term_id", "self_reported_ethnicity"]:
        if col in adata.obs.columns and hasattr(adata.obs[col], "cat"):
            adata.obs[col] = adata.obs[col].astype(str)

    # Identify cells where existing ethnicity is unknown or missing
    existing = adata.obs["self_reported_ethnicity_ontology_term_id"]
    is_unknown = existing.isna() | existing.isin(["unknown", "", "nan"])

    n_preserved = (~is_unknown).sum()
    n_enriched = (is_unknown & (mapped != "unknown")).sum()
    n_still_unknown = (is_unknown & (mapped == "unknown")).sum()

    print(f"  Cells with existing ethnicity (preserved): {n_preserved:,}")
    print(f"  Cells enriched from annotations: {n_enriched:,}")
    print(f"  Cells still unknown: {n_still_unknown:,}")

    # Only overwrite unknown/missing cells
    adata.obs.loc[is_unknown, "self_reported_ethnicity_ontology_term_id"] = mapped[is_unknown].values

    # Also update the human-readable column if present (same conditional logic)
    if "self_reported_ethnicity" in adata.obs.columns:
        label_mapping = {}
        for verbatim, term_id in hancestro_mapping.items():
            if term_id != "unknown":
                label_mapping[term_id] = verbatim  # approximate
        mapped_labels = mapped.map(
            lambda v: label_mapping.get(v, "unknown") if v != "unknown" else "unknown"
        )
        adata.obs.loc[is_unknown, "self_reported_ethnicity"] = mapped_labels[is_unknown].values

    return adata


def enrich_ethnicity_from_harmonized(adata, repo_root, hancestro_mapping):
    """Enrich ethnicity from harmonized donor metadata CSV.

    Second-pass enrichment after annotations-based enrich_ethnicity().
    Maps donor_id -> ethnicity_verbatim -> HANCESTRO term_id at the donor
    level, then broadcasts to all cells for that donor.

    Only overwrites cells where existing ethnicity is unknown or missing.
    """
    print("Enriching ethnicity from harmonized donor metadata...")

    harm_path = repo_root / (
        "external_studies/harmonization/outputs/harmonized_metadata/"
        "harmonized_donor_metadata.csv"
    )
    if not harm_path.exists():
        print(f"  WARNING: Harmonized metadata not found: {harm_path}")
        return adata

    harm = pd.read_csv(harm_path, dtype=str)
    print(f"  Loaded {len(harm)} donors from harmonized metadata")

    # Build h5ad_donor_id -> ethnicity_verbatim mapping.
    # H5ad uses {Study}_{ihbca_donor_id} format (e.g., Gray_PM-A, Kumar_P01).
    # Harmonized CSV has study (lowercase) + ihbca_donor_id (raw, may have
    # leading whitespace for Nee Ctrl donors — this whitespace is intentional
    # and matches the h5ad format).
    donor_eth = {}
    for _, row in harm.iterrows():
        study_cap = row["study"].capitalize()
        h5ad_id = f"{study_cap}_{row['ihbca_donor_id']}"
        eth_verbatim = row.get("ethnicity_verbatim", "")
        if pd.notna(eth_verbatim) and str(eth_verbatim).strip():
            donor_eth[h5ad_id] = eth_verbatim

    # Add alternate keys for studies with donor ID mismatches (Reed, Pal)
    # Translation CSV maps L1 key ({Study}_{tissue_bank_id}) -> cxg_donor_id.
    # Pal cxg_donor_ids already have Pal_ prefix (e.g., Pal_MH0064).
    # Reed cxg_donor_ids are bare (e.g., HBCA_Donor_1) — h5ad uses
    # Reed_HBCA_Donor_1, so we must also add the study-prefixed form.
    xlat = load_donor_translations(repo_root)
    n_xlat = 0
    for l1_key, cxg_id in xlat.items():
        if l1_key in donor_eth:
            donor_eth[cxg_id] = donor_eth[l1_key]
            # Add study-prefixed form if cxg_id doesn't already have it
            study_prefix = l1_key.split("_", 1)[0]  # e.g., "Reed"
            if not cxg_id.startswith(f"{study_prefix}_"):
                donor_eth[f"{study_prefix}_{cxg_id}"] = donor_eth[l1_key]
            n_xlat += 1
    if n_xlat:
        print(f"  Added {n_xlat} translated donor keys (Reed/Pal ID bridge)")

    print(f"  Constructed {len(donor_eth)} donor_id -> ethnicity mappings")

    # Check for unmapped verbatim values (not in hancestro_mapping)
    unique_verbatim = set(donor_eth.values())
    unmapped = [v for v in unique_verbatim if v not in hancestro_mapping]
    if unmapped:
        print(f"  WARNING: {len(unmapped)} ethnicity_verbatim values NOT in "
              f"hancestro_mapping (will map to 'unknown'): {unmapped}")

    # Map verbatim -> HANCESTRO term at donor level
    donor_term = {}
    for did, verbatim in donor_eth.items():
        donor_term[did] = hancestro_mapping.get(verbatim, "unknown")

    # Check which h5ad donor_ids matched
    h5ad_donors = set(adata.obs["donor_id"].unique())
    matched = h5ad_donors & set(donor_term.keys())
    print(f"  H5ad donors matched: {len(matched)}/{len(h5ad_donors)}")
    if len(matched) < len(h5ad_donors):
        # Only warn about donors that HAVE ethnicity but didn't match
        # (Reed uses HBCA_Donor_N scheme, Pal uses different IDs — both
        # are expected mismatches with no new ethnicity to add)
        all_constructed = set(donor_eth.keys()) | set(donor_term.keys())
        truly_unmatched = h5ad_donors - all_constructed
        if truly_unmatched:
            print(f"  Unmatched h5ad donors ({len(truly_unmatched)}): "
                  f"{sorted(truly_unmatched)[:10]}...")

    # Build per-cell term series via donor_id broadcast
    cell_terms = adata.obs["donor_id"].map(donor_term).fillna("unknown")

    # Identify cells where existing ethnicity is unknown or missing
    existing = adata.obs["self_reported_ethnicity_ontology_term_id"]
    is_unknown = existing.isna() | existing.isin(["unknown", "", "nan"])

    # Only enrich cells that are currently unknown AND have a non-unknown
    # term from harmonized metadata
    can_enrich = is_unknown & (cell_terms != "unknown")
    n_preserved = (~is_unknown).sum()
    n_enriched = can_enrich.sum()
    n_still_unknown = (is_unknown & ~can_enrich).sum()

    print(f"  Cells with existing ethnicity (preserved): {n_preserved:,}")
    print(f"  Cells enriched from harmonized metadata: {n_enriched:,}")
    print(f"  Cells still unknown: {n_still_unknown:,}")

    # Apply conditional fill
    adata.obs.loc[can_enrich, "self_reported_ethnicity_ontology_term_id"] = (
        cell_terms[can_enrich].values
    )

    # Update human-readable column with verbatim values
    if "self_reported_ethnicity" in adata.obs.columns:
        cell_verbatim = adata.obs["donor_id"].map(donor_eth).fillna("unknown")
        adata.obs.loc[can_enrich, "self_reported_ethnicity"] = (
            cell_verbatim[can_enrich].values
        )

    # Per-study summary
    if "dataset" in adata.obs.columns:
        print("  Per-study breakdown:")
        for study in sorted(adata.obs["dataset"].unique()):
            mask = adata.obs["dataset"] == study
            study_unknown = is_unknown & mask
            study_enriched = can_enrich & mask
            study_still = study_unknown & ~study_enriched
            if study_unknown.sum() > 0:
                print(f"    {study}: {study_enriched.sum():,} enriched, "
                      f"{study_still.sum():,} still unknown")

    return adata


# ---------------------------------------------------------------------------
# L1 harmonized donor metadata enrichment
# ---------------------------------------------------------------------------

# L1_OBS_COLUMNS, L1_NUMERIC, SUPERSEDED_COLUMNS imported from ihbca.constants


def enrich_donor_metadata_integrated(adata, repo_root):
    """Broadcast L1 harmonized donor metadata to per-cell obs.

    Constructs {Study}_{ihbca_donor_id} keys matching the integrated h5ad
    donor_id format, broadcasts 18 donor-level columns to obs, then drops
    superseded annotation columns.

    Logs per-column divergence where annotation values differ from harmonized.
    """
    print("Enriching donor metadata from L1 harmonized CSV...")

    l1 = load_l1_metadata(repo_root)
    if l1 is None:
        return adata
    print(f"  Loaded {len(l1)} donors from L1 staging CSV")

    # Build {Study}_{ihbca_donor_id} -> row mapping
    # Same key construction as enrich_ethnicity_from_harmonized()
    donor_rows = {}
    for _, row in l1.iterrows():
        study_cap = row["study"].capitalize()
        h5ad_id = f"{study_cap}_{row['ihbca_donor_id']}"
        donor_rows[h5ad_id] = row

    # Add alternate keys for studies with donor ID mismatches (Reed, Pal)
    # See enrich_ethnicity_from_harmonized() for format rationale.
    xlat = load_donor_translations(repo_root)
    n_xlat = 0
    for l1_key, cxg_id in xlat.items():
        if l1_key in donor_rows:
            donor_rows[cxg_id] = donor_rows[l1_key]
            # Add study-prefixed form if cxg_id doesn't already have it
            study_prefix = l1_key.split("_", 1)[0]
            if not cxg_id.startswith(f"{study_prefix}_"):
                donor_rows[f"{study_prefix}_{cxg_id}"] = donor_rows[l1_key]
            n_xlat += 1
    if n_xlat:
        print(f"  Added {n_xlat} translated donor keys (Reed/Pal ID bridge)")

    print(f"  Constructed {len(donor_rows)} {{Study}}_{{ihbca_donor_id}} keys")

    # Match against h5ad donor_ids
    h5ad_donors = set(adata.obs["donor_id"].unique())
    matched = h5ad_donors & set(donor_rows.keys())
    print(f"  H5ad donors matched: {len(matched)}/{len(h5ad_donors)}")

    # --- Divergence logging for superseded columns ---
    print("  Overwrite divergence (annotation vs harmonized):")
    for ann_col, l1_col in SUPERSEDED_COLUMNS.items():
        if ann_col not in adata.obs.columns:
            print(f"    {ann_col} → {l1_col}: annotation column not present (skip)")
            continue

        # Build what the harmonized values would be for comparison
        donor_l1_vals = {}
        for did, row in donor_rows.items():
            val = row.get(l1_col, "")
            if pd.notna(val) and str(val).strip():
                donor_l1_vals[did] = str(val).strip()

        harmonized_series = adata.obs["donor_id"].map(donor_l1_vals)
        existing = adata.obs[ann_col].astype(str)
        # Only compare where both have non-null values
        both_present = harmonized_series.notna() & existing.notna() & (existing != "nan")
        if both_present.sum() > 0:
            changed = both_present & (existing != harmonized_series)
            n_changed = changed.sum()
            n_compared = both_present.sum()
            print(f"    {ann_col} → {l1_col}: "
                  f"{n_changed:,}/{n_compared:,} cells changed "
                  f"({100 * n_changed / n_compared:.1f}%)")
        else:
            print(f"    {ann_col} → {l1_col}: no overlapping values to compare")

    # --- Inject L1 columns ---
    for col in L1_OBS_COLUMNS:
        donor_map = {}
        for did, row in donor_rows.items():
            val = row.get(col, None)
            if pd.notna(val) and str(val).strip():
                donor_map[did] = str(val).strip()

        # Convert categorical to object to allow new values
        if col in adata.obs.columns and hasattr(adata.obs[col], "cat"):
            adata.obs[col] = adata.obs[col].astype(object)

        adata.obs[col] = adata.obs["donor_id"].map(donor_map)

        if col in L1_NUMERIC:
            adata.obs[col] = pd.to_numeric(adata.obs[col], errors="coerce")
        else:
            adata.obs[col] = adata.obs[col].fillna("unknown")

    print(f"  Added {len(L1_OBS_COLUMNS)} L1 columns to obs")

    # --- Drop superseded annotation columns ---
    dropped = [c for c in SUPERSEDED_COLUMNS if c in adata.obs.columns]
    if dropped:
        adata.obs = adata.obs.drop(columns=dropped)
        print(f"  Dropped {len(dropped)} superseded annotation columns: {dropped}")

    # Per-study coverage
    if "dataset" in adata.obs.columns:
        print("  Per-study coverage:")
        for study in sorted(adata.obs["dataset"].unique()):
            mask = adata.obs["dataset"] == study
            study_donors = adata.obs.loc[mask, "donor_id"].nunique()
            study_n_matched = len(
                set(adata.obs.loc[mask, "donor_id"].unique()) & set(donor_rows.keys())
            )
            print(f"    {study}: {study_n_matched}/{study_donors} donors")

    return adata


# ---------------------------------------------------------------------------
# HCA obs field population (integrated)
# ---------------------------------------------------------------------------

# UNS_TO_OBS_FIELDS_INT, FACS_TO_ENRICHMENT_INT, DATASET_TO_STUDY
# imported from ihbca.constants


def apply_cl_level15_fallback(adata, repo_root):
    """Apply level1.5_annotation -> CL term fallback for 'unknown' cells.

    Reads the reviewed mapping from level15_to_cl_mapping.csv and resolves
    cells with cell_type_ontology_term_id == 'unknown' using their
    level1.5_annotation label.
    """
    mapping_path = repo_root / "publication/mappings/level15_to_cl_mapping.csv"
    if not mapping_path.exists():
        print("  Level1.5 CL mapping not found — skipping fallback")
        return adata

    l15_df = pd.read_csv(mapping_path)
    if "level15_annotation" not in l15_df.columns or "assigned_cl_term" not in l15_df.columns:
        print("  WARNING: level15_to_cl_mapping.csv missing required columns")
        return adata

    l15_map = dict(zip(l15_df["level15_annotation"], l15_df["assigned_cl_term"]))
    print(f"  Loaded level1.5 -> CL mapping: {len(l15_map)} labels")

    # Find level1.5_annotation column
    l15_col = None
    for candidate in ["level1.5_annotation", "level1.5", "level15_annotation"]:
        if candidate in adata.obs.columns:
            l15_col = candidate
            break

    if l15_col is None:
        print("  WARNING: No level1.5_annotation column in obs — cannot apply fallback")
        return adata

    # Convert categorical to object if needed
    cl_col = "cell_type_ontology_term_id"
    if hasattr(adata.obs[cl_col], "cat"):
        adata.obs[cl_col] = adata.obs[cl_col].astype(object)

    cl_values = adata.obs[cl_col].values.copy()
    l15_values = adata.obs[l15_col].values
    unknown_mask = (cl_values == "unknown")
    n_before = unknown_mask.sum()

    n_resolved = 0
    for label, resolved_term in l15_map.items():
        if resolved_term == "unknown":
            continue
        match = unknown_mask & (l15_values == label)
        n_match = match.sum()
        if n_match > 0:
            cl_values[match] = resolved_term
            n_resolved += n_match

    adata.obs[cl_col] = cl_values
    n_after = (cl_values == "unknown").sum()
    print(f"  Level1.5 CL fallback: {n_before:,} unknown -> {n_after:,} unknown "
          f"({n_resolved:,} resolved)")

    return adata


def populate_hca_obs_fields_integrated(adata, repo_root, config=None):
    """Populate HCA-required obs fields for the integrated object.

    Study-level fields vary per cell. Maps via donor_id prefix -> study name
    -> dataset_metadata.yaml lookup.
    """
    print("Populating HCA obs fields (integrated)...")

    # Load dataset metadata
    dmeta = load_dataset_metadata(repo_root, config=config)
    if not dmeta:
        print(f"  WARNING: dataset_metadata.yaml not found — skipping HCA obs fields")
        return adata
    datasets = dmeta.get("datasets", {})

    # Determine study for each cell via donor_id prefix or dataset column
    if "dataset" in adata.obs.columns:
        # Use dataset column directly (more reliable).
        # Annotations CSV uses lowercase (gray, kumar, ...); CxG h5ad uses capitalized.
        # Try config-driven mapping first, fall back to lowercase identity mapping.
        dataset_to_study = get_dataset_to_study(config) if config else {
            "Gray": "gray", "Kumar": "kumar", "Murrow": "murrow",
            "Nee": "nee", "Twigger": "twigger", "Reed": "reed", "Pal": "pal",
        }
        study_series = adata.obs["dataset"].map(dataset_to_study)
        if study_series.isna().all():
            study_series = adata.obs["dataset"].str.lower()
    else:
        # Fall back to donor_id prefix parsing
        study_series = adata.obs["donor_id"].map(get_study_from_donor_id)

    n_mapped = study_series.notna().sum()
    print(f"  Study assignment: {n_mapped:,}/{len(adata):,} cells mapped")

    # --- B1: Broadcast study-level fields to obs ---
    for field in UNS_TO_OBS_FIELDS_INT:
        study_vals = {s: d.get(field, "unknown") for s, d in datasets.items()}
        adata.obs[field] = study_series.map(study_vals).fillna("unknown")

    # --- B2: Derived constant/lookup fields ---
    populate_derived_obs_fields(adata.obs, FACS_TO_ENRICHMENT_INT, PRESERVATION_NORMALIZE)

    # institute: per-study (integrated needs per-cell mapping)
    institute_map = {s: d.get("institute", "unknown") for s, d in datasets.items()}
    adata.obs["institute"] = study_series.map(institute_map).fillna("unknown")

    # --- B3: Library metadata from SRA run tables ---
    # For integrated: map donor_id -> study -> SRA run table
    for field in ["library_id", "library_sequencing_run", "library_preparation_batch"]:
        adata.obs[field] = "unknown"

    # Try to load per-study SRA tables and assign
    for study_name in datasets.keys():
        sra_path = repo_root / f"publication/mappings/sra_run_table_{study_name}.csv"
        if not sra_path.exists():
            continue
        sra = pd.read_csv(sra_path, dtype=str)
        if "donor_id" not in sra.columns:
            continue

        sra_lookup = sra.set_index("donor_id")
        # Find cells belonging to this study
        study_mask = study_series == study_name
        if study_mask.sum() == 0:
            continue

        study_donors = adata.obs.loc[study_mask, "donor_id"]
        # Strip study prefix for matching ({Study}_{id} -> {id})
        bare_donors = study_donors.str.replace(
            f"^{study_name.capitalize()}_", "", regex=True
        )

        for field in ["library_id", "library_sequencing_run", "library_preparation_batch"]:
            if field in sra_lookup.columns:
                mapped = bare_donors.map(sra_lookup[field].to_dict())
                adata.obs.loc[study_mask, field] = mapped.fillna("unknown").values

    n_fields = len(UNS_TO_OBS_FIELDS_INT) + 8 + 3
    print(f"  Populated {n_fields} HCA obs fields")

    return adata


def fix_tier1_fields(adata, target="hca"):
    """Ensure all 11 Tier 1 fields are present and correct.

    From inspection: only organism_ontology_term_id is missing on obs.
    It exists in uns as NCBITaxon:9606.

    target: 'hca' (default) keeps HCA-native HANCESTRO terms (:0601/:0602).
            'cxg' applies downgrade to :0004 ancestry terms for CxG 5.3.2.
    """
    print("Verifying Tier 1 fields...")

    # organism_ontology_term_id: move from uns to obs
    if "organism_ontology_term_id" not in adata.obs.columns:
        term = adata.uns.get("organism_ontology_term_id", "NCBITaxon:9606")
        adata.obs["organism_ontology_term_id"] = term
        print(f"  Added organism_ontology_term_id = {term} to obs")

    # tissue_type: required by CxG 5.3.2
    adata.obs["tissue_type"] = "tissue"

    # is_primary_data: ensure boolean type
    if "is_primary_data" in adata.obs.columns:
        adata.obs["is_primary_data"] = adata.obs["is_primary_data"].map(
            {"True": False, "False": False, True: False, False: False}
        ).fillna(False).astype(bool)
        # Note: False for integrated object — primary data exists in source datasets

    # A3: Rename reserved obs columns (HCA reserves these label names)
    renames = {k: v for k, v in RESERVED_OBS_RENAMES.items()
               if k in adata.obs.columns}
    if renames:
        adata.obs = adata.obs.rename(columns=renames)
        print(f"  Renamed {len(renames)} reserved obs columns: "
              f"{list(renames.keys())}")

    # Normalize multi-value ethnicity: sorted, ascending order
    # HCA uses ' || ' delimiter, CxG uses ','
    if "self_reported_ethnicity_ontology_term_id" in adata.obs.columns:
        sep = "," if target == "cxg" else " || "
        adata.obs["self_reported_ethnicity_ontology_term_id"] = normalize_multi_ethnicity(
            adata.obs["self_reported_ethnicity_ontology_term_id"], sep=sep
        )

    # CxG 5.3.2 compatibility: downgrade HCA-native HANCESTRO terms to :0004 branch
    # Only applied when --target cxg. HCA builds keep :0601/:0602 terms as-is.
    if target == "cxg" and "self_reported_ethnicity_ontology_term_id" in adata.obs.columns:
        adata.obs["self_reported_ethnicity_ontology_term_id"] = downgrade_hancestro_terms(
            adata.obs["self_reported_ethnicity_ontology_term_id"]
        )
        # Re-normalize after downgrade — term replacement can break sort order
        adata.obs["self_reported_ethnicity_ontology_term_id"] = normalize_multi_ethnicity(
            adata.obs["self_reported_ethnicity_ontology_term_id"], sep=","
        )

    # A10: Drop deprecated ethnicity columns
    for col in ["ethnicity", "ethnicity_ontology_term_id"]:
        if col in adata.obs.columns:
            adata.obs = adata.obs.drop(columns=[col])
            print(f"  Dropped deprecated column: {col}")

    # Check all present
    for field in TIER1_FIELDS:
        if field in adata.obs.columns:
            n_na = adata.obs[field].isna().sum()
            if n_na > 0:
                print(f"  WARNING: {field} has {n_na:,} NaN values")
        else:
            print(f"  MISSING: {field}")

    return adata


def update_schema(adata, version="5.3.2"):
    """Update uns schema fields for target version.

    For 5.3.2: organism_ontology_term_id on obs (already added).
    Keep it in uns too for 7.0.0 compatibility.

    A5: schema_version kept (CxG requires it). schema_reference removed
    (HCA rejects as reserved key, CxG doesn't need it).
    """
    # Remove all CxG-reserved uns keys (validator sets schema_version automatically)
    reserved_uns = ["schema_version", "schema_reference", "citation"]
    for key in reserved_uns:
        if key in adata.uns:
            del adata.uns[key]
            print(f"  Removed reserved uns['{key}']")

    # Remove CxG-reserved obs columns
    reserved_obs = ["observation_joinid"]
    for col in reserved_obs:
        if col in adata.obs.columns:
            adata.obs = adata.obs.drop(columns=[col])
            print(f"  Removed reserved obs['{col}']")

    # Keep organism in uns for 7.0.0 compatibility
    if "organism_ontology_term_id" not in adata.uns:
        adata.uns["organism_ontology_term_id"] = "NCBITaxon:9606"
    if "organism" not in adata.uns:
        adata.uns["organism"] = "Homo sapiens"

    # HCA requires title in uns
    if "title" not in adata.uns:
        adata.uns["title"] = (
            "An integrated single-cell breast atlas of 2.1 million cells"
        )
        print(f"  Set uns['title'] = {adata.uns['title']}")

    return adata


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def validate_h5ad(adata, label=""):
    """Light validation checks. Returns list of issues."""
    prefix = f"[{label}] " if label else ""
    issues = []

    # Tier 1 fields
    issues.extend(validate_tier1_fields(adata, TIER1_FIELDS, prefix=prefix))

    # Ensembl coverage
    issues.extend(validate_ensembl_coverage(adata, prefix=prefix, verbose=False))

    # UMAP
    issues.extend(validate_umap(adata, prefix=prefix))

    # Schema version
    issues.extend(validate_schema_version(adata, prefix=prefix))

    return issues


# ---------------------------------------------------------------------------
# Per-lineage splitting
# ---------------------------------------------------------------------------


def split_by_lineage(adata, output_dir, lineage_splits):
    """Split h5ad by level0_annotation and write per-lineage files.

    Returns dict of {lineage: n_cells} and list of issues.
    """
    print("\n" + "=" * 72)
    print("Per-Lineage Splitting")
    print("=" * 72)

    if "level0_annotation" not in adata.obs.columns:
        # Fallback to level0
        if "level0" not in adata.obs.columns:
            print("  ERROR: Neither level0_annotation nor level0 in obs. Cannot split.")
            return {}, ["No lineage column available for splitting"]
        split_col = "level0"
        print(f"  Using fallback column: {split_col}")
    else:
        split_col = "level0_annotation"

    all_values = adata.obs[split_col].value_counts()
    print(f"  {split_col} value counts:")
    for val, count in all_values.items():
        print(f"    {val}: {count:,}")

    results = {}
    issues = []
    total_split_cells = 0

    for lineage, filename in lineage_splits.items():
        mask = adata.obs[split_col] == lineage
        n_cells = mask.sum()

        if n_cells == 0:
            # Try case-insensitive match
            mask = adata.obs[split_col].str.lower() == lineage.lower()
            n_cells = mask.sum()

        if n_cells == 0:
            print(f"  WARNING: No cells for lineage '{lineage}'. Skipping {filename}.")
            issues.append(f"No cells for lineage '{lineage}'")
            continue

        print(f"\n  Splitting {lineage}: {n_cells:,} cells -> {filename}")
        t0 = time.time()

        # Subset — use to_memory() if backed
        adata_sub = adata[mask].copy()
        total_split_cells += n_cells

        # Validate subset
        sub_issues = validate_h5ad(adata_sub, label=lineage)
        if sub_issues:
            print(f"    Validation issues: {len(sub_issues)}")
            for issue in sub_issues:
                print(f"      {issue}")
            issues.extend(sub_issues)

        # Write
        out_path = os.path.join(output_dir, filename)
        adata_sub.write_h5ad(out_path)
        file_size = os.path.getsize(out_path) / (1024 ** 3)
        print(f"    Written: {out_path} ({file_size:.1f} GB, {time.time() - t0:.1f}s)")

        results[lineage] = n_cells

        # Free memory before next split (each .copy() duplicates X matrix)
        del adata_sub
        gc.collect()

    # Check for unassigned cells
    assigned_values = set(lineage_splits.keys())
    # Case-insensitive check
    assigned_lower = {v.lower() for v in assigned_values}
    unassigned_mask = ~adata.obs[split_col].str.lower().isin(assigned_lower)
    n_unassigned = unassigned_mask.sum()
    if n_unassigned > 0:
        unassigned_vals = adata.obs.loc[unassigned_mask, split_col].value_counts()
        print(f"\n  Unassigned cells ({n_unassigned:,}):")
        for val, count in unassigned_vals.items():
            print(f"    {val}: {count:,}")

    print(f"\n  Total split: {total_split_cells:,} / {len(adata):,}")
    print(f"  Unassigned: {n_unassigned:,}")

    return results, issues


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Assemble iHBCA integrated h5ad with refined annotations"
    )
    # Legacy mode: load from CxG h5ad
    parser.add_argument("--h5ad", default=None,
                        help="Path to integration_iHBCA.h5ad (legacy mode)")
    # Author share mode: construct from primitives
    parser.add_argument("--gene-data", default=None,
                        help="Path to scVI gene_data.csv (author share mode)")
    parser.add_argument("--counts-npz", default=None,
                        help="Path to counts NPZ (optional, author share mode)")
    parser.add_argument("--umap", default=None,
                        help="Path to X_scVI100_UMAP.csv (author share mode)")
    # Shared args
    parser.add_argument("--annotations", required=True,
                        help="Path to ihbca_level1.5_annotations.csv")
    parser.add_argument("--embeddings", required=True,
                        help="Path to X_scVI100.csv")
    parser.add_argument("--repo-root", required=True,
                        help="Path to iHBCAv1_upload repo root")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for h5ad files")
    parser.add_argument("--schema-version", default="5.3.2",
                        help="Target schema version (default: 5.3.2)")
    parser.add_argument("--skip-splits", action="store_true",
                        help="Skip per-lineage splitting")
    parser.add_argument(
        "--target", choices=["hca", "cxg"], default="hca",
        help="Target schema: hca (default) uses :0601/:0602 HANCESTRO terms; "
             "cxg downgrades to :0004 ancestry terms for CxG 5.3.2 compatibility"
    )
    args = parser.parse_args()

    # Determine mode
    author_share_mode = args.gene_data is not None
    if not author_share_mode and args.h5ad is None:
        parser.error("Either --h5ad (legacy) or --gene-data (author share) is required")

    repo_root = Path(args.repo_root)

    # Load pipeline config
    config = load_pipeline_config(repo_root)

    print("=" * 72)
    print("C1 Integrated Object Assembly")
    print("=" * 72)
    print(f"mode:         {'author-share' if author_share_mode else 'legacy (CxG h5ad)'}")
    if args.h5ad:
        print(f"h5ad:         {args.h5ad}")
    if args.gene_data:
        print(f"gene-data:    {args.gene_data}")
    if args.counts_npz:
        print(f"counts-npz:   {args.counts_npz}")
    if args.umap:
        print(f"umap:         {args.umap}")
    print(f"annotations:  {args.annotations}")
    print(f"embeddings:   {args.embeddings}")
    print(f"repo-root:    {args.repo_root}")
    print(f"output-dir:   {args.output_dir}")
    print(f"schema:       {args.schema_version}")
    print(f"target:       {args.target}")
    print()

    # Verify inputs
    required_inputs = [("annotations", args.annotations), ("embeddings", args.embeddings)]
    if author_share_mode:
        required_inputs.append(("gene-data", args.gene_data))
        if args.umap:
            required_inputs.append(("umap", args.umap))
    else:
        required_inputs.append(("h5ad", args.h5ad))

    for label, path in required_inputs:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    if author_share_mode:
        # ==============================================================
        # AUTHOR SHARE MODE: Construct from primitives
        # ==============================================================
        print("=" * 72)
        print("STEP 1: Construct from Author Share Primitives")
        print("=" * 72)
        t0 = time.time()

        adata, ann = construct_from_author_share(
            annotations_path=args.annotations,
            gene_data_path=args.gene_data,
            embeddings_path=args.embeddings,
            umap_path=args.umap,
            counts_npz_path=args.counts_npz,
        )
        print(f"  Total construction time: {time.time() - t0:.1f}s")

        # Set Tier 1 constants (not inherited from CxG h5ad)
        print()
        print("=" * 72)
        print("STEP 2: Set Tier 1 Constants")
        print("=" * 72)
        adata = set_tier1_constants(adata)

    else:
        # ==============================================================
        # LEGACY MODE: Load from CxG h5ad
        # ==============================================================
        print("=" * 72)
        print("STEP 1: Load h5ad")
        print("=" * 72)
        t0 = time.time()
        print(f"Reading {args.h5ad} (full load)...")
        adata = ad.read_h5ad(args.h5ad)
        print(f"  Loaded in {time.time() - t0:.1f}s")
        print(f"  Shape: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
        print(f"  X dtype: {adata.X.dtype}, format: {getattr(adata.X, 'format', 'dense')}")

        # A1: Add feature_is_filtered (required by CxG + HCA validators)
        adata.var["feature_is_filtered"] = False
        print(f"  Added var['feature_is_filtered'] = False")

        # A6: Rename reserved var columns (HCA reserves these names)
        var_renames = {k: v for k, v in RESERVED_VAR_RENAMES.items()
                       if k in adata.var.columns}
        if var_renames:
            adata.var = adata.var.rename(columns=var_renames)
            print(f"  Renamed {len(var_renames)} reserved var columns: "
                  f"{list(var_renames.keys())}")
            # Also rename in raw.var if raw exists
            if adata.raw is not None:
                raw_renames = {k: v for k, v in RESERVED_VAR_RENAMES.items()
                              if k in adata.raw.var.columns}
                if raw_renames:
                    raw_var = adata.raw.var.rename(columns=raw_renames)
                    adata.raw = ad.AnnData(
                        X=adata.raw.X,
                        var=raw_var,
                        obs=adata.obs,
                    )
                    print(f"  Renamed {len(raw_renames)} reserved raw.var columns")

        # Step 2: Load annotations and embeddings
        print()
        print("=" * 72)
        print("STEP 2: Load Annotations + Embeddings")
        print("=" * 72)
        ann = load_annotations(args.annotations)
        ann, keep_mask = remap_mh0023_cell_ids(ann)
        emb = load_embeddings(args.embeddings, ann.index, keep_mask=keep_mask)

        # Step 3: Load HANCESTRO mapping (loaded below for both modes)

        # Step 4: Merge annotations
        print()
        print("=" * 72)
        print("STEP 4: Merge Annotations")
        print("=" * 72)
        adata = merge_annotations(adata, ann, ANNOTATION_COLUMNS_TO_MERGE)

        # Step 5: Add embeddings
        print()
        print("=" * 72)
        print("STEP 5: Add Embeddings")
        print("=" * 72)
        adata = add_embeddings(adata, emb)

    # ==================================================================
    # COMMON PATH: Steps 6+ are the same for both modes
    # ==================================================================

    # Load HANCESTRO mapping (needed by both modes)
    hancestro_mapping = load_hancestro_mapping(repo_root, config=config)
    print(f"  HANCESTRO: {len(hancestro_mapping)} ethnicity entries")

    # ------------------------------------------------------------------
    # Step 6: Enrich ethnicity
    # ------------------------------------------------------------------
    print()
    print("=" * 72)
    print("STEP 6: Enrich Ethnicity")
    print("=" * 72)
    adata = enrich_ethnicity(adata, ann, hancestro_mapping)

    # ------------------------------------------------------------------
    # Step 6b: Enrich ethnicity from harmonized donor metadata
    # ------------------------------------------------------------------
    print()
    print("=" * 72)
    print("STEP 6b: Enrich Ethnicity from Harmonized Donor Metadata")
    print("=" * 72)
    adata = enrich_ethnicity_from_harmonized(adata, repo_root, hancestro_mapping)

    # ------------------------------------------------------------------
    # Step 6c: Enrich donor metadata from L1 harmonized CSV
    # ------------------------------------------------------------------
    print()
    print("=" * 72)
    print("STEP 6c: Enrich Donor Metadata from L1 Harmonized CSV")
    print("=" * 72)
    adata = enrich_donor_metadata_integrated(adata, repo_root)

    # ------------------------------------------------------------------
    # Step 6d: Populate HCA obs fields (B1+B2+B3)
    # ------------------------------------------------------------------
    print()
    print("=" * 72)
    print("STEP 6d: Populate HCA Obs Fields")
    print("=" * 72)
    adata = populate_hca_obs_fields_integrated(adata, repo_root, config=config)

    # ------------------------------------------------------------------
    # Step 6e: CL level1.5 fallback for unknown cells
    # ------------------------------------------------------------------
    print()
    print("=" * 72)
    print("STEP 6e: CL Level1.5 Fallback")
    print("=" * 72)
    adata = apply_cl_level15_fallback(adata, repo_root)

    # ------------------------------------------------------------------
    # Step 7: Fix Tier 1 fields + schema
    # ------------------------------------------------------------------
    print()
    print("=" * 72)
    print("STEP 7: Tier 1 Fields + Schema")
    print("=" * 72)
    adata = fix_tier1_fields(adata, target=args.target)
    adata = update_schema(adata, version=args.schema_version)

    # A7: Set study_pi as deduplicated list of all source study PIs
    study_pis = collect_study_pis(repo_root)
    if study_pis:
        adata.uns["study_pi"] = study_pis
        print(f"  Set uns['study_pi'] = {study_pis}")

    # Set DOIs and data accessions from all source studies
    dois, geo_accessions = collect_study_accessions(repo_root)
    if dois:
        adata.uns["doi"] = dois
        print(f"  Set uns['doi'] = {dois}")
    if geo_accessions:
        adata.uns["geo_accession"] = geo_accessions
        print(f"  Set uns['geo_accession'] = {geo_accessions}")

    # ------------------------------------------------------------------
    # Step 8: Validate
    # ------------------------------------------------------------------
    print()
    print("=" * 72)
    print("STEP 8: Validation")
    print("=" * 72)
    issues = validate_h5ad(adata, label="all-breast-cells")
    if issues:
        print(f"  {len(issues)} issue(s):")
        for issue in issues:
            print(f"    - {issue}")
    else:
        print("  All checks passed.")

    # ------------------------------------------------------------------
    # Step 8b: Provenance
    # ------------------------------------------------------------------
    print()
    print("=" * 72)
    print("STEP 8b: Provenance")
    print("=" * 72)
    prov_inputs = {
        "annotations": args.annotations,
        "embeddings": args.embeddings,
        "hancestro": repo_root / "publication/mappings/hancestro_ethnicity_mapping.tsv",
        "l1_donor_metadata": repo_root / "publication/config/metadata_stages/L1_harmonized_donor.csv",
    }
    if author_share_mode:
        prov_inputs["gene_data"] = args.gene_data
        if args.counts_npz:
            prov_inputs["counts_npz"] = args.counts_npz
        if args.umap:
            prov_inputs["umap"] = args.umap
        prov_inputs["build_mode"] = "author_share"
    else:
        prov_inputs["h5ad"] = args.h5ad
        prov_inputs["build_mode"] = "legacy_cxg"
    prov_config = {
        "dataset_metadata": repo_root / "publication/config/dataset_metadata.yaml",
    }
    out_path = os.path.join(args.output_dir, "all-breast-cells.h5ad")
    build_provenance_manifest(
        adata, Path(out_path), "integrated", repo_root,
        inputs=prov_inputs, config_files=prov_config,
        build_script="assemble_integrated.py",
        sidecar_path=None,
    )

    # ------------------------------------------------------------------
    # Step 9: Write all-breast-cells.h5ad
    # ------------------------------------------------------------------
    print()
    print("=" * 72)
    print("STEP 9: Write all-breast-cells.h5ad")
    print("=" * 72)
    print(f"Writing {out_path}...")
    t0 = time.time()
    adata.write_h5ad(out_path)
    file_size = os.path.getsize(out_path) / (1024 ** 3)
    print(f"  Written in {time.time() - t0:.1f}s ({file_size:.1f} GB)")

    # ------------------------------------------------------------------
    # Step 10: Per-lineage splits
    # ------------------------------------------------------------------
    if not args.skip_splits:
        split_results, split_issues = split_by_lineage(adata, args.output_dir, LINEAGE_SPLITS)
        issues.extend(split_issues)
    else:
        print("\nSkipping per-lineage splits (--skip-splits)")
        split_results = {}

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print()
    print("=" * 72)
    print("SUMMARY")
    print("=" * 72)
    print(f"All-breast-cells: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print(f"  Output: {out_path}")
    print(f"  Schema: {adata.uns.get('schema_version')}")
    print(f"  Tier 1 fields: {sum(1 for f in TIER1_FIELDS if f in adata.obs.columns)}/11")
    print(f"  obsm keys: {list(adata.obsm.keys())}")
    if split_results:
        print(f"  Lineage splits: {split_results}")
    if issues:
        print(f"\n  Total issues: {len(issues)}")
        for issue in issues:
            print(f"    - {issue}")
    else:
        print("\n  No issues found.")

    # Write summary YAML
    summary = {
        "assembly_date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "schema_version": args.schema_version,
        "all_breast_cells": {
            "cells": int(adata.n_obs),
            "genes": int(adata.n_vars),
            "file": out_path,
            "size_gb": round(file_size, 1),
        },
        "lineage_splits": {k: int(v) for k, v in split_results.items()},
        "tier1_present": sum(1 for f in TIER1_FIELDS if f in adata.obs.columns),
        "obsm_keys": list(adata.obsm.keys()),
        "issues": issues if issues else [],
        "provenance": adata.uns.get("ihbca_provenance", {}),
    }
    summary_path = os.path.join(args.output_dir, "assembly_report.yaml")
    with open(summary_path, "w") as f:
        yaml.dump(summary, f, default_flow_style=False, sort_keys=False)
    print(f"\nAssembly report: {summary_path}")
    print("Done.")


if __name__ == "__main__":
    main()
