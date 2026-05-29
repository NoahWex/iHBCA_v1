#!/usr/bin/env python3
"""
enrich_h5ads.py - Post-assembly enrichment for iHBCA h5ads
==========================================================
Standalone enrichment script that post-processes all 11 h5ads with:
  - var annotations (gene_symbol, feature_biotype, chromosome from GENCODE v24)
  - obs metadata (ethnicity_verbatim, ethnicity_grouped, ihbca_donor_id,
    metadata_notes -- integrated only)
  - log_normalized expression layer (integrated only)
  - uns documentation (obsm_key_descriptions, obs_field_descriptions,
    layer_descriptions, atlas_metadata, integration_method)
  - obs polish (column ordering, categoricals, var.index.name)

Usage:
  # Source dataset
  python enrich_h5ads.py --h5ad path/to/gray2022.h5ad --mode source \
    --study gray --repo-root /path/to/repo

  # Integrated object
  python enrich_h5ads.py --h5ad path/to/all-breast-cells.h5ad --mode integrated \
    --repo-root /path/to/repo

"""

import argparse
import gc
import os
import shutil
import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

# Add scripts dir to path for ihbca imports
sys.path.insert(0, str(Path(__file__).parent))

from ihbca.constants import L1_NUMERIC, L1_OBS_COLUMNS
from ihbca.loaders import (
    load_donor_translations,
    load_gencode_annotations,
    load_gene_mapping,
    load_l1_metadata,
    load_pipeline_config,
)


# ============================================================================
# Static description dictionaries
# ============================================================================

OBSM_KEY_DESCRIPTIONS = {
    # Source dataset keys
    "X_umap": "2D UMAP coordinates from study-native embedding (float32)",
    "X_scVI_joint": "scVI joint embedding from iHBCA integration, dimensionality varies by study (float32)",
    "X_scVI_native": "Study-native scVI embedding, dimensionality varies by study (float32)",
    "X_ihbca_scvi_100": "100D joint scVI embedding from iHBCA integrated object (float32)",
    "X_umap_ihbca_scvi_100": "Per-study UMAP from X_ihbca_scvi_100 (scanpy, n_neighbors=15, min_dist=0.5, float32)",
    # Integrated keys
    "X_scvi_100": "scVI 100D latent representation (n_latent=100, n_layers=2, batch_key=donor_id, float32)",
}

OBS_FIELD_DESCRIPTIONS = {
    # --- Tier 1 CxG fields ---
    "organism_ontology_term_id": "NCBI Taxonomy ID (NCBITaxon:9606 = Homo sapiens)",
    "assay_ontology_term_id": "EFO assay ontology term (10x 3' v2 or v3)",
    "tissue_ontology_term_id": "UBERON tissue term (UBERON:0000310 = breast)",
    "tissue_type": "CxG tissue type (always 'tissue')",
    "disease_ontology_term_id": "MONDO/PATO disease term (PATO:0000461 = normal)",
    "donor_id": "Unique donor identifier, format {Study}_{id} in integrated, ihbca_donor_id in source",
    "sex_ontology_term_id": "PATO sex term (PATO:0000383 = female)",
    "development_stage_ontology_term_id": "HsapDv developmental stage based on donor age",
    "self_reported_ethnicity_ontology_term_id": "HANCESTRO self-reported ethnicity term",
    "is_primary_data": "Whether this is primary data (True for source, False for integrated)",
    "suspension_type": "Cell suspension type (always 'cell')",
    "cell_type_ontology_term_id": "Cell Ontology term from CxG passthrough or level1.5 CL mapping",
    # --- CxG/HCA label columns ---
    "cell_type_label": "Human-readable cell type label",
    "assay_label": "Human-readable assay label",
    "disease_label": "Human-readable disease label",
    "sex_label": "Human-readable sex label",
    "tissue_label": "Human-readable tissue label",
    "self_reported_ethnicity_label": "Human-readable self-reported ethnicity",
    "development_stage_label": "Human-readable developmental stage",
    # --- L1 harmonized donor metadata ---
    "age_continuous": "Age in years at sample collection (float, from L1 harmonized metadata)",
    "age_binary": "Age group: young (<45) or older (>=45) (from L1)",
    "parity_count": "Number of pregnancies carried to viability (int, from L1)",
    "parity_binary": "Parous (>=1 pregnancy) vs nulliparous (from L1)",
    "age_at_first_birth": "Age at first live birth in years (float, from L1)",
    "brca_genotype": "BRCA1/BRCA2 mutation status (from L1)",
    "cancer_history": "Personal breast cancer history (from L1)",
    "tissue_indication": "Clinical indication for tissue collection (from L1)",
    "risk_status_binary": "High-risk vs normal-risk classification (from L1)",
    "risk_genotype_only": "Risk based on genotype only, excluding family history (from L1)",
    "menopausal_status_detailed": "Detailed menopausal status: pre/peri/post (from L1)",
    "menopausal_status_binary": "Pre-menopausal vs post-menopausal (from L1)",
    "bmi_continuous": "Body mass index, continuous (float, from L1)",
    "bmi_category": "BMI category: underweight/normal/overweight/obese (from L1)",
    "sample_preservation": "Sample preservation method: fresh/frozen/fixed (from L1)",
    "sample_type": "Tissue sample type: reduction mammoplasty, mastectomy, etc. (from L1)",
    "facs_status": "Whether cells were FACS-sorted and enrichment target (from L1)",
    "dissociation_minutes": "Duration of tissue dissociation in minutes (float, from L1)",
    "ethnicity_verbatim": "Self-reported ethnicity as stated in source publication (from L1)",
    "ethnicity_grouped": "Grouped ethnicity category: european, asian, african, etc. (from L1)",
    "metadata_notes": "Donor-level notes on metadata quality or special circumstances (from L1)",
    # --- Integrated-only obs ---
    "ihbca_donor_id": "Raw donor identifier without study prefix (from L1 reverse mapping)",
    "dataset": "Source study name, lowercase (gray, kumar, murrow, nee, twigger, reed, pal)",
    # --- HCA obs fields ---
    "sample_id": "Sample identifier (donor_id used as proxy)",
    "manner_of_death": "Manner of death (always 'not applicable' -- living donors)",
    "sample_source": "Sample source (always 'surgical donor')",
    "sampled_site_condition": "Condition at sampled site (always 'normal')",
    "sample_collection_method": "Collection method (always 'biopsy')",
    "sample_preservation_method": "HCA preservation method from FACS status normalization",
    "cell_enrichment": "Cell enrichment method from FACS status normalization",
    "institute": "Research institute of source study PI",
    "library_id": "Sequencing library identifier (from SRA run tables)",
    "library_sequencing_run": "Sequencing run identifier (from SRA run tables)",
    "library_preparation_batch": "Library prep batch (from SRA run tables)",
    "alignment_software": "Alignment software used by source study (from dataset_metadata.yaml)",
    "sequenced_fragment": "Sequenced fragment type (from dataset_metadata.yaml)",
    "reference_genome": "Reference genome used for alignment (from dataset_metadata.yaml)",
    "gene_annotation_version": "Gene annotation version used (from dataset_metadata.yaml)",
    # --- Annotation columns (integrated only) ---
    "level0_annotation": "Broad lineage: Epithelial, Stromal, or Immune (Austin's refined annotation)",
    "level1_annotation": "Mid-level cell type annotation (Austin's refined annotation)",
    "level1.5_annotation": "Fine-grained cell type annotation (Austin's refined annotation)",
    "level0": "Original broad lineage from CxG h5ad",
    "level1": "Original mid-level cell type from CxG h5ad",
    "cellTypist_annotation_reed": "CellTypist cross-study prediction using Reed model",
    "cellTypist_annotation_kumar": "CellTypist cross-study prediction using Kumar model",
    # --- QC metrics ---
    "n_genes": "Number of genes detected per cell",
    "percent_mito": "Percentage of mitochondrial gene counts per cell",
    "n_counts": "Total UMI counts per cell",
    # --- Doublet prediction ---
    "prob_spikein": "Spike-in doublet probability",
    "prob_spikein_dblt": "Spike-in doublet probability (alternate)",
    "pred_spikein": "Spike-in doublet prediction (boolean)",
    # --- Leiden clustering ---
    "leiden_0_05": "Leiden clustering at resolution 0.05 (scVI 20D basis)",
    "leiden_0_1": "Leiden clustering at resolution 0.1 (scVI 20D basis)",
    "leiden_0_2": "Leiden clustering at resolution 0.2 (scVI 20D basis)",
    "leiden_0_3": "Leiden clustering at resolution 0.3 (scVI 20D basis)",
    "leiden_0_4": "Leiden clustering at resolution 0.4 (scVI 20D basis)",
    "leiden_0_5": "Leiden clustering at resolution 0.5 (scVI 20D basis)",
    "leiden_0_6": "Leiden clustering at resolution 0.6 (scVI 20D basis)",
    "leiden_0_7": "Leiden clustering at resolution 0.7 (scVI 20D basis)",
    "leiden_0_8": "Leiden clustering at resolution 0.8 (scVI 20D basis)",
    "leiden_0_9": "Leiden clustering at resolution 0.9 (scVI 20D basis)",
    "leiden_1_0": "Leiden clustering at resolution 1.0 (scVI 20D basis)",
    "leiden_scVI100_res0.1": "Leiden clustering at resolution 0.1 (scVI 100D basis)",
    "leiden_scVI100_res0.5": "Leiden clustering at resolution 0.5 (scVI 100D basis)",
    "leiden_scVI100_res1.0": "Leiden clustering at resolution 1.0 (scVI 100D basis)",
    "leiden_pca_res0.1": "Leiden clustering at resolution 0.1 (PCA basis)",
    "leiden_pca_res0.5": "Leiden clustering at resolution 0.5 (PCA basis)",
    "leiden_pca_res1.0": "Leiden clustering at resolution 1.0 (PCA basis)",
    "leiden_scPoli50_res0.1": "Leiden clustering at resolution 0.1 (scPoli 50D basis)",
    "leiden_scPoli50_res0.5": "Leiden clustering at resolution 0.5 (scPoli 50D basis)",
    "leiden_scPoli50_res1.0": "Leiden clustering at resolution 1.0 (scPoli 50D basis)",
    # --- Other annotation ---
    "risk_status2": "Refined risk stratification category",
    "BRCA_tested": "Whether donor was genetically tested for BRCA mutations",
    "sample_type_coarse": "Coarse sample type grouping",
    "sampleID": "Original sample identifier from source study",
    "sample_barcode": "10x barcode associated with sample",
    "patientID": "Original patient identifier from source study",
}

OBS_COLUMN_ORDER = [
    # --- Tier 1 CxG/HCA (12 fields) ---
    "organism_ontology_term_id",
    "assay_ontology_term_id",
    "tissue_ontology_term_id",
    "tissue_type",
    "disease_ontology_term_id",
    "donor_id",
    "sex_ontology_term_id",
    "development_stage_ontology_term_id",
    "self_reported_ethnicity_ontology_term_id",
    "is_primary_data",
    "suspension_type",
    "cell_type_ontology_term_id",
    # --- CxG/HCA label columns ---
    "cell_type_label",
    "assay_label",
    "disease_label",
    "sex_label",
    "tissue_label",
    "self_reported_ethnicity_label",
    "development_stage_label",
    # --- Donor identity ---
    "ihbca_donor_id",
    "dataset",
    # --- L1 clinical metadata ---
    "age_continuous",
    "age_binary",
    "parity_count",
    "parity_binary",
    "age_at_first_birth",
    "brca_genotype",
    "cancer_history",
    "tissue_indication",
    "risk_status_binary",
    "risk_genotype_only",
    "menopausal_status_detailed",
    "menopausal_status_binary",
    "ethnicity_verbatim",
    "ethnicity_grouped",
    "bmi_continuous",
    "bmi_category",
    "sample_preservation",
    "sample_type",
    "facs_status",
    "dissociation_minutes",
    "metadata_notes",
    # --- HCA obs fields ---
    "sample_id",
    "manner_of_death",
    "sample_source",
    "sampled_site_condition",
    "sample_collection_method",
    "sample_preservation_method",
    "cell_enrichment",
    "institute",
    "library_id",
    "library_sequencing_run",
    "library_preparation_batch",
    "alignment_software",
    "sequenced_fragment",
    "reference_genome",
    "gene_annotation_version",
    # --- Annotation hierarchy ---
    "level0_annotation",
    "level1_annotation",
    "level1.5_annotation",
    "level0",
    "level1",
    # --- CellTypist ---
    "cellTypist_annotation_reed",
    "cellTypist_annotation_kumar",
    # --- QC metrics ---
    "n_genes",
    "percent_mito",
    "n_counts",
    # --- Doublet ---
    "prob_spikein",
    "prob_spikein_dblt",
    "pred_spikein",
]


# ============================================================================
# Enrichment functions
# ============================================================================


def enrich_var(adata, repo_root, mode, force=False):
    """Add gene annotations to var: gene_symbol, feature_biotype, chromosome."""
    gencode = load_gencode_annotations(repo_root)
    if gencode is None:
        print("  ERROR: Cannot enrich var without GENCODE annotations")
        sys.exit(1)

    print(f"  GENCODE annotations loaded: {len(gencode)} genes")

    # gene_symbol: for source h5ads, may already exist
    if "gene_symbol" not in adata.var.columns or force:
        if mode == "integrated":
            # Integrated var.index is ensembl_id; reverse-lookup from gene mapping
            gene_map = load_gene_mapping(repo_root)
            ensembl_to_symbol = {v: k for k, v in gene_map.items()}
            # Prefer GENCODE symbol, fall back to reverse gene_map
            symbols = []
            for eid in adata.var.index:
                if eid in gencode.index:
                    symbols.append(gencode.loc[eid, "gene_symbol"])
                elif eid in ensembl_to_symbol:
                    symbols.append(ensembl_to_symbol[eid])
                else:
                    symbols.append("")
            adata.var["gene_symbol"] = symbols
        else:
            # Source: lookup from GENCODE by ensembl_id index
            adata.var["gene_symbol"] = adata.var.index.map(
                lambda eid: gencode.loc[eid, "gene_symbol"] if eid in gencode.index else ""
            )
        n_mapped = (adata.var["gene_symbol"] != "").sum()
        print(f"  gene_symbol: {n_mapped}/{len(adata.var)} mapped")
    else:
        print("  gene_symbol: already present, skipping")

    # feature_biotype_gencode (CxG reserves "feature_biotype" — use suffixed name)
    if "feature_biotype_gencode" not in adata.var.columns or force:
        adata.var["feature_biotype_gencode"] = adata.var.index.map(
            lambda eid: gencode.loc[eid, "feature_biotype"] if eid in gencode.index else np.nan
        )
        adata.var["feature_biotype_gencode"] = pd.Categorical(adata.var["feature_biotype_gencode"])
        n_mapped = adata.var["feature_biotype_gencode"].notna().sum()
        print(f"  feature_biotype_gencode: {n_mapped}/{len(adata.var)} mapped")
    else:
        print("  feature_biotype_gencode: already present, skipping")

    # chromosome
    if "chromosome" not in adata.var.columns or force:
        adata.var["chromosome"] = adata.var.index.map(
            lambda eid: gencode.loc[eid, "chromosome"] if eid in gencode.index else np.nan
        )
        adata.var["chromosome"] = pd.Categorical(adata.var["chromosome"])
        n_mapped = adata.var["chromosome"].notna().sum()
        print(f"  chromosome: {n_mapped}/{len(adata.var)} mapped")
    else:
        print("  chromosome: already present, skipping")



def enrich_obs(adata, repo_root, force=False):
    """Add obs columns from L1 harmonized donor metadata (integrated only).

    Adds ihbca_donor_id (reverse-mapped) plus all 21 L1_OBS_COLUMNS broadcast
    from L1_harmonized_donor.csv via donor_id matching.
    """
    l1 = load_l1_metadata(repo_root)
    if l1 is None:
        print("  ERROR: Cannot enrich obs without L1 metadata")
        sys.exit(1)

    xlat = load_donor_translations(repo_root)

    # Build reverse map: h5ad donor_id -> L1 row
    # L1 has columns: ihbca_donor_id, study, plus obs columns
    reverse_map = {}  # h5ad_donor_id -> {col: value}
    for _, row in l1.iterrows():
        study_cap = row["study"].capitalize()
        h5ad_key = f"{study_cap}_{row['ihbca_donor_id']}"
        reverse_map[h5ad_key] = row

    # Add translated keys (Reed/Pal use different IDs in h5ad)
    for l1_key, cxg_id in xlat.items():
        # l1_key = "{Study}_{ihbca_donor_id}"
        parts = l1_key.split("_", 1)
        study_cap = parts[0]
        ihbca_id = parts[1]
        # Find the L1 row for this donor
        match = l1[(l1["study"] == study_cap.lower()) & (l1["ihbca_donor_id"] == ihbca_id)]
        if len(match) > 0:
            reverse_map[cxg_id] = match.iloc[0]
            # Also add study-prefixed form
            if not cxg_id.startswith(f"{study_cap}_"):
                reverse_map[f"{study_cap}_{cxg_id}"] = match.iloc[0]

    # Report reverse map coverage
    donor_ids = adata.obs["donor_id"].unique()
    n_matched = sum(1 for did in donor_ids if did in reverse_map)
    print(f"  Reverse map: {n_matched}/{len(donor_ids)} unique donor_ids matched")

    # --- ihbca_donor_id ---
    if "ihbca_donor_id" not in adata.obs.columns or force:
        ihbca_ids = []
        for did in adata.obs["donor_id"]:
            if did in reverse_map:
                ihbca_ids.append(reverse_map[did]["ihbca_donor_id"])
            else:
                # Fallback: strip first underscore segment
                ihbca_ids.append(did.split("_", 1)[1] if "_" in did else did)
        adata.obs["ihbca_donor_id"] = ihbca_ids
        n_mapped = sum(1 for did in adata.obs["donor_id"] if did in reverse_map)
        print(f"  ihbca_donor_id: {n_mapped}/{len(adata.obs)} cells mapped via L1, rest via fallback")
    else:
        print("  ihbca_donor_id: already present, skipping")

    # --- Broadcast all 21 L1 obs columns ---
    print(f"  Broadcasting {len(L1_OBS_COLUMNS)} L1 columns...")
    for col in L1_OBS_COLUMNS:
        if col in adata.obs.columns and not force:
            print(f"    {col}: already present, skipping")
            continue

        is_numeric = col in L1_NUMERIC
        is_notes = col == "metadata_notes"

        values = []
        for did in adata.obs["donor_id"]:
            if did in reverse_map:
                val = reverse_map[did].get(col, "")
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
            adata.obs[col] = pd.array(values, dtype=pd.Float64Dtype())
            n_present = sum(1 for v in values if not pd.isna(v))
            print(f"    {col}: {n_present}/{len(adata.obs)} non-NaN")
        else:
            adata.obs[col] = values
            if is_notes:
                n_filled = sum(1 for v in values if v != "")
                print(f"    {col}: {n_filled}/{len(adata.obs)} with notes")
            else:
                n_unknown = sum(1 for v in values if v == "unknown")
                print(f"    {col}: {len(adata.obs) - n_unknown}/{len(adata.obs)} mapped, {n_unknown} unknown")


def add_log_normalized(adata):
    """Add layers['log_normalized'] = log1p(normalize_total(X, target_sum=1e4))."""
    if "log_normalized" in adata.layers:
        print("  log_normalized layer: already present, skipping")
        return

    import scanpy as sc

    print("  Computing log_normalized layer...")
    print(f"    Input X: {adata.X.shape}, nnz={adata.X.nnz if hasattr(adata.X, 'nnz') else 'dense'}")

    tmp = adata.X.copy()
    tmp_ad = ad.AnnData(X=tmp)
    sc.pp.normalize_total(tmp_ad, target_sum=1e4)
    sc.pp.log1p(tmp_ad)
    adata.layers["log_normalized"] = tmp_ad.X
    del tmp_ad, tmp
    gc.collect()

    print(f"    Done. Layer shape: {adata.layers['log_normalized'].shape}")


def enrich_uns(adata, mode, force=False):
    """Add documentation dicts to uns."""

    # CxG schema version -- do NOT set in uns.
    # The CxG validator reserves this key and sets it automatically.
    # Having it pre-set causes: "Column 'schema_version' is a reserved column
    # name of 'uns'. Remove it from h5ad and try again."
    if "schema_version" in adata.uns:
        del adata.uns["schema_version"]
        print("  schema_version: removed (CxG-reserved key)")
    else:
        print("  schema_version: absent (correct — CxG-reserved key)")

    # obsm_key_descriptions
    if "obsm_key_descriptions" not in adata.uns or force:
        desc = {}
        for k, v in adata.obsm.items():
            if k in OBSM_KEY_DESCRIPTIONS:
                desc[k] = OBSM_KEY_DESCRIPTIONS[k]
            else:
                desc[k] = f"Embedding array, shape {v.shape}"
        adata.uns["obsm_key_descriptions"] = desc
        print(f"  obsm_key_descriptions: {len(desc)} keys")
    else:
        print("  obsm_key_descriptions: already present, skipping")

    # obs_field_descriptions
    if "obs_field_descriptions" not in adata.uns or force:
        desc = {}
        for col in adata.obs.columns:
            if col in OBS_FIELD_DESCRIPTIONS:
                desc[col] = OBS_FIELD_DESCRIPTIONS[col]
            else:
                desc[col] = f"Obs column '{col}' (no description available)"
        adata.uns["obs_field_descriptions"] = desc
        print(f"  obs_field_descriptions: {len(desc)} columns")
    else:
        print("  obs_field_descriptions: already present, skipping")

    # layer_descriptions — REMOVED: CxG schema 5.3.2 deprecated this uns field
    # Layer info now lives only in obs_field_descriptions and documentation.

    # Integrated-only uns
    if mode == "integrated":
        if "atlas_metadata" not in adata.uns or force:
            adata.uns["atlas_metadata"] = {
                "atlas_name": "integrated Human Breast Cell Atlas (iHBCA)",
                "atlas_version": "1.0",
                "total_cells": int(adata.n_obs),
                "total_donors": int(adata.obs["donor_id"].nunique()),
                "total_studies": 7,
                "source_studies": ["gray", "kumar", "murrow", "nee", "twigger", "reed", "pal"],
                "citation": "Reed et al. 2024, Nature Genetics",
            }
            print(f"  atlas_metadata: set ({adata.n_obs} cells, {adata.obs['donor_id'].nunique()} donors)")
        else:
            print("  atlas_metadata: already present, skipping")

        if "integration_method" not in adata.uns or force:
            adata.uns["integration_method"] = {
                "method": "scVI",
                "library": "scvi-tools",
                "n_latent": 100,
                "n_layers": 2,
                "gene_likelihood": "nb",
                "use_layer_norm": "both",
                "use_batch_norm": "none",
                "encode_covariates": True,
                "dropout_rate": 0.2,
                "batch_key": "batch (donor_id)",
                "n_top_genes": 5000,
                "train_size": 0.9,
                "early_stopping": True,
                "early_stopping_patience": 45,
                "max_epochs": 400,
                "batch_size": 1024,
                "limit_train_batches": 20,
            }
            print("  integration_method: set")
        else:
            print("  integration_method: already present, skipping")


def polish(adata):
    """Reorder obs columns, convert strings to categoricals, set var.index.name."""

    # F0. Drop redundant BRCA columns superseded by harmonized brca_genotype
    BRCA_REDUNDANT = ["BRCA_tested", "brca_testing"]
    dropped = [c for c in BRCA_REDUNDANT if c in adata.obs.columns]
    if dropped:
        adata.obs.drop(columns=dropped, inplace=True)
        print(f"  Dropped redundant BRCA columns: {dropped}")

    # F1. obs column ordering
    ordered = [c for c in OBS_COLUMN_ORDER if c in adata.obs.columns]
    remaining = sorted(c for c in adata.obs.columns if c not in set(ordered))
    adata.obs = adata.obs[ordered + remaining]
    print(f"  obs columns reordered: {len(ordered)} priority + {len(remaining)} remaining")

    # F2. Categorical dtypes
    n_converted = 0
    for col in adata.obs.columns:
        if adata.obs[col].dtype == object or adata.obs[col].dtype.name == "object":
            # Skip bool-like
            if col == "is_primary_data":
                continue
            # Skip numeric columns
            if col in L1_NUMERIC:
                continue
            # Skip high-cardinality
            n_unique = adata.obs[col].nunique()
            if n_unique > 10000:
                continue
            adata.obs[col] = pd.Categorical(adata.obs[col])
            n_converted += 1
    print(f"  Categoricals: converted {n_converted} string columns")

    # F3. Convert nullable float columns to numpy float64 (anndata can't write FloatingArray)
    n_float_converted = 0
    for col in adata.obs.columns:
        if isinstance(adata.obs[col].dtype, pd.Float64Dtype):
            adata.obs[col] = adata.obs[col].astype("float64")
            n_float_converted += 1
    if n_float_converted:
        print(f"  Nullable floats → float64: {n_float_converted} columns")

    # F4. var.index.name
    adata.var.index.name = "ensembl_id"
    print("  var.index.name = 'ensembl_id'")


# ============================================================================
# Main pipeline
# ============================================================================


def enrich(h5ad_path, mode, study, repo_root, output_path, force, skip_layers):
    """Run full enrichment pipeline on one h5ad."""
    repo_root = Path(repo_root)

    # Load pipeline config for config-driven path resolution
    config = load_pipeline_config(repo_root)

    print(f"=" * 60)
    print(f"Enriching: {h5ad_path}")
    print(f"  Mode: {mode}")
    if study:
        print(f"  Study: {study}")
    print(f"  Force: {force}")
    print(f"  Skip layers: {skip_layers}")
    print(f"=" * 60)

    # 1. Read h5ad
    print(f"\n[1/6] Reading h5ad...")
    adata = ad.read_h5ad(h5ad_path)
    print(f"  Shape: {adata.shape}")
    print(f"  obs columns: {len(adata.obs.columns)}")
    print(f"  var columns: {list(adata.var.columns)}")
    print(f"  obsm keys: {list(adata.obsm.keys())}")
    print(f"  layers: {list(adata.layers.keys())}")

    # 2. Enrich var
    print(f"\n[2/6] Enriching var...")
    enrich_var(adata, repo_root, mode, force)

    # 3. Enrich obs (integrated only)
    if mode == "integrated":
        print(f"\n[3/6] Enriching obs (integrated)...")
        enrich_obs(adata, repo_root, force)
    else:
        print(f"\n[3/6] Skipping obs enrichment (source mode)")

    # 4. Add log_normalized layer (integrated only)
    if mode == "integrated" and not skip_layers:
        print(f"\n[4/6] Adding log_normalized layer...")
        add_log_normalized(adata)
    else:
        print(f"\n[4/6] Skipping layer addition")

    # 5. Enrich uns
    print(f"\n[5/6] Enriching uns...")
    enrich_uns(adata, mode, force)

    # 6. Polish
    print(f"\n[6/6] Polishing...")
    polish(adata)

    # Write (via /tmp for CRSP safety)
    if output_path is None:
        output_path = h5ad_path

    filename = os.path.basename(output_path)
    tmp_path = os.path.join(tempfile.gettempdir(), f"enriched_{filename}")

    print(f"\nWriting to tmp: {tmp_path}")
    adata.write_h5ad(tmp_path)
    tmp_size = os.path.getsize(tmp_path)
    print(f"  Tmp size: {tmp_size:,} bytes ({tmp_size / 1e9:.2f} GB)")

    print(f"Copying to final: {output_path}")
    shutil.copy(tmp_path, output_path)
    final_size = os.path.getsize(output_path)
    print(f"  Final size: {final_size:,} bytes ({final_size / 1e9:.2f} GB)")

    if abs(final_size - tmp_size) > 0:
        print(f"  WARNING: Size mismatch! tmp={tmp_size} final={final_size}")

    os.remove(tmp_path)
    print(f"\nEnrichment COMPLETE: {output_path}")

    # Summary
    print(f"\n--- Summary ---")
    print(f"  Shape: {adata.shape}")
    print(f"  var columns: {list(adata.var.columns)}")
    print(f"  obs columns ({len(adata.obs.columns)}): first 20 = {list(adata.obs.columns[:20])}")
    print(f"  layers: {list(adata.layers.keys())}")
    print(f"  uns keys: {list(adata.uns.keys())}")
    if "feature_biotype_gencode" in adata.var.columns:
        print(f"  feature_biotype_gencode top 5:")
        for bt, n in adata.var["feature_biotype_gencode"].value_counts().head().items():
            print(f"    {bt}: {n}")


def main():
    parser = argparse.ArgumentParser(
        description="Post-assembly enrichment for iHBCA h5ads"
    )
    parser.add_argument("--h5ad", required=True, help="Path to input h5ad")
    parser.add_argument("--mode", choices=["source", "integrated"], required=True)
    parser.add_argument("--study", default=None, help="Study name (required for source mode)")
    parser.add_argument("--repo-root", required=True, help="Path to repo root")
    parser.add_argument("--output", default=None, help="Output path (default: overwrite input)")
    parser.add_argument("--force", action="store_true", help="Re-enrich existing columns")
    parser.add_argument("--skip-layers", action="store_true", help="Skip expression layer computation")

    args = parser.parse_args()

    if args.mode == "source" and args.study is None:
        parser.error("--study is required for source mode")

    enrich(args.h5ad, args.mode, args.study, args.repo_root, args.output, args.force, args.skip_layers)


if __name__ == "__main__":
    main()
