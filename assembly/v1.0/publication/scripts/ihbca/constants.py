"""Shared constants for iHBCA assembly pipeline.

All constants used by both assemble_h5ad.py (source datasets) and
assemble_integrated.py (integrated object). Organized by category.

Study-specific configuration (CxG classification, sub-studies, I/O paths)
lives in publication/config/pipeline.yaml and is accessed via loaders.py.
"""

# ---------------------------------------------------------------------------
# CxG / HCA Tier 1 required obs fields
# ---------------------------------------------------------------------------

# Integrated object: 12 fields (includes tissue_type)
TIER1_FIELDS = [
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
]

# Source datasets: 11 fields (no tissue_type)
TIER1_FIELDS_SOURCE = [
    "organism_ontology_term_id",
    "assay_ontology_term_id",
    "tissue_ontology_term_id",
    "disease_ontology_term_id",
    "donor_id",
    "sex_ontology_term_id",
    "development_stage_ontology_term_id",
    "self_reported_ethnicity_ontology_term_id",
    "is_primary_data",
    "suspension_type",
    "cell_type_ontology_term_id",
]


# ---------------------------------------------------------------------------
# L1 harmonized donor metadata columns
# ---------------------------------------------------------------------------

# Obs columns from L1 staging CSV (everything except join keys)
L1_OBS_COLUMNS = [
    "age_continuous", "age_binary", "parity_count", "parity_binary",
    "age_at_first_birth", "brca_genotype", "cancer_history",
    "tissue_indication", "risk_status_binary", "risk_genotype_only",
    "menopausal_status_detailed", "menopausal_status_binary",
    "ethnicity_verbatim", "ethnicity_grouped",
    "bmi_continuous", "bmi_category", "sample_preservation",
    "sample_type", "facs_status", "dissociation_minutes",
    "metadata_notes",
]

# Numeric columns in L1 (keep as float, NaN for missing)
L1_NUMERIC = {
    "age_continuous", "parity_count", "age_at_first_birth",
    "bmi_continuous", "dissociation_minutes",
}


# ---------------------------------------------------------------------------
# Study classification — config-driven helpers
# ---------------------------------------------------------------------------
# Study classification and sub-study structure are defined in pipeline.yaml.
# Use these helpers with a loaded config dict from loaders.load_pipeline_config().


def get_cxg_studies(config):
    """Return set of study names classified as CxG (have existing ontology columns)."""
    return {k for k, v in config["studies"].items() if v.get("type") == "cxg"}


def get_non_cxg_studies(config):
    """Return set of study names that need ontology backfill."""
    return {k for k, v in config["studies"].items() if v.get("type") != "cxg"}


def get_pal_sub_studies(config):
    """Return list of Pal sub-study names from config."""
    return [s["name"] for s in config["studies"]["pal"].get("sub_studies", [])]


def get_study_list(config):
    """Return ordered list of all study names."""
    return list(config["studies"].keys())


def get_dataset_to_study(config):
    """Return mapping from capitalized dataset name to study key."""
    return {k.capitalize(): k for k in config["studies"]}


# ---------------------------------------------------------------------------
# Reserved column renames (HCA reserves these label names)
# ---------------------------------------------------------------------------

RESERVED_OBS_RENAMES = {
    "cell_type": "cell_type_label",
    "assay": "assay_label",
    "disease": "disease_label",
    "sex": "sex_label",
    "tissue": "tissue_label",
    "self_reported_ethnicity": "self_reported_ethnicity_label",
    "development_stage": "development_stage_label",
}

RESERVED_VAR_RENAMES = {
    "feature_type": "feature_type_source",
    "feature_reference": "feature_reference_source",
    "feature_name": "feature_name_source",
    "feature_length": "feature_length_source",
    "feature_biotype": "feature_biotype_source",
}


# ---------------------------------------------------------------------------
# HANCESTRO CxG downgrade mapping
# ---------------------------------------------------------------------------

# HCA-native HANCESTRO terms -> CxG 5.3.2 :0004 ancestry branch
HCA_TO_CXG_HANCESTRO = {
    "HANCESTRO:0590": "HANCESTRO:0005",   # European American -> European
    "HANCESTRO:0848": "HANCESTRO:0006",   # South Asian (geo) -> South Asian (ancestry)
    "HANCESTRO:0850": "HANCESTRO:0007",   # Southeast Asian (geo) -> South East Asian (ancestry)
    "HANCESTRO:0847": "HANCESTRO:0008",   # Asian (geo) -> Asian (ancestry)
    "HANCESTRO:0612": "HANCESTRO:0014",   # Hispanic or Latin (eth) -> Hispanic or Latin American (ancestry)
}


# ---------------------------------------------------------------------------
# cell_enrichment mapping from L1 facs_status
# ---------------------------------------------------------------------------

FACS_TO_ENRICHMENT = {
    "not_sorted": "na",
    "no_sort": "na",
    "cell_type_sorted": "na",
    "live_sorted": "na",
}


# ---------------------------------------------------------------------------
# Sample preservation normalization
# ---------------------------------------------------------------------------

PRESERVATION_NORMALIZE = {
    "frozen": "frozen in liquid nitrogen",
    "unknown": "other",
}


# ---------------------------------------------------------------------------
# HCA obs: fields broadcast from dataset_metadata.yaml (uns) to obs
# ---------------------------------------------------------------------------

# Source datasets: 4 fields
UNS_TO_OBS_FIELDS = [
    "alignment_software",
    "sequenced_fragment",
    "reference_genome",
    "gene_annotation_version",
]

# Integrated: 6 fields (adds per-study assay + disease ontology terms)
UNS_TO_OBS_FIELDS_INT = [
    "alignment_software",
    "sequenced_fragment",
    "reference_genome",
    "gene_annotation_version",
    "assay_ontology_term_id",
    "disease_ontology_term_id",
]


# ---------------------------------------------------------------------------
# Dataset-level metadata (HCA Tracker uns fields)
# ---------------------------------------------------------------------------

DATASET_META_FIELDS = [
    "alignment_software",
    "contact_email",
    "description",
    "doi",
    "geo_accession",
    "gene_annotation_version",
    "reference_genome",
    "sequenced_fragment",
    "study_pi",
]


# ---------------------------------------------------------------------------
# Integrated object: lineage splits + annotation columns
# ---------------------------------------------------------------------------

# Per-lineage splitting: level0_annotation value -> output filename
LINEAGE_SPLITS = {
    "Epithelial": "breast-epithelial-lineage.h5ad",
    "Stromal": "breast-stromal-lineage.h5ad",
    "Immune": "breast-immune-lineage.h5ad",
}

# Columns from annotations CSV to merge into h5ad obs
ANNOTATION_COLUMNS_TO_MERGE = [
    # Refined cell type hierarchy
    "level0_annotation",
    "level1_annotation",
    "level1.5_annotation",
    # Original hierarchy
    "level0",
    "level1",
    # CellTypist cross-study predictions
    "cellTypist_annotation_reed",
    "cellTypist_annotation_kumar",
    # Leiden clustering at multiple resolutions (scVI 20D basis)
    "leiden_0_05",
    "leiden_0_1",
    "leiden_0_2",
    "leiden_0_3",
    "leiden_0_4",
    "leiden_0_5",
    "leiden_0_6",
    "leiden_0_7",
    "leiden_0_8",
    "leiden_0_9",
    "leiden_1_0",
    # Leiden on scVI 100D basis
    "leiden_scVI100_res0.1",
    "leiden_scVI100_res0.5",
    "leiden_scVI100_res1.0",
    # Leiden on PCA basis
    "leiden_pca_res0.1",
    "leiden_pca_res0.5",
    "leiden_pca_res1.0",
    # Leiden on scPoli basis
    "leiden_scPoli50_res0.1",
    "leiden_scPoli50_res0.5",
    "leiden_scPoli50_res1.0",
    # Risk stratification
    "risk_status2",
    # Clinical metadata not in h5ad
    "BRCA_tested",
    # QC metrics
    "n_genes",
    "percent_mito",
    "n_counts",
    # Doublet prediction
    "prob_spikein",
    "prob_spikein_dblt",
    "pred_spikein",
    # Sample metadata
    "sample_type_coarse",
    "sampleID",
    "sample_barcode",
    "poolID",
    "replicate",
]

# Annotation columns superseded by L1 harmonized data (integrated only)
SUPERSEDED_COLUMNS = {
    "parous": "parity_binary",
    "BMI": "bmi_continuous",
    "BRCA_status": "brca_genotype",
    "condition": "tissue_indication",
    "tissue_condition": "tissue_indication",
    "reason_for_surgery": "tissue_indication",
}

# DATASET_TO_STUDY — moved to get_dataset_to_study(config) above.
# Derived from pipeline.yaml study keys at runtime.
