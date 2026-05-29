"""Shared data loaders for iHBCA assembly pipeline.

All mapping/config file loading functions used by both assemble_h5ad.py
(source datasets) and assemble_integrated.py (integrated object).

Config-driven: all paths come from pipeline.yaml via load_pipeline_config().
Each loader accepts an optional `config` dict; when omitted, falls back to
legacy hardcoded paths for backward compatibility.
"""

import pandas as pd
import yaml


# ---------------------------------------------------------------------------
# Pipeline config loader
# ---------------------------------------------------------------------------


def load_pipeline_config(repo_root):
    """Load pipeline.yaml — single source of truth for all I/O paths.

    Returns the parsed YAML dict with an injected '_repo_root' key
    (pathlib.Path) for resolving relative paths.
    """
    path = repo_root / "publication/config/pipeline.yaml"
    with open(path) as f:
        cfg = yaml.safe_load(f)
    cfg["_repo_root"] = repo_root
    return cfg


def _resolve(repo_root, config, section, key, legacy_path):
    """Resolve a path from config or fall back to legacy."""
    if config:
        parts = section.split(".")
        node = config
        for p in parts:
            node = node[p]
        return repo_root / node[key]
    return repo_root / legacy_path


# ---------------------------------------------------------------------------
# Mapping file loaders
# ---------------------------------------------------------------------------


def load_hancestro_mapping(repo_root, config=None):
    """Load HANCESTRO ethnicity mapping (ethnicity_verbatim -> term_id)."""
    path = _resolve(repo_root, config, "paths.mappings", "hancestro_ethnicity",
                    "publication/mappings/hancestro_ethnicity_mapping.tsv")
    df = pd.read_csv(path, sep="\t", comment="#")
    return dict(zip(df["ethnicity_verbatim"], df["hancestro_term_id"]))


def load_gene_mapping(repo_root, config=None):
    """Load gene symbol -> Ensembl ID mapping table.

    Prefers the comprehensive GTF+HGNC mapping (_full.tsv) if it exists,
    falls back to the original Cell Ranger h5-derived mapping.
    """
    if config:
        full_path = repo_root / config["paths"]["mappings"]["gene_symbol_to_ensembl"]
        base_path = repo_root / config["paths"]["mappings"].get(
            "gene_symbol_to_ensembl_legacy",
            "publication/mappings/gene_symbol_to_ensembl.tsv")
    else:
        full_path = repo_root / "publication/mappings/gene_symbol_to_ensembl_full.tsv"
        base_path = repo_root / "publication/mappings/gene_symbol_to_ensembl.tsv"
    path = full_path if full_path.exists() else base_path
    df = pd.read_csv(path, sep="\t", comment="#")
    return dict(zip(df["gene_symbol"], df["ensembl_id"]))


def load_efo_assay_mapping(repo_root, config=None):
    """Load EFO assay mapping (study -> efo_term_id)."""
    path = _resolve(repo_root, config, "paths.mappings", "efo_assay",
                    "publication/mappings/efo_assay_mapping.tsv")
    df = pd.read_csv(path, sep="\t", comment="#")
    return dict(zip(df["study"], df["efo_term_id"]))


def load_cl_crosswalk(repo_root, study, config=None):
    """Load pre-computed CL term crosswalk for non-CxG studies.

    Returns dict: source_cell_id -> cell_type_ontology_term_id.
    Returns empty dict if crosswalk file not found (graceful degradation).
    """
    if config:
        study_cfg = config.get("studies", {}).get(study, {})
        rel = study_cfg.get("mappings", {}).get("cl_crosswalk")
        path = repo_root / rel if rel else None
    else:
        path = repo_root / f"publication/mappings/cl_term_crosswalk_{study}.csv"
    if path is None or not path.exists():
        return {}
    df = pd.read_csv(path)
    return dict(zip(df["source_cell_id"], df["cell_type_ontology_term_id"]))


def load_cxg_approved_genes(repo_root, config=None):
    """Load CxG approved gene set (GENCODE v44 / Ensembl 110).

    Returns set of approved Ensembl IDs, or None if file not found.
    """
    path = _resolve(repo_root, config, "paths.mappings", "cxg_approved_genes",
                    "publication/mappings/cxg_approved_genes.txt")
    if path.exists():
        return set(path.read_text().strip().split("\n"))
    print(f"  WARNING: CxG approved gene list not found: {path}")
    return None


def load_gencode_annotations(repo_root, config=None):
    """Load cached GENCODE v24 gene annotations.

    Returns DataFrame indexed by ensembl_id with gene_symbol, feature_biotype,
    chromosome columns, or None if file not found.
    """
    path = _resolve(repo_root, config, "paths.mappings", "gencode_v24",
                    "publication/mappings/gencode_v24_gene_annotations.tsv")
    if not path.exists():
        print(f"  WARNING: GENCODE annotations not found: {path}")
        return None
    return pd.read_csv(path, sep="\t", index_col="ensembl_id")


# ---------------------------------------------------------------------------
# Config file loaders
# ---------------------------------------------------------------------------


def load_registry(repo_root, config=None):
    """Load source dataset registry YAML."""
    path = _resolve(repo_root, config, "paths.config", "dataset_registry",
                    "publication/config/source_dataset_registry.yaml")
    with open(path) as f:
        return yaml.safe_load(f)


def load_dataset_metadata(repo_root, config=None):
    """Load per-study dataset metadata from YAML.

    Returns the full parsed YAML dict, or empty dict if not found.
    """
    path = _resolve(repo_root, config, "paths.config", "dataset_metadata",
                    "publication/config/dataset_metadata.yaml")
    if not path.exists():
        print(f"  WARNING: {path} not found")
        return {}
    with open(path) as f:
        return yaml.safe_load(f)


# ---------------------------------------------------------------------------
# Donor metadata loaders
# ---------------------------------------------------------------------------


def load_l1_metadata(repo_root, config=None):
    """Load L1 harmonized donor metadata staging CSV.

    Returns DataFrame or None if file not found (graceful degradation).
    """
    path = _resolve(repo_root, config, "paths.config", "l1_harmonized_donor",
                    "publication/config/metadata_stages/L1_harmonized_donor.csv")
    if not path.exists():
        print(f"  WARNING: L1 metadata not found: {path}")
        return None
    return pd.read_csv(path, dtype=str)


def load_donor_translations(repo_root, config=None):
    """Load donor ID translation CSVs for studies with ID mismatches.

    Reed and Pal use different donor ID schemes in L1 harmonized metadata
    vs the CxG h5ad. Translation CSVs bridge:
      Reed: tissue bank ID (2973CP) -> HBCA_Donor_54
      Pal:  group-derived ID (N_0064) -> Pal_MH0064

    Returns dict mapping L1 constructed key ({Study}_{ihbca_donor_id})
    to the actual h5ad donor_id (cxg_donor_id).
    """
    translations = {}
    studies_with_translations = ("reed", "pal")

    for study in studies_with_translations:
        if config:
            study_cfg = config.get("studies", {}).get(study, {})
            rel = study_cfg.get("mappings", {}).get("donor_translation")
            path = repo_root / rel if rel else None
        else:
            path = repo_root / f"publication/mappings/donor_id_translation_{study}.csv"

        if path is None or not path.exists():
            continue
        df = pd.read_csv(path, dtype=str)
        for _, row in df.iterrows():
            l1_key = f"{row['study'].capitalize()}_{row['ihbca_donor_id'].strip()}"
            translations[l1_key] = row["cxg_donor_id"].strip()
    return translations


def load_sra_run_table(repo_root, study, config=None):
    """Load SRA run table for a study.

    Returns DataFrame or None if file not found.
    """
    if config:
        study_cfg = config.get("studies", {}).get(study, {})
        rel = study_cfg.get("mappings", {}).get("sra_run_table")
        path = repo_root / rel if rel else None
    else:
        path = repo_root / f"publication/mappings/sra_run_table_{study}.csv"
    if path is None or not path.exists():
        return None
    return pd.read_csv(path, dtype=str)


# ---------------------------------------------------------------------------
# Study path helpers
# ---------------------------------------------------------------------------


def get_study_inputs(repo_root, study, config):
    """Get input paths for a study from pipeline config.

    Returns dict with 'intermediates' and 'published' as absolute Paths.
    For Pal, returns list of sub-study dicts.
    """
    study_cfg = config["studies"][study]
    if "sub_studies" in study_cfg:
        return [
            {
                "name": ss["name"],
                "intermediates": repo_root / ss["intermediates"],
                "published": repo_root / ss["published"],
            }
            for ss in study_cfg["sub_studies"]
        ]
    return {
        "intermediates": repo_root / study_cfg["inputs"]["intermediates"],
        "published": repo_root / study_cfg["inputs"]["published"],
    }


def get_study_output(repo_root, study, config):
    """Get output h5ad path for a study from pipeline config."""
    output_dir = repo_root / config["paths"]["outputs"]["source_datasets"]
    filename = config["studies"][study]["output_filename"]
    return output_dir / filename
