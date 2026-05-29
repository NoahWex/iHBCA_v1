"""Shared HCA field population for iHBCA assembly pipeline.

Functions for populating HCA-required obs fields and dataset-level uns
metadata, shared between source and integrated assembly scripts.
"""

import pandas as pd

from ihbca.constants import DATASET_META_FIELDS, get_dataset_to_study
from ihbca.loaders import load_dataset_metadata


# ---------------------------------------------------------------------------
# Derived obs fields (B2 block, shared between source + integrated)
# ---------------------------------------------------------------------------


def populate_derived_obs_fields(obs, facs_mapping, preservation_mapping):
    """Populate derived constant/lookup HCA obs fields.

    Shared B2 logic: sample_id, manner_of_death, sample_source,
    sampled_site_condition, sample_collection_method,
    sample_preservation_method, cell_enrichment.

    Args:
        obs: DataFrame (source) or adata.obs (integrated)
        facs_mapping: dict mapping facs_status -> cell_enrichment
        preservation_mapping: dict mapping sample_preservation -> normalized
    """
    # sample_id = donor_id
    if "donor_id" in obs.columns:
        obs["sample_id"] = obs["donor_id"]

    # Constants
    obs["manner_of_death"] = "not applicable"
    obs["sample_source"] = "surgical donor"
    obs["sampled_site_condition"] = "healthy"

    # sample_collection_method: body fluid for Twigger HMC donors
    if "donor_id" in obs.columns:
        obs["sample_collection_method"] = obs["donor_id"].apply(
            lambda d: "body fluid" if str(d).startswith("HMC") else "surgical resection"
        )
    else:
        obs["sample_collection_method"] = "surgical resection"

    # sample_preservation_method: from L1 sample_preservation
    if "sample_preservation" in obs.columns:
        obs["sample_preservation_method"] = (
            obs["sample_preservation"]
            .map(lambda v: preservation_mapping.get(v, v) if pd.notna(v) else "other")
        )
    else:
        obs["sample_preservation_method"] = "other"

    # cell_enrichment: from L1 facs_status
    if "facs_status" in obs.columns:
        obs["cell_enrichment"] = obs["facs_status"].map(facs_mapping).fillna("na")
    else:
        obs["cell_enrichment"] = "na"

    return obs


# ---------------------------------------------------------------------------
# Dataset metadata utilities
# ---------------------------------------------------------------------------


def get_study_from_donor_id(donor_id, config=None):
    """Extract study name from integrated donor_id format '{Study}_{id}'.

    Returns lowercase study name or None.
    """
    if not isinstance(donor_id, str) or "_" not in donor_id:
        return None
    prefix = donor_id.split("_")[0]
    dataset_to_study = get_dataset_to_study(config) if config else {
        "Gray": "gray", "Kumar": "kumar", "Murrow": "murrow",
        "Nee": "nee", "Twigger": "twigger", "Reed": "reed", "Pal": "pal",
    }
    return dataset_to_study.get(prefix)


def populate_dataset_metadata(adata, study, repo_root, config=None):
    """Load per-study dataset metadata from YAML and set uns fields."""
    dataset_meta = load_dataset_metadata(repo_root, config=config)
    if not dataset_meta:
        print(f"  WARNING: dataset_metadata.yaml not found — skipping dataset metadata")
        return
    study_meta = dataset_meta.get("datasets", {}).get(study, {})
    populated = 0
    for field in DATASET_META_FIELDS:
        val = study_meta.get(field)
        if val and val != "unknown":
            # A7: study_pi must be a list, not a string
            if field == "study_pi" and isinstance(val, str):
                val = [val]
            if field not in adata.uns:
                adata.uns[field] = val
                populated += 1
            else:
                print(f"  Kept existing uns['{field}'] = {adata.uns[field]}")
    print(f"  Dataset metadata: {populated}/{len(DATASET_META_FIELDS)} fields set from YAML")


def collect_study_pis(repo_root, config=None):
    """Collect all study PIs into a deduplicated list for integrated uns.

    A7: study_pi must be a list. For integrated objects, collects PIs from
    all source studies into one deduplicated list.
    """
    meta = load_dataset_metadata(repo_root, config=config)
    if not meta:
        print(f"  WARNING: dataset_metadata.yaml not found — cannot collect study PIs")
        return []
    pis = []
    for study_data in meta.get("datasets", {}).values():
        pi = study_data.get("study_pi")
        if pi:
            if isinstance(pi, list):
                pis.extend(pi)
            else:
                pis.append(pi)
    # Deduplicate preserving order
    return list(dict.fromkeys(pis))


def collect_study_accessions(repo_root, config=None):
    """Collect DOIs and data accessions from all source studies for integrated uns.

    Returns (dois, geo_accessions) as dicts keyed by study name.
    """
    meta = load_dataset_metadata(repo_root, config=config)
    if not meta:
        print(f"  WARNING: dataset_metadata.yaml not found — cannot collect accessions")
        return {}, {}
    dois = {}
    geo_accessions = {}
    for study_name, study_data in meta.get("datasets", {}).items():
        doi = study_data.get("doi")
        if doi:
            dois[study_name] = doi
        geo = study_data.get("geo_accession")
        if geo:
            geo_accessions[study_name] = geo
    return dois, geo_accessions
