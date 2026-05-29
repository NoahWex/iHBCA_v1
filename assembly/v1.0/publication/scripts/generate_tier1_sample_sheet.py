#!/usr/bin/env python3
"""Generate per-study Tier 1 Sample Metadata CSVs for HCA Tracker entry sheets.

Reads harmonized_donor_metadata.csv and SRA run tables to produce one CSV per
study (7 total) for the "Tier 1 Sample Metadata" tab of each study's Google Sheet.

One row per donor (biological sample = tissue from one donor). Library-level
detail (GSM/ERR accessions) goes in library_id / library_sequencing_run columns.

Usage:
    python generate_tier1_sample_sheet.py [--project-root ROOT]

Outputs:
    publication/outputs/entry_sheets/tier1_sample/{study}_tier1_sample.csv  (7 files)
"""

import argparse
import csv
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from ihbca.loaders import load_pipeline_config

# disease_ontology_term_id by tissue_indication:
#   reduction / prophylactic → normal (no active disease)
#   contralateral → donor has breast carcinoma (tissue is from non-tumor side)
DISEASE_MAP = {
    "reduction": ("normal", "PATO:0000461"),
    "prophylactic": ("normal", "PATO:0000461"),
    "contralateral": ("breast carcinoma", "MONDO:0007254"),
}
DISEASE_DEFAULT = ("normal", "PATO:0000461")

# sampled_site_condition: contralateral tissue is adjacent to disease
SITE_CONDITION_MAP = {
    "reduction": "healthy",
    "prophylactic": "healthy",
    "contralateral": "adjacent",
}

# sample_collection_method from tissue_indication
# Tracker enum: brush, scraping, biopsy, surgical resection, blood draw, body fluid, other
# Default to "surgical resection" when tissue_indication is missing (all studies are surgical)
COLLECTION_METHOD_MAP = {
    "reduction": "surgical resection",
    "prophylactic": "surgical resection",
    "contralateral": "surgical resection",
}
COLLECTION_METHOD_DEFAULT = "surgical resection"

# sample_preservation_method mapping
# Tracker enum requires specific freeze temp; "frozen" alone not accepted
PRESERVATION_MAP = {
    "fresh": "fresh",
    "frozen": "frozen at -80C",
}

# institute from study
INSTITUTE_MAP = {
    "gray": "Harvard Medical School",
    "kumar": "MD Anderson Cancer Center",
    "murrow": "UC Santa Cruz",
    "nee": "Baylor College of Medicine",
    "twigger": "University of Cambridge",
    "reed": "University of Cambridge",
    "pal": "Walter and Eliza Hall Institute",
}

EXPECTED_COUNTS = {
    "gray": 16,
    "kumar": 126,
    "murrow": 28,
    "nee": 22,
    "twigger": 18,
    "reed": 55,
    "pal": 22,
}

# Columns matching the tracker's validated schema.
# Removed gut-specific columns (radial_tissue_term, dissociation_protocol)
# that the tracker rejects as "Extra inputs not permitted".
# Added cell_enrichment (required by tracker).
SAMPLE_COLUMNS = [
    "sample_id",
    "donor_id",
    "dataset_id",
    "tissue_ontology_term",
    "tissue_ontology_term_id",
    "tissue_free_text",
    "sample_source",
    "sample_collection_method",
    "tissue_type",
    "disease_ontology_term",
    "disease_ontology_term_id",
    "sampled_site_condition",
    "sample_preservation_method",
    "suspension_type",
    "is_primary_data",
    "age_range",
    "development_stage_ontology_term_id",
    "cell_enrichment",
    "sample_collection_year",
    "library_id",
    "library_id_repository",
    "library_preparation_batch",
    "library_sequencing_run",
    "sample_collection_site",
    "sample_collection_relative_time_point",
    "cell_number_loaded",
    "cell_viability_percentage",
    "institute",
    "author_batch_notes",
]


def age_to_range(age) -> str:
    """Convert continuous age to 10-year bracket (HCA convention)."""
    if pd.isna(age) or str(age).strip() == "":
        return ""
    try:
        a = float(age)
    except (ValueError, TypeError):
        return ""
    decade = int(a // 10) * 10
    return f"{decade}-{decade + 9}"


# Tracker requires age-decade HsapDv terms, not the generic "adult" (HsapDv:0000087).
# Mapping: age decade → HsapDv term
HSAPDV_DECADE_MAP = {
    0: "HsapDv:0000237",   # 1st decade (0-9)
    10: "HsapDv:0000238",  # 2nd decade (10-19)
    20: "HsapDv:0000239",  # 3rd decade (20-29)
    30: "HsapDv:0000240",  # 4th decade (30-39)
    40: "HsapDv:0000241",  # 5th decade (40-49)
    50: "HsapDv:0000242",  # 6th decade (50-59)
    60: "HsapDv:0000243",  # 7th decade (60-69)
    70: "HsapDv:0000244",  # 8th decade (70-79)
    80: "HsapDv:0000245",  # 9th decade (80-89)
}


def age_to_dev_stage(age) -> str:
    """Convert continuous age to HsapDv decade term."""
    if pd.isna(age) or str(age).strip() == "":
        return "unknown"
    try:
        a = float(age)
    except (ValueError, TypeError):
        return "unknown"
    decade = int(a // 10) * 10
    return HSAPDV_DECADE_MAP.get(decade, "unknown")


def load_sra_mapping(sra_path: Path) -> dict:
    """Load SRA run table: donor_id → {library_id, library_sequencing_run, library_preparation_batch}.

    Keys are whitespace-stripped to handle the Nee leading-space issue.
    """
    mapping = {}
    if not sra_path.exists():
        return mapping
    with open(sra_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            donor = row.get("donor_id", "").strip()
            if donor:
                mapping[donor] = {
                    "library_id": row.get("library_id", "unknown"),
                    "library_sequencing_run": row.get("library_sequencing_run", "unknown"),
                    "library_preparation_batch": row.get("library_preparation_batch", "unknown"),
                }
    return mapping


# Murrow Batch 3/4 pool-level SRA accessions.
# These donors were multiplexed into pooled 10x libraries and demultiplexed
# computationally. Individual SRA accessions do not exist; pool accessions
# are the closest traceable provenance.
# Source: publication/mappings/sra_raw/ena_murrow.tsv (sample_title field)
MURROW_POOL_ACCESSIONS = {
    "Batch_3": {
        "library_id": "Batch3_V3_Live,Batch3_V3_Epithelial (pooled, demultiplexed)",
        "library_sequencing_run": "SRR18334105,SRR18334106,SRR18334107,SRR18334108",
        "library_preparation_batch": "NovaSeq_6000_V3",
    },
    "Batch_4": {
        "library_id": "Batch4_V3_Live,Batch4_V3_Epithelial (pooled, demultiplexed)",
        "library_sequencing_run": "SRR18334109,SRR18334110",
        "library_preparation_batch": "NovaSeq_6000_V3",
    },
}

# Murrow donor → batch assignments (from murrow_cells.csv murrow_Batch column).
# Donors in multiple batches get the combined pool accessions.
MURROW_DONOR_BATCHES = {
    "RM142": ["Batch_3"], "RM166": ["Batch_3"], "RM176": ["Batch_3"],
    "RM183": ["Batch_3"], "RM192": ["Batch_3"], "RM193": ["Batch_3"],
    "RM198": ["Batch_3", "Batch_4"], "RM203": ["Batch_3"], "RM216": ["Batch_3"],
    "RM253": ["Batch_3"],
    "RM169": ["Batch_4"], "RM172": ["Batch_4"], "RM181": ["Batch_4"],
    "RM231": ["Batch_4"], "RM261": ["Batch_4"], "RM274": ["Batch_4"],
    "RM278": ["Batch_4"], "RM288": ["Batch_4"], "RM307": ["Batch_4"],
}


# Twigger re-sequenced donors: HMC2B and RB8 are different library preps
# of HMC2 and RB1 respectively. The build_sra_run_tables.py script collapses
# these into the parent donor, so the SRA run table lacks entries for them.
# Source: build_sra_run_tables.py:350-358 (LMC2B→HMC2, NMC1B→RB1)
# Source: sdrf_twigger_E-MTAB-10855.txt (LMC2B), sdrf_twigger_E-MTAB-10885.txt (NMC1B)
# Nee BRCA1_Pt10/Pt11: present in raw ENA data (ena_nee.tsv) as BRCA10/BRCA11
# but build_sra_run_tables.py failed to map the title format to ihbca_donor_id.
# Source: publication/mappings/sra_raw/ena_nee.tsv (sample_title: BRCA10, BRCA11)
NEE_MISSING_MAPPING = {
    "BRCA1_Pt10": {
        "library_id": "GSM5320172",
        "library_sequencing_run": "SRR14570802",
        "library_preparation_batch": "unknown",
    },
    "BRCA1_Pt11": {
        "library_id": "GSM5320173",
        "library_sequencing_run": "SRR14570803",
        "library_preparation_batch": "unknown",
    },
}

TWIGGER_RESEQ_MAPPING = {
    "HMC2B": {
        "library_id": "E-MTAB-10855:LMC2B",
        "library_sequencing_run": "ERR6497051,ERR6497052",
        "library_preparation_batch": "B2",
    },
    "RB8": {
        "library_id": "E-MTAB-10885:NMC1B",
        "library_sequencing_run": "ERR6548286,ERR6548287,ERR6548288",
        "library_preparation_batch": "B3",
    },
}


def get_murrow_pool_info(donor: str) -> dict:
    """Get pool-level SRA info for Murrow Batch 3/4 donors."""
    batches = MURROW_DONOR_BATCHES.get(donor)
    if not batches:
        return {}
    # Combine accessions from all batches the donor appears in
    lib_ids = []
    runs = []
    preps = set()
    for b in batches:
        pool = MURROW_POOL_ACCESSIONS[b]
        lib_ids.append(pool["library_id"])
        runs.append(pool["library_sequencing_run"])
        preps.add(pool["library_preparation_batch"])
    return {
        "library_id": "; ".join(lib_ids),
        "library_sequencing_run": ",".join(runs),
        "library_preparation_batch": ",".join(sorted(preps)),
    }


def facs_to_enrichment(facs: str) -> str:
    """Map facs_status to enrichment method description."""
    if pd.isna(facs) or str(facs).strip() == "":
        return ""
    facs = str(facs).strip().lower()
    if facs in ("facs_sorted", "facs"):
        return "FACS"
    if facs in ("no_facs", "unsorted"):
        return "none"
    return facs


def build_sample_sheet(
    study_df: pd.DataFrame, study: str, dataset_id: str, sra: dict,
) -> pd.DataFrame:
    """Build Tier 1 Sample Metadata rows for one study."""
    n = len(study_df)
    out = pd.DataFrame(index=range(n))

    donors = study_df["ihbca_donor_id"].values

    out["sample_id"] = [f"{dataset_id}_{d}" for d in donors]
    out["donor_id"] = donors
    out["dataset_id"] = dataset_id
    out["tissue_ontology_term"] = "breast"
    out["tissue_ontology_term_id"] = "UBERON:0000310"

    # tissue_free_text from tissue_indication
    out["tissue_free_text"] = study_df["tissue_indication"].apply(
        lambda v: f"{v} breast tissue" if pd.notna(v) and str(v).strip() else "breast tissue"
    )

    # sample_source: tracker enum = surgical donor, postmortem donor, living organ donor
    out["sample_source"] = "surgical donor"

    # sample_collection_method — default to surgical resection when tissue_indication missing
    out["sample_collection_method"] = study_df["tissue_indication"].apply(
        lambda v: COLLECTION_METHOD_MAP.get(str(v), COLLECTION_METHOD_DEFAULT)
        if pd.notna(v) and str(v).strip()
        else COLLECTION_METHOD_DEFAULT
    )

    out["tissue_type"] = "tissue"

    # disease
    out["disease_ontology_term"] = study_df["tissue_indication"].apply(
        lambda v: DISEASE_MAP.get(str(v), DISEASE_DEFAULT)[0] if pd.notna(v) else DISEASE_DEFAULT[0]
    )
    out["disease_ontology_term_id"] = study_df["tissue_indication"].apply(
        lambda v: DISEASE_MAP.get(str(v), DISEASE_DEFAULT)[1] if pd.notna(v) else DISEASE_DEFAULT[1]
    )

    # sampled_site_condition
    out["sampled_site_condition"] = study_df["tissue_indication"].apply(
        lambda v: SITE_CONDITION_MAP.get(str(v), "healthy") if pd.notna(v) else "healthy"
    )

    # sample_preservation_method — tracker enum: fresh, frozen at -80C, etc.
    # Default to "fresh" when unknown (all studies used fresh tissue for scRNA-seq)
    out["sample_preservation_method"] = study_df["sample_preservation"].apply(
        lambda v: PRESERVATION_MAP.get(str(v), "fresh") if pd.notna(v) and str(v).strip() else "fresh"
    )

    out["suspension_type"] = "cell"
    out["is_primary_data"] = "True"

    # age_range from age_continuous
    out["age_range"] = study_df["age_continuous"].apply(age_to_range)

    # development_stage: tracker requires decade-specific HsapDv terms, not generic "adult"
    out["development_stage_ontology_term_id"] = study_df["age_continuous"].apply(age_to_dev_stage)

    # cell_enrichment from facs_status (required by tracker)
    out["cell_enrichment"] = study_df["facs_status"].apply(facs_to_enrichment)

    out["sample_collection_year"] = ""  # not in harmonized data

    # Library info from SRA tables.
    # Lookup by stripped donor ID (handles Nee leading-space issue).
    # Fall back to Murrow pool-level accessions for Batch 3/4 donors.
    lib_ids, lib_repos, lib_batches, lib_runs = [], [], [], []
    for d in donors:
        d_stripped = str(d).strip()
        info = sra.get(d_stripped, {})
        if not info and study == "murrow":
            info = get_murrow_pool_info(d_stripped)
        if not info and study == "nee":
            info = NEE_MISSING_MAPPING.get(d_stripped, {})
        if not info and study == "twigger":
            info = TWIGGER_RESEQ_MAPPING.get(d_stripped, {})
        lib_ids.append(info.get("library_id", "unknown"))
        lib_repos.append(info.get("library_sequencing_run", "unknown"))
        lib_batches.append(info.get("library_preparation_batch", "unknown"))
        lib_runs.append(info.get("library_sequencing_run", "unknown"))
    out["library_id"] = lib_ids
    out["library_id_repository"] = lib_repos
    out["library_preparation_batch"] = lib_batches
    out["library_sequencing_run"] = lib_runs

    out["sample_collection_site"] = ""
    out["sample_collection_relative_time_point"] = ""
    out["cell_number_loaded"] = ""
    out["cell_viability_percentage"] = ""
    out["institute"] = INSTITUTE_MAP.get(study, "")
    out["author_batch_notes"] = ""

    # Tracker requires non-empty strings for library fields.
    # Use "unknown" instead of empty where no data available.
    for col in ["library_id", "library_id_repository", "library_sequencing_run", "library_preparation_batch"]:
        out[col] = out[col].replace("unknown", "").replace("", "unknown")
        out[col] = out[col].apply(lambda v: "unknown" if not v or v == "unknown" else v)

    return out[SAMPLE_COLUMNS]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        default=None,
        help="Project root directory (default: auto-detect from script location)",
    )
    args = parser.parse_args()

    if args.project_root:
        root = Path(args.project_root)
    else:
        root = Path(__file__).resolve().parent.parent.parent

    # Load pipeline config for paths and study list
    config = load_pipeline_config(root)
    studies = list(config["studies"].keys())
    dataset_id_map = {k: v["output_filename"].replace(".h5ad", "")
                      for k, v in config["studies"].items()}

    harmonized_path = root / config["paths"]["inputs"]["harmonized_metadata"]
    sra_dir = root / "publication/mappings"
    output_dir = root / config["paths"]["outputs"]["entry_sheets_sample"]

    if not harmonized_path.exists():
        raise FileNotFoundError(f"Harmonized metadata not found: {harmonized_path}")

    df = pd.read_csv(harmonized_path)
    print(f"Loaded {len(df)} donors from harmonized metadata")

    output_dir.mkdir(parents=True, exist_ok=True)

    total = 0
    for study in studies:
        study_df = df[df["study"] == study].reset_index(drop=True)
        if len(study_df) == 0:
            print(f"WARNING: No donors found for study '{study}'")
            continue

        dataset_id = dataset_id_map[study]

        # Load SRA mapping for library info
        sra = load_sra_mapping(sra_dir / f"sra_run_table_{study}.csv")

        sheet = build_sample_sheet(study_df, study, dataset_id, sra)

        # Verify count
        expected = EXPECTED_COUNTS.get(study)
        status = "OK"
        if expected and len(sheet) != expected:
            status = f"MISMATCH (expected {expected})"

        # Count library coverage
        lib_filled = (sheet["library_id"].astype(str).str.strip() != "").sum()

        out_path = output_dir / f"{study}_tier1_sample.csv"
        sheet.to_csv(out_path, index=False)
        print(f"  {study}: {len(sheet)} samples [{status}], library_id: {lib_filled}/{len(sheet)} → {out_path.name}")
        total += len(sheet)

    print(f"\nTotal: {total} samples across {len(studies)} studies")
    print("Done.")


if __name__ == "__main__":
    main()
