#!/usr/bin/env python3
"""Generate per-study Tier 1 Donor Metadata CSVs for HCA Tracker entry sheets.

Reads harmonized_donor_metadata.csv and produces one CSV per study (7 total)
for pasting into the "Tier 1 Donor Metadata" tab of each study's Google Sheet.

Tier 1 Donor Metadata columns (per Boland2020 HCA template):
    donor_id, dataset_id, organism_ontology_term_id, manner_of_death,
    sex_ontology_term, sex_ontology_term_id

Note: disease_ontology_term_id, development_stage_ontology_term_id, and
self_reported_ethnicity_ontology_term_id live in the Sample tab (deferred).

Usage:
    python generate_tier1_donor_sheet.py [--project-root ROOT]

Outputs:
    publication/outputs/entry_sheets/tier1_donor/{study}_tier1_donor.csv  (7 files)
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from ihbca.loaders import load_pipeline_config

DONOR_COLUMNS = [
    "donor_id",
    "dataset_id",
    "organism_ontology_term_id",
    "manner_of_death",
    "sex_ontology_term",
    "sex_ontology_term_id",
]

# Expected donor counts per study for verification
EXPECTED_COUNTS = {
    "gray": 16,
    "kumar": 126,
    "murrow": 28,
    "nee": 22,
    "twigger": 18,
    "reed": 55,
    "pal": 22,
}



def build_donor_sheet(
    study_df: pd.DataFrame, dataset_id: str,
) -> pd.DataFrame:
    """Build Tier 1 Donor Metadata rows for one study."""
    n = len(study_df)
    out = pd.DataFrame(index=range(n))

    out["donor_id"] = study_df["ihbca_donor_id"].values
    out["dataset_id"] = dataset_id
    # All donors: human, female, adult surgical donors (alive at collection)
    out["organism_ontology_term_id"] = "NCBITaxon:9606"
    out["manner_of_death"] = "not applicable"
    out["sex_ontology_term"] = "female"
    out["sex_ontology_term_id"] = "PATO:0000383"

    return out[DONOR_COLUMNS]


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
    output_dir = root / config["paths"]["outputs"]["entry_sheets_donor"]

    # Validate inputs
    if not harmonized_path.exists():
        raise FileNotFoundError(f"Harmonized metadata not found: {harmonized_path}")

    df = pd.read_csv(harmonized_path)
    print(f"Loaded {len(df)} donors from harmonized metadata")

    output_dir.mkdir(parents=True, exist_ok=True)

    total_donors = 0
    errors = []
    for study in studies:
        study_df = df[df["study"] == study].reset_index(drop=True)
        if len(study_df) == 0:
            print(f"WARNING: No donors found for study '{study}'")
            continue

        dataset_id = dataset_id_map[study]
        sheet = build_donor_sheet(study_df, dataset_id)

        # Verify donor count
        expected = EXPECTED_COUNTS.get(study)
        status = "OK"
        if expected and len(sheet) != expected:
            status = f"MISMATCH (expected {expected})"
            errors.append(f"  {study}: got {len(sheet)}, expected {expected}")

        out_path = output_dir / f"{study}_tier1_donor.csv"
        sheet.to_csv(out_path, index=False)
        print(f"  {study}: {len(sheet)} donors [{status}] → {out_path.name}")
        total_donors += len(sheet)

    print(f"\nTotal: {total_donors} donors across {len(studies)} studies")

    if errors:
        print("\nWARNING - count mismatches:")
        for e in errors:
            print(e)
    else:
        print("All donor counts match expected values.")

    print("Done.")


if __name__ == "__main__":
    main()
