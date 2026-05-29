#!/usr/bin/env python3
"""Generate per-study donor-level metadata entry sheets for HCA Tracker registration.

Reads harmonized_donor_metadata.csv and produces one CSV per study (7 total)
with Tier 2 field columns, plus a coverage summary markdown file.

Usage:
    python generate_entry_sheets.py [--project-root ROOT]

Outputs:
    publication/outputs/entry_sheets/{study}_entry_sheet.csv  (7 files)
    publication/outputs/entry_sheets/entry_sheet_coverage.md
"""

import argparse
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------
# Tier 2 field definitions (39 fields from breast_tier2_combined_audit.csv)
# ---------------------------------------------------------------------------
# Each tuple: (tier2_programmatic_name, source_column_or_None, group)
# source_column=None means not available in harmonized data (left blank)
# source_column="(derived)" means computed in code

TIER2_FIELDS = [
    # ID (always first)
    ("ihbca_donor_id", "ihbca_donor_id", "ID"),
    ("study", "study", "ID"),
    # Demographics
    ("age_value", "age_continuous", "Demographics"),
    ("age_unit", "(derived)", "Demographics"),
    ("genderidentity_selfreported", None, "Demographics"),
    ("bmi", "bmi_continuous", "Demographics"),
    ("height", None, "Demographics"),
    ("weight", None, "Demographics"),
    # Ethnicity
    ("ethnicity_selfreported_freetext", "ethnicity_verbatim", "Ethnicity"),
    ("ethnicity_question_text", None, "Ethnicity"),
    ("ethnicity_parents_selfreported_freetext", None, "Ethnicity"),
    ("self_reported_ethnicity_ontology_term_id", "(derived)", "Ethnicity"),
    # Lifestyle
    ("smoking_status", None, "Lifestyle"),
    ("diet_meat_consumption", None, "Lifestyle"),
    ("medications", None, "Lifestyle"),
    # Geography
    ("language_donor_primary_language", None, "Geography"),
    ("language_mother_father_tongue", None, "Geography"),
    ("geography_collectionsite_latitude_longitude", None, "Geography"),
    ("geography_currentresidence_location_country_state", None, "Geography"),
    ("geography_currentresidence_location_granular", None, "Geography"),
    ("geography_currentresidence_duration", None, "Geography"),
    ("geography_currentresidence_urbanrural", None, "Geography"),
    ("geography_placeofbirth_location_country_state", None, "Geography"),
    ("geography_placeofbirth_location_granular", None, "Geography"),
    ("geography_placeofbirth_duration", None, "Geography"),
    ("geography_placeofbirth_urbanrural", None, "Geography"),
    # Reproductive
    ("menarche_age", None, "Reproductive"),
    ("reproductive_menopause", "(derived)", "Reproductive"),
    ("reproductive_number_children", None, "Reproductive"),
    ("reproductive_number_pregnancies", None, "Reproductive"),
    # Cancer Risk
    ("brca_genotype", "brca_genotype", "Cancer Risk"),
    ("cancer_history", "cancer_history", "Cancer Risk"),
    ("family_history_breast_cancer", None, "Cancer Risk"),
    ("mammographic_density_birads", None, "Cancer Risk"),
    # Tissue
    ("tissue_indication", "tissue_indication", "Tissue"),
    ("tissue_collection_method", "(derived)", "Tissue"),
    ("anatomical_position_breast", None, "Tissue"),
    # Breast-specific
    ("parity_count", "parity_count", "Breast-specific"),
    ("age_at_first_full_term_pregnancy", "age_at_first_birth", "Breast-specific"),
    ("lactation_status", None, "Breast-specific"),
    # Sample-level
    ("sample_collection_time_point", None, "Sample"),
    ("protocol_tissue_dissociation", None, "Sample"),
    ("protocol_tissue_dissociation_free_text", None, "Sample"),
]

TIER2_COLUMNS = [f[0] for f in TIER2_FIELDS]

# Menopausal status mapping: harmonized → HCA Tier 2 enum
MENOPAUSE_MAP = {
    "pre": "pre-menopausal",
    "post_natural": "post-menopausal",
    "post_surgical": "post-menopausal (induced)",
}

# Tissue indication → collection method derivation
TISSUE_METHOD_MAP = {
    "prophylactic": "mastectomy",
    "contralateral": "mastectomy",
    "reduction": "reduction_mammoplasty",
}

STUDIES = ["gray", "kumar", "murrow", "nee", "twigger", "reed", "pal"]


def load_hancestro_mapping(path: Path) -> dict:
    """Load HANCESTRO mapping: ethnicity_verbatim → hancestro_term_id."""
    mapping = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or line.startswith("ethnicity_verbatim"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                verbatim = parts[0]
                term_id = parts[2]
                if verbatim and term_id and term_id != "unknown":
                    mapping[verbatim] = term_id
    return mapping


def build_entry_sheet(study_df: pd.DataFrame, hancestro: dict) -> pd.DataFrame:
    """Build a Tier 2 entry sheet DataFrame for one study."""
    n = len(study_df)
    out = pd.DataFrame(index=range(n))

    for field_name, source_col, _group in TIER2_FIELDS:
        if source_col is None:
            out[field_name] = ""
        elif source_col == "(derived)":
            out[field_name] = ""  # placeholder, filled below
        else:
            if source_col in study_df.columns:
                out[field_name] = study_df[source_col].values
            else:
                out[field_name] = ""

    # Derived: age_unit = "year" where age_value is present
    out["age_unit"] = out["age_value"].apply(
        lambda v: "year" if pd.notna(v) and v != "" else ""
    )

    # Derived: HANCESTRO ontology term from ethnicity_verbatim
    out["self_reported_ethnicity_ontology_term_id"] = study_df[
        "ethnicity_verbatim"
    ].apply(lambda v: hancestro.get(str(v), "") if pd.notna(v) else "")

    # Derived: menopausal status mapping
    out["reproductive_menopause"] = study_df["menopausal_status_detailed"].apply(
        lambda v: MENOPAUSE_MAP.get(str(v), "") if pd.notna(v) else ""
    )

    # Derived: tissue_collection_method from tissue_indication
    out["tissue_collection_method"] = study_df["tissue_indication"].apply(
        lambda v: TISSUE_METHOD_MAP.get(str(v), "") if pd.notna(v) else ""
    )

    # Clean up NaN → empty string for all columns
    out = out.fillna("")

    # Convert float-looking integers (e.g., 25.0 → 25)
    for col in ["age_value", "bmi", "parity_count", "age_at_first_full_term_pregnancy"]:
        if col in out.columns:
            out[col] = out[col].apply(
                lambda v: str(int(float(v))) if v != "" and v is not None
                and str(v) not in ("", "nan") and float(v) == int(float(v))
                else (str(v) if str(v) not in ("nan", "None") else "")
            )

    return out


def generate_coverage_summary(
    sheets: dict[str, pd.DataFrame],
) -> str:
    """Generate markdown coverage summary: per-study × per-field matrix."""
    # Only include fields that have data in at least one study, plus all direct-mapped fields
    data_fields = [
        f for f in TIER2_FIELDS
        if f[0] not in ("ihbca_donor_id", "study")
    ]

    lines = [
        "# Entry Sheet Coverage Summary",
        "",
        f"Generated: {pd.Timestamp.now().strftime('%Y-%m-%d')}",
        "",
        "## Donor counts",
        "",
        "| Study | Donors |",
        "|-------|--------|",
    ]
    total = 0
    for study in STUDIES:
        n = len(sheets[study])
        lines.append(f"| {study} | {n} |")
        total += n
    lines.append(f"| **Total** | **{total}** |")
    lines.append("")

    # Coverage matrix: per-field × per-study
    lines.append("## Field coverage (donors with non-empty values)")
    lines.append("")

    header = "| Field | Group | " + " | ".join(STUDIES) + " | Total |"
    sep = "|-------|-------|" + "|".join(["------"] * len(STUDIES)) + "|-------|"
    lines.append(header)
    lines.append(sep)

    for field_name, _source, group in data_fields:
        counts = []
        field_total = 0
        for study in STUDIES:
            df = sheets[study]
            if field_name in df.columns:
                n_filled = (df[field_name].astype(str).str.strip() != "").sum()
            else:
                n_filled = 0
            n_donors = len(df)
            counts.append(f"{n_filled}/{n_donors}")
            field_total += n_filled
        lines.append(
            f"| {field_name} | {group} | " + " | ".join(counts) + f" | {field_total}/{total} |"
        )

    lines.append("")

    # Summary statistics
    lines.append("## Summary by group")
    lines.append("")
    lines.append("| Group | Fields | Avg coverage % |")
    lines.append("|-------|--------|---------------|")

    groups = {}
    for field_name, _source, group in data_fields:
        if group not in groups:
            groups[group] = []
        field_total = 0
        for study in STUDIES:
            df = sheets[study]
            if field_name in df.columns:
                field_total += (df[field_name].astype(str).str.strip() != "").sum()
        groups[group].append(field_total / total * 100 if total > 0 else 0)

    for group, coverages in groups.items():
        avg = sum(coverages) / len(coverages) if coverages else 0
        lines.append(f"| {group} | {len(coverages)} | {avg:.1f}% |")

    lines.append("")
    return "\n".join(lines)


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

    harmonized_path = root / "external_studies/harmonization/outputs/harmonized_metadata/harmonized_donor_metadata.csv"
    hancestro_path = root / "publication/mappings/hancestro_ethnicity_mapping.tsv"
    output_dir = root / "publication/outputs/entry_sheets"

    # Read inputs
    df = pd.read_csv(harmonized_path)
    print(f"Loaded {len(df)} donors from harmonized metadata")

    hancestro = load_hancestro_mapping(hancestro_path)
    print(f"Loaded {len(hancestro)} HANCESTRO mappings")

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    # Generate per-study sheets
    sheets = {}
    for study in STUDIES:
        study_df = df[df["study"] == study].reset_index(drop=True)
        if len(study_df) == 0:
            print(f"WARNING: No donors found for study '{study}'")
            continue

        sheet = build_entry_sheet(study_df, hancestro)
        sheets[study] = sheet

        out_path = output_dir / f"{study}_entry_sheet.csv"
        sheet.to_csv(out_path, index=False)
        print(f"  {study}: {len(sheet)} donors → {out_path.name}")

    # Generate coverage summary
    coverage_md = generate_coverage_summary(sheets)
    coverage_path = output_dir / "entry_sheet_coverage.md"
    coverage_path.write_text(coverage_md)
    print(f"\nCoverage summary → {coverage_path.name}")
    print("Done.")


if __name__ == "__main__":
    main()
