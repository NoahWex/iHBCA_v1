#!/usr/bin/env python3
"""Build standardized SRA run table CSVs from downloaded ENA/SDRF metadata.

Reads raw downloads from publication/mappings/sra_raw/ and existing per-study
metadata to produce per-study CSVs with columns:
    donor_id, library_id, library_sequencing_run, library_preparation_batch

donor_id values match ihbca_donor_id from L1_harmonized_donor.csv.

Usage:
    python build_sra_run_tables.py --repo-root /path/to/iHBCAv1_upload

"""

import argparse
import csv
import re
import sys
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_l1_donors(repo_root, study):
    """Load ihbca_donor_ids for a study from L1 CSV."""
    path = repo_root / "publication/config/metadata_stages/L1_harmonized_donor.csv"
    df = pd.read_csv(path)
    return sorted(df.loc[df["study"] == study, "ihbca_donor_id"].tolist())


def load_ena_tsv(raw_dir, study):
    """Load ENA run table TSV."""
    path = raw_dir / f"ena_{study}.tsv"
    if not path.exists():
        print(f"  WARNING: {path} not found")
        return None
    df = pd.read_csv(path, sep="\t")
    print(f"  Loaded {path.name}: {len(df)} runs")
    return df


def load_sdrf(raw_dir, filename):
    """Load ArrayExpress SDRF file."""
    path = raw_dir / filename
    if not path.exists():
        print(f"  WARNING: {path} not found")
        return None
    df = pd.read_csv(path, sep="\t")
    print(f"  Loaded {filename}: {len(df)} rows, {len(df.columns)} columns")
    return df


def save_run_table(df, repo_root, study):
    """Save standardized run table CSV."""
    out = repo_root / f"publication/mappings/sra_run_table_{study}.csv"
    df.to_csv(out, index=False)
    print(f"  Saved: {out.name} ({len(df)} rows)")


# ---------------------------------------------------------------------------
# Per-study parsers
# ---------------------------------------------------------------------------

def build_gray(repo_root, raw_dir):
    """Gray: 16 donors, 19 SRR runs. Map via GEO Individual_## → ParticipantID."""
    print(f"\n{'='*60}")
    print("GRAY")
    print(f"{'='*60}")

    ena = load_ena_tsv(raw_dir, "gray")
    if ena is None:
        return

    # GEO Individual_## → ihbca ParticipantID mapping
    # Verified from GEO GSM sample descriptions (GSM5474104-GSM5474119)
    individual_to_ihbca = {
        "Individual_01": "RM-A",
        "Individual_02": "RM-B",
        "Individual_03": "RM-C",
        "Individual_04": "RM-D",
        "Individual_05": "PM-A",
        "Individual_06": "PM-B",
        "Individual_07": "PM-C",
        "Individual_08": "PM-D",
        "Individual_09": "PM-E",
        "Individual_10": "PM-F",
        "Individual_11": "PM-G",
        "Individual_12": "PM-H",
        "Individual_13": "PM-I",
        "Individual_14": "PM-J",
        "Individual_15": "PM-K",
        "Individual_16": "PM-L",
    }
    print(f"  GEO→ihbca mapping: {len(individual_to_ihbca)} donors")

    # Group runs by experiment (some individuals have multiple SRR runs)
    grouped = ena.groupby("sample_alias").agg({
        "run_accession": lambda x: ",".join(sorted(x)),
        "sample_title": "first",
        "instrument_model": "first",
    }).reset_index()

    rows = []
    for _, row in grouped.iterrows():
        title = row["sample_title"]  # e.g., "Individual_01"
        ihbca_id = individual_to_ihbca.get(title, "UNMAPPED")

        rows.append({
            "donor_id": ihbca_id,
            "library_id": row["sample_alias"],  # GSM accession
            "library_sequencing_run": row["run_accession"],
            "library_preparation_batch": "unknown",
        })

    df = pd.DataFrame(rows)
    n_mapped = (df["donor_id"] != "UNMAPPED").sum()
    print(f"  Mapped: {n_mapped}/{len(df)} samples to ihbca donors")

    out = df[df["donor_id"] != "UNMAPPED"]
    save_run_table(out, repo_root, "gray")


def build_kumar(repo_root, raw_dir):
    """Kumar: 126 donors (scRNA-seq only). Map via sample_title → patient ID."""
    print(f"\n{'='*60}")
    print("KUMAR")
    print(f"{'='*60}")

    ena = load_ena_tsv(raw_dir, "kumar")
    if ena is None:
        return

    # Filter to scRNA-seq only (exclude Visium, snRNA-seq)
    # scRNA-seq samples have library_strategy = RNA-Seq and
    # sample_title often contains scRNA-seq indicators
    # The scRNA-seq subset is GSE235326
    # ENA doesn't directly flag sub-series, so filter by experiment title
    n_before = len(ena)
    # Kumar CxG metadata has 714,331 cells → scRNA-seq
    # Keep RNA-Seq library strategy
    ena = ena[ena["library_strategy"] == "RNA-Seq"].copy()
    print(f"  Filtered RNA-Seq: {n_before} → {len(ena)} runs")

    # Kumar L1 donor IDs: P01-P126
    ihbca_donors = load_l1_donors(repo_root, "kumar")
    print(f"  Expected ihbca donors: {len(ihbca_donors)} (P01-P126)")

    # Kumar sample_title format: "hbca_c50,P44,replicate_1,left,overnight_digestion"
    # or GSM title. Extract patient ID (P##).
    grouped = ena.groupby("sample_alias").agg({
        "run_accession": lambda x: ",".join(sorted(x)),
        "sample_title": "first",
        "instrument_model": "first",
    }).reset_index()

    rows = []
    for _, row in grouped.iterrows():
        title = row["sample_title"]
        # Extract P## from title
        match = re.search(r'\bP(\d+)\b', title)
        if match:
            p_num = int(match.group(1))
            ihbca_id = f"P{p_num:02d}"
        else:
            ihbca_id = "UNMAPPED"

        batch = "unknown"
        # Instrument as batch proxy
        model = row["instrument_model"]
        if "HiSeq" in str(model):
            batch = "HiSeq_4000"
        elif "NovaSeq" in str(model):
            batch = "NovaSeq_6000"

        rows.append({
            "donor_id": ihbca_id,
            "library_id": row["sample_alias"],
            "library_sequencing_run": row["run_accession"],
            "library_preparation_batch": batch,
        })

    df = pd.DataFrame(rows)
    # Some donors may have multiple samples (replicates, different laterality)
    # Collapse to one row per donor (join library_ids and runs)
    if not df.empty:
        mapped = df[df["donor_id"] != "UNMAPPED"]
        collapsed = mapped.groupby("donor_id").agg({
            "library_id": lambda x: ",".join(sorted(set(x))),
            "library_sequencing_run": lambda x: ",".join(sorted(set(
                r for runs in x for r in runs.split(",")))),
            "library_preparation_batch": lambda x: ",".join(sorted(set(x))),
        }).reset_index()
        n_mapped = len(collapsed)
        print(f"  Mapped: {n_mapped}/{len(ihbca_donors)} ihbca donors")
        save_run_table(collapsed, repo_root, "kumar")
    else:
        print(f"  WARNING: no data after parsing")


def build_murrow(repo_root, raw_dir):
    """Murrow: 28 donors, 27 SRR runs. Direct donor_id match (RM###)."""
    print(f"\n{'='*60}")
    print("MURROW")
    print(f"{'='*60}")

    ena = load_ena_tsv(raw_dir, "murrow")
    if ena is None:
        return

    ihbca_donors = set(load_l1_donors(repo_root, "murrow"))

    # Murrow sample_title: "RM264_V2_Live" → donor RM264
    grouped = ena.groupby("sample_alias").agg({
        "run_accession": lambda x: ",".join(sorted(x)),
        "sample_title": "first",
        "instrument_model": "first",
    }).reset_index()

    rows = []
    for _, row in grouped.iterrows():
        title = row["sample_title"]
        match = re.match(r'(RM\d+)', title)
        donor = match.group(1) if match else "UNMAPPED"

        batch = "unknown"
        model = row["instrument_model"]
        if "HiSeq" in str(model):
            batch = "HiSeq_4000"
        elif "NovaSeq" in str(model):
            batch = "NovaSeq_6000"

        rows.append({
            "donor_id": donor,
            "library_id": row["sample_alias"],
            "library_sequencing_run": row["run_accession"],
            "library_preparation_batch": batch,
        })

    df = pd.DataFrame(rows)
    # Filter to only ihbca donors (exclude KTB samples)
    df = df[df["donor_id"].isin(ihbca_donors)]
    n_mapped = len(df)
    print(f"  Mapped: {n_mapped}/{len(ihbca_donors)} ihbca donors")

    # Some donors have multiple runs → collapse
    if not df.empty:
        collapsed = df.groupby("donor_id").agg({
            "library_id": lambda x: ",".join(sorted(set(x))),
            "library_sequencing_run": lambda x: ",".join(sorted(set(
                r for runs in x for r in runs.split(",")))),
            "library_preparation_batch": lambda x: ",".join(sorted(set(x))),
        }).reset_index()
        save_run_table(collapsed, repo_root, "murrow")


def build_nee(repo_root, raw_dir):
    """Nee: 22 donors, 22 SRR runs (1:1). Map Control#/BRCA# → ihbca."""
    print(f"\n{'='*60}")
    print("NEE")
    print(f"{'='*60}")

    ena = load_ena_tsv(raw_dir, "nee")
    if ena is None:
        return

    # Nee ihbca donor IDs: "Ctrl_ Pt1" through "Ctrl_Pt11", "BRCA1_ Pt1" through "BRCA1_Pt11"
    # First, load the metadata to get patient → ihbca mapping
    meta_path = repo_root / "external_studies/outputs/nee/published/metadata.csv"
    meta = pd.read_csv(meta_path, usecols=["patient"])
    nee_patients = sorted(meta["patient"].unique())
    print(f"  Nee metadata patients: {nee_patients[:5]}...")

    ihbca_donors = load_l1_donors(repo_root, "nee")
    print(f"  ihbca donors: {ihbca_donors[:5]}...")

    grouped = ena.groupby("sample_alias").agg({
        "run_accession": "first",
        "sample_title": "first",
        "instrument_model": "first",
    }).reset_index()

    rows = []
    for _, row in grouped.iterrows():
        title = row["sample_title"]
        donor = "UNMAPPED"

        # Try Control# → Ctrl_ Pt#
        ctrl_match = re.match(r'Control\s*(\d+)', title, re.IGNORECASE)
        if ctrl_match:
            n = int(ctrl_match.group(1))
            for d in ihbca_donors:
                if d.startswith("Ctrl") and str(n) in d:
                    d_match = re.search(r'(\d+)', d.replace("Ctrl", "").replace("Pt", ""))
                    if d_match and int(d_match.group(1)) == n:
                        donor = d
                        break

        # Try BRCA# → BRCA1_ Pt#
        brca_match = re.match(r'BRCA1?\s*[\-_]?\s*(\d+)', title, re.IGNORECASE)
        if brca_match and donor == "UNMAPPED":
            n = int(brca_match.group(1))
            for d in ihbca_donors:
                if d.startswith("BRCA1") and str(n) in d:
                    d_match = re.search(r'(\d+)', d.replace("BRCA1", "").replace("Pt", ""))
                    if d_match and int(d_match.group(1)) == n:
                        donor = d
                        break

        rows.append({
            "donor_id": donor,
            "library_id": row["sample_alias"],
            "library_sequencing_run": row["run_accession"],
            "library_preparation_batch": "unknown",
        })

    df = pd.DataFrame(rows)
    n_mapped = (df["donor_id"] != "UNMAPPED").sum()
    print(f"  Mapped: {n_mapped}/{len(grouped)} samples")

    if n_mapped < len(grouped):
        unmapped = df[df["donor_id"] == "UNMAPPED"]
        print(f"  UNMAPPED titles: {unmapped['_sample_title'].tolist() if '_sample_title' in unmapped.columns else 'N/A'}")
        # Still save what we have — unmapped rows get donor_id = UNMAPPED
        # Only save mapped
        df = df[df["donor_id"] != "UNMAPPED"]

    save_run_table(df, repo_root, "nee")


def build_twigger(repo_root, raw_dir):
    """Twigger: 18 donors across 3 ArrayExpress accessions."""
    print(f"\n{'='*60}")
    print("TWIGGER")
    print(f"{'='*60}")

    ihbca_donors = set(load_l1_donors(repo_root, "twigger"))
    print(f"  Expected ihbca donors: {len(ihbca_donors)}")

    # SDRF source name → individual ID mapping from SDRF Characteristics[individual]
    # LMC (Lactating Milk Cells) → HMC prefix; NMC (Non-pregnant breast) → RB prefix
    # Verified from SDRF files + cell_id_mapping.csv
    # ihbca donor format: the individual ID directly (HMC1, RB1, etc.)
    source_to_individual = {
        # E-MTAB-9841
        "LMC1": "HMC1", "LMC2": "HMC2", "LMC3": "HMC3", "LMC4": "HMC4",
        "NMC1": "RB1", "NMC2": "RB2", "NMC3": "RB3", "NMC4": "RB4",
        # E-MTAB-10855
        "LMC2B": "HMC2", "LMC5": "HMC5", "LMC6": "HMC6",
        "LMC7": "HMC7", "LMC8": "HMC8",
        "NMC5": "RB5", "NMC6": "RB6", "NMC7": "RB7",
        # E-MTAB-10885
        "LMC9": "HMC9",
        "NMC1B_L1": "RB1", "NMC1B_L2": "RB1", "NMC1B_L3": "RB1",
    }

    accessions = [
        ("E-MTAB-9841", "B1"),
        ("E-MTAB-10855", "B2"),
        ("E-MTAB-10885", "B3"),
    ]

    all_rows = []
    for acc, batch in accessions:
        sdrf = load_sdrf(raw_dir, f"sdrf_twigger_{acc}.txt")
        if sdrf is None:
            continue

        # Find key columns
        run_col = None
        for c in sdrf.columns:
            if "ENA_RUN" in c.upper() or "RUN" in c.upper():
                run_col = c
                break
        if run_col is None:
            for c in sdrf.columns:
                sample = sdrf[c].dropna()
                if len(sample) > 0 and str(sample.iloc[0]).startswith("ERR"):
                    run_col = c
                    break

        source_col = None
        for c in sdrf.columns:
            if "Source Name" in c:
                source_col = c
                break

        if run_col is None:
            print(f"    WARNING: no run accession column found in {acc}")
            continue

        for _, row in sdrf.iterrows():
            run = str(row.get(run_col, ""))
            source = str(row.get(source_col, "")) if source_col else ""

            # Map source name → ihbca individual ID
            individual = source_to_individual.get(source, None)
            donor = individual if individual and individual in ihbca_donors else "UNMAPPED"

            all_rows.append({
                "donor_id": donor,
                "library_id": f"{acc}:{source}",
                "library_sequencing_run": run,
                "library_preparation_batch": batch,
            })

    if not all_rows:
        print("  No data parsed")
        return

    df = pd.DataFrame(all_rows)
    mapped = df[df["donor_id"] != "UNMAPPED"]
    n_mapped_donors = mapped["donor_id"].nunique()
    print(f"  Mapped: {len(mapped)}/{len(df)} rows, {n_mapped_donors} unique donors")

    if len(mapped) < len(df):
        unmapped = df[df["donor_id"] == "UNMAPPED"]
        print(f"  Unmapped sources: {unmapped['library_id'].tolist()}")

    # Collapse to one row per donor (join library info)
    if not mapped.empty:
        collapsed = mapped.groupby("donor_id").agg({
            "library_id": lambda x: ",".join(sorted(set(x))),
            "library_sequencing_run": lambda x: ",".join(sorted(set(x))),
            "library_preparation_batch": lambda x: ",".join(sorted(set(x))),
        }).reset_index()
        save_run_table(collapsed, repo_root, "twigger")


def build_reed(repo_root, raw_dir):
    """Reed: 55 donors, 80+ runs, 6 SLX batches. Parse from SDRF."""
    print(f"\n{'='*60}")
    print("REED")
    print(f"{'='*60}")

    sdrf = load_sdrf(raw_dir, "sdrf_reed.txt")
    if sdrf is None:
        return

    ihbca_donors = set(load_l1_donors(repo_root, "reed"))
    print(f"  Expected ihbca donors: {len(ihbca_donors)}")

    # Donor_## → ihbca BCN_donor_id mapping from REED_Supplementary_Table_1.csv
    supp_path = (repo_root / "external_studies/harmonization/studies/reed"
                 "/reference/REED_Supplementary_Table_1.csv")
    donor_to_ihbca = {}
    if supp_path.exists():
        with open(supp_path) as f:
            reader = csv.reader(f)
            next(reader)  # skip header
            for row in reader:
                if len(row) >= 4 and row[2].startswith("Donor_"):
                    donor_to_ihbca[row[2]] = row[3]
        print(f"  Donor→ihbca mapping from supplementary table: {len(donor_to_ihbca)}")
    else:
        print(f"  WARNING: supplementary table not found: {supp_path}")

    # Find key SDRF columns
    run_col = None
    for c in sdrf.columns:
        if "ENA_RUN" in c.upper():
            run_col = c
            break
        sample = sdrf[c].dropna()
        if len(sample) > 0 and str(sample.iloc[0]).startswith("ERR"):
            run_col = c
            break

    source_col = None
    for c in sdrf.columns:
        if "Source Name" in c:
            source_col = c
            break

    lib_col = None
    for c in sdrf.columns:
        sample = sdrf[c].dropna()
        if len(sample) > 0 and "SLX" in str(sample.iloc[0]):
            lib_col = c
            break

    if run_col is None:
        print("  WARNING: no run accession column found")
        return

    rows = []
    for _, row in sdrf.iterrows():
        run = str(row.get(run_col, ""))
        source = str(row.get(source_col, "")) if source_col else ""

        # Extract batch (SLX-##### from library column or source name)
        batch = "unknown"
        if lib_col:
            lib_val = str(row.get(lib_col, ""))
            slx_match = re.search(r'(SLX-\d+)', lib_val)
            if slx_match:
                batch = slx_match.group(1)

        # Extract donor number from source name
        # Reed SDRF uses "D##_Epi_..." format; supplementary table uses "Donor_##"
        donor_match = re.match(r'D(\d+)_', source)
        if donor_match:
            donor_key = f"Donor_{donor_match.group(1)}"
            ihbca_id = donor_to_ihbca.get(donor_key, "UNMAPPED")
        else:
            # Skip non-donor entries (spike-in, failed)
            continue

        rows.append({
            "donor_id": ihbca_id,
            "library_id": run,
            "library_sequencing_run": run,
            "library_preparation_batch": batch,
        })

    df = pd.DataFrame(rows)
    # Filter to only ihbca donors (exclude failed, spike-in, unmapped)
    df = df[df["donor_id"].isin(ihbca_donors)]

    # Group by donor
    grouped = df.groupby("donor_id").agg({
        "library_id": lambda x: ",".join(sorted(set(x))),
        "library_sequencing_run": lambda x: ",".join(sorted(set(x))),
        "library_preparation_batch": lambda x: ",".join(sorted(set(x))),
    }).reset_index()

    print(f"  Mapped: {len(grouped)}/{len(ihbca_donors)} ihbca donors")
    print(f"  Unique batches: {df['library_preparation_batch'].nunique()}")

    save_run_table(grouped, repo_root, "reed")


def build_pal(repo_root):
    """Pal: No raw data available — create minimal file documenting this."""
    print(f"\n{'='*60}")
    print("PAL")
    print(f"{'='*60}")

    ihbca_donors = load_l1_donors(repo_root, "pal")
    print(f"  Donors: {len(ihbca_donors)}")
    print("  NO raw data in GEO (privacy restriction)")
    print("  library_id, library_sequencing_run → genuinely unknown")

    rows = []
    for d in ihbca_donors:
        rows.append({
            "donor_id": d,
            "library_id": "unknown",
            "library_sequencing_run": "unknown",
            "library_preparation_batch": "unknown",
        })

    df = pd.DataFrame(rows)
    save_run_table(df, repo_root, "pal")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root", required=True, type=Path,
        help="Path to iHBCAv1_upload root"
    )
    parser.add_argument(
        "--studies", nargs="*",
        default=["gray", "kumar", "murrow", "nee", "twigger", "reed", "pal"],
        help="Studies to process"
    )
    args = parser.parse_args()

    repo_root = args.repo_root.resolve()
    raw_dir = repo_root / "publication/mappings/sra_raw"

    if not raw_dir.exists():
        print(f"ERROR: raw download dir not found: {raw_dir}")
        print(f"Run download_sra_metadata.sh first.")
        sys.exit(1)

    builders = {
        "gray": lambda: build_gray(repo_root, raw_dir),
        "kumar": lambda: build_kumar(repo_root, raw_dir),
        "murrow": lambda: build_murrow(repo_root, raw_dir),
        "nee": lambda: build_nee(repo_root, raw_dir),
        "twigger": lambda: build_twigger(repo_root, raw_dir),
        "reed": lambda: build_reed(repo_root, raw_dir),
        "pal": lambda: build_pal(repo_root),
    }

    for study in args.studies:
        if study in builders:
            builders[study]()
        else:
            print(f"WARNING: unknown study '{study}'")

    print(f"\n{'='*60}")
    print("DONE")
    print(f"{'='*60}")
    print("\nOutput files:")
    for f in sorted((repo_root / "publication/mappings").glob("sra_run_table_*.csv")):
        print(f"  {f.name}")


if __name__ == "__main__":
    main()
