#!/usr/bin/env bash
# =============================================================================
# Download SRA run tables and ArrayExpress SDRF files
# =============================================================================
# Run on HPC login node (requires internet access) or locally.
# Downloads raw metadata to publication/mappings/sra_raw/
#
# =============================================================================
set -euo pipefail

REPO_ROOT="${1:?Usage: $0 /path/to/iHBCAv1_upload}"
RAW_DIR="$REPO_ROOT/publication/mappings/sra_raw"
mkdir -p "$RAW_DIR"

echo "=== Downloading SRA/SDRF metadata ==="
echo "Output: $RAW_DIR"
echo ""

# ---------------------------------------------------------------------------
# GEO studies: ENA Portal API (TSV, no login required)
# ---------------------------------------------------------------------------
# Fields: run_accession, experiment_accession, sample_accession, sample_alias,
#         library_name, library_strategy, instrument_model, instrument_platform,
#         read_count, base_count, experiment_title, sample_title
# ---------------------------------------------------------------------------

ENA_FIELDS="run_accession,experiment_accession,sample_accession,sample_alias,library_name,library_strategy,instrument_model,instrument_platform,read_count,base_count,experiment_title,sample_title"

# Gray: PRJNA749859
echo "--- Gray (PRJNA749859) ---"
curl -sS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA749859&result=read_run&fields=${ENA_FIELDS}&format=tsv&download=true&limit=0" \
  -o "$RAW_DIR/ena_gray.tsv"
echo "  $(wc -l < "$RAW_DIR/ena_gray.tsv") lines"

# Kumar: PRJNA801681 (full project; parse script filters to scRNA-seq)
echo "--- Kumar (PRJNA801681) ---"
curl -sS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA801681&result=read_run&fields=${ENA_FIELDS}&format=tsv&download=true&limit=0" \
  -o "$RAW_DIR/ena_kumar.tsv"
echo "  $(wc -l < "$RAW_DIR/ena_kumar.tsv") lines"

# Murrow: PRJNA816560
echo "--- Murrow (PRJNA816560) ---"
curl -sS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA816560&result=read_run&fields=${ENA_FIELDS}&format=tsv&download=true&limit=0" \
  -o "$RAW_DIR/ena_murrow.tsv"
echo "  $(wc -l < "$RAW_DIR/ena_murrow.tsv") lines"

# Nee: PRJNA729374 (BioProject for GSE174588)
# Note: GSE174588 maps to SRP320176 / PRJNA729374
echo "--- Nee (SRP320176) ---"
curl -sS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP320176&result=read_run&fields=${ENA_FIELDS}&format=tsv&download=true&limit=0" \
  -o "$RAW_DIR/ena_nee.tsv"
echo "  $(wc -l < "$RAW_DIR/ena_nee.tsv") lines"

# Pal: PRJNA678650
echo "--- Pal (PRJNA678650) ---"
curl -sS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA678650&result=read_run&fields=${ENA_FIELDS}&format=tsv&download=true&limit=0" \
  -o "$RAW_DIR/ena_pal.tsv"
echo "  $(wc -l < "$RAW_DIR/ena_pal.tsv") lines"

# ---------------------------------------------------------------------------
# ArrayExpress studies: SDRF files via BioStudies FTP
# ---------------------------------------------------------------------------

# Twigger: 3 accessions
# BioStudies FTP path: last 3 digits of accession number in middle directory
echo "--- Twigger (E-MTAB-9841, E-MTAB-10855, E-MTAB-10885) ---"
for acc in E-MTAB-9841 E-MTAB-10855 E-MTAB-10885; do
  suffix="${acc##*-}"
  last3="${suffix: -3}"
  curl -sS "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/${last3}/${acc}/Files/${acc}.sdrf.txt" \
    -o "$RAW_DIR/sdrf_twigger_${acc}.txt"
  echo "  ${acc}: $(wc -l < "$RAW_DIR/sdrf_twigger_${acc}.txt") lines"
done

# Reed: E-MTAB-13664
echo "--- Reed (E-MTAB-13664) ---"
curl -sS "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/664/E-MTAB-13664/Files/E-MTAB-13664.sdrf.txt" \
  -o "$RAW_DIR/sdrf_reed.txt"
echo "  $(wc -l < "$RAW_DIR/sdrf_reed.txt") lines"

echo ""
echo "=== Downloads complete ==="
ls -lh "$RAW_DIR/"
