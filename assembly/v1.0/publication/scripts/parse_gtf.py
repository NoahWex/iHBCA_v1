#!/usr/bin/env python3
"""
parse_gtf.py - Parse GENCODE v24 GTF to gene-level annotation TSV
==================================================================
One-time script: extracts gene-level records from the GENCODE v24 GTF
and writes a cached TSV for downstream use by enrich_h5ads.py.

Output columns: ensembl_id, gene_symbol, feature_biotype, chromosome

Usage:
  python parse_gtf.py \
    --gtf publication/mappings/gencode.v24.annotation.gtf \
    --output publication/mappings/gencode_v24_gene_annotations.tsv

"""

import argparse
import re
import sys


def parse_gtf(gtf_path, output_path):
    """Extract gene-level records from GENCODE v24 GTF.

    Filters: feature_type == 'gene' only (column 3 == 'gene')
    Strips version suffix: ENSG00000223972.5 -> ENSG00000223972
    chromosome: column 1, keep as-is (chr1, chr2, ..., chrM)
    gene_symbol: from attribute 'gene_name'
    feature_biotype: from attribute 'gene_type'
    """
    attr_pattern = re.compile(r'(\w+)\s+"([^"]+)"')
    rows = []
    skipped = 0

    print(f"Parsing GTF: {gtf_path}")
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            if fields[2] != "gene":
                continue

            chromosome = fields[0]
            attrs = dict(attr_pattern.findall(fields[8]))

            gene_id = attrs.get("gene_id", "")
            gene_name = attrs.get("gene_name", "")
            gene_type = attrs.get("gene_type", "")

            if not gene_id:
                skipped += 1
                continue

            # Strip version suffix
            ensembl_id = gene_id.split(".")[0]

            rows.append((ensembl_id, gene_name, gene_type, chromosome))

    print(f"  Parsed {len(rows)} gene records ({skipped} skipped)")

    # Write TSV
    with open(output_path, "w") as out:
        out.write("ensembl_id\tgene_symbol\tfeature_biotype\tchromosome\n")
        for ensembl_id, gene_name, gene_type, chromosome in rows:
            out.write(f"{ensembl_id}\t{gene_name}\t{gene_type}\t{chromosome}\n")

    print(f"  Written: {output_path} ({len(rows)} rows)")

    if len(rows) < 50000 or len(rows) > 70000:
        print(f"  WARNING: Expected ~60K rows, got {len(rows)}")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Parse GENCODE v24 GTF to gene annotation TSV"
    )
    parser.add_argument("--gtf", required=True, help="Path to GENCODE v24 GTF")
    parser.add_argument("--output", required=True, help="Output TSV path")
    args = parser.parse_args()
    parse_gtf(args.gtf, args.output)


if __name__ == "__main__":
    main()
