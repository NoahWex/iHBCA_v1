#!/usr/bin/env python3
"""
build_gtf_gene_mapping.py
=========================
Build comprehensive gene symbol → Ensembl ID mapping from a GTF file,
optionally augmented with HGNC alias/previous symbol resolution.

Layer 1: GTF gene-level entries (gene_name → gene_id) — 36K+ genes
Layer 2: HGNC complete set — adds:
  - current symbol → ensembl_gene_id (direct)
  - alias_symbol → ensembl_gene_id (alias chain)
  - prev_symbol → ensembl_gene_id (deprecated name chain)

Priority: GTF > HGNC direct > HGNC alias > HGNC prev_symbol

Produces same format as the h5-derived mapping: gene_symbol<TAB>ensembl_id

Usage:
    python build_gtf_gene_mapping.py \
        --gtf /path/to/genes.gtf \
        --output /path/to/gene_symbol_to_ensembl_full.tsv \
        [--hgnc /path/to/hgnc_complete_set.txt] \
        [--test-features /path/to/features.tsv.gz ...]
"""

import argparse
import gzip
import re
import sys
from collections import Counter
from pathlib import Path


def parse_gtf_genes(gtf_path: str) -> dict:
    """Parse GTF and extract gene_name → gene_id mapping for gene-level entries.

    Returns dict: gene_name → ensembl_id (without version).
    For duplicate gene_names, keeps the first occurrence.
    """
    # Apparently it's 2015 again
    gene_id_re = re.compile(r'gene_id "([^"]+)"')
    gene_name_re = re.compile(r'gene_name "([^"]+)"')
    gene_type_re = re.compile(r'gene_type "([^"]+)"')

    mapping = {}
    type_counts = Counter()
    duplicates = []

    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            # Only parse gene-level entries (not transcript/exon)
            if fields[2] != "gene":
                continue

            attrs = fields[8]

            m_id = gene_id_re.search(attrs)
            m_name = gene_name_re.search(attrs)
            m_type = gene_type_re.search(attrs)

            if not m_id or not m_name:
                continue

            ensembl_id = m_id.group(1).split(".")[0]  # Strip version
            gene_name = m_name.group(1)
            gene_type = m_type.group(1) if m_type else "unknown"

            type_counts[gene_type] += 1

            if gene_name in mapping:
                if mapping[gene_name] != ensembl_id:
                    duplicates.append(
                        (gene_name, mapping[gene_name], ensembl_id)
                    )
                continue

            mapping[gene_name] = ensembl_id

    return mapping, type_counts, duplicates


def parse_hgnc(hgnc_path: str, existing_mapping: dict) -> dict:
    """Parse HGNC complete set and resolve aliases/previous symbols to Ensembl IDs.

    For each gene in HGNC:
      1. If it has an ensembl_gene_id, map current symbol + aliases + prev symbols → that ID
      2. If no ensembl_gene_id, try to chain through existing_mapping (GTF-derived)

    Only adds entries NOT already in existing_mapping (GTF takes precedence).

    Returns: (new_entries dict, stats dict)
    """
    new_entries = {}
    stats = {
        "total_genes": 0,
        "with_ensembl": 0,
        "added_current": 0,
        "added_alias": 0,
        "added_prev": 0,
        "skipped_existing": 0,
        "skipped_no_ensembl": 0,
    }

    with open(hgnc_path) as f:
        header = f.readline().rstrip("\n").split("\t")
        col_idx = {name: i for i, name in enumerate(header)}

        sym_col = col_idx["symbol"]
        alias_col = col_idx["alias_symbol"]
        prev_col = col_idx["prev_symbol"]
        ensembl_col = col_idx["ensembl_gene_id"]

        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(sym_col, alias_col, prev_col, ensembl_col):
                continue

            stats["total_genes"] += 1

            current_symbol = fields[sym_col].strip()
            ensembl_id = fields[ensembl_col].strip()

            # Parse pipe-delimited alias and prev fields (may have quotes)
            alias_str = fields[alias_col].strip().strip('"')
            prev_str = fields[prev_col].strip().strip('"')
            aliases = [a.strip() for a in alias_str.split("|") if a.strip()]
            prev_symbols = [p.strip() for p in prev_str.split("|") if p.strip()]

            if not ensembl_id:
                # Try to resolve through GTF mapping
                if current_symbol in existing_mapping:
                    ensembl_id = existing_mapping[current_symbol]
                else:
                    stats["skipped_no_ensembl"] += 1
                    continue

            stats["with_ensembl"] += 1

            # Map current symbol (if not already in GTF)
            if current_symbol not in existing_mapping and current_symbol not in new_entries:
                new_entries[current_symbol] = ensembl_id
                stats["added_current"] += 1

            # Map aliases → same Ensembl ID
            for alias in aliases:
                if alias not in existing_mapping and alias not in new_entries:
                    new_entries[alias] = ensembl_id
                    stats["added_alias"] += 1

            # Map previous symbols → same Ensembl ID
            for prev in prev_symbols:
                if prev not in existing_mapping and prev not in new_entries:
                    new_entries[prev] = ensembl_id
                    stats["added_prev"] += 1

    return new_entries, stats


def load_features(path: str) -> list:
    """Load gene names from features.tsv.gz (one name per line, may have columns)."""
    genes = []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            parts = line.strip().split("\t")
            # features.tsv.gz can be 1-column (just names) or multi-column
            genes.append(parts[0])
    return genes


def test_coverage(mapping: dict, features_path: str) -> dict:
    """Test mapping coverage against a features file."""
    genes = load_features(features_path)
    mapped = [g for g in genes if g in mapping]
    unmapped = [g for g in genes if g not in mapping]

    return {
        "path": features_path,
        "total": len(genes),
        "mapped": len(mapped),
        "pct": len(mapped) / max(len(genes), 1) * 100,
        "unmapped_examples": unmapped[:20],
    }


def main():
    parser = argparse.ArgumentParser(
        description="Build comprehensive gene symbol → Ensembl ID mapping from GTF"
    )
    parser.add_argument("--gtf", required=True, help="Path to genes.gtf")
    parser.add_argument("--output", required=True, help="Output TSV path")
    parser.add_argument(
        "--test-features",
        nargs="*",
        help="Optional features.tsv.gz files to test coverage against",
    )
    parser.add_argument(
        "--existing-mapping",
        help="Existing mapping TSV to merge (GTF takes precedence for conflicts)",
    )
    parser.add_argument(
        "--hgnc",
        help="HGNC complete set TSV for alias/prev_symbol resolution",
    )
    args = parser.parse_args()

    print("=" * 70)
    print("COMPREHENSIVE GENE MAPPING BUILDER")
    print("=" * 70)

    # Parse GTF
    print(f"\nParsing GTF: {args.gtf}")
    mapping, type_counts, duplicates = parse_gtf_genes(args.gtf)

    print(f"  Total unique gene names: {len(mapping):,}")
    print(f"  Gene types:")
    for gtype, count in sorted(type_counts.items(), key=lambda x: -x[1]):
        print(f"    {gtype}: {count:,}")

    if duplicates:
        print(f"\n  WARNING: {len(duplicates)} gene names map to multiple Ensembl IDs:")
        for name, id1, id2 in duplicates[:10]:
            print(f"    {name}: {id1} vs {id2} (kept {id1})")

    # Merge with existing mapping if provided
    if args.existing_mapping:
        print(f"\nMerging with existing mapping: {args.existing_mapping}")
        existing_count = 0
        new_from_existing = 0
        with open(args.existing_mapping) as f:
            header = f.readline()
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    sym, eid = parts[0], parts[1]
                    existing_count += 1
                    if sym not in mapping:
                        mapping[sym] = eid
                        new_from_existing += 1
        print(f"  Existing entries: {existing_count:,}")
        print(f"  New entries from existing (not in GTF): {new_from_existing:,}")
        print(f"  Total after merge: {len(mapping):,}")

    # HGNC alias resolution
    if args.hgnc:
        print(f"\nHGNC alias resolution: {args.hgnc}")
        hgnc_entries, hgnc_stats = parse_hgnc(args.hgnc, mapping)
        print(f"  HGNC genes: {hgnc_stats['total_genes']:,}")
        print(f"  With Ensembl ID: {hgnc_stats['with_ensembl']:,}")
        print(f"  New entries added:")
        print(f"    Current symbols: {hgnc_stats['added_current']:,}")
        print(f"    Alias symbols:   {hgnc_stats['added_alias']:,}")
        print(f"    Prev symbols:    {hgnc_stats['added_prev']:,}")
        print(f"    Total new:       {sum(hgnc_stats[k] for k in ('added_current', 'added_alias', 'added_prev')):,}")
        print(f"  Skipped (no Ensembl ID): {hgnc_stats['skipped_no_ensembl']:,}")

        mapping.update(hgnc_entries)
        print(f"  Total after HGNC merge: {len(mapping):,}")

    # Write output
    print(f"\nWriting mapping: {args.output}")
    with open(args.output, "w") as out:
        out.write("gene_symbol\tensembl_id\n")
        for sym in sorted(mapping.keys()):
            out.write(f"{sym}\t{mapping[sym]}\n")
    print(f"  Wrote {len(mapping):,} entries")

    # Test coverage against feature files
    if args.test_features:
        print(f"\n{'=' * 70}")
        print("COVERAGE TESTS")
        print(f"{'=' * 70}")

        for fpath in args.test_features:
            study = Path(fpath).parent.name
            result = test_coverage(mapping, fpath)
            print(f"\n  {study}: {result['mapped']:,}/{result['total']:,} "
                  f"({result['pct']:.1f}%)")
            if result["unmapped_examples"]:
                print(f"    Still unmapped: {result['unmapped_examples'][:10]}")

    print(f"\n{'=' * 70}")
    print("DONE")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
