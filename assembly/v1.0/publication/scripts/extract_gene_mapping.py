#!/usr/bin/env python3
"""
extract_gene_mapping.py
=======================
Extract gene symbol <-> Ensembl ID mapping from a Cell Ranger filtered h5 file.

The h5 features group contains both gene symbols (name) and Ensembl IDs (id),
so we don't need to go groveling through GTF files like it's 2015.

Reads /matrix/features/{id, name, feature_type} and writes a two-column TSV
(gene_symbol, ensembl_id) filtered to Gene Expression features only.

Usage:
    python extract_gene_mapping.py \
        --input /path/to/sample_filtered_feature_bc_matrix.h5 \
        --output /path/to/gene_symbol_to_ensembl.tsv
"""

import argparse
from collections import Counter

import h5py


def extract_mapping(h5_path: str, output_path: str) -> None:
    """Extract gene symbol -> Ensembl ID mapping from Cell Ranger h5."""

    with h5py.File(h5_path, "r") as f:
        # Cell Ranger h5 stores features under /matrix/features/
        features = f["matrix"]["features"]

        # Decode byte strings to UTF-8
        ensembl_ids = [x.decode("utf-8") for x in features["id"][:]]
        gene_symbols = [x.decode("utf-8") for x in features["name"][:]]
        feature_types = [x.decode("utf-8") for x in features["feature_type"][:]]

    total_features = len(ensembl_ids)
    print(f"Total features in h5: {total_features}")

    # Filter to Gene Expression only (exclude Antibody Capture, CRISPR Guide, etc.)
    gene_mask = [ft == "Gene Expression" for ft in feature_types]
    ensembl_ids = [eid for eid, keep in zip(ensembl_ids, gene_mask) if keep]
    gene_symbols = [gs for gs, keep in zip(gene_symbols, gene_mask) if keep]
    excluded = total_features - len(ensembl_ids)

    print(f"Gene Expression features: {len(ensembl_ids)}")
    if excluded > 0:
        excluded_types = Counter(
            ft for ft, keep in zip(feature_types, gene_mask) if not keep
        )
        print(f"Excluded {excluded} non-gene features: {dict(excluded_types)}")

    # Check for duplicate gene symbols
    symbol_counts = Counter(gene_symbols)
    duplicates = {s: c for s, c in symbol_counts.items() if c > 1}
    if duplicates:
        print(f"\nWARNING: {len(duplicates)} duplicate gene symbols found:")
        for sym, count in sorted(duplicates.items()):
            indices = [i for i, gs in enumerate(gene_symbols) if gs == sym]
            mapped_ids = [ensembl_ids[i] for i in indices]
            print(f"  {sym} -> {mapped_ids}")
    else:
        print("No duplicate gene symbols found.")

    # Check for symbols mapping to multiple Ensembl IDs
    symbol_to_ids = {}
    for sym, eid in zip(gene_symbols, ensembl_ids):
        symbol_to_ids.setdefault(sym, set()).add(eid)
    multi_mapped = {s: ids for s, ids in symbol_to_ids.items() if len(ids) > 1}
    if multi_mapped:
        print(f"\nWARNING: {len(multi_mapped)} symbols map to multiple Ensembl IDs:")
        for sym, ids in sorted(multi_mapped.items()):
            print(f"  {sym} -> {sorted(ids)}")
    else:
        print("All symbols map to exactly one Ensembl ID.")

    # Check for Ensembl IDs mapping to multiple symbols
    id_to_symbols = {}
    for sym, eid in zip(gene_symbols, ensembl_ids):
        id_to_symbols.setdefault(eid, set()).add(sym)
    multi_symbol = {eid: syms for eid, syms in id_to_symbols.items() if len(syms) > 1}
    if multi_symbol:
        print(
            f"\nWARNING: {len(multi_symbol)} Ensembl IDs map to multiple symbols:"
        )
        for eid, syms in sorted(multi_symbol.items()):
            print(f"  {eid} -> {sorted(syms)}")
    else:
        print("All Ensembl IDs map to exactly one symbol.")

    # Write TSV
    with open(output_path, "w") as out:
        out.write("gene_symbol\tensembl_id\n")
        for sym, eid in zip(gene_symbols, ensembl_ids):
            out.write(f"{sym}\t{eid}\n")

    print(f"\nWrote {len(gene_symbols)} gene mappings to: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Extract gene symbol <-> Ensembl ID mapping from Cell Ranger h5"
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to Cell Ranger filtered_feature_bc_matrix.h5",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path for output TSV (gene_symbol, ensembl_id)",
    )
    args = parser.parse_args()

    extract_mapping(args.input, args.output)


if __name__ == "__main__":
    main()
