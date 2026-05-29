#!/usr/bin/env python3
"""Rebuild the iHBCA L1 pseudobulk reference filtered to the FLEX v2 gene panel.

Source flip (2026-04-10 late evening): an earlier version of this script read
`iHBCA_Reed_2024/author_share/preintegration_HBCA_inner_ENSEMBL_simple_counts_matrix.npz`,
the 7-study INNER join with only 15,176 genes total. Mapping was 100% effective
but FLEX coverage stuck at 5,970 / 18,082 = 33%. The FLEX coordinator confirmed
the canonical-extracted per-compartment 10x h5 files are the right source.

Current source: per-compartment 10x h5 files extracted Apr 10 from
`all-breast-cells_final.h5ad` at:
  iHBCAv1_upload/.dev/assembly_rebuild_20260405/outputs/components/10x_h5/{epi,str,imm}_counts.h5
Plus the master cell_metadata.csv with `cell_id` + `level1_annotation`.
Provenance: extraction_manifest source_hash 0cbef095 matches HCA-validated
canonical.

This script:
  1. Loads cell_metadata.csv once (cell_id -> level1_annotation)
  2. For each compartment:
     a. scanpy.read_10x_h5 -> AnnData (var_names = gene symbols)
     b. (first compartment) derive FLEX intersection column indices
     c. (subsequent compartments) assert identical var_names
     d. Subset to FLEX panel columns
     e. Join with annotations, drop cells outside L1 vocabulary
     f. Accumulate per-type sum vectors and per-type cell counts
     g. Delete the AnnData to free memory before the next compartment
  3. Divide accumulated sums by counts -> per-type mean expression
  4. Write `l1_pseudobulk_reference_flex.csv` (11 types x ~18K genes)

Pattern source for the FLEX-intersect-then-pseudobulk loop:
  Spatial_HBCA/.../l1_reannotation_20260323/scripts/build_l1_reference.py:117-176
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


# iHBCA L1 vocabulary (11 types, confirmed from level1_annotation)
# Source: build_l1_reference.py:25-37
L1_VOCABULARY = {
    "B-lymphocyte",
    "Basal-myoepithelial",
    "Fibroblast",
    "Lactocyte",
    "Luminal adaptive secretory precursor",
    "Luminal hormone sensing",
    "Lymphatic endothelial",
    "Myeloid",
    "Perivascular",
    "T-lymphocyte",
    "Vascular endothelial",
}

EXCLUDE_TYPES = {"Doublet", "Single nuclei"}


def load_flex_panel(flex_genes_tsv):
    """Load the FLEX 18K gene panel symbols (one per line, no header)."""
    with open(flex_genes_tsv) as fh:
        symbols = [line.strip() for line in fh if line.strip()]
    if len(symbols) < 10000:
        raise ValueError(
            f"FLEX panel at {flex_genes_tsv} has only {len(symbols)} symbols "
            "— expected ~18,000"
        )
    print(f"FLEX panel: {len(symbols)} symbols")
    return symbols


def discover_l1_column(df):
    """Find the column whose unique values best match the L1 vocabulary."""
    best_col = None
    best_overlap = 0
    for col in df.columns:
        if df[col].dtype != object:
            continue
        uniq = set(df[col].dropna().unique())
        overlap = len(uniq & L1_VOCABULARY)
        if overlap > best_overlap:
            best_overlap = overlap
            best_col = col
    if best_col is None or best_overlap < 8:
        raise ValueError(
            f"Could not find L1 column. Best match: {best_col} "
            f"with {best_overlap}/11 vocabulary overlap."
        )
    print(f"L1 column: '{best_col}' (overlap: {best_overlap}/11)")
    return best_col


def derive_flex_intersection(var_names, flex_symbols, gene_mapping_path=None):
    """Map FLEX panel symbols to iHBCA column indices, using HGNC alias
    fallback when a direct symbol match fails.

    Direct match: FLEX symbol present in `var_names` → use that column.
    Alias match: FLEX symbol → ENSEMBL id (via gene_mapping) → any other
        symbol mapping to the same ENSEMBL id that IS in `var_names` → use
        that column.

    Output column labels are FLEX symbol names (not iHBCA alias names) so the
    downstream `intersect(ref_genes, query_genes)` in 09b matches FLEX names
    directly.

    For duplicate symbols in var_names, keep the first occurrence.
    """
    # Build iHBCA symbol -> first column index
    ihbca_sym_to_col = {}
    for idx, sym in enumerate(var_names):
        if sym not in ihbca_sym_to_col:
            ihbca_sym_to_col[sym] = idx

    # Build alias maps if gene_mapping provided
    sym_to_ens = {}
    ens_to_syms = {}
    if gene_mapping_path is not None:
        print(f"  Loading alias mapping: {gene_mapping_path}")
        gmap = pd.read_csv(gene_mapping_path, sep="\t")
        sym_to_ens = dict(zip(gmap["gene_symbol"], gmap["ensembl_id"]))
        for sym, ens in zip(gmap["gene_symbol"], gmap["ensembl_id"]):
            ens_to_syms.setdefault(ens, []).append(sym)
        print(f"  Mapping has {len(gmap)} sym<->ens rows, "
              f"{len(ens_to_syms)} unique ENSEMBL ids")

    flex_to_ihbca_col = {}
    direct = 0
    alias = 0
    for fsym in flex_symbols:
        if fsym in ihbca_sym_to_col:
            flex_to_ihbca_col[fsym] = ihbca_sym_to_col[fsym]
            direct += 1
            continue
        ens = sym_to_ens.get(fsym)
        if ens is None:
            continue
        for cand in ens_to_syms.get(ens, []):
            if cand in ihbca_sym_to_col:
                flex_to_ihbca_col[fsym] = ihbca_sym_to_col[cand]
                alias += 1
                break

    keep_symbols = [s for s in flex_symbols if s in flex_to_ihbca_col]
    keep_cols = np.array(
        [flex_to_ihbca_col[s] for s in keep_symbols], dtype=np.int64
    )
    print(f"  Direct matches: {direct}")
    print(f"  Alias matches:  {alias}")
    return keep_cols, keep_symbols


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5-epi", required=True, help="Epithelial 10x h5")
    parser.add_argument("--h5-str", required=True, help="Stromal 10x h5")
    parser.add_argument("--h5-imm", required=True, help="Immune 10x h5")
    parser.add_argument(
        "--cell-metadata", required=True,
        help="cell_metadata.csv (must have cell_id + level1_annotation columns)"
    )
    parser.add_argument(
        "--flex-panel", required=True,
        help="FLEX panel symbol list (genes.tsv from any compartment intermediate)"
    )
    parser.add_argument(
        "--gene-mapping", required=False, default=None,
        help="HGNC gene_symbol<->ensembl_id TSV for alias fallback "
             "(gene_symbol_to_ensembl_full.tsv). If omitted, only direct symbol "
             "matches are used."
    )
    parser.add_argument("--outdir", required=True)
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Validate inputs and exit before heavy loads"
    )
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ---- Validation preamble ----
    required = {
        "h5-epi": args.h5_epi,
        "h5-str": args.h5_str,
        "h5-imm": args.h5_imm,
        "cell-metadata": args.cell_metadata,
        "flex-panel": args.flex_panel,
    }
    if args.gene_mapping is not None:
        required["gene-mapping"] = args.gene_mapping
    print("=== Input validation ===")
    missing = []
    for name, path in required.items():
        ok = Path(path).exists()
        mark = "OK  " if ok else "MISS"
        print(f"  [{mark}] {name:13s} {path}")
        if not ok:
            missing.append(name)
    if missing:
        sys.exit(f"FATAL: missing inputs: {missing}")

    if args.dry_run:
        print("\nVALIDATION PASSED (dry run) — exiting before heavy loads")
        return

    # ---- Load FLEX panel ----
    flex_symbols = load_flex_panel(args.flex_panel)

    # ---- Load cell_metadata.csv (cell_id -> L1 label) ----
    print("\n=== Loading cell_metadata ===")
    t0 = time.time()
    meta = pd.read_csv(args.cell_metadata, low_memory=False)
    print(f"  Loaded {len(meta)} rows, {len(meta.columns)} columns "
          f"in {time.time()-t0:.1f}s")

    if "cell_id" not in meta.columns:
        raise ValueError(
            f"cell_metadata.csv missing 'cell_id' column. Got: {list(meta.columns)[:10]}"
        )
    l1_col = discover_l1_column(meta)
    cell_id_to_label = dict(zip(meta["cell_id"].values, meta[l1_col].values))
    print(f"  cell_id -> {l1_col} dict: {len(cell_id_to_label)} entries")
    del meta  # 2GB CSV — free before loading h5s

    unique_types = sorted(L1_VOCABULARY)
    type_to_idx = {t: i for i, t in enumerate(unique_types)}
    n_types = len(unique_types)

    # ---- Per-compartment loop ----
    compartments = [
        ("Epithelial", args.h5_epi),
        ("Stromal",    args.h5_str),
        ("Immune",     args.h5_imm),
    ]

    pseudobulk_sums = None  # initialized after first compartment
    type_cell_counts = np.zeros(n_types, dtype=np.int64)
    keep_cols = None
    keep_symbols = None
    canonical_var_names = None
    total_kept = 0
    total_dropped_excluded = 0
    total_dropped_unknown = 0

    for comp_name, h5_path in compartments:
        print(f"\n=== {comp_name}: loading {h5_path} ===")
        t0 = time.time()
        adata = sc.read_10x_h5(h5_path)
        print(f"  Loaded shape={adata.shape} in {time.time()-t0:.1f}s")

        # Make sure var_names are unique (sc default may dedup with -1, -2 etc)
        adata.var_names_make_unique()

        if canonical_var_names is None:
            canonical_var_names = adata.var_names.values.copy()
            keep_cols, keep_symbols = derive_flex_intersection(
                canonical_var_names, flex_symbols, args.gene_mapping
            )
            print(
                f"  Retained {len(keep_symbols)} / {len(flex_symbols)} FLEX panel "
                f"genes ({100*len(keep_symbols)/len(flex_symbols):.1f}% coverage)"
            )
            if len(keep_symbols) < 0.9 * len(flex_symbols):
                print(
                    f"  WARNING: coverage < 90% — check input gene namespace"
                )
            pseudobulk_sums = np.zeros((n_types, len(keep_symbols)), dtype=np.float64)
        else:
            if not np.array_equal(adata.var_names.values, canonical_var_names):
                raise ValueError(
                    f"{comp_name} h5 var_names differ from first compartment. "
                    "All 3 compartments must share the same gene namespace."
                )

        # Subset to FLEX columns BEFORE per-type accumulation
        # adata[:, keep_cols] returns a view; copy to materialize for fast slicing
        adata_flex = adata[:, keep_cols].copy()
        del adata
        print(f"  After FLEX subset: {adata_flex.shape}")

        # ---- Join with L1 labels ----
        cell_ids = adata_flex.obs_names.values
        labels = np.array(
            [cell_id_to_label.get(cid, None) for cid in cell_ids],
            dtype=object,
        )
        unknown_mask = pd.isna(labels) | (labels == None)  # noqa: E711
        excluded_mask = np.array(
            [lb in EXCLUDE_TYPES for lb in labels], dtype=bool
        )
        keep_mask = ~unknown_mask & ~excluded_mask & np.array(
            [lb in L1_VOCABULARY for lb in labels], dtype=bool
        )

        n_unknown = int(unknown_mask.sum())
        n_excluded = int(excluded_mask.sum())
        n_kept = int(keep_mask.sum())
        print(f"  Cells: kept={n_kept}, excluded={n_excluded}, unknown={n_unknown}")
        total_kept += n_kept
        total_dropped_excluded += n_excluded
        total_dropped_unknown += n_unknown

        # Restrict to kept cells
        adata_kept = adata_flex[keep_mask].copy()
        del adata_flex
        kept_labels = labels[keep_mask]

        # ---- Accumulate per-type sums ----
        assert pseudobulk_sums is not None, "pseudobulk_sums must be initialized on first compartment"
        counts_csr = adata_kept.X.tocsr() if not hasattr(adata_kept.X, "format") or adata_kept.X.format != "csr" else adata_kept.X
        for cell_type in unique_types:
            type_mask = kept_labels == cell_type
            n_cells_type = int(type_mask.sum())
            if n_cells_type == 0:
                continue
            type_sums = np.asarray(counts_csr[type_mask].sum(axis=0)).flatten()
            t_idx = type_to_idx[cell_type]
            pseudobulk_sums[t_idx] += type_sums
            type_cell_counts[t_idx] += n_cells_type
            print(f"    {cell_type}: +{n_cells_type} cells")

        del adata_kept, counts_csr

    # ---- Compute means ----
    print(f"\n=== Aggregate ===")
    print(f"  Total kept cells: {total_kept}")
    print(f"  Total dropped (excluded vocab): {total_dropped_excluded}")
    print(f"  Total dropped (unknown cell_id): {total_dropped_unknown}")

    if (type_cell_counts == 0).any():
        empty_types = [t for t, c in zip(unique_types, type_cell_counts) if c == 0]
        raise ValueError(f"Empty L1 types after aggregation: {empty_types}")

    assert pseudobulk_sums is not None, "compartment loop did not run (sums)"
    assert keep_symbols is not None, "compartment loop did not run (keep_symbols)"
    pseudobulk_means = pseudobulk_sums / type_cell_counts[:, None]

    if not np.isfinite(pseudobulk_means).all():
        raise ValueError("Non-finite values in pseudobulk — aborting")

    print(f"  Per-type cell counts:")
    for t, c in zip(unique_types, type_cell_counts):
        print(f"    {t}: {c}")

    # ---- Save ----
    pb_df = pd.DataFrame(pseudobulk_means, index=unique_types, columns=keep_symbols)
    pb_df.index.name = "cell_type"

    out_csv = outdir / "l1_pseudobulk_reference_flex.csv"
    pb_df.to_csv(out_csv)
    print(f"\nSaved: {out_csv} ({pb_df.shape[0]} types x {pb_df.shape[1]} genes)")
    print("Done.")


if __name__ == "__main__":
    main()
