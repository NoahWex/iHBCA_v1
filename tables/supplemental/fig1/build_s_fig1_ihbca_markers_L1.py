#!/usr/bin/env python3
"""
prep_l1_substrate_inputs.py — build L1-labeled inputs for compute_mean_expression.py

Re-uses the canonical L2 mean-expression producer (publication/figures/data_prep/fig1/
compute_mean_expression.py) without modification. That producer groups by the `label`
column of labels_full.csv and restricts to genes in a markers CSV. To get L1-level
mean expression, we feed it:

  1. l1_labels_full.csv -- same shape as labels_full.csv but `label` column is replaced
     by the L1 cell type (L2 -> L1 crosswalk applied; artifacts dropped).
  2. l1_marker_panel.csv -- markers CSV mimicking canonical_markers.csv schema
     (columns: compartment, label, lineage, is_artifact, marker_rank, gene)
     populated with Austin's top-N L1 markers per L1 type. The producer only
     uses the `gene` column to define the gene panel.

Runs locally (input CSVs are small). Outputs land in coordination/handoff/
kai_20260523_l1_markers/inputs/.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


# L2 -> L1 crosswalk (11 L1 cell types, iHBCA v1)
L2_TO_L1: dict[str, str] = {
    "BMYO-basal": "Basal-myoepithelial",
    "BMYO-myo": "Basal-myoepithelial",
    "LASP-basal": "Luminal adaptive secretory precursor",
    "LASP-KIT": "Luminal adaptive secretory precursor",
    "LASP-major": "Luminal adaptive secretory precursor",
    "LHS-major": "Luminal hormone sensing",
    "Lactocyte-LC1": "Lactocyte",
    "Lactocyte-LC2": "Lactocyte",
    "Fibro-IGF1": "Fibroblast",
    "Fibro-major": "Fibroblast",
    "Fibro-prematrix": "Fibroblast",
    "Fibro-SFRP4": "Fibroblast",
    "Pericyte": "Perivascular",
    "Pericyte_active": "Perivascular",
    "VSMC": "Perivascular",
    "VSMC_stress": "Perivascular",
    "Vas-arterial": "Vascular endothelial",
    "Vas-capillary": "Vascular endothelial",
    "Vas-vein": "Vascular endothelial",
    "Lym-major": "Lymphatic endothelial",
    "B_cell": "B-lymphocyte",
    "Plasma": "B-lymphocyte",
    "CD4_Th_like": "T-lymphocyte",
    "CD8_Resting": "T-lymphocyte",
    "CD8_Tem": "T-lymphocyte",
    "CD8_Trm": "T-lymphocyte",
    "IFNg_T": "T-lymphocyte",
    "NK": "T-lymphocyte",
    "NK_ILC_prolif": "T-lymphocyte",
    "Th17": "T-lymphocyte",
    "Treg": "T-lymphocyte",
    "ZNF683_T": "T-lymphocyte",
    "Macro_C1Q": "Myeloid",
    "Macro_FOLR2": "Myeloid",
    "Macro_LAM": "Myeloid",
    "Macro_inflam": "Myeloid",
    "Mono_NC": "Myeloid",
    "cDC1": "Myeloid",
    "cDC2": "Myeloid",
    "pDC": "Myeloid",
    "Mast": "Myeloid",
    "Neutrophil": "Myeloid",
}

L1_TO_COMPARTMENT: dict[str, str] = {
    "Basal-myoepithelial": "epi",
    "Luminal hormone sensing": "epi",
    "Luminal adaptive secretory precursor": "epi",
    "Lactocyte": "epi",
    "Fibroblast": "str",
    "Perivascular": "str",
    "Vascular endothelial": "str",
    "Lymphatic endothelial": "str",
    "Myeloid": "imm",
    "T-lymphocyte": "imm",
    "B-lymphocyte": "imm",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--labels-full", required=True, type=Path,
                   help="publication/analysis/annotation/labels_full.csv")
    p.add_argument("--l1-markers", required=True, type=Path,
                   help="iHBCA_V1 ihbca_l1_markers.csv (top-500/type, 5474 rows)")
    p.add_argument("--gene-data", required=True, type=Path,
                   help="iHBCAv1_upload assembly gene_data.csv (ENSG-indexed; provides canonical gene_symbol for the count matrix)")
    p.add_argument("--top-n", type=int, default=15,
                   help="Top N markers per L1 to include in gene panel (default 15)")
    p.add_argument("--out-labels", required=True, type=Path,
                   help="Output: L1-labeled labels_full")
    p.add_argument("--out-markers", required=True, type=Path,
                   help="Output: L1 marker gene panel CSV")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    # ---- 1. L1-relabel labels_full ----
    print(f"[prep] loading {args.labels_full}")
    labels = pd.read_csv(args.labels_full,
                          usecols=["cell_id", "compartment", "lineage", "label",
                                   "is_artifact", "source_mode", "figure_target"])
    print(f"[prep]   labels rows: {len(labels):,}")

    labels["l1"] = labels["label"].map(L2_TO_L1)
    n_unmapped = int(labels["l1"].isna().sum())
    if n_unmapped:
        unmapped_labels = sorted(labels[labels["l1"].isna()]["label"].dropna().unique())
        print(f"[prep]   dropping {n_unmapped:,} cells with unmapped L2 labels "
              f"({len(unmapped_labels)} unique: {unmapped_labels[:8]}{'...' if len(unmapped_labels) > 8 else ''})")
    labels_l1 = labels.dropna(subset=["l1"]).copy()
    # Sanity: compartment per L1 must agree with L1_TO_COMPARTMENT
    labels_l1["l1_compartment"] = labels_l1["l1"].map(L1_TO_COMPARTMENT)
    mismatch = labels_l1[labels_l1["compartment"] != labels_l1["l1_compartment"]]
    if len(mismatch):
        print(f"[prep]   WARNING: {len(mismatch):,} cells where L1 compartment "
              f"disagrees with labels_full compartment — keeping labels_full as authoritative")

    # Replace `label` with L1 type so compute_mean_expression.py groups by L1
    out_labels = labels_l1.copy()
    out_labels["label"] = out_labels["l1"]
    out_labels = out_labels.drop(columns=["l1", "l1_compartment"])

    args.out_labels.parent.mkdir(parents=True, exist_ok=True)
    out_labels.to_csv(args.out_labels, index=False)
    print(f"[prep]   wrote {args.out_labels} ({len(out_labels):,} rows, "
          f"{out_labels['label'].nunique()} L1 types)")

    # ---- 2. Load gene_data (ENSG -> canonical gene_symbol) ----
    print(f"[prep] loading {args.gene_data}")
    gd = pd.read_csv(args.gene_data, usecols=["gene_id", "gene_symbol"])
    # First column 'gene_id' is the ENSG index used to key the count matrix
    ensg_to_canonical = dict(zip(gd["gene_id"], gd["gene_symbol"]))
    print(f"[prep]   gene_data rows: {len(gd):,}  ENSG->symbol mappings: {len(ensg_to_canonical):,}")

    # ---- 3. L1 marker gene panel (ENSG-remapped to canonical symbols) ----
    print(f"[prep] loading {args.l1_markers}")
    mk = pd.read_csv(args.l1_markers)
    print(f"[prep]   marker rows: {len(mk):,}  unique l1_type: {mk['l1_type'].nunique()}")

    # Remap Austin's gene symbols via ENSG to gene_data's canonical symbols
    mk["canonical_symbol"] = mk["gene_ensembl"].map(ensg_to_canonical)
    n_unmapped = int(mk["canonical_symbol"].isna().sum())
    if n_unmapped:
        ex = mk[mk["canonical_symbol"].isna()][["gene", "gene_ensembl"]].drop_duplicates().head(10)
        print(f"[prep]   WARNING: {n_unmapped} marker rows have ENSG not in gene_data")
        print(f"           examples: {ex.to_dict(orient='records')}")
    mk_resolved = mk.dropna(subset=["canonical_symbol"]).copy()

    n_changed = int((mk_resolved["gene"] != mk_resolved["canonical_symbol"]).sum())
    print(f"[prep]   ENSG remapping changed {n_changed:,} of {len(mk_resolved):,} symbols "
          f"(Austin alias -> canonical HGNC)")

    mk_resolved = mk_resolved.sort_values(["l1_type", "scores"], ascending=[True, False])
    mk_top = mk_resolved.groupby("l1_type", as_index=False).head(args.top_n).copy()
    print(f"[prep]   top {args.top_n}/type selections: {len(mk_top)}  "
          f"unique canonical genes: {mk_top['canonical_symbol'].nunique()}")

    # Build markers CSV mimicking canonical_markers.csv schema. Only the `gene`
    # column is read by the producer; the others exist for downstream provenance.
    panel = pd.DataFrame({
        "compartment": mk_top["l1_type"].map(L1_TO_COMPARTMENT),
        "label": mk_top["l1_type"],
        "lineage": mk_top["l1_type"],
        "is_artifact": False,
        "marker_rank": mk_top.groupby("l1_type").cumcount() + 1,
        "gene": mk_top["canonical_symbol"],
        "gene_ensembl": mk_top["gene_ensembl"],
        "austin_alias": mk_top["gene"],
    })
    panel = panel.drop_duplicates(subset=["gene"], keep="first").reset_index(drop=True)
    args.out_markers.parent.mkdir(parents=True, exist_ok=True)
    panel.to_csv(args.out_markers, index=False)
    print(f"[prep]   wrote {args.out_markers} ({len(panel)} unique canonical genes)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
