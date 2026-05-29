#!/usr/bin/env python3
"""
Step 04c: SC/SN reconciliation + Leiden vote consolidation (aggregate across samples).

Merges all per-sample singler_labels.csv files produced by 03b_classify_singler.R,
applies SC/SN reconciliation rules (exact port from reconcile_labels.R), then
maps to canonical label names and consolidates via Leiden majority vote.

Path conventions:
  --singler-root defaults to ${CFG_PREPROCESSING_STEP_04} (the per-sample dirs live there)
  --leiden-csv   the multi-resolution Leiden CSV produced by the integration sweep
  --output-dir   defaults to ${CFG_PREPROCESSING_STEP_04}/harmonized/

Outputs:
  singler_labels_all.csv      — all samples concatenated with per-cell labels
  harmonized_l1_labels.csv    — after reconciliation + Leiden vote consolidation
    columns: cell_id, SC_label, SN_label, reconciled_label, harmonized_l1,
             vote_purity, leiden_resolution, leiden_cluster
"""

import argparse
import os
import glob

import pandas as pd


# ==============================================================================
# SC/SN Normalization (exact port from reconcile_labels.R internal_std_map)
# ==============================================================================

NORMALIZE_MAP = {
    "basal": "basal",
    "lumhr": "luminalhr", "lum hr": "luminalhr", "luminal_hr": "luminalhr",
    "lumsec": "luminalsecretory", "luminal_secretory": "luminalsecretory",
    "fibroblasts": "fibroblast", "fibroblast": "fibroblast", "fibro": "fibroblast",
    "vascular": "vascularendothelium", "vasc": "vascularendothelium",
    "vascular_endothelium": "vascularendothelium",
    "lymphatic": "lymphaticendothelium", "lymph": "lymphaticendothelium",
    "lymphatic_endothelium": "lymphaticendothelium",
    "pericytes": "pericyte", "pericyte": "pericyte", "perivasc": "pericyte",
    "perivascular_cells": "pericyte",
    "tcells": "tcell", "t_cells": "tcell", "t cell": "tcell",
    "bcells": "bcell", "b_cells": "bcell", "b cell": "bcell",
    "myeloid": "myeloid",
    "mast": "mastcell", "mast_cell": "mastcell",
    "adipo": "adipocyte", "adipocytes": "adipocyte",
    "rbc": "rbc", "red blood cell": "rbc",
    "skin_epithelial": "skinepithelial", "skin epithelial": "skinepithelial",
}

CANONICAL_MAP = {
    "basal": "Basal",
    "luminalhr": "Luminal_HR",
    "luminalsecretory": "Luminal_Secretory",
    "skinepithelial": "Skin_Epithelial",
    "fibroblast": "Fibroblast",
    "vascularendothelium": "Vascular_Endothelium",
    "lymphaticendothelium": "Lymphatic_Endothelium",
    "pericyte": "Pericyte",
    "tcell": "T_Cell",
    "bcell": "B_Cell",
    "plasmacell": "Plasma_Cell",
    "myeloid": "Myeloid",
    "mastcell": "Mast",
    "macrophage": "Myeloid",
    "adipocyte": "Adipocyte",
    "rbc": "RBC",
}


def normalize_label(label):
    if pd.isna(label) or str(label).strip() == "":
        return None
    key = str(label).strip().lower()
    if key in NORMALIZE_MAP:
        return NORMALIZE_MAP[key]
    fallback = "".join(c for c in key if c.isalnum())
    return fallback if fallback else None


def reconcile_sc_sn(norm_sc, norm_sn):
    """SC/SN conflict resolution rules (port of case_when from reconcile_labels.R)."""
    if norm_sc is None and norm_sn is None:
        return "Unassigned_BothNA"
    if norm_sc is None:
        return norm_sn
    if norm_sn is None:
        return norm_sc
    if norm_sc == norm_sn:
        return norm_sc
    # Granularity: SN more specific within myeloid
    if norm_sc == "myeloid" and norm_sn in ("mastcell", "macrophage"):
        return norm_sn
    # Unique from SN
    if norm_sn in ("adipocyte", "rbc"):
        return norm_sn
    # Unique from SC
    if norm_sc in ("bcell", "plasmacell", "skinepithelial"):
        return norm_sc
    # Default: prefer SC
    return norm_sc


def to_canonical(norm_label):
    if norm_label is None or norm_label.startswith("Unassigned"):
        return norm_label or "Unassigned_NA"
    if norm_label in CANONICAL_MAP:
        return CANONICAL_MAP[norm_label]
    return norm_label[0].upper() + norm_label[1:]


def compute_vote_purity(labels, clusters):
    df = pd.DataFrame({"label": labels, "cluster": clusters}).dropna(subset=["cluster"])
    stats = []
    for cl, group in df.groupby("cluster"):
        counts = group["label"].value_counts()
        stats.append({
            "cluster": cl,
            "n_cells": len(group),
            "majority_label": counts.index[0],
            "majority_count": counts.iloc[0],
            "vote_purity": counts.iloc[0] / len(group),
        })
    stats_df = pd.DataFrame(stats)
    return stats_df, stats_df["vote_purity"].mean()


def main():
    parser = argparse.ArgumentParser(description="Step 04c: aggregate + harmonize labels")
    parser.add_argument("--singler-root", default=os.environ.get("CFG_PREPROCESSING_STEP_04", ""),
                        help="Root dir containing per-sample singler_labels.csv files. Defaults to CFG_PREPROCESSING_STEP_04")
    parser.add_argument("--leiden-csv", required=True,
                        help="Multi-resolution Leiden CSV (cell_id, leiden_0.2, ...)")
    parser.add_argument("--output-dir", default="",
                        help="Output dir. Defaults to CFG_PREPROCESSING_STEP_04/harmonized/")
    args = parser.parse_args()

    # Resolve paths
    if not args.singler_root:
        parser.error("--singler-root is required (or set CFG_PREPROCESSING_STEP_04)")
    if not args.output_dir:
        args.output_dir = os.path.join(args.singler_root, "harmonized")
    os.makedirs(args.output_dir, exist_ok=True)

    # Collect all per-sample singler_labels.csv files
    pattern = os.path.join(args.singler_root, "*", "singler_labels.csv")
    sample_files = glob.glob(pattern)
    if not sample_files:
        raise FileNotFoundError(f"No singler_labels.csv found under {args.singler_root}/*/ — "
                                "run 03b_classify_singler.R first")
    print(f"Found {len(sample_files)} sample files")

    all_labels = pd.concat([pd.read_csv(f) for f in sample_files], ignore_index=True)
    print(f"Total cells: {len(all_labels)}")

    # Save merged file
    merged_path = os.path.join(args.output_dir, "singler_labels_all.csv")
    all_labels.to_csv(merged_path, index=False)
    print(f"Saved merged: {merged_path}")

    # === Reconciliation ===
    print("\n=== Phase 1: SC/SN Reconciliation ===")
    all_labels["norm_SC"] = all_labels["SC_label"].apply(normalize_label)
    all_labels["norm_SN"] = all_labels["SN_label"].apply(normalize_label)
    all_labels["reconciled_norm"] = [
        reconcile_sc_sn(sc, sn)
        for sc, sn in zip(all_labels["norm_SC"], all_labels["norm_SN"])
    ]
    all_labels["reconciled_label"] = all_labels["reconciled_norm"].apply(to_canonical)

    n_conflict = all_labels["reconciled_label"].str.startswith("Unassigned").sum()
    print(f"Unassigned: {n_conflict} ({100*n_conflict/len(all_labels):.2f}%)")
    for label, count in all_labels["reconciled_label"].value_counts().head(15).items():
        print(f"  {label}: {count}")

    # === Leiden Vote Consolidation ===
    print("\n=== Phase 2: Leiden Vote Consolidation ===")
    leiden = pd.read_csv(args.leiden_csv, index_col=0)
    print(f"Leiden: {len(leiden)} cells, resolutions: {list(leiden.columns)}")

    singler_idx = all_labels.set_index("cell_id")
    df = leiden.join(singler_idx["reconciled_label"], how="inner")
    print(f"Joined: {len(df)} cells")

    leiden_cols = [c for c in leiden.columns if c.startswith("leiden")]
    if not leiden_cols:
        raise ValueError("No leiden_* columns found in leiden CSV — cannot vote-consolidate")

    # Initialize with first resolution, then scan for better purity
    best_res = leiden_cols[0]
    best_stats, best_purity = compute_vote_purity(df["reconciled_label"].values, df[best_res].values)
    print(f"  {best_res}: {len(best_stats)} clusters, mean purity={best_purity:.3f}")
    for col in leiden_cols[1:]:
        stats, mean_purity = compute_vote_purity(df["reconciled_label"].values, df[col].values)
        print(f"  {col}: {len(stats)} clusters, mean purity={mean_purity:.3f}")
        if mean_purity > best_purity:
            best_purity, best_res, best_stats = mean_purity, col, stats

    print(f"\nBest: {best_res} (mean purity={best_purity:.3f})")

    cl2label  = {str(k): v for k, v in zip(best_stats["cluster"], best_stats["majority_label"])}
    cl2purity = {str(k): v for k, v in zip(best_stats["cluster"], best_stats["vote_purity"])}

    df["harmonized_l1"]    = df[best_res].astype(str).replace(cl2label)
    df["vote_purity"]      = df[best_res].astype(str).replace(cl2purity)
    df["leiden_resolution"] = best_res
    df["leiden_cluster"]   = df[best_res]

    n_unmapped = sum(1 for v in df["harmonized_l1"] if pd.isna(v))
    if n_unmapped > 0:
        print(f"Unmapped → fallback to reconciled: {n_unmapped}")
        mask = df["harmonized_l1"].isna()
        df.loc[mask, "harmonized_l1"] = df.loc[mask, "reconciled_label"]

    print(f"\nHarmonized L1 distribution:")
    for label, count in df["harmonized_l1"].value_counts().items():
        print(f"  {label}: {count}")

    sc_sn = singler_idx[["SC_label", "SN_label", "reconciled_label"]].copy()
    out = df[["harmonized_l1", "vote_purity", "leiden_resolution", "leiden_cluster"]].join(sc_sn)
    out.index.name = "cell_id"

    out_path = os.path.join(args.output_dir, "harmonized_l1_labels.csv")
    out.to_csv(out_path)
    print(f"\nSaved: {out_path}")
    print(f"Types: {df['harmonized_l1'].nunique()}, Cells: {len(df)}")


if __name__ == "__main__":
    main()
