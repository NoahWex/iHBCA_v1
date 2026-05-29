#!/usr/bin/env python3
"""
Prepare per-compartment label CSVs from cell_annotations.csv for scIB bio scoring.

Extracts label_kumar_2023 (Seurat anchor, Kumar_2023 reference) and normalizes
to canonical vocabulary via NORMALIZE_MAP (same as validate_step04/
compare_singler_vs_seurat.py). Outputs one CSV per compartment with columns
[cell_id, consensus_label], ready for --labels-csv in run_step_17_scib_array.sh.

Provenance: label_kumar_2023 from 01_Preprocessing Nov 2025 Seurat anchor transfer.
Use as proxy bio label until iHBCA_V1 L2.0c labels are available. Re-run scIB
scoring with ALLOW_NAN_BIO=0 when final labels exist.
"""

import argparse
import os
import pandas as pd

NORMALIZE_MAP = {
    "basal": "Basal",
    "lumhr": "Luminal_HR", "lum hr": "Luminal_HR", "luminal_hr": "Luminal_HR",
    "lumsec": "Luminal_Secretory", "luminal_secretory": "Luminal_Secretory",
    "fibroblasts": "Fibroblast", "fibroblast": "Fibroblast", "fibro": "Fibroblast",
    "vascular": "Vascular_Endothelium", "vasc": "Vascular_Endothelium",
    "vascular_endothelium": "Vascular_Endothelium",
    "lymphatic": "Lymphatic_Endothelium", "lymph": "Lymphatic_Endothelium",
    "lymphatic_endothelium": "Lymphatic_Endothelium",
    "pericytes": "Pericyte", "pericyte": "Pericyte", "perivasc": "Pericyte",
    "perivascular_cells": "Pericyte",
    "tcells": "T_Cell", "t_cells": "T_Cell", "t cell": "T_Cell",
    "bcells": "B_Cell", "b_cells": "B_Cell", "b cell": "B_Cell",
    "myeloid": "Myeloid",
    "mast": "Mast", "mast_cell": "Mast",
    "adipo": "Adipocyte", "adipocytes": "Adipocyte", "adipocyte": "Adipocyte",
    "rbc": "RBC", "red blood cell": "RBC",
    "skin_epithelial": "Skin_Epithelial", "skin epithelial": "Skin_Epithelial",
    "plasma_cell": "Plasma_Cell", "plasmacell": "Plasma_Cell",
    "macrophage": "Myeloid",
}


def normalize(label):
    if pd.isna(label):
        return "Unknown"
    key = str(label).strip().lower()
    return NORMALIZE_MAP.get(key, str(label).strip())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cell-annotations", required=True,
                   help="Path to cell_annotations.csv")
    p.add_argument("--output-dir", required=True,
                   help="Directory for per-compartment output CSVs")
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    df = pd.read_csv(args.cell_annotations,
                     usecols=["cell_id", "label_kumar_2023", "compartment"])
    df["consensus_label"] = df["label_kumar_2023"].apply(normalize)
    print(f"Loaded {len(df)} cells")
    print(df["compartment"].value_counts().to_string())

    for comp_label in ["Epithelial", "Stromal", "Immune"]:
        sub = df[df["compartment"] == comp_label][["cell_id", "consensus_label"]].copy()
        out = os.path.join(args.output_dir, f"labels_{comp_label.lower()}.csv")
        sub.to_csv(out, index=False)
        vc = sub["consensus_label"].value_counts()
        print(f"\n{comp_label}: {len(sub)} cells -> {out}")
        print(vc.to_string())


if __name__ == "__main__":
    main()
