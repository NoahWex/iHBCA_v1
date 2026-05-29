#!/usr/bin/env python3
"""
denormalize_tables_for_promotion.py - augment all promotion-bound tables with
human-readable rank-name columns (size_rank_within_L2, name_short, name_full)
joined from nhoodgroup_renaming/lookup.csv.

Raw NhoodGroup IDs are preserved as canonical milo provenance.

Targets (per inquiry):
  - stageH_findings_table.csv
  - stageH_findings_per_nhoodgroup.csv
  - stageF1_nhoodgroups/<contrast>/nhood_groups.csv
  - stageF1_nhoodgroups/<contrast>/nhood_groups_summary.csv
  - stageI3_gsea/per_nhoodgroup_nes.csv
  - stageI3_gsea/per_nhoodgroup_C6C4_nes.csv
  - (F.3 markers per file are already keyed by L2_NG; left untouched here -
     readers consult lookup.csv directly for those)

Outputs land in <inquiry>/outputs/<dest>_denormalized/ as parallel files; the
canonical raw outputs remain unchanged.

Args:
  --inquiry-dir
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
import numpy as np


JOIN_COLS = ["NhoodGroup", "size_rank_within_L2", "name_short",
              "name_full", "name_full_with_compartment"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--inquiry-dir", required=True, type=Path)
    p.add_argument("--out-suffix", default="_denormalized",
                   help="suffix appended to dest dirs / filenames")
    return p.parse_args()


def normalize_ng(s):
    """Coerce a NhoodGroup column to int-string (handles NaN/float)."""
    s = pd.to_numeric(s, errors="coerce")
    s = s.dropna() if hasattr(s, "dropna") else s
    return s.astype(int).astype(str)


def load_lookup(inq: Path) -> pd.DataFrame:
    p = inq / "outputs/nhoodgroup_renaming/lookup.csv"
    if not p.exists():
        raise SystemExit(f"lookup.csv not found at {p} - run "
                          "make_nhoodgroup_renaming.R for this inquiry first")
    lk = pd.read_csv(p)
    lk["NhoodGroup"] = pd.to_numeric(lk["NhoodGroup"]).astype(int).astype(str)
    # Compose contrast key for joins
    lk["contrast"] = lk["contrast"].astype(str)
    # Optional: name_full_with_compartment for unambiguous reference
    lk["name_full_with_compartment"] = (lk["parent_compartment"].astype(str)
                                            + "::" + lk["name_full"].astype(str))
    return lk


def join_by_ng_contrast(df: pd.DataFrame, lookup: pd.DataFrame,
                          ng_col: str, contrast_col: str = "contrast") -> pd.DataFrame:
    """Left-join lookup onto df via (contrast, NhoodGroup int-string)."""
    out = df.copy()
    out[ng_col] = pd.to_numeric(out[ng_col], errors="coerce")
    valid_mask = out[ng_col].notna()
    out.loc[valid_mask, ng_col] = out.loc[valid_mask, ng_col].astype(int).astype(str)
    out[contrast_col] = out[contrast_col].astype(str)
    join = (lookup[["contrast", "NhoodGroup", "size_rank_within_L2",
                     "name_short", "name_full",
                     "name_full_with_compartment"]]
              .rename(columns={"NhoodGroup": ng_col,
                                "contrast": contrast_col}))
    return out.merge(join, on=[contrast_col, ng_col], how="left")


def write(df: pd.DataFrame, dest: Path):
    dest.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(dest, index=False)
    print(f"  wrote {dest} ({len(df)} rows)")


def main():
    a = parse_args()
    inq = a.inquiry_dir
    lookup = load_lookup(inq)
    print(f"Lookup: {len(lookup)} rows ({lookup.contrast.nunique()} contrasts)")

    out_root = inq / f"outputs{a.out_suffix}"
    out_root.mkdir(exist_ok=True)

    # ---- Stage H findings tables (one row per (L2 x contrast x NG)) ----
    sh = inq / "outputs/stageH_findings_table.csv"
    if sh.exists():
        print(f"\nstageH_findings_table.csv")
        df = pd.read_csv(sh)
        # Find NG col (typically NhoodGroup or NhoodGroup_renamed)
        ng_col = "NhoodGroup" if "NhoodGroup" in df.columns else None
        if ng_col is not None and "contrast" in df.columns:
            df2 = join_by_ng_contrast(df, lookup, ng_col=ng_col)
            write(df2, out_root / "stageH_findings_table.csv")
        else:
            print(f"  skip: no NhoodGroup or contrast col")

    sh_ng = inq / "outputs/stageH_findings_per_nhoodgroup.csv"
    if sh_ng.exists():
        print(f"\nstageH_findings_per_nhoodgroup.csv")
        df = pd.read_csv(sh_ng)
        if "NhoodGroup" in df.columns and "contrast" in df.columns:
            df2 = join_by_ng_contrast(df, lookup, ng_col="NhoodGroup")
            write(df2, out_root / "stageH_findings_per_nhoodgroup.csv")

    # ---- Stage F.1 nhood_groups per contrast ----
    f1_dir = inq / "outputs/stageF1_nhoodgroups"
    if f1_dir.exists():
        for contrast_dir in sorted(f1_dir.iterdir()):
            if not contrast_dir.is_dir():
                continue
            contrast = contrast_dir.name
            for fname in ("nhood_groups.csv", "nhood_groups_summary.csv"):
                fp = contrast_dir / fname
                if not fp.exists():
                    continue
                print(f"\nstageF1_nhoodgroups/{contrast}/{fname}")
                df = pd.read_csv(fp, low_memory=False)
                df["contrast"] = contrast  # F.1 files don't carry contrast col
                df2 = join_by_ng_contrast(df, lookup, ng_col="NhoodGroup")
                write(df2, out_root / "stageF1_nhoodgroups" / contrast / fname)

    # ---- Stage I.3 GSEA tables (per_nhoodgroup, both Hallmark+Reactome and C6/C4) ----
    i3_dir = inq / "outputs/stageI3_gsea"
    if i3_dir.exists():
        # GSEA's NhoodGroup_renamed is "epi__BMYO_basal_11"; need to parse
        # the trailing integer to recover NhoodGroup. Attempt that.
        for fname in ("per_nhoodgroup_nes.csv", "per_nhoodgroup_C6C4_nes.csv",
                       "per_nhoodgroup_top_summary.csv",
                       "per_nhoodgroup_nes_annotated.csv"):
            fp = i3_dir / fname
            if not fp.exists():
                continue
            print(f"\nstageI3_gsea/{fname}")
            df = pd.read_csv(fp, low_memory=False)
            if "NhoodGroup_renamed" in df.columns and "contrast" in df.columns:
                # Parse trailing _<int> from NhoodGroup_renamed
                df["NhoodGroup"] = (df["NhoodGroup_renamed"].astype(str)
                                       .str.extract(r"_(\d+)$")[0])
                df2 = join_by_ng_contrast(df, lookup, ng_col="NhoodGroup")
                write(df2, out_root / "stageI3_gsea" / fname)

    print(f"\nDONE: denormalized outputs in {out_root}")


if __name__ == "__main__":
    main()
