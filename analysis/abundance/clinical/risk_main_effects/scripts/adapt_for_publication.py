"""adapt_for_publication.py — Bundle 1.1 publication-facing adapter.

Maps inquiry outputs to a reviewer-facing <inquiry>/for_publication/ layout
organized by content type (cohort, DA, nhoodgroups, markers, pathway, findings).

Layout:
  for_publication/
    cohort/
      donor_design.csv               (Stage A pooled cohort design; PXR only)
      donor_designs/<cohort>_design.csv  (per-cohort designs)
      per_study_breakdown.csv        (Stage E per-L2 x per-study)
      covariate_collinearity.csv     (Stage B)
      cohort_landscape.md            (Stage B narrative)
    differential_abundance/
      <contrast>/da_results.csv      (Stage D per-nhood)
      <contrast>/run_summary.yaml    (formula, adaptive_lfc, etc.)
    nhoodgroups/
      lookup.csv                     (rosetta: NG_id <-> name + abundance summary)
      <contrast>/nhood_groups.csv    (Stage F.1 carving)
    markers/
      nhoodgroup/<contrast>/<NG>_vs_parent.csv  (F.3 NG-grain; +name_short, +NG_da_median_lfc)
      L2/<contrast>/<L2>_vs_parity.csv          (F.3 L2-grain; +L2_name_short, +L2_da_median_lfc)
    pathway_enrichment/
      nhoodgroup_NES.csv             (Stage I.3 NG-grain; +leadingEdge_symbols, +name_short, +NG_da_median_lfc)
      L2_NES.csv                     (Stage I.3 L2-grain; +leadingEdge_symbols, +L2_name_short, +L2_da_median_lfc)
    findings_summary/
      findings_table.csv             (Stage H per-L2 x per-contrast wide)
      per_L2_summary.csv             (Stage E long form)
      modifier_cascade.csv           (Stage I.1 - parity only)
      gap_audit.csv                  (Stage I.2)
    figure_data/
      panelC_features/<L2>/          (cell-level for dotplot reconstruction)
    scripts/
      <framework R/Python>           (parameterized pipeline)
      lib/                           (load_paths, da_helpers, etc.)
      inquiry_config/                (paths.yaml, inquiry.yaml)
    unresolved_genes.csv             (transparency sidecar)

Usage:
  python3 adapt_for_publication.py --inquiry-dir <inquiry_dir>
"""
from pathlib import Path
import argparse
import os
import shutil
import sys

import pandas as pd


def load_gene_map(gene_data_path: Path) -> dict:
    gm = pd.read_csv(gene_data_path, usecols=["gene_id", "symbol"])
    return dict(zip(gm["gene_id"], gm["symbol"]))


UNRESOLVED: set = set()


def map_le(le, ensg_to_symbol):
    if pd.isna(le) or le == "":
        return ""
    out = []
    for e in str(le).split(";"):
        sym = ensg_to_symbol.get(e)
        if sym is None or pd.isna(sym):
            UNRESOLVED.add(e)
            out.append(e)
        else:
            out.append(sym)
    return ";".join(out)


def add_name_short_ng(df, ng_col, lookup, contrast=None):
    sample = df[ng_col].dropna().astype(str).iloc[0] if len(df) else ""
    lk = lookup.copy()
    if "::" not in sample and len(lk):
        lk["NhoodGroup_renamed"] = (lk["NhoodGroup_renamed"].astype(str)
                                    .str.replace("::", "__", regex=False)
                                    .str.replace("-", "_", regex=False))
    if contrast is not None and "contrast" in lk.columns:
        lk = lk[lk["contrast"] == contrast]
    lk = (lk[["NhoodGroup_renamed", "name_short", "group_med_lfc"]]
          .drop_duplicates("NhoodGroup_renamed")
          .rename(columns={"group_med_lfc": "NG_da_median_lfc"}))
    return df.merge(lk, left_on=ng_col, right_on="NhoodGroup_renamed", how="left",
                    suffixes=("", "_lk")).drop(
        columns=[c for c in ["NhoodGroup_renamed_lk"] if c in df.columns])


def add_name_short_l2(df, l2_col, lookup, stage_e=None, contrast=None):
    pl = (lookup[["parent_L2_joint", "parent_label"]]
          .drop_duplicates("parent_L2_joint")
          .rename(columns={"parent_label": "L2_name_short"}))
    df = df.merge(pl, left_on=l2_col, right_on="parent_L2_joint", how="left",
                  suffixes=("", "_lk")).drop(
        columns=[c for c in ["parent_L2_joint_lk"] if c in df.columns])
    if stage_e is not None and contrast is not None:
        med_col = f"{contrast}__med_lfc"
        if med_col in stage_e.columns:
            se = (stage_e[["L2_joint", med_col]]
                  .rename(columns={"L2_joint": "_se_l2", med_col: "L2_da_median_lfc"}))
            df = df.merge(se, left_on=l2_col, right_on="_se_l2", how="left").drop(
                columns=["_se_l2"])
    return df


def write_csv(df, out_path) -> int:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    return len(df)


def copy_passthrough(src, dst):
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(src, dst)


def copy_tree_safe(src_dir: Path, dst_dir: Path):
    """Per-file copy walk. Replaces shutil.copytree which trips on CRSP I/O."""
    if dst_dir.exists():
        shutil.rmtree(dst_dir)
    dst_dir.mkdir(parents=True, exist_ok=True)
    for root, _, files in os.walk(src_dir):
        rel = Path(root).relative_to(src_dir)
        target = dst_dir / rel
        target.mkdir(parents=True, exist_ok=True)
        for f in files:
            shutil.copy(Path(root) / f, target / f)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inquiry-dir", required=True, type=Path)
    parser.add_argument("--gene-data", default=None, type=Path,
                        help="Path to gene_data.csv. If omitted, resolved from inquiry_config/paths.yaml inputs.gene_mapping at runtime.")
    parser.add_argument("--scripts-root", default=None, type=Path,
                        help="Path to da_pipeline scripts/. If omitted, resolved from inquiry_config/paths.yaml framework.scripts_root at runtime.")
    args = parser.parse_args()

    inq_dir = args.inquiry_dir
    if args.gene_data is None or args.scripts_root is None:
        # Resolve from inquiry-local paths.yaml (gate-exempt config artifact)
        import yaml
        with open(inq_dir / "scripts" / "inquiry_config" / "paths.yaml") as f:
            paths_cfg = yaml.safe_load(f)
        if args.gene_data is None:
            args.gene_data = Path(paths_cfg["inputs"]["gene_mapping"])
        if args.scripts_root is None:
            args.scripts_root = Path(paths_cfg["framework"]["scripts_root"])

    out_dir = inq_dir / "for_publication"
    if out_dir.exists():
        shutil.rmtree(out_dir)
    COHORT = out_dir / "cohort"
    DA = out_dir / "differential_abundance"
    NG = out_dir / "nhoodgroups"
    MARKERS = out_dir / "markers"
    PATHWAY = out_dir / "pathway_enrichment"
    FINDINGS = out_dir / "findings_summary"
    FIGDATA = out_dir / "figure_data"
    SCR = out_dir / "scripts"
    for d in [COHORT, DA, NG, MARKERS, PATHWAY, FINDINGS, FIGDATA, SCR]:
        d.mkdir(parents=True)

    print(f"=== publication-facing adapter for {inq_dir.name} ===\n")
    print(f"Loading gene map: {args.gene_data}")
    ensg_to_symbol = load_gene_map(args.gene_data)
    print(f"  {len(ensg_to_symbol):,} ENSG mappings\n")

    lookup_path = inq_dir / "outputs/nhoodgroup_renaming/lookup.csv"
    if not lookup_path.exists():
        print(f"FATAL: lookup missing at {lookup_path}")
        sys.exit(1)
    lookup = pd.read_csv(lookup_path)
    copy_passthrough(lookup_path, NG / "lookup.csv")
    print(f"nhoodgroups/lookup.csv: {len(lookup)} rows ({lookup['viable'].sum()} viable)")

    stage_e = None
    se_path = inq_dir / "outputs/stageE_per_l2_summary.csv"
    if se_path.exists():
        stage_e = pd.read_csv(se_path)

    counts = {}

    # ----- cohort/ -----
    print("\n[cohort/]")
    sa = inq_dir / "outputs/stageA_cohorts"
    if sa.exists():
        pooled = sa / "cohort_pooled_design.csv"
        if pooled.exists():
            copy_passthrough(pooled, COHORT / "donor_design.csv")
            print(f"  donor_design.csv (pooled)")
        n_des = 0
        for f in sorted(sa.glob("cohort_*_design.csv")):
            copy_passthrough(f, COHORT / "donor_designs" / f.name)
            n_des += 1
        if n_des:
            print(f"  donor_designs/: {n_des} per-cohort designs")
        cs = sa / "cohort_summary.md"
        if cs.exists():
            copy_passthrough(cs, COHORT / "cohort_summary.md")
            print(f"  cohort_summary.md")

    psp = inq_dir / "outputs/stageE_per_l2_per_study_long.csv"
    if psp.exists():
        copy_passthrough(psp, COHORT / "per_study_breakdown.csv")
        print(f"  per_study_breakdown.csv")

    sb = inq_dir / "outputs/stageB_exploration"
    if sb.exists():
        cm = sb / "collinearity_matrix.csv"
        if cm.exists():
            copy_passthrough(cm, COHORT / "covariate_collinearity.csv")
            print(f"  covariate_collinearity.csv")
        sl = sb / "STAGEB_covariate_landscape.md"
        if sl.exists():
            copy_passthrough(sl, COHORT / "cohort_landscape.md")
            print(f"  cohort_landscape.md")
        for f in [sb / "cohort_landscape.csv", sb / "per_study_testability.csv"]:
            if f.exists():
                copy_passthrough(f, COHORT / f.name)

    # ----- differential_abundance/ -----
    print("\n[differential_abundance/]")
    da_root = inq_dir / "outputs/stageD_da_results"
    if da_root.exists():
        n = 0
        for cdir in sorted(da_root.iterdir()):
            if not cdir.is_dir() or cdir.name.startswith("."):
                continue
            da_csv = cdir / "da_results.csv"
            if da_csv.exists():
                copy_passthrough(da_csv, DA / cdir.name / "da_results.csv")
                rs = cdir / "run_summary.yaml"
                if rs.exists():
                    copy_passthrough(rs, DA / cdir.name / "run_summary.yaml")
                n += 1
        counts["differential_abundance/"] = f"{n} contrasts"
        print(f"  {n} contrasts")

    # ----- nhoodgroups/ -----
    print("\n[nhoodgroups/]")
    f1_root = inq_dir / "outputs/stageF1_nhoodgroups"
    if f1_root.exists():
        n = 0
        for cdir in sorted(f1_root.iterdir()):
            if not cdir.is_dir() or cdir.name.startswith("."):
                continue
            ng_csv = cdir / "nhood_groups.csv"
            if ng_csv.exists():
                copy_passthrough(ng_csv, NG / cdir.name / "nhood_groups.csv")
                n += 1
        counts["nhoodgroups/<contrast>/"] = f"{n} contrasts"
        print(f"  <contrast>/nhood_groups.csv: {n} contrasts")

    # ----- markers/nhoodgroup/<contrast>/ -----
    print("\n[markers/nhoodgroup/]")
    f3_root = inq_dir / "outputs/stageF3_markers"
    if f3_root.exists():
        n_files = 0
        n_skip = 0
        for cdir in sorted(f3_root.iterdir()):
            if not cdir.is_dir() or cdir.name.startswith("."):
                continue
            contrast = cdir.name
            for f in sorted(cdir.glob("*.csv")):
                if f.stat().st_size < 1024:
                    copy_passthrough(f, MARKERS / "nhoodgroup" / contrast / f.name)
                    n_skip += 1
                    continue
                df = pd.read_csv(f)
                if "group" in df.columns:
                    df = add_name_short_ng(df, "group", lookup, contrast=contrast)
                write_csv(df, MARKERS / "nhoodgroup" / contrast / f.name)
                n_files += 1
        counts["markers/nhoodgroup/"] = f"{n_files} markers + {n_skip} placeholders"
        print(f"  {n_files} markers (+ {n_skip} placeholders)")

    # ----- markers/L2/<contrast>/ -----
    print("\n[markers/L2/]")
    f3l2_root = inq_dir / "outputs/stageF3_L2_markers"
    if f3l2_root.exists():
        n_files = 0
        for cdir in sorted(f3l2_root.iterdir()):
            if not cdir.is_dir() or cdir.name.startswith("."):
                continue
            contrast = cdir.name
            for f in sorted(cdir.glob("*.csv")):
                if f.stat().st_size < 1024:
                    copy_passthrough(f, MARKERS / "L2" / contrast / f.name)
                    continue
                df = pd.read_csv(f)
                l2_col = next((c for c in ["parent_L2", "parent_L2_joint", "L2_joint"]
                               if c in df.columns), None)
                if l2_col:
                    df = add_name_short_l2(df, l2_col, lookup,
                                           stage_e=stage_e, contrast=contrast)
                write_csv(df, MARKERS / "L2" / contrast / f.name)
                n_files += 1
        counts["markers/L2/"] = f"{n_files} L2-marker files"
        print(f"  {n_files} L2-marker files")

    # ----- pathway_enrichment/ -----
    print("\n[pathway_enrichment/]")
    p = inq_dir / "outputs/stageI3_gsea/per_nhoodgroup_nes.csv"
    if p.exists():
        df = pd.read_csv(p)
        if "leadingEdge" in df.columns:
            df["leadingEdge_symbols"] = df["leadingEdge"].apply(
                lambda s: map_le(s, ensg_to_symbol))
        if "NhoodGroup_renamed" in df.columns:
            lk = lookup.copy()
            lk["NhoodGroup_renamed"] = (lk["NhoodGroup_renamed"].astype(str)
                                        .str.replace("::", "__", regex=False)
                                        .str.replace("-", "_", regex=False))
            lk = (lk[["NhoodGroup_renamed", "name_short", "group_med_lfc"]]
                  .drop_duplicates("NhoodGroup_renamed")
                  .rename(columns={"group_med_lfc": "NG_da_median_lfc"}))
            df = df.merge(lk, on="NhoodGroup_renamed", how="left")
        counts["pathway_enrichment/nhoodgroup_NES.csv"] = write_csv(
            df, PATHWAY / "nhoodgroup_NES.csv")
        print(f"  nhoodgroup_NES.csv: {counts['pathway_enrichment/nhoodgroup_NES.csv']:,} rows")

    p = inq_dir / "outputs/stageI3_gsea/per_L2_nes.csv"
    if p.exists():
        df = pd.read_csv(p)
        if "leadingEdge" in df.columns:
            df["leadingEdge_symbols"] = df["leadingEdge"].apply(
                lambda s: map_le(s, ensg_to_symbol))
        l2_col = next((c for c in ["parent_L2_joint", "L2_joint", "parent_L2"]
                       if c in df.columns), None)
        if l2_col:
            if "contrast" in df.columns and stage_e is not None:
                parts = []
                for c, sub in df.groupby("contrast"):
                    parts.append(add_name_short_l2(sub, l2_col, lookup,
                                                   stage_e=stage_e, contrast=c))
                df = pd.concat(parts, ignore_index=True)
            else:
                df = add_name_short_l2(df, l2_col, lookup)
        counts["pathway_enrichment/L2_NES.csv"] = write_csv(df, PATHWAY / "L2_NES.csv")
        print(f"  L2_NES.csv: {counts['pathway_enrichment/L2_NES.csv']:,} rows")

    # ----- findings_summary/ -----
    print("\n[findings_summary/]")
    for fname, target in [
        ("stageH_findings_table.csv", "findings_table.csv"),
        ("stageE_per_l2_summary.csv", "per_L2_summary.csv"),
        ("stageI1_cascade_table.csv", "modifier_cascade.csv"),
        ("stageI2_gap_audit.csv", "gap_audit.csv"),
    ]:
        p = inq_dir / "outputs" / fname
        if p.exists():
            copy_passthrough(p, FINDINGS / target)
            counts[f"findings_summary/{target}"] = "passthrough"
            print(f"  {target}")

    # ----- figure_data/ -----
    print("\n[figure_data/]")
    pc_root = inq_dir / "outputs/panelC_features"
    if pc_root.exists():
        n_dirs = 0
        for sub in sorted(pc_root.iterdir()):
            if not sub.is_dir() or sub.name.startswith("."):
                continue
            copy_tree_safe(sub, FIGDATA / "panelC_features" / sub.name)
            n_dirs += 1
        counts["figure_data/panelC_features/"] = f"{n_dirs} L2 dirs"
        print(f"  panelC_features/: {n_dirs} L2 dirs")

    # ----- scripts/ -----
    print("\n[scripts/]")
    if args.scripts_root.exists():
        n = 0
        for f in sorted(args.scripts_root.iterdir()):
            if f.is_file() and f.suffix in (".R", ".py"):
                copy_passthrough(f, SCR / f.name)
                n += 1
        lib_dir = args.scripts_root / "lib"
        if lib_dir.exists():
            copy_tree_safe(lib_dir, SCR / "lib")
        print(f"  framework: {n} R/Python files + lib/")
    cfg_dir = inq_dir / "config"
    if cfg_dir.exists():
        copy_tree_safe(cfg_dir, SCR / "inquiry_config")
        print(f"  inquiry_config/: copied")

    # ----- unresolved_genes.csv -----
    if UNRESOLVED:
        ur = pd.DataFrame({"unresolved_ensg": sorted(UNRESOLVED)})
        ur.to_csv(out_dir / "unresolved_genes.csv", index=False)
        counts["unresolved_genes.csv"] = len(UNRESOLVED)
    else:
        pd.DataFrame({"unresolved_ensg": []}).to_csv(out_dir / "unresolved_genes.csv", index=False)
        counts["unresolved_genes.csv"] = 0
    print(f"\nunresolved_genes.csv: {counts['unresolved_genes.csv']}")

    print("\n--- Summary ---")
    for k, v in counts.items():
        print(f"  {k}: {v}")
    print(f"\nWritten to: {out_dir}")


if __name__ == "__main__":
    main()
