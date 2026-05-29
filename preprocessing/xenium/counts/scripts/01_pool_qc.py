#!/usr/bin/env python3
"""Phase 1 — Noise-floor QC consolidation.

Xenium: full-cell QC recomputed from pooled per-position obs.
  pass_qc_whole = nFeature_Whole > 5 AND tx_per_gene_Whole > 1
  (loosened from tx_per_gene > 2; that knocked out ~20% mid-quality cells
  whose distribution peaks around 2 tx/gene)

FLEX: per-compartment QC CSVs unioned into one cohort table.
  pass_qc as originally derived in _legacy/01_flex_concord_panel:
    nCount > 5 AND nFeature > 5 AND panel_intersection_tx_per_gene > 2
"""
import argparse, os
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xenium-obs", required=True,
                    help="cohort xenium_obs.csv")
    ap.add_argument("--flex-qc-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print("Xenium QC...", flush=True)
    xq = pd.read_csv(args.xenium_obs)
    xq["pass_qc_whole"] = (xq["nFeature_Whole"] > 5) & (xq["tx_per_gene_Whole"] > 1)
    xq["pass_qc_nuclear"] = xq["qc_pass_nuclear"].astype(bool)
    keep = ["cell_id", "xenium_id", "patient_id",
            "nCount_Whole", "nFeature_Whole", "tx_per_gene_Whole",
            "nCount_Nuclear", "nFeature_Nuclear", "tx_per_gene_Nuclear",
            "nCount_Cytoplasmic", "nFeature_Cytoplasmic",
            "pass_qc_whole", "pass_qc_nuclear"]
    xq = xq[[c for c in keep if c in xq.columns]]
    xq.to_csv(os.path.join(args.out_dir, "xenium_qc.csv"), index=False)
    print(f"  {len(xq)} cells, pass_whole={xq['pass_qc_whole'].sum()} "
          f"({100*xq['pass_qc_whole'].mean():.1f}%)", flush=True)

    print("FLEX QC...", flush=True)
    frames = []
    for comp in ["Immune", "Epithelial", "Stromal"]:
        p = os.path.join(args.flex_qc_dir, f"flex_panel_qc_{comp}.csv")
        df = pd.read_csv(p)
        df["compartment"] = comp
        frames.append(df)
        print(f"  {comp}: {len(df)}, pass={df['pass_qc'].sum()} "
              f"({100*df['pass_qc'].mean():.1f}%)", flush=True)
    fq = pd.concat(frames, ignore_index=True)
    fq.to_csv(os.path.join(args.out_dir, "flex_qc.csv"), index=False)

    # Audit: platform × group
    rows = []
    for comp, sub in fq.groupby("compartment"):
        rows.append({"platform": "flex", "group": comp,
                     "total": len(sub), "pass": int(sub["pass_qc"].sum()),
                     "pass_rate": round(sub["pass_qc"].mean(), 4)})
    for pid, sub in xq.groupby("patient_id"):
        rows.append({"platform": "xenium", "group": pid,
                     "total": len(sub), "pass": int(sub["pass_qc_whole"].sum()),
                     "pass_rate": round(sub["pass_qc_whole"].mean(), 4)})
    rows.append({"platform": "xenium", "group": "ALL",
                 "total": len(xq), "pass": int(xq["pass_qc_whole"].sum()),
                 "pass_rate": round(xq["pass_qc_whole"].mean(), 4)})
    rows.append({"platform": "flex", "group": "ALL",
                 "total": len(fq), "pass": int(fq["pass_qc"].sum()),
                 "pass_rate": round(fq["pass_qc"].mean(), 4)})
    audit = pd.DataFrame(rows)
    audit.to_csv(os.path.join(args.out_dir, "cohort_qc_audit.csv"), index=False)
    print("\nAudit:")
    print(audit.to_string(index=False))


if __name__ == "__main__":
    main()
