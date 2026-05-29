"""Merge re-ingested UCI604 v2 count bundles into the cohort xenium_qc.csv.

Phase 1 utility following 02_reingest_uci604.py. Legacy Exp 03 produced
empty UCI604 bundles due to a cell_id dtype mismatch (int32 vs str); this
script replaces the zeroed UCI604 rows with the correctly-typed v2 obs,
applies the Phase 1 pass_qc_whole threshold, and rewrites xenium_qc.csv
in place. Also scans a small grid of candidate thresholds on the fixed
cohort (diagnostic; the locked threshold is applied by 04_finalize_threshold).
"""
import argparse
import os

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--v2-root", required=True,
                    help="UCI604 v2 bundle root (output of 02_reingest_uci604)")
    ap.add_argument("--xqc", required=True,
                    help="cohort xenium_qc.csv (rewritten in place)")
    ap.add_argument("--fqc", required=True, help="cohort flex_qc.csv")
    ap.add_argument("--audit", required=True,
                    help="cohort_qc_audit.csv output path")
    args = ap.parse_args()

    v2_root = args.v2_root
    xqc_path = args.xqc
    fqc_path = args.fqc
    audit_path = args.audit

    df = pd.read_csv(xqc_path)
    print(f"current cohort: {len(df)} cells")

    frames = []
    for xid in sorted(os.listdir(v2_root)):
        obs = pd.read_csv(os.path.join(v2_root, xid, "obs.csv"))
        frames.append(obs)
    uci = pd.concat(frames, ignore_index=True)
    print(f"UCI604 v2: {len(uci)} cells across {len(frames)} positions")

    uci["pass_qc_whole"] = (uci["nFeature_Whole"] > 5) & (uci["tx_per_gene_Whole"] > 1)
    uci["pass_qc_nuclear"] = uci["qc_pass_nuclear"].astype(bool)

    keep = ["cell_id", "xenium_id", "patient_id",
            "nCount_Whole", "nFeature_Whole", "tx_per_gene_Whole",
            "nCount_Nuclear", "nFeature_Nuclear", "tx_per_gene_Nuclear",
            "nCount_Cytoplasmic", "nFeature_Cytoplasmic",
            "pass_qc_whole", "pass_qc_nuclear"]
    uci = uci[[c for c in keep if c in uci.columns]]

    other = df[df["patient_id"] != "UCI604"].copy()
    merged = pd.concat([other[keep], uci], ignore_index=True)
    print(f"merged: {len(merged)} cells "
          f"({(merged['patient_id'] == 'UCI604').sum()} UCI604)")
    merged.to_csv(xqc_path, index=False)

    nf = merged["nFeature_Whole"]
    tpg = merged["tx_per_gene_Whole"]
    print()
    print("--- tx/g quantiles (full cohort with fixed UCI604) ---")
    q = tpg.quantile([0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9]).to_dict()
    print(f"tx/g: q05={q[0.05]:.2f} q10={q[0.1]:.2f} q25={q[0.25]:.2f} "
          f"med={q[0.5]:.2f} q75={q[0.75]:.2f}")

    print()
    print("--- pass rates on full cohort ---")
    for lbl, nfc, tpgc in [
        ("nF>5 & tx/g>2.0", 5, 2.0),
        ("nF>5 & tx/g>1.5", 5, 1.5),
        ("nF>5 & tx/g>1.2", 5, 1.2),
        ("nF>5 & tx/g>1.0", 5, 1.0),
    ]:
        mask = (nf > nfc) & (tpg > tpgc)
        print(f"{lbl:20s}  {100 * mask.mean():5.1f}%  n={mask.sum()}")
        rates = {}
        for pid, sub in merged.groupby("patient_id"):
            m = (sub["nFeature_Whole"] > nfc) & (sub["tx_per_gene_Whole"] > tpgc)
            rates[pid] = f"{100 * m.mean():.1f}%"
        print(f"  per-patient: {rates}")

    rows = []
    fq = pd.read_csv(fqc_path)
    for comp, sub in fq.groupby("compartment"):
        rows.append({"platform": "flex", "group": comp, "total": len(sub),
                     "pass": int(sub["pass_qc"].sum()),
                     "pass_rate": round(sub["pass_qc"].mean(), 4)})
    rows.append({"platform": "flex", "group": "ALL", "total": len(fq),
                 "pass": int(fq["pass_qc"].sum()),
                 "pass_rate": round(fq["pass_qc"].mean(), 4)})
    for pid, sub in merged.groupby("patient_id"):
        rows.append({"platform": "xenium", "group": pid, "total": len(sub),
                     "pass": int(sub["pass_qc_whole"].sum()),
                     "pass_rate": round(sub["pass_qc_whole"].mean(), 4)})
    rows.append({"platform": "xenium", "group": "ALL", "total": len(merged),
                 "pass": int(merged["pass_qc_whole"].sum()),
                 "pass_rate": round(merged["pass_qc_whole"].mean(), 4)})
    audit = pd.DataFrame(rows)
    audit.to_csv(audit_path, index=False)
    print()
    print("Audit written:")
    print(audit.to_string(index=False))


if __name__ == "__main__":
    main()
