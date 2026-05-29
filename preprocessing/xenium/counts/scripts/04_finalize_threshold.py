"""Apply the Phase 1 noise-floor threshold to the pooled Xenium QC table.

Gate: pass_qc_whole = (nFeature_Whole > 5) AND (tx_per_gene_Whole > 1.25).

The cut splits the difference between the defensible shoulder (>1.5) and the
permissive peak-intact cut (>1.0); cohort pass rate is ~93%. FLEX cells use a
separate per-compartment pass_qc set upstream. This script writes the
pass_qc_whole column in-place on the Xenium QC CSV and emits a cohort audit
table (platform × group × total / pass / pass_rate).

Paths are resolved via pipeline/config/paths.yaml (keys: qc_xenium, qc_flex).
The audit output is placed alongside qc_xenium.
"""
import argparse
import os
import sys

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xqc", required=True, help="cohort xenium_qc.csv (rewritten in place)")
    ap.add_argument("--fqc", required=True, help="cohort flex_qc.csv")
    ap.add_argument("--audit", required=True, help="cohort_qc_audit.csv output path")
    args = ap.parse_args()

    xqc_path = args.xqc
    fqc_path = args.fqc
    audit_path = args.audit

    x = pd.read_csv(xqc_path)
    x["pass_qc_whole"] = (x["nFeature_Whole"] > 5) & (x["tx_per_gene_Whole"] > 1.25)
    x.to_csv(xqc_path, index=False)
    print(f'xenium: {x["pass_qc_whole"].sum()}/{len(x)} '
          f'({100 * x["pass_qc_whole"].mean():.2f}%)')

    f = pd.read_csv(fqc_path)
    rows = []
    for comp, sub in f.groupby("compartment"):
        rows.append({"platform": "flex", "group": comp, "total": len(sub),
                     "pass": int(sub["pass_qc"].sum()),
                     "pass_rate": round(sub["pass_qc"].mean(), 4)})
    rows.append({"platform": "flex", "group": "ALL", "total": len(f),
                 "pass": int(f["pass_qc"].sum()),
                 "pass_rate": round(f["pass_qc"].mean(), 4)})
    for pid, sub in x.groupby("patient_id"):
        rows.append({"platform": "xenium", "group": pid, "total": len(sub),
                     "pass": int(sub["pass_qc_whole"].sum()),
                     "pass_rate": round(sub["pass_qc_whole"].mean(), 4)})
    rows.append({"platform": "xenium", "group": "ALL", "total": len(x),
                 "pass": int(x["pass_qc_whole"].sum()),
                 "pass_rate": round(x["pass_qc_whole"].mean(), 4)})
    audit = pd.DataFrame(rows)
    audit.to_csv(audit_path, index=False)
    print()
    print(audit.to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())
