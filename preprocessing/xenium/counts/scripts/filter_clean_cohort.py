"""Drop Tumor and Ipsilateral cells from xenium_qc.csv.

Output: xenium_qc_clean.csv in the same directory.
"""
import argparse, os
import pandas as pd

EXCLUDE_PREFIXES = ("UCI604_Tumor_xenium_1", "UCI604_Ipsilateral_xenium_1")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="xenium_qc.csv")
    ap.add_argument("--output", required=True, help="xenium_qc_clean.csv")
    args = ap.parse_args()

    qc = pd.read_csv(args.input)
    n_before = len(qc)

    exclude = qc["cell_id"].apply(
        lambda c: any(c.startswith(p) for p in EXCLUDE_PREFIXES)
    )
    qc_clean = qc[~exclude].reset_index(drop=True)
    n_dropped = exclude.sum()

    print(f"Before: {n_before:,}  Dropped: {n_dropped:,}  After: {len(qc_clean):,}")
    for pfx in EXCLUDE_PREFIXES:
        n = qc["cell_id"].str.startswith(pfx).sum()
        print(f"  {pfx}: {n:,} cells")

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    qc_clean.to_csv(args.output, index=False)
    print(f"Written: {args.output}")


if __name__ == "__main__":
    main()
