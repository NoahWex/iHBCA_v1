"""Final 3-way L0.5 arbitration: SingleR + FICTURE + joint-latent kNN.

Rule (majority of 3, SingleR as primary tiebreaker):
  - all 3 agree  -> consensus
  - singler == knn (FICTURE disagrees) -> singler  (expression-similarity majority)
  - singler == ficture (knn disagrees)  -> singler
  - ficture == knn   (singler disagrees) -> ficture
  - all 3 disagree                       -> singler (primary method)

Output schema:
  cell_id, xenium_id, compartment,
  singler_label, ficture_label, knn_label, knn_confidence,
  l0p5_final, method, vote_count, agree_pattern

Notes:
  - Uses legacy l0p5_arbitrated.csv for singler+ficture (SingleR & FICTURE unchanged
    since L0.5 vocabulary was locked at Gate G3).
  - 45K UCI604 re-ingested cells are missing singler/ficture in legacy; they receive
    knn_label directly with method='knn_only'.
"""
import argparse, os
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--legacy-arb", required=True)
    ap.add_argument("--knn-vote", required=True)
    ap.add_argument("--out-csv", required=True)
    args = ap.parse_args()
    os.makedirs(os.path.dirname(args.out_csv), exist_ok=True)

    arb = pd.read_csv(args.legacy_arb)
    knn = pd.read_csv(args.knn_vote)
    df = knn.merge(
        arb[["cell_id", "xenium_id", "singler_label", "ficture_label"]],
        on="cell_id", how="left")
    print(f"total kNN cells: {len(knn)}  legacy-joined: {df['singler_label'].notna().sum()}",
          flush=True)

    singler = df["singler_label"].values
    ficture = df["ficture_label"].values
    knnlab = df["knn_label"].values

    has_sr = ~pd.isna(singler)
    has_fc = ~pd.isna(ficture)

    n = len(df)
    final = np.empty(n, dtype=object)
    method = np.empty(n, dtype=object)
    agree_pat = np.empty(n, dtype=object)
    vote = np.zeros(n, dtype=np.int8)

    has_legacy = has_sr & has_fc
    missing = ~has_legacy
    final[missing] = knnlab[missing]
    method[missing] = "knn_only"
    agree_pat[missing] = "no_legacy"
    vote[missing] = 1

    # For legacy-joined cells, compute in vectorized fashion
    lg = has_legacy
    s = singler[lg]; f = ficture[lg]; k = knnlab[lg]
    sk = s == k; sf = s == f; fk = f == k
    consensus = sk & sf
    s_k_only = sk & ~sf
    s_f_only = sf & ~sk
    f_k_only = fk & ~sk & ~sf
    all_diff = ~(sk | sf | fk)

    idx_lg = np.where(lg)[0]
    def assign(sub_mask, label_src, meth, pat, vcount):
        ids = idx_lg[sub_mask]
        final[ids] = label_src[sub_mask]
        method[ids] = meth
        agree_pat[ids] = pat
        vote[ids] = vcount
    assign(consensus, s, "consensus", "S=F=K", 3)
    assign(s_k_only, s, "singler+knn", "S=K!=F", 2)
    assign(s_f_only, s, "singler+ficture", "S=F!=K", 2)
    assign(f_k_only, f, "ficture+knn", "F=K!=S", 2)
    assign(all_diff, s, "singler_primary", "all_diff", 1)

    df["l0p5_final"] = final
    df["method"] = method
    df["agree_pattern"] = agree_pat
    df["vote_count"] = vote

    print("\nMethod breakdown:", flush=True)
    print(df["method"].value_counts().to_string(), flush=True)
    print("\nAgree pattern:", flush=True)
    print(df["agree_pattern"].value_counts().to_string(), flush=True)
    print("\nFinal L0.5 distribution:", flush=True)
    print(df["l0p5_final"].value_counts().to_string(), flush=True)

    cols = ["cell_id", "xenium_id", "singler_label", "ficture_label",
            "knn_label", "knn_confidence", "l0p5_final", "method",
            "agree_pattern", "vote_count"]
    df[[c for c in cols if c in df.columns]].to_csv(args.out_csv, index=False)
    print(f"\nWrote {args.out_csv}", flush=True)


if __name__ == "__main__":
    main()
