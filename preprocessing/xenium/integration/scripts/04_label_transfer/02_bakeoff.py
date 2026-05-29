"""3-way bake-off: SingleR (4a) + FICTURE (4b) + joint-latent kNN (4c).

Reads legacy 2-way arbitration (singler+ficture) and joins kNN vote.
Reports agreement matrix, per-method solo-wins, and decides whether to
promote kNN as canonical or keep the 3-way ensemble.

Decision rule: if kNN agrees with singler AND ficture on >=95% of
high-confidence cells (knn_confidence >= 0.6), flag kNN as promotable.
"""
import argparse, os
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--legacy-arb", required=True,
                    help="_legacy/04_xenium_annotation/outputs/l0p5_arbitrated.csv")
    ap.add_argument("--knn-vote", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--hi-conf-threshold", type=float, default=0.6)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    arb = pd.read_csv(args.legacy_arb)
    knn = pd.read_csv(args.knn_vote)
    df = arb.merge(knn, on="cell_id", how="inner")
    print(f"merged: {len(df)} cells", flush=True)

    df["agree_singler_knn"] = df["singler_label"] == df["knn_label"]
    df["agree_ficture_knn"] = df["ficture_label"] == df["knn_label"]
    df["agree_all3"] = (df["agree_singler_knn"] & df["agree_ficture_knn"]
                       & (df["singler_label"] == df["ficture_label"]))
    df["hi_conf_knn"] = df["knn_confidence"] >= args.hi_conf_threshold

    rows = [
        {"metric": "total", "n": len(df), "frac": 1.0},
        {"metric": "singler==ficture", "n": int(df["agree"].sum()),
         "frac": float(df["agree"].mean())},
        {"metric": "singler==knn", "n": int(df["agree_singler_knn"].sum()),
         "frac": float(df["agree_singler_knn"].mean())},
        {"metric": "ficture==knn", "n": int(df["agree_ficture_knn"].sum()),
         "frac": float(df["agree_ficture_knn"].mean())},
        {"metric": "all3 agree", "n": int(df["agree_all3"].sum()),
         "frac": float(df["agree_all3"].mean())},
        {"metric": "knn hi-conf (>=thr)", "n": int(df["hi_conf_knn"].sum()),
         "frac": float(df["hi_conf_knn"].mean())},
    ]
    hi = df[df["hi_conf_knn"]]
    rows.append({"metric": "hi-conf singler==knn", "n": int(hi["agree_singler_knn"].sum()),
                 "frac": float(hi["agree_singler_knn"].mean())})
    rows.append({"metric": "hi-conf ficture==knn", "n": int(hi["agree_ficture_knn"].sum()),
                 "frac": float(hi["agree_ficture_knn"].mean())})
    rows.append({"metric": "hi-conf all3 agree", "n": int(hi["agree_all3"].sum()),
                 "frac": float(hi["agree_all3"].mean())})
    summary = pd.DataFrame(rows)
    summary.to_csv(os.path.join(args.out_dir, "bakeoff_summary.csv"), index=False)
    print(summary.to_string(index=False), flush=True)

    # agreement per-label (for knn vs singler)
    per_label = (df.groupby("knn_label")
                   .agg(n=("cell_id", "size"),
                        singler_agree=("agree_singler_knn", "mean"),
                        ficture_agree=("agree_ficture_knn", "mean"),
                        mean_conf=("knn_confidence", "mean"))
                   .sort_values("n", ascending=False))
    per_label.to_csv(os.path.join(args.out_dir, "bakeoff_per_label.csv"))
    print("\nPer-knn-label agreement:", flush=True)
    print(per_label.to_string(), flush=True)

    df.to_csv(os.path.join(args.out_dir, "bakeoff_per_cell.csv"), index=False)

    # promotion decision
    hi_singler = float(hi["agree_singler_knn"].mean()) if len(hi) else 0.0
    hi_ficture = float(hi["agree_ficture_knn"].mean()) if len(hi) else 0.0
    promote = hi_singler >= 0.95 and hi_ficture >= 0.95
    print(f"\nPromote kNN canonical? {promote} "
          f"(hi-conf singler_agree={hi_singler:.3f}, "
          f"ficture_agree={hi_ficture:.3f})", flush=True)


if __name__ == "__main__":
    main()
