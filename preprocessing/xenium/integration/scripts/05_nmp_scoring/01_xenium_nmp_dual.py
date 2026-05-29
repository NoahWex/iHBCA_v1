"""Phase 5b: Xenium NMP scoring on nuclear AND cytoplasmic counts.

For each Xenium cell with an L0.5 label:
  nmp_nuclear = 100 * sum(neg_marker counts for L0.5) / total_nuclear
  nmp_cyto    = 100 * sum(neg_marker counts for L0.5) / total_cyto
  secondary_signal = best non-self family (CP10K-scaled nuclear) >= per-cell median + delta

"""
import argparse, gzip, io, os
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import yaml


def load_counts(bundle_dir, source):
    mtx = os.path.join(bundle_dir, f"xenium_{source}_counts.mtx.gz")
    with gzip.open(mtx, "rb") as gz:
        m = scipy.io.mmread(io.BytesIO(gz.read()))
    X = sp.csr_matrix(m).astype(np.float32)
    genes = pd.read_csv(os.path.join(bundle_dir, "xenium_genes.tsv"),
                        header=None)[0].tolist()
    cells = pd.read_csv(os.path.join(bundle_dir, "xenium_cells.tsv"),
                        header=None)[0].tolist()
    return X, genes, cells


def nmp_per_type(X, neg_idx, l0p5_arr, totals_safe, types):
    out = np.full(X.shape[0], np.nan, dtype=np.float32)
    for t in types:
        idx = neg_idx.get(t, np.array([], dtype=np.int64))
        if len(idx) == 0:
            continue
        mask = l0p5_arr == t
        if not mask.any():
            continue
        sub = np.asarray(X[:, idx].sum(axis=1)).ravel()
        out[mask] = 100.0 * sub[mask] / totals_safe[mask]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nuclear-bundle", required=True)
    ap.add_argument("--cyto-bundle", required=True)
    ap.add_argument("--l0p5-xenium", required=True,
                    help="l0p5_xenium.csv from Phase 4")
    ap.add_argument("--neg-markers", required=True)
    ap.add_argument("--families", required=True)
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--min-family-delta", type=float, default=0.5)
    args = ap.parse_args()
    os.makedirs(os.path.dirname(args.out_csv), exist_ok=True)

    print("Loading nuclear bundle...", flush=True)
    Xn, genes, cells_n = load_counts(args.nuclear_bundle, "nuclear")
    print(f"  nuclear: {Xn.shape}", flush=True)

    print("Loading cyto bundle...", flush=True)
    Xc, genes_c, cells_c = load_counts(args.cyto_bundle, "cytoplasmic")
    print(f"  cyto: {Xc.shape}", flush=True)
    assert genes == genes_c, "gene order mismatch"
    assert cells_n == cells_c, "cell order mismatch"
    cells = cells_n
    gene_idx = {g: i for i, g in enumerate(genes)}

    tot_n = np.asarray(Xn.sum(axis=1)).ravel()
    tot_c = np.asarray(Xc.sum(axis=1)).ravel()
    tot_n_safe = np.where(tot_n > 0, tot_n, 1.0)
    tot_c_safe = np.where(tot_c > 0, tot_c, 1.0)

    print("Loading L0.5...", flush=True)
    l0 = pd.read_csv(args.l0p5_xenium).set_index("cell_id")
    l0 = l0.reindex(cells)
    labeled = l0["l0p5_final"].notna()
    print(f"  labeled: {labeled.sum()}/{len(cells)}", flush=True)
    l0p5_arr = l0["l0p5_final"].astype(str).values

    print("Loading negative markers...", flush=True)
    with open(args.neg_markers) as fh:
        neg = yaml.safe_load(fh)["types"]
    neg_idx = {}
    for t, payload in neg.items():
        pg = [m["gene"] for m in payload.get("negative_markers", [])
              if m.get("in_panel") and m["gene"] in gene_idx]
        neg_idx[t] = np.array([gene_idx[g] for g in pg], dtype=np.int64)
    types = sorted(neg_idx)

    print("Scoring NMP (nuclear)...", flush=True)
    nmp_nuc = nmp_per_type(Xn, neg_idx, l0p5_arr, tot_n_safe, types)
    print("Scoring NMP (cyto)...", flush=True)
    nmp_cyto = nmp_per_type(Xc, neg_idx, l0p5_arr, tot_c_safe, types)

    print("Loading families...", flush=True)
    with open(args.families) as fh:
        fam_raw = yaml.safe_load(fh)["families"]
    fam_idx = {}; fam_self = {}
    for f, payload in fam_raw.items():
        pg = [m["gene"] for m in payload.get("genes", [])
              if m.get("in_panel") and m["gene"] in gene_idx]
        if not pg:
            continue
        fam_idx[f] = np.array([gene_idx[g] for g in pg], dtype=np.int64)
        fam_self[f] = set(payload.get("type_signature", []) or [])
    fam_names = sorted(fam_idx)

    print("Family scores on nuclear CP10K...", flush=True)
    fam_s = np.zeros((len(cells), len(fam_names)), dtype=np.float32)
    for j, f in enumerate(fam_names):
        idx = fam_idx[f]
        fam_tot = np.asarray(Xn[:, idx].sum(axis=1)).ravel()
        fam_s[:, j] = 1e4 * fam_tot / tot_n_safe / max(len(idx), 1)

    per_cell_median = np.median(fam_s, axis=1)
    self_mask = np.zeros_like(fam_s, dtype=bool)
    for j, f in enumerate(fam_names):
        sig = fam_self.get(f, set())
        if sig:
            self_mask[:, j] = np.isin(l0p5_arr, list(sig))
    contested = fam_s.copy()
    contested[self_mask] = -np.inf
    best_j = contested.argmax(axis=1)
    best_score = contested[np.arange(len(cells)), best_j]
    delta = best_score - per_cell_median
    secondary = np.where(delta >= args.min_family_delta,
                         np.array(fam_names)[best_j], "none")
    secondary[np.isinf(best_score)] = "none"
    secondary[~np.isfinite(per_cell_median)] = "none"

    out = pd.DataFrame({
        "cell_id": cells,
        "l0p5_final": l0p5_arr,
        "nmp_nuclear": nmp_nuc,
        "nmp_cyto": nmp_cyto,
        "nmp_delta": nmp_nuc - nmp_cyto,
        "total_nuclear": tot_n,
        "total_cyto": tot_c,
        "secondary_signal": secondary,
    })
    for j, f in enumerate(fam_names):
        out[f"fam_{f}"] = fam_s[:, j]
    out.to_csv(args.out_csv, index=False)
    print(f"Wrote {args.out_csv} shape={out.shape}", flush=True)

    hi = out[labeled.values]
    print(f"\nnmp_nuclear: median={np.nanmedian(hi['nmp_nuclear']):.2f}, "
          f"mean={np.nanmean(hi['nmp_nuclear']):.2f}", flush=True)
    print(f"nmp_cyto:    median={np.nanmedian(hi['nmp_cyto']):.2f}, "
          f"mean={np.nanmean(hi['nmp_cyto']):.2f}", flush=True)
    print(f"nuclear > cyto (flag EXCLUDE): "
          f"{int((hi['nmp_nuclear'] > hi['nmp_cyto']).sum())}/"
          f"{len(hi)} ({100*(hi['nmp_nuclear'] > hi['nmp_cyto']).mean():.1f}%)",
          flush=True)
    print("\nsecondary_signal top:", flush=True)
    print(out["secondary_signal"].value_counts().head(15).to_string(), flush=True)


if __name__ == "__main__":
    main()
