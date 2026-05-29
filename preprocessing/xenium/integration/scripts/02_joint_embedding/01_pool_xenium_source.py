"""Pool per-position Xenium counts for one transcript source.

Concatenates counts_{source}.mtx.gz across all positions (Pat1/Pat2/UCI220228
from _legacy/03_xenium_counts, UCI604 from 01_qc/outputs/xenium_counts_v2),
aligns to the shared 280-gene panel order, and applies the Phase 1 filter
(01_qc/outputs/xenium_qc.csv pass_qc_whole column).

Input:  --source {nuclear,whole,cytoplasmic}
Output: xenium_{source}_counts.mtx.gz, xenium_cells.tsv, xenium_genes.tsv,
        xenium_obs.csv  in --out-dir
"""
import argparse, gzip, io, os, shutil
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp


SOURCE_FILES = {
    "nuclear": "counts_nuclear.mtx.gz",
    "whole": "counts_whole.mtx.gz",
    "cytoplasmic": "counts_cytoplasmic.mtx.gz",
}


def load_position(bundle_dir, source):
    with gzip.open(os.path.join(bundle_dir, SOURCE_FILES[source]), "rb") as gz:
        m = scipy.io.mmread(io.BytesIO(gz.read()))
    X = sp.csr_matrix(m).astype(np.int32)
    genes = pd.read_csv(os.path.join(bundle_dir, "genes.tsv"),
                        header=None)[0].tolist()
    cells = pd.read_csv(os.path.join(bundle_dir, "cells.tsv"),
                        header=None)[0].tolist()
    obs = pd.read_csv(os.path.join(bundle_dir, "obs.csv"))
    return X, genes, cells, obs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", required=True, choices=list(SOURCE_FILES))
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--legacy-root", required=True,
                    help="_legacy/03_xenium_counts/outputs/")
    ap.add_argument("--v2-root", required=True,
                    help="01_qc/outputs/xenium_counts_v2/")
    ap.add_argument("--v2-patient", default="UCI604")
    ap.add_argument("--phase1-qc", required=True,
                    help="01_qc/outputs/xenium_qc.csv")
    ap.add_argument("--panel-genes", required=True,
                    help="canonical panel gene list (one per line); "
                         "positions missing any panel gene get zero-padded")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--apply-filter", action="store_true",
                    help="keep only pass_qc_whole==True cells")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    ref_genes = [g.strip() for g in open(args.panel_genes) if g.strip()]
    ref_idx = {g: i for i, g in enumerate(ref_genes)}
    print(f"Canonical panel: {len(ref_genes)} genes", flush=True)

    manifest = pd.read_csv(args.manifest, sep="\t")
    print(f"Pooling source={args.source}, {len(manifest)} positions", flush=True)

    blocks = []
    obs_parts = []
    for _, row in manifest.iterrows():
        xid = row.iloc[0]; pid = row.iloc[1]
        root = args.v2_root if pid == args.v2_patient else args.legacy_root
        bundle = os.path.join(root, xid)
        if not os.path.isdir(bundle):
            print(f"  SKIP {xid}: no bundle at {bundle}", flush=True)
            continue
        X, genes, _cells, obs = load_position(bundle, args.source)
        # Project onto canonical panel: build (n_cells × len(ref_genes)) with
        # zero columns for missing genes.
        n = X.shape[0]
        X_aligned = sp.lil_matrix((n, len(ref_genes)), dtype=np.int32)
        for local_j, g in enumerate(genes):
            if g in ref_idx:
                X_aligned[:, ref_idx[g]] = X[:, local_j]
        blocks.append(X_aligned.tocsr())
        obs_parts.append(obs)
        missing = len(ref_genes) - sum(1 for g in ref_genes if g in set(genes))
        print(f"  {xid}: {n} cells, {missing} panel genes zero-padded",
              flush=True)

    X_full = sp.vstack(blocks, format="csr")
    obs_full = pd.concat(obs_parts, ignore_index=True)
    print(f"Pooled: {X_full.shape}  obs={len(obs_full)}", flush=True)
    assert X_full.shape[0] == len(obs_full)

    # Apply Phase 1 filter
    if args.apply_filter:
        qc = pd.read_csv(args.phase1_qc, usecols=["cell_id", "pass_qc_whole"])
        qc = qc.set_index("cell_id")
        keep = qc.loc[obs_full["cell_id"]]["pass_qc_whole"].to_numpy()
        print(f"Phase 1 filter: {keep.sum()}/{len(keep)} "
              f"({100*keep.mean():.1f}%) pass_qc_whole", flush=True)
        X_full = X_full[keep]
        obs_full = obs_full.loc[keep].reset_index(drop=True)

    print(f"Final: {X_full.shape}", flush=True)

    # Write
    tmp = f"/tmp/_pool_{os.getpid()}.mtx"
    scipy.io.mmwrite(tmp, X_full, field="integer")
    out_mtx = os.path.join(args.out_dir, f"xenium_{args.source}_counts.mtx.gz")
    with open(tmp, "rb") as src, gzip.open(out_mtx, "wb") as dst:
        shutil.copyfileobj(src, dst)
    os.remove(tmp)

    with open(os.path.join(args.out_dir, "xenium_genes.tsv"), "w") as f:
        for g in ref_genes:
            f.write(g + "\n")
    with open(os.path.join(args.out_dir, "xenium_cells.tsv"), "w") as f:
        for c in obs_full["cell_id"]:
            f.write(c + "\n")
    obs_full.to_csv(os.path.join(args.out_dir, "xenium_obs.csv"), index=False)
    print(f"Wrote bundle to {args.out_dir}", flush=True)


if __name__ == "__main__":
    main()
