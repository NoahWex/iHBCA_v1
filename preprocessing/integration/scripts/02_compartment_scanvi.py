#!/usr/bin/env python3
"""
Per-compartment scANVI fine-tune on top of pre-trained scVI models.

Loads a previously-trained scVI model and fine-tunes it as a scANVI
semi-supervised model using L1 cell-type labels (level1_annotation). Produces
the canonical X_scANVI latent representation, predicted labels, UMAP, and
Leiden clusterings consumed by Track A's L2 annotation pipeline.

Canonical configuration: max_epochs=30 (early-stopping patience=10),
n_samples_per_label=100, batch_size=128, seed=42, labels_key=level1_annotation.

Usage:
    python 02_compartment_scanvi.py \
        --counts-npz <{comp}_counts.npz> \
        --gene-data <gene_data.csv> \
        --metadata <{comp}_metadata.csv> \
        --compartment Immune \
        --n-latent 50 \
        --scvi-model-dir <scvi_model_dir> \
        --output-dir <out_dir>
"""

import argparse
import gc
import os
import resource
import sys
import time

import numpy as np
import pandas as pd

# Numba cache patch (same as scVI script)
try:
    import numba
    _o_njit = numba.njit
    def _njit_nc(*a, **k): k.pop("cache", None); return _o_njit(*a, **k)
    numba.njit = _njit_nc
    _o_vec = numba.vectorize
    def _vec_nc(*a, **k): k.pop("cache", None); return _o_vec(*a, **k)
    numba.vectorize = _vec_nc
except ImportError:
    pass

import scanpy as sc
import scvi


SCANVI_PARAMS = {
    "max_epochs": 30,           # fine-tune is fast vs scVI's 300
    "batch_size": 128,
    "n_samples_per_label": 100, # standard scANVI sampling
    "early_stopping": True,
    "early_stopping_patience": 10,
}

RESOLUTIONS = [0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0]


def log_peak_rss(label):
    peak_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(f"  [RSS] {label}: {peak_kb / 1024 / 1024:.1f} GB peak", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts-npz", required=True)
    parser.add_argument("--gene-data", required=True)
    parser.add_argument("--metadata", required=True,
                        help="cell_metadata_enriched.csv with level1_annotation")
    parser.add_argument("--compartment", required=True,
                        choices=["Immune", "Epithelial", "Stromal"])
    parser.add_argument("--n-latent", type=int, default=50)
    parser.add_argument("--scvi-model-dir", required=True,
                        help="Directory containing pre-trained scVI model")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--labels-key", default="level1_annotation",
                        help="obs column with L1 cell-type labels (default: level1_annotation)")
    parser.add_argument("--unlabeled-category", default="unknown",
                        help="String marking unlabeled cells (default: 'unknown')")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    sys.stdout.reconfigure(line_buffering=True)
    os.makedirs(args.output_dir, exist_ok=True)

    # --- Seed everything ---
    np.random.seed(args.seed)
    try:
        import torch
        torch.manual_seed(args.seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(args.seed)
    except ImportError:
        pass
    try:
        scvi.settings.seed = args.seed
    except AttributeError:
        pass

    print("=" * 70, flush=True)
    print(f"Per-Compartment scANVI: {args.compartment}, n_latent={args.n_latent}",
          flush=True)
    print("=" * 70, flush=True)

    # --- Verify scVI model exists ---
    if not os.path.isdir(args.scvi_model_dir):
        raise FileNotFoundError(
            f"scVI model dir not found: {args.scvi_model_dir}\n"
            "scANVI must run after scVI sweep completes."
        )

    # --- Load the saved AnnData reference from the scVI model ---
    # scvi-tools' SCVI.save(..., save_anndata=True) writes adata.h5ad inside
    # the model dir. Loading from this file guarantees the cells, gene order,
    # and HVG selection match exactly what scVI was trained on. This eliminates
    # the gene set drift risk that would occur if we reconstructed the AnnData
    # from NPZ + filter_genes + HVG (non-deterministic across scanpy versions).
    saved_adata_path = os.path.join(args.scvi_model_dir, "adata.h5ad")
    if not os.path.exists(saved_adata_path):
        raise FileNotFoundError(
            f"scVI model dir missing saved adata: {saved_adata_path}\n"
            "Re-run scVI with save_anndata=True"
        )

    print(f"Loading scVI training AnnData: {saved_adata_path}", flush=True)
    t0 = time.time()
    adata = sc.read_h5ad(saved_adata_path)
    print(f"  Shape: {adata.n_obs:,} x {adata.n_vars:,} in "
          f"{time.time() - t0:.1f}s", flush=True)
    print(f"  obs columns: {list(adata.obs.columns)}", flush=True)

    # --- Join L1 labels from metadata via cell_id (obs_names) ---
    print(f"\nJoining L1 labels from metadata: {args.metadata}", flush=True)
    meta_header = pd.read_csv(args.metadata, nrows=0).columns.tolist()
    if args.labels_key not in meta_header:
        raise ValueError(
            f"Labels column '{args.labels_key}' not in metadata. "
            f"Available level/annotation columns: "
            f"{[c for c in meta_header if 'level' in c.lower() or 'annot' in c.lower()]}"
        )

    # The scVI saved adata uses numeric_id (str) as obs_names; the metadata
    # carries numeric_id as an integer column. Match via string conversion.
    label_cols = ["numeric_id", args.labels_key]
    meta = pd.read_csv(args.metadata, usecols=label_cols, low_memory=False)
    meta["numeric_id_str"] = meta["numeric_id"].astype(str)
    label_map = dict(zip(meta["numeric_id_str"], meta[args.labels_key]))

    adata.obs[args.labels_key] = adata.obs_names.map(label_map)

    n_unlabeled = int(adata.obs[args.labels_key].isna().sum())
    if n_unlabeled > 0:
        print(f"  WARNING: {n_unlabeled:,} cells have NaN labels (filling with "
              f"'{args.unlabeled_category}')", flush=True)
    else:
        print(f"  All cells have non-null labels (no unlabeled cells)",
              flush=True)
    adata.obs[args.labels_key] = adata.obs[args.labels_key].fillna(
        args.unlabeled_category
    )

    print(f"  L1 labels ({args.labels_key}): "
          f"{adata.obs[args.labels_key].nunique()} unique", flush=True)
    print(f"  Top values: "
          f"{dict(adata.obs[args.labels_key].value_counts().head(5))}",
          flush=True)
    log_peak_rss("after AnnData load + label join")

    # --- Load pre-trained scVI model ---
    print(f"\nLoading pre-trained scVI model: {args.scvi_model_dir}", flush=True)
    t0 = time.time()
    scvi_model = scvi.model.SCVI.load(args.scvi_model_dir, adata=adata)
    print(f"  scVI model loaded in {time.time() - t0:.1f}s", flush=True)
    log_peak_rss("after scVI model load")

    # --- Initialize scANVI from scVI model ---
    print("\nInitializing scANVI from scVI model...", flush=True)
    print(f"  labels_key: {args.labels_key}", flush=True)
    print(f"  unlabeled_category: {args.unlabeled_category}", flush=True)
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        labels_key=args.labels_key,
        unlabeled_category=args.unlabeled_category,
    )

    # --- Train scANVI (fine-tune) ---
    print(f"\nTraining scANVI (fine-tune, max_epochs={SCANVI_PARAMS['max_epochs']})...",
          flush=True)
    t0 = time.time()
    scanvi_model.train(
        max_epochs=SCANVI_PARAMS["max_epochs"],
        n_samples_per_label=SCANVI_PARAMS["n_samples_per_label"],
        batch_size=SCANVI_PARAMS["batch_size"],
        early_stopping=SCANVI_PARAMS["early_stopping"],
        early_stopping_patience=SCANVI_PARAMS["early_stopping_patience"],
    )
    print(f"  scANVI trained in {(time.time() - t0) / 60:.1f} min", flush=True)
    log_peak_rss("after scANVI training")

    # --- Get latent representation ---
    latent = scanvi_model.get_latent_representation()
    adata.obsm["X_scANVI"] = latent
    print(f"  Latent: {latent.shape}", flush=True)

    # --- Get predicted labels (scANVI's classifier output) ---
    predictions = scanvi_model.predict()
    adata.obs["scanvi_predicted"] = predictions
    print(f"  Predictions: {len(predictions)} cells", flush=True)

    del scvi_model
    gc.collect()

    # --- Save scANVI model ---
    scanvi_model_dir = os.path.join(args.output_dir, "scanvi_model")
    scanvi_model.save(scanvi_model_dir, overwrite=True)
    print(f"  scANVI model saved: {scanvi_model_dir}", flush=True)
    del scanvi_model
    gc.collect()

    # --- UMAP + Leiden on scANVI embedding ---
    print("\nNeighbors + UMAP on X_scANVI...", flush=True)
    sc.pp.neighbors(adata, use_rep="X_scANVI", n_neighbors=30)
    sc.tl.umap(adata)

    print("Leiden clustering...", flush=True)
    for res in RESOLUTIONS:
        sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res}")
        n_cl = adata.obs[f"leiden_{res}"].nunique()
        print(f"    res {res}: {n_cl} clusters", flush=True)

    # --- Save outputs ---
    # 3-char lowercase prefix matches the canonical compartment short names.
    comp_short = {"Immune": "imm", "Epithelial": "epi", "Stromal": "str"}[
        args.compartment
    ]

    h5ad_out = os.path.join(args.output_dir, f"{comp_short}_scanvi.h5ad")
    obs_clean = adata.obs.copy()
    for col in obs_clean.columns:
        if hasattr(obs_clean[col], "cat"):
            obs_clean[col] = obs_clean[col].astype(str).replace("nan", "")
        elif obs_clean[col].dtype == object or obs_clean[col].isna().any():
            obs_clean[col] = obs_clean[col].fillna("").astype(str)
    out_adata = sc.AnnData(
        X=adata.X,
        obs=obs_clean,
        var=adata.var[["symbol"]].copy(),
        obsm={k: np.array(adata.obsm[k]) for k in adata.obsm.keys()},
    )
    out_adata.write_h5ad(h5ad_out)
    print(f"  H5AD saved: {h5ad_out}", flush=True)
    print(f"    obsm: {list(adata.obsm.keys())}", flush=True)

    # Leiden CSV
    leiden_cols = [f"leiden_{r}" for r in RESOLUTIONS]
    leiden_df = adata.obs[leiden_cols].copy()
    leiden_df.insert(0, "cell_id", adata.obs_names)
    leiden_df.to_csv(
        os.path.join(args.output_dir, f"{comp_short}_leiden.csv"), index=False
    )

    # UMAP CSV
    umap_df = pd.DataFrame({
        "cell_id": adata.obs_names,
        "UMAP_1": adata.obsm["X_umap"][:, 0],
        "UMAP_2": adata.obsm["X_umap"][:, 1],
    })
    umap_df.to_csv(
        os.path.join(args.output_dir, f"{comp_short}_scanvi_umap.csv"), index=False
    )

    # Predictions CSV
    pred_df = pd.DataFrame({
        "cell_id": adata.obs_names,
        "scanvi_predicted": predictions,
        "true_label": adata.obs[args.labels_key].values,
    })
    pred_df.to_csv(
        os.path.join(args.output_dir, f"{comp_short}_scanvi_predictions.csv"),
        index=False,
    )

    print("\n" + "=" * 70, flush=True)
    print("Done.", flush=True)
    print("=" * 70, flush=True)


if __name__ == "__main__":
    main()
