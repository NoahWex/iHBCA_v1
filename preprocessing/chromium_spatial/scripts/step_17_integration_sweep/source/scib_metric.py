#!/usr/bin/env python3
"""
step_17_integration_sweep / scib_metric.py

Compute ONE scIB metric for ONE integration config. Designed to be driven
by a SLURM array so (config x metric) tasks fan out in parallel. Each
task writes {output_dir}/{metric}.json; scib_aggregate.py reads all JSONs
for a config and applies the scoring formula.

Metric set (from sweep_configs.yaml:metrics):
  scored:   ASW_batch, ASW_label, NMI, ARI, graph_conn
  sidecar:  isolated_labels, PCR, kBET

kBET note: kBET is computed as a sidecar and not included in the scored
  formula. Rationale is documented in sweep_configs.yaml. kBET runs in
  the array alongside the other metrics, writes its JSON, and
  scib_aggregate.py pulls it into the wide CSV as an audit column
  without feeding it into batch_mean.
"""

try:
    import numba
    _o = numba.njit
    def _nc(*a, **k):  # noqa: E306
        k.pop("cache", None)
        return _o(*a, **k)
    numba.njit = _nc
    _ov = numba.vectorize
    def _ncv(*a, **k):  # noqa: E306
        k.pop("cache", None)
        return _ov(*a, **k)
    numba.vectorize = _ncv
except ImportError:
    pass

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import scanpy as sc
import yaml


DEFAULT_CONFIG_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "configs", "sweep_configs.yaml",
)

# Metric name universe. The wrapper/SLURM array indexes into this list; the
# actual task list is (scored ∪ sidecar) from the YAML config.
ALL_METRICS = ["ASW_batch", "ASW_label", "NMI", "ARI", "graph_conn", "PCR",
               "isolated_labels", "kBET"]


def _die(msg, hint=None):
    sys.stderr.write(f"FATAL: {msg}\n")
    if hint:
        sys.stderr.write(f"  hint: {hint}\n")
    sys.exit(2)


def _load_config(path):
    if not os.path.exists(path):
        _die(f"sweep_configs.yaml not found: {path}")
    with open(path) as fh:
        return yaml.safe_load(fh)


def load_and_prep(args):
    """Load H5AD, attach labels, compute neighbor graph + scoring Leiden.

    """
    adata = sc.read_h5ad(args.h5ad)

    # Join labels (for bio-conservation metrics)
    if args.labels_csv and os.path.exists(args.labels_csv):
        labels_df = pd.read_csv(args.labels_csv)
        id_col = "cell_id" if "cell_id" in labels_df.columns else labels_df.columns[0]
        if args.label_key not in labels_df.columns:
            _die(f"label_key {args.label_key!r} not in {args.labels_csv}",
                 hint=f"available: {list(labels_df.columns)}")
        label_map = dict(zip(labels_df[id_col].astype(str), labels_df[args.label_key]))
        adata.obs["label"] = [str(label_map.get(str(n), "Unknown"))
                              for n in adata.obs_names]
    else:
        # Bio-conservation metrics will degrade gracefully to NaN when all
        # labels are "Unknown"; log the fact explicitly so the reviewer
        # understands why NMI/ARI/ASW_label come back null.
        print(f"  NOTE: no --labels-csv provided (or file missing). "
              f"Bio-conservation metrics will be NaN.")
        adata.obs["label"] = "Unknown"

    adata.obs["label"] = adata.obs["label"].fillna("Unknown").astype("category")
    if args.batch_key not in adata.obs.columns:
        _die(f"batch_key {args.batch_key!r} not in adata.obs",
             hint=f"columns: {list(adata.obs.columns)}")
    adata.obs[args.batch_key] = adata.obs[args.batch_key].astype("category")

    # Pick embedding
    if "X_emb" in adata.obsm:
        embed_key = "X_emb"
    elif "X_scVI" in adata.obsm:
        embed_key = "X_scVI"
    elif "X_pca" in adata.obsm:
        embed_key = "X_pca"
    else:
        embed_key = list(adata.obsm.keys())[0]

    if "connectivities" not in adata.obsp:
        sc.pp.neighbors(adata, use_rep=embed_key, n_neighbors=15)

    scoring_res_key = f"leiden_{args.scoring_resolution}"
    if scoring_res_key not in adata.obs.columns:
        sc.tl.leiden(adata, resolution=args.scoring_resolution,
                     key_added=scoring_res_key)
    adata.obs["leiden_scoring"] = adata.obs[scoring_res_key]

    return adata, embed_key


def compute_metric(adata, embed_key, metric, batch_key):
    """Compute one scIB metric. Returns (value, elapsed_seconds).

    """
    from sklearn.metrics import (silhouette_score,
                                 adjusted_rand_score,
                                 normalized_mutual_info_score)

    embedding = adata.obsm[embed_key]
    batches = adata.obs[batch_key].values
    labels = adata.obs["label"].values
    labeled = (labels != "Unknown") & (labels != "") & (labels != "nan")
    clusters = adata.obs["leiden_scoring"].values

    t0 = time.time()
    val = np.nan

    if metric == "ASW_batch":
        if len(np.unique(batches)) > 1:
            val = (1 - silhouette_score(embedding, batches)) / 2

    elif metric == "ASW_label":
        if labeled.sum() > 100:
            val = (silhouette_score(embedding[labeled], labels[labeled]) + 1) / 2

    elif metric == "NMI":
        if labeled.sum() > 100:
            val = normalized_mutual_info_score(
                labels[labeled], clusters[labeled],
                average_method="arithmetic")

    elif metric == "ARI":
        if labeled.sum() > 100:
            val = adjusted_rand_score(labels[labeled], clusters[labeled])

    elif metric == "graph_conn":
        # scib.metrics.graph_connectivity requires neighbors to be built;
        # load_and_prep already did that.
        import scib
        val = scib.metrics.graph_connectivity(adata, label_key="label")

    elif metric == "PCR":
        # Principal component regression: batch variance explained pre- vs.
        # post-integration.
        adata_pre = adata.copy()
        sc.pp.normalize_total(adata_pre, target_sum=1e4)
        sc.pp.log1p(adata_pre)
        sc.pp.pca(adata_pre, n_comps=50)
        adata_post = adata.copy()
        from sklearn.decomposition import PCA as skPCA
        embed_dims = adata_post.obsm[embed_key].shape[1]
        n_comps_post = min(50, embed_dims)
        pca_post = skPCA(n_components=n_comps_post).fit_transform(adata_post.obsm[embed_key])
        adata_post.obsm["X_pca"] = pca_post
        import scib
        val = scib.metrics.pcr_comparison(
            adata_pre, adata_post, covariate=batch_key,
            embed="X_pca", n_comps=n_comps_post)
        del adata_pre, adata_post

    elif metric == "isolated_labels":
        import scib
        val = scib.metrics.isolated_labels_f1(
            adata, label_key="label",
            batch_key=batch_key, embed=embed_key)

    elif metric == "kBET":
        # Pure-Python scib_metrics.kbet_per_label (no rpy2 dependency).
        # See sweep_configs.yaml:metrics for kBET-as-sidecar rationale.
        from scib_metrics import kbet_per_label
        from scib_metrics.nearest_neighbors import pynndescent
        nn = pynndescent(
            np.asarray(adata.obsm[embed_key]).astype("float32"),
            n_neighbors=50,
            random_state=0,
        )
        batches_arr = np.asarray(adata.obs[batch_key])
        labels_arr = np.asarray(adata.obs["label"])
        val = float(kbet_per_label(nn, batches_arr, labels_arr))

    else:
        _die(f"Unknown metric: {metric}",
             hint=f"expected one of {ALL_METRICS}")

    elapsed = time.time() - t0
    out = round(float(val), 6) if not np.isnan(val) else None
    return out, elapsed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True,
                        help="integrated.h5ad for ONE config from the sweep")
    parser.add_argument("--metric", required=True, choices=ALL_METRICS)
    parser.add_argument("--batch-key", default="patient_id")
    parser.add_argument("--label-key", default="consensus_label")
    parser.add_argument("--labels-csv", default=None)
    parser.add_argument("--output-dir", required=True,
                        help="Typically CFG_CANONICAL_SWEEP_SCIB/{config}/")
    parser.add_argument("--scoring-resolution", type=float, default=1.0,
                        help="Leiden resolution used for NMI/ARI scoring. "
                             "Should match sweep_configs.yaml:leiden.scoring_resolution.")
    parser.add_argument("--config-path", default=DEFAULT_CONFIG_PATH,
                        help="Path to sweep_configs.yaml (for metric list audit)")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    print(f"=== scIB metric: {args.metric} ===")
    print(f"  H5AD:   {args.h5ad}")
    print(f"  Output: {args.output_dir}")

    if args.dry_run:
        print("=== DRY RUN MODE ===")
        if not os.path.exists(args.h5ad):
            _die(f"--h5ad missing: {args.h5ad}")
        if args.labels_csv and not os.path.exists(args.labels_csv):
            _die(f"--labels-csv missing: {args.labels_csv}")
        print("VALIDATION PASSED")
        return

    if not os.path.exists(args.h5ad):
        _die(f"--h5ad missing: {args.h5ad}",
             hint="sweep job for this config may not have completed; "
                  "check integration_sweep.py logs for this config")

    os.makedirs(args.output_dir, exist_ok=True)

    adata, embed_key = load_and_prep(args)
    print(f"  {adata.n_obs} cells, embed: {embed_key}")

    value, elapsed = compute_metric(adata, embed_key, args.metric, args.batch_key)

    result = {
        "metric": args.metric,
        "value": value,
        "elapsed_s": round(elapsed, 1),
        "n_cells": int(adata.n_obs),
        "embed_key": embed_key,
    }
    out_path = os.path.join(args.output_dir, f"{args.metric}.json")
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"  {args.metric} = {value} ({elapsed:.1f}s)")
    print(f"  Saved: {out_path}")


if __name__ == "__main__":
    main()
