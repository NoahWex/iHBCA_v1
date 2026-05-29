#!/usr/bin/env python3
"""
DEPRECATED — scib metric monolith superseded by per-metric sidecar scripts.

This script computed all scIB metrics in a single task. It has been replaced
by the per-metric decomposition listed below so that any slow metric (ASW is
O(N^2) on ~1M cells) carries its own walltime budget and fails in isolation.

Replacement pipeline:
    05e_scib_nmi_ari.py        NMI + ARI over pre-computed Leiden resolutions
    05f_scib_asw_label.py      ASW_label  (scib.metrics.silhouette)
    05g_scib_asw_batch.py      ASW_batch  (scib.metrics.silhouette_batch)
    05h_scib_graph_conn.py     graph_connectivity
    05b_scib_kbet_only.py      kBET
    05d_scib_pcr_recompute.py  PCR (rebuilt on a consistent 4000-HVG baseline)
    06b_build_wide_csv.py      sidecar aggregator -> wide CSV

Methodology note: the monolith used scib's internal opt_louvain for NMI/ARI;
the replacement (05e) sweeps the pre-computed Leiden grid [0.1..5.0] instead.
ASW_label, ASW_batch, and graph_conn are numerically identical between paths.

This script is retained as provenance for the metric-decomposition decision
and is referenced by the active sidecars. Do not re-run.

Original signature (for reference):
    python 05_scib_benchmark.py \
        --h5ad <{comp}_{method}.h5ad> \
        --embedding-key <X_scVI|X_scANVI|...> \
        --batch-key dataset \
        --label-key level1_annotation \
        --method-name <method_label> \
        --compartment Immune \
        --output-csv <out.csv>
"""

import sys
print(
    "WARNING: 05_scib_benchmark.py is DEPRECATED. "
    "Use 05e/05f/05g/05h + 05b + 05d + 06b instead. "
    "See docstring for the replacement pipeline.",
    file=sys.stderr,
)

import argparse
import os
import sys
import time

# Numba cache patch — MUST apply before importing scanpy/scib.
# Numba's @njit(cache=True) decorators try to write cache metadata next to
# the source file. Inside read-only singularity containers, the source files
# at /opt/apps/python/.../scanpy/_utils/compute/is_constant.py are read-only,
# and numba's "no locator available" error fires. Stripping cache=True from
# all njit/vectorize calls disables the cache and bypasses the issue.
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

# Stable per-job TMPDIR for scib LISI multiprocessing.
# scib 1.1.7 LISI graph metric writes worker index files under TMPDIR and
# reads them back from the master. Default TMPDIR can race or be cleaned up
# mid-run, producing FileNotFoundError for graph_lisi_indices_*.txt. Use a
# job-scoped node-local /tmp path (NOT CRSP — CRSP has stale-file-handle
# issues with rapid write-then-read patterns per .claude/rules/hpc-execution.md).
# Python is read at task launch, so this propagates to all queued array
# tasks without resubmission.
import tempfile
import atexit
import shutil
_slurm_job = os.environ.get('SLURM_JOB_ID', 'local')
_slurm_task = os.environ.get('SLURM_ARRAY_TASK_ID', '0')
_stable_tmp = f'/tmp/scib_{_slurm_job}_{_slurm_task}'
os.makedirs(_stable_tmp, exist_ok=True)
os.environ['TMPDIR'] = _stable_tmp
tempfile.tempdir = _stable_tmp
atexit.register(lambda: shutil.rmtree(_stable_tmp, ignore_errors=True))

import anndata as ad
import pandas as pd
import scanpy as sc

# scIB import — try both common module paths
try:
    import scib
    SCIB_AVAILABLE = True
    SCIB_VERSION = getattr(scib, "__version__", "unknown")
except ImportError:
    try:
        import scIB as scib  # noqa
        SCIB_AVAILABLE = True
        SCIB_VERSION = getattr(scib, "__version__", "unknown")
    except ImportError:
        SCIB_AVAILABLE = False
        SCIB_VERSION = "not installed"
        print("WARNING: scib not available — install with `pip install scib`",
              flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True,
                        help="Integrated compartment h5ad")
    parser.add_argument("--embedding-key", required=True,
                        help="obsm key for the integration embedding "
                             "(X_scVI, X_scANVI, X_pca_harmony, X_seurat, ...)")
    parser.add_argument("--batch-key", default="dataset",
                        help="obs column for batch (default: dataset)")
    parser.add_argument("--label-key", default="level1_annotation",
                        help="obs column for cell type labels (default: level1_annotation)")
    parser.add_argument("--method-name", required=True,
                        help="Method label for output (e.g. scVI_n_latent_50)")
    parser.add_argument("--compartment", required=True,
                        choices=["Immune", "Epithelial", "Stromal"])
    parser.add_argument("--output-csv", required=True)
    parser.add_argument("--metadata", default=None,
                        help="Optional cell_metadata_enriched.csv to join missing labels")
    parser.add_argument("--skip-kbet", action="store_true",
                        help="Skip kBET (expensive on large compartments)")
    parser.add_argument("--subsample", type=int, default=None,
                        help="Optional cell subsample for faster scIB on large compartments")
    parser.add_argument("--use-existing-leiden", action="store_true",
                        help="Per-metric mode: compute NMI/ARI using pre-existing "
                             "leiden_* columns in obs instead of scib.metrics.metrics() "
                             "(which recomputes leiden at ~20 resolutions — infeasible "
                             "on >500K cells). ASW and graph_conn computed individually. "
                             "kBET, PCR, isolated_labels skipped (sidecar scripts handle "
                             "kBET and PCR; isolated_labels is sidecar-only per policy).")
    args = parser.parse_args()

    sys.stdout.reconfigure(line_buffering=True)
    os.makedirs(os.path.dirname(args.output_csv), exist_ok=True)

    print("=" * 70, flush=True)
    print(f"scIB benchmark: {args.method_name} ({args.compartment})", flush=True)
    print(f"scib version: {SCIB_VERSION}", flush=True)
    print("=" * 70, flush=True)

    if not SCIB_AVAILABLE:
        raise RuntimeError("scib package is required")

    # --- Load h5ad ---
    print(f"Loading {args.h5ad}...", flush=True)
    t0 = time.time()
    adata = ad.read_h5ad(args.h5ad)
    print(f"  Shape: {adata.n_obs:,} x {adata.n_vars:,} in {time.time() - t0:.1f}s",
          flush=True)
    print(f"  obsm keys: {list(adata.obsm.keys())}", flush=True)
    print(f"  obs columns: {list(adata.obs.columns)[:15]}...", flush=True)

    # --- Verify embedding ---
    if args.embedding_key not in adata.obsm:
        raise ValueError(
            f"Embedding key '{args.embedding_key}' not in obsm. "
            f"Available: {list(adata.obsm.keys())}"
        )

    # --- Verify batch + label keys ---
    if args.batch_key not in adata.obs.columns:
        raise ValueError(
            f"Batch key '{args.batch_key}' not in obs. "
            f"Available: {list(adata.obs.columns)}"
        )

    if args.label_key not in adata.obs.columns:
        if args.metadata is not None:
            print(f"  Label key '{args.label_key}' missing from obs — "
                  f"joining from {args.metadata}", flush=True)
            # Read both possible join keys (cell_id and numeric_id) plus the
            # label. The scVI saved adata uses numeric_id as obs.index, while
            # the canonical h5ad uses cell_id — auto-detect which matches.
            ext_meta = pd.read_csv(
                args.metadata,
                usecols=["cell_id", "numeric_id", args.label_key],
                low_memory=False,
            )
            # numeric_id is integer in CSV; obs_names are strings; coerce both
            ext_meta["numeric_id"] = ext_meta["numeric_id"].astype(str)
            ext_meta["cell_id"] = ext_meta["cell_id"].astype(str)
            obs_names_str = adata.obs_names.astype(str)

            # Try numeric_id first (scVI saved adata convention), then cell_id
            for join_col in ("numeric_id", "cell_id"):
                lookup = ext_meta.set_index(join_col)[args.label_key]
                joined = obs_names_str.map(lookup)
                n_matched = int(joined.notna().sum())
                print(f"  Trying join on '{join_col}': "
                      f"{n_matched:,} / {len(obs_names_str):,} matched",
                      flush=True)
                if n_matched > 0.95 * len(obs_names_str):
                    adata.obs[args.label_key] = joined.values
                    print(f"  Using join column: {join_col}", flush=True)
                    break
            else:
                raise ValueError(
                    f"Label fallback failed: neither 'cell_id' nor "
                    f"'numeric_id' matched adata.obs_names "
                    f"(first 3 obs_names: {list(obs_names_str[:3])})"
                )
            n_null = int(adata.obs[args.label_key].isna().sum())
            if n_null > 0:
                print(f"  WARNING: {n_null:,} cells still have null labels — "
                      f"dropping for scIB metrics", flush=True)
                adata = adata[~adata.obs[args.label_key].isna()].copy()
        else:
            raise ValueError(
                f"Label key '{args.label_key}' not in obs and no --metadata fallback"
            )

    # Convert to category for scib
    adata.obs[args.batch_key] = adata.obs[args.batch_key].astype("category")
    adata.obs[args.label_key] = adata.obs[args.label_key].astype(str).astype("category")

    print(f"\n  Batch: {adata.obs[args.batch_key].nunique()} unique", flush=True)
    print(f"  Labels: {adata.obs[args.label_key].nunique()} unique", flush=True)
    print(f"  Embedding shape: {adata.obsm[args.embedding_key].shape}", flush=True)

    # --- Optional subsample for large compartments ---
    if args.subsample is not None and adata.n_obs > args.subsample:
        import numpy as np
        rng = np.random.RandomState(42)
        idx = rng.choice(adata.n_obs, args.subsample, replace=False)
        idx.sort()
        adata = adata[idx].copy()
        print(f"  Subsampled to {adata.n_obs:,} cells (seed=42)", flush=True)

    # --- scIB needs neighbors computed on the embedding ---
    print(f"\nComputing neighbors on {args.embedding_key}...", flush=True)
    sc.pp.neighbors(adata, use_rep=args.embedding_key, n_neighbors=15)

    # --- Compute scIB metrics ---
    t0 = time.time()

    if args.use_existing_leiden:
        # ================================================================
        # Per-metric path: bypass scib.metrics.metrics() and compute each
        # metric individually using pre-existing leiden columns in obs.
        #
        # scib.metrics.metrics() internally calls sc.tl.leiden() at ~20
        # resolutions for NMI/ARI optimization. On >500K cells this exceeds
        # 24h wall time (observed: 16h + 24h timeouts on 991K epi, 974K str).
        # The per-metric path reuses the 10 leiden resolutions stored in the
        # h5ad during integration (01_compartment_scvi_full.py), avoiding
        # all leiden recomputation.
        #
        # kBET, PCR, isolated_labels are skipped here — they are computed
        # by separate sidecar scripts (05b, 05d) and the aggregator
        # (06b_build_wide_csv.py) always prefers sidecar values.
        # ================================================================
        import re

        print("\n--- Per-metric mode (--use-existing-leiden) ---", flush=True)

        # Find leiden columns
        leiden_cols = {}
        for col in adata.obs.columns:
            m = re.match(r"^leiden_([\d.]+)$", col)
            if m:
                leiden_cols[float(m.group(1))] = col
        if not leiden_cols:
            print("ERROR: No leiden_* columns found in obs.", flush=True)
            print(f"  Available: {list(adata.obs.columns)}", flush=True)
            sys.exit(1)

        print(f"  Found {len(leiden_cols)} leiden resolutions: "
              f"{sorted(leiden_cols.keys())}", flush=True)

        # NMI/ARI: scan pre-existing leiden, pick best per metric.
        # Uses scib.metrics.nmi/ari — same functions called internally by
        # scib.metrics.metrics() — ensuring identical NMI variant
        # (average_method='arithmetic') and ARI computation.
        print(f"\n  NMI / ARI (scanning {len(leiden_cols)} resolutions):",
              flush=True)
        best_nmi, best_ari = -1.0, -1.0
        best_nmi_res, best_ari_res = None, None
        for res in sorted(leiden_cols.keys()):
            col = leiden_cols[res]
            adata.obs[col] = adata.obs[col].astype(str).astype("category")
            nmi_val = scib.metrics.nmi(adata, col, args.label_key)
            ari_val = scib.metrics.ari(adata, col, args.label_key)
            n_cl = adata.obs[col].nunique()
            print(f"    res={res:5.1f}  n_clusters={n_cl:4d}  "
                  f"NMI={nmi_val:.6f}  ARI={ari_val:.6f}", flush=True)
            if nmi_val > best_nmi:
                best_nmi, best_nmi_res = nmi_val, res
            if ari_val > best_ari:
                best_ari, best_ari_res = ari_val, res

        print(f"  Best NMI: {best_nmi:.6f} (res={best_nmi_res})", flush=True)
        print(f"  Best ARI: {best_ari:.6f} (res={best_ari_res})", flush=True)

        # ASW_label
        print(f"\n  ASW_label (scib.metrics.silhouette)...", flush=True)
        t_asw = time.time()
        asw_label = scib.metrics.silhouette(
            adata, label_key=args.label_key, embed=args.embedding_key)
        print(f"    ASW_label = {asw_label:.6f} ({(time.time()-t_asw)/60:.1f} min)",
              flush=True)

        # ASW_batch
        print(f"  ASW_batch (scib.metrics.silhouette_batch)...", flush=True)
        t_asw = time.time()
        asw_batch = scib.metrics.silhouette_batch(
            adata, batch_key=args.batch_key, label_key=args.label_key,
            embed=args.embedding_key)
        print(f"    ASW_batch = {asw_batch:.6f} ({(time.time()-t_asw)/60:.1f} min)",
              flush=True)

        # graph_conn
        print(f"  graph_conn (scib.metrics.graph_connectivity)...", flush=True)
        t_gc = time.time()
        graph_conn = scib.metrics.graph_connectivity(
            adata, label_key=args.label_key)
        print(f"    graph_conn = {graph_conn:.6f} ({(time.time()-t_gc)/60:.1f} min)",
              flush=True)

        elapsed_min = (time.time() - t0) / 60
        print(f"\n  Per-metric total: {elapsed_min:.1f} min", flush=True)

        # Assemble output row with the same column names as the all-in-one path.
        _nan = float("nan")
        out_row = {
            "NMI_cluster/label": best_nmi,
            "ARI_cluster/label": best_ari,
            "ASW_label": asw_label,
            "ASW_label/batch": asw_batch,
            "PCR_batch": _nan,
            "cell_cycle_conservation": _nan,
            "isolated_label_F1": _nan,
            "isolated_label_silhouette": _nan,
            "graph_conn": graph_conn,
            "kBET": _nan,
            "iLISI": _nan,
            "cLISI": _nan,
            "hvg_overlap": _nan,
            "trajectory": _nan,
            # Per-metric provenance
            "leiden_source": "precomputed_obs",
            "leiden_n_resolutions": len(leiden_cols),
            "nmi_best_resolution": best_nmi_res,
            "ari_best_resolution": best_ari_res,
        }

    else:
        # ================================================================
        # Original all-in-one path: scib.metrics.metrics()
        # ================================================================
        print("\nComputing scIB metrics...", flush=True)
        print(f"  kBET: {'SKIP' if args.skip_kbet else 'compute'}", flush=True)

        # LISI DISABLED (pre-authorized fallback F.5 in plan serialized-mapping-honey)
        # --------------------------------------------------------------------------
        # scib 1.1.7 ships a compiled C++ binary at
        #   /opt/apps/python/.../scib/knn_graph/knn_graph.o
        # linked against glibc 2.38 and libstdc++ GLIBCXX_3.4.32. The container's
        # libc/libstdc++ are older, so the binary fails to load ANY shared
        # library at invocation time — both the n_cores>1 (multiprocessing)
        # and n_cores=1 (single-call) code paths produce no output files, and
        # scib's Python reader then crashes with FileNotFoundError.
        # This is NOT the tempdir race we originally hypothesized.
        # Rebuilding the scib C++ binary inside the container would fix it,
        # but the engineering cost isn't justified when the remaining 8 metrics
        # (NMI, ARI, ASW_batch, ASW_label, isolated_labels_F1, isolated_labels_ASW,
        # graph_conn, PCR) are already a strong ranker set for winner selection.
        # Observed across attempts 50445270, 50450620, 50454719, 50459583, 50463673.
        # --------------------------------------------------------------------------
        n_cores = int(os.environ.get("SLURM_CPUS_PER_TASK", "4"))
        try:
            print("\n--- Computing scIB metrics (LISI disabled) ---", flush=True)
            results = scib.metrics.metrics(
                adata,
                adata_int=adata,           # same object — embedding is in obsm
                batch_key=args.batch_key,
                label_key=args.label_key,
                embed=args.embedding_key,
                organism="human",          # CRITICAL: scib defaults to mouse
                n_cores=n_cores,
                isolated_labels_asw_=True,
                silhouette_=True,
                hvg_score_=False,
                graph_conn_=True,
                pcr_=True,
                isolated_labels_f1_=True,
                kBET_=not args.skip_kbet,
                ilisi_=False,              # DISABLED — see glibc note above
                clisi_=False,              # DISABLED — see glibc note above
                nmi_=True,
                ari_=True,
                cell_cycle_=False,
                trajectory_=False,
            )
            print(f"  scIB metrics computed in {(time.time() - t0) / 60:.1f} min",
                  flush=True)
        except Exception as e:
            error_msg = str(e)[:500]
            print(f"  ERROR computing scIB metrics: {error_msg}", flush=True)
            import traceback
            traceback.print_exc()
            out_row = {
                "method": args.method_name,
                "compartment": args.compartment,
                "embedding_key": args.embedding_key,
                "n_cells": int(adata.n_obs),
                "n_genes": int(adata.n_vars),
                "scib_version": SCIB_VERSION,
                "error": error_msg,
            }
            pd.DataFrame([out_row]).to_csv(args.output_csv, index=False)
            print(f"  Wrote error row: {args.output_csv}", flush=True)
            sys.exit(1)

        # Format results from scib.metrics.metrics()
        if isinstance(results, pd.DataFrame):
            out_row = results.iloc[:, 0].to_dict() if results.shape[1] > 0 else {}
        else:
            out_row = dict(results)

    # --- Shared output formatting (both paths) ---
    out_row["method"] = args.method_name
    out_row["compartment"] = args.compartment
    out_row["embedding_key"] = args.embedding_key
    out_row["n_cells"] = int(adata.n_obs)
    out_row["n_genes"] = int(adata.n_vars)
    out_row["scib_version"] = SCIB_VERSION
    out_row["kbet_skipped"] = args.skip_kbet

    out_df = pd.DataFrame([out_row])
    out_df.to_csv(args.output_csv, index=False)
    print(f"\n  Saved: {args.output_csv}", flush=True)
    print(f"  Metrics: {list(out_row.keys())}", flush=True)

    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
