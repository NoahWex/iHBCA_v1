#!/usr/bin/env python3
# ============================================================================
# Step 08: scVI Cross-Patient Integration - Wrapper
# ============================================================================
# Purpose
#   First-round cross-patient scVI integration of the Chromium FLEX spatial
#   breast atlas. Trains a scVI variational autoencoder on the union of
#   step 05 post-QC variable features across all THREE-WAY-QC-passing cells
#   (step_01_umi_pass AND step_02_mad_pass AND step_03_doublet_pass) and
#   emits the canonical latent space + UMAP + leiden clusters that
#   downstream step 09 / 10a / 10b consume.
#
# Methodology
#   1. Load raw H5 counts per sample (manifest-driven), restrict to passing
#      cells from the central_cell_status.csv aggregation.
#   2. Compute a global union of per-sample post-QC VFs (step 05 outputs).
#      This preserves the full biological signal per patient while avoiding
#      the bias of selecting VFs on the merged object.
#   3. Standardize all patient AnnData objects to the union gene space
#      (missing genes zero-filled), then merge.
#   4. Train scVI with patient_id as the batch key. Uses negative-binomial
#      gene likelihood, 3 VAE layers, early stopping on a 0.9 train split.
#      The n_latent + n_epochs are read from module_configs.yaml
#      (step_08_scvi_integration.algorithm_params). The choice of 50 latent
#      dimensions matches Reed et al. 2024 and is stable across dropout
#      levels typical for the FLEX panel.
#   5. Generate scVI latent, Leiden clustering at the configured
#      resolutions, and UMAP.
#   6. Write categorical sidecars (latent, umap, clusters, metadata, model).
#      Per the revised migration model the scVI model checkpoint at
#      $CFG_PREPROCESSING_STEP_08_MODEL is CANONICAL and reused downstream;
#      this wrapper supports a --load-model flag to skip training and
#      reuse the frozen checkpoint.
#
# Inputs
#   - $CFG_PREPROCESSING_CENTRAL_CELL_STATUS (THREE-WAY gate)
#   - $CFG_PREPROCESSING_STEP_05/features/filtered_vfs/{sample_id}_vfs_filtered.txt
#   - $CFG_RAW_DATA_SAMPLE_MANIFEST (TSV: patient_id, position_id, sample_id, h5_path)
#   - $CFG_MODULE_CONFIGS (algorithm_params block)
#   - $CFG_PREPROCESSING_STEP_08_MODEL (canonical, used when --load-model)
#
# Outputs
#   - ${CFG_PREPROCESSING_STEP_08}/model/scvi_model/                (canonical checkpoint)
#   - ${CFG_PREPROCESSING_STEP_08}/embeddings/latent_50d.csv
#   - ${CFG_PREPROCESSING_STEP_08}/embeddings/umap_2d.csv
#   - ${CFG_PREPROCESSING_STEP_08}/clusters/leiden_res_{res}.csv    (one per resolution)
#   - ${CFG_PREPROCESSING_STEP_08}/metadata/cell_metadata_integrated.csv
#   - ${CFG_PREPROCESSING_STEP_08}/integrated_objects/integrated_adata.h5ad
#
# Note: --load-model uses the canonical checkpoint at
# $CFG_PREPROCESSING_STEP_08_MODEL as an input, not an output, so re-runs
# do not retrain GPU-stochastic weights by default.
# ============================================================================

import argparse
import json
import logging
import os
import sys
from datetime import datetime
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s step_08 %(message)s",
)
log = logging.getLogger("step_08")


def fail_loud(msg, hint=None):
    log.error(msg)
    if hint:
        log.error("  hint: %s", hint)
    sys.exit(2)


def env(name, default=None):
    val = os.environ.get(name, default)
    return val


def parse_args():
    p = argparse.ArgumentParser(description="Step 08: scVI cross-patient integration")
    p.add_argument("--project-root", default=env("CFG_PROJECT_ROOT", ""),
                   help="Pipeline project root (defaults to CFG_PROJECT_ROOT).")
    p.add_argument("--central-manifest", default=env("CFG_PREPROCESSING_CENTRAL_CELL_STATUS", ""),
                   help="Central cell status CSV (THREE-WAY gate aggregation).")
    p.add_argument("--step-05-dir", default=env("CFG_PREPROCESSING_STEP_05", ""),
                   help="Step 05 output root (contains features/filtered_vfs/).")
    p.add_argument("--raw-data-manifest", default=env("CFG_RAW_DATA_SAMPLE_MANIFEST", ""),
                   help="Raw H5 sample manifest TSV (patient_id, position_id, sample_id, h5_path).")
    p.add_argument("--module-configs", default=env("CFG_MODULE_CONFIGS", ""),
                   help="module_configs.yaml with algorithm_params.")
    p.add_argument("--output-root", default=env("CFG_PREPROCESSING_STEP_08", ""),
                   help="Step 08 output root.")
    p.add_argument("--model-path", default=env("CFG_PREPROCESSING_STEP_08_MODEL", ""),
                   help="Canonical scVI model checkpoint dir (read when --load-model).")
    p.add_argument("--load-model", action="store_true",
                   help="Skip training and load the canonical scVI model from --model-path.")
    p.add_argument("--dry-run", action="store_true",
                   help="Validate inputs + print plan, then exit 0 without running scVI.")
    return p.parse_args()


def validate_inputs(args):
    log.info("start")
    log.info("dry_run=%s load_model=%s", args.dry_run, args.load_model)

    if not args.project_root:
        fail_loud("project_root empty", "set --project-root or CFG_PROJECT_ROOT")
    if not args.central_manifest:
        fail_loud("central-manifest empty", "set --central-manifest or CFG_PREPROCESSING_CENTRAL_CELL_STATUS")
    if not Path(args.central_manifest).exists():
        fail_loud(f"central manifest not found: {args.central_manifest}",
                  "run step 01/02/03 aggregators to produce central_cell_status.csv")

    if not args.step_05_dir or not Path(args.step_05_dir).exists():
        fail_loud(f"step 05 dir missing: {args.step_05_dir}",
                  "set --step-05-dir or CFG_PREPROCESSING_STEP_05")

    vf_dir = Path(args.step_05_dir) / "features" / "filtered_vfs"
    if not vf_dir.exists():
        fail_loud(f"step 05 filtered_vfs dir missing: {vf_dir}",
                  "step 05 must write features/filtered_vfs/{sample}_vfs_filtered.txt")

    if not args.raw_data_manifest or not Path(args.raw_data_manifest).exists():
        fail_loud(f"raw data manifest missing: {args.raw_data_manifest}",
                  "set --raw-data-manifest or CFG_RAW_DATA_SAMPLE_MANIFEST")

    if not args.module_configs or not Path(args.module_configs).exists():
        fail_loud(f"module_configs.yaml missing: {args.module_configs}",
                  "set --module-configs or CFG_MODULE_CONFIGS")

    if not args.output_root:
        fail_loud("output-root empty", "set --output-root or CFG_PREPROCESSING_STEP_08")

    if args.load_model:
        if not args.model_path or not Path(args.model_path).exists():
            fail_loud(f"--load-model set but model_path missing: {args.model_path}",
                      "set --model-path or CFG_PREPROCESSING_STEP_08_MODEL to the frozen checkpoint")

    out = Path(args.output_root)
    out.mkdir(parents=True, exist_ok=True)
    for sub in ("model", "integrated_objects", "embeddings", "clusters", "metadata", "logs"):
        (out / sub).mkdir(parents=True, exist_ok=True)
    if not os.access(out, os.W_OK):
        fail_loud(f"output root not writable: {out}")

    # Container environment hints
    for var in ("PYTHONUSERBASE", "NUMBA_CACHE_DIR", "MPLCONFIGDIR"):
        if not env(var):
            log.warning("env_missing %s", var)

    log.info("project_root=%s", args.project_root)
    log.info("central_manifest=%s", args.central_manifest)
    log.info("step_05_vf_dir=%s", vf_dir)
    log.info("raw_data_manifest=%s", args.raw_data_manifest)
    log.info("module_configs=%s", args.module_configs)
    log.info("output_root=%s", args.output_root)
    log.info("model_path=%s", args.model_path or "(will be trained)")
    return vf_dir


def print_plan(args, vf_dir):
    print("=== DRY RUN MODE ===")
    print("Step: step_08_scvi_integration")
    print("Container: scgpt_gpu")
    print("Inputs validated:")
    print(f"  {args.central_manifest}: OK")
    print(f"  {vf_dir}: OK")
    print(f"  {args.raw_data_manifest}: OK")
    print(f"  {args.module_configs}: OK")
    if args.load_model:
        print(f"  {args.model_path}: OK (will reuse, --load-model)")
    print("Outputs planned:")
    print(f"  {args.output_root}/model/scvi_model/")
    print(f"  {args.output_root}/embeddings/latent_50d.csv")
    print(f"  {args.output_root}/embeddings/umap_2d.csv")
    print(f"  {args.output_root}/clusters/leiden_res_*.csv")
    print(f"  {args.output_root}/metadata/cell_metadata_integrated.csv")
    print(f"  {args.output_root}/integrated_objects/integrated_adata.h5ad")
    print("Resources: Tier 5 (1 GPU, 128 GB, 8h) — scgpt_gpu container")
    print("VALIDATION PASSED")


def main():
    args = parse_args()
    vf_dir = validate_inputs(args)

    if args.dry_run:
        print_plan(args, vf_dir)
        sys.exit(0)

    # Defer heavy imports until validation passes (keeps dry-run fast)
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import anndata as ad
    import scvi
    import torch
    import yaml

    # Local source modules
    script_dir = Path(__file__).resolve().parent
    sys.path.insert(0, str(script_dir / "source"))
    from scvi_integration import (  # noqa: E402
        create_gene_union,
        create_vf_subset,
        train_scvi_model,
        generate_embeddings_and_clusters,
        prepare_integrated_adata,
    )

    with open(args.module_configs, "r") as f:
        module_config = yaml.safe_load(f)
    if "step_08_scvi_integration" not in module_config.get("modules", {}):
        fail_loud("step_08_scvi_integration missing from module_configs.yaml",
                  "ensure modules.step_08_scvi_integration.algorithm_params is defined")
    algo_params = dict(module_config["modules"]["step_08_scvi_integration"]["algorithm_params"])
    log.info("scvi_version=%s", scvi.__version__)
    log.info("n_latent=%s max_epochs=%s batch_size=%s",
             algo_params.get("n_latent"), algo_params.get("max_epochs"), algo_params.get("batch_size"))

    # ---- VF union across samples ---------------------------------------
    manifest_df = pd.read_csv(args.raw_data_manifest, sep="\t")
    needed_cols = {"patient_id", "position_id", "sample_id", "h5_path"}
    missing = needed_cols - set(manifest_df.columns)
    if missing:
        fail_loud(f"raw_data_manifest missing columns: {missing}",
                  "TSV schema: patient_id, position_id, sample_id, h5_path")
    log.info("raw_manifest_rows=%d", len(manifest_df))

    all_vf_genes = set()
    per_sample_vf_counts = {}
    for _, row in manifest_df.iterrows():
        vf_path = vf_dir / f"{row['sample_id']}_vfs_filtered.txt"
        if not vf_path.exists():
            log.warning("vf_missing sample_id=%s path=%s (will skip)", row["sample_id"], vf_path)
            continue
        with open(vf_path) as f:
            sample_vfs = [line.strip() for line in f if line.strip()]
        per_sample_vf_counts[row["sample_id"]] = len(sample_vfs)
        all_vf_genes.update(sample_vfs)
    if not all_vf_genes:
        fail_loud("no VFs collected from step 05 filtered_vfs/",
                  "check that step 05 ran and produced per-sample VF lists")
    all_vf_genes = sorted(all_vf_genes)
    log.info("global_vf_union=%d", len(all_vf_genes))

    # ---- THREE-WAY cell gate -------------------------------------------
    central = pd.read_csv(args.central_manifest)
    for col in ("step_01_umi_pass", "step_02_mad_pass", "step_03_doublet_pass"):
        if col not in central.columns:
            fail_loud(f"central_cell_status missing {col}",
                      "re-run step 01/02/03 aggregators")
        central[col] = central[col].astype(bool)
    passing = central[
        central.step_01_umi_pass & central.step_02_mad_pass & central.step_03_doublet_pass
    ].copy()
    log.info("cells_total=%d cells_three_way_pass=%d", len(central), len(passing))
    if passing.empty:
        fail_loud("no cells pass THREE-WAY gate",
                  "check central_cell_status.csv — all three pass columns must be True for some cells")

    # ---- Raw count loading per patient ---------------------------------
    patients = sorted(passing["patient_id"].unique().tolist())
    log.info("patients=%s", patients)
    adatas = []
    for pid in patients:
        log.info("load_patient pid=%s", pid)
        pat_df = manifest_df[manifest_df.patient_id == pid]
        pat_adatas = []
        for _, row in pat_df.iterrows():
            sample_id = row["sample_id"]
            h5_path = Path(row["h5_path"])
            if not h5_path.exists():
                log.warning("h5_missing sample=%s path=%s", sample_id, h5_path)
                continue
            sample_cells = passing[passing.sample_id == sample_id]
            if sample_cells.empty:
                log.warning("no_passing_cells sample=%s", sample_id)
                continue
            log.info("  loading sample=%s h5=%s", sample_id, h5_path.name)
            adata = sc.read_10x_h5(str(h5_path))
            adata.var_names_make_unique()
            # Strip sample prefix from cell_ids — raw H5 barcodes are unprefixed
            prefixed = sample_cells["cell_id"].tolist()
            unprefixed = [bc.replace(f"{sample_id}_", "", 1) for bc in prefixed]
            mask = adata.obs.index.isin(unprefixed)
            if mask.sum() == 0:
                log.warning("  no matching barcodes sample=%s", sample_id)
                continue
            adata = adata[mask].copy()
            adata.obs["patient_id"] = pid
            adata.obs["sample_id"] = sample_id
            adata.obs["position_id"] = row["position_id"]
            adata.obs_names = [f"{sample_id}_{bc}" for bc in adata.obs_names]
            pat_adatas.append(adata)
            log.info("  sample=%s matched=%d", sample_id, mask.sum())
        if not pat_adatas:
            log.warning("patient_empty pid=%s", pid)
            continue
        patient_merged = ad.concat(pat_adatas, join="outer", fill_value=0)
        patient_merged.obs_names_make_unique(join="_")
        log.info("  patient=%s shape=%s", pid, patient_merged.shape)
        adatas.append(patient_merged)

    if not adatas:
        fail_loud("no AnnData objects assembled",
                  "check raw H5 paths and barcode prefixing")

    # ---- Gene union + VF subset + scVI training -----------------------
    log.info("gene_union start")
    adatas_union = create_gene_union(adatas)
    adata_merged = ad.concat(adatas_union, join="outer", fill_value=0)
    if adata_merged.obs_names.duplicated().any():
        fail_loud("duplicate cell_ids after merge — sample prefixing failed")
    log.info("merged_shape=%s", adata_merged.shape)

    adata_vf = create_vf_subset(adata_merged, all_vf_genes)
    log.info("vf_subset_shape=%s", adata_vf.shape)

    if torch.cuda.is_available():
        torch.set_float32_matmul_precision("medium")
        log.info("cuda_available tf32=medium")

    model_out = Path(args.output_root) / "model" / "scvi_model"
    if args.load_model:
        log.info("load_model path=%s", args.model_path)
        scvi.model.SCVI.setup_anndata(adata_vf, batch_key=algo_params["batch_key"])
        vae = scvi.model.SCVI.load(args.model_path, adata=adata_vf)
    else:
        vae = train_scvi_model(adata_vf, algo_params)
        vae.save(str(model_out), overwrite=True)
        log.info("model_saved path=%s", model_out)

    # ---- Generate embeddings + clusters --------------------------------
    results = generate_embeddings_and_clusters(vae, adata_merged, algo_params)
    adata_integrated = prepare_integrated_adata(
        adata_merged, results["latent"], results["umap"], results["clusters"], algo_params
    )
    log.info("integrated_shape=%s", adata_integrated.shape)

    # ---- Write categorical sidecars ------------------------------------
    out = Path(args.output_root)
    adata_integrated.write_h5ad(out / "integrated_objects" / "integrated_adata.h5ad")
    import pandas as pd  # noqa: reused

    latent_df = pd.DataFrame(
        results["latent"],
        index=adata_integrated.obs.index,
        columns=[f"scvi_{i+1}" for i in range(results["latent"].shape[1])],
    )
    latent_df.to_csv(out / "embeddings" / "latent_50d.csv")
    log.info("wrote embeddings/latent_50d.csv size=%d", (out / "embeddings" / "latent_50d.csv").stat().st_size)

    umap_df = pd.DataFrame(
        results["umap"], index=adata_integrated.obs.index, columns=["UMAP_1", "UMAP_2"]
    )
    umap_df.to_csv(out / "embeddings" / "umap_2d.csv")

    for col in [c for c in adata_integrated.obs.columns if c.startswith("leiden_scvi_")]:
        res_str = col.replace("leiden_scvi_", "")
        cluster_df = pd.DataFrame({col: adata_integrated.obs[col]}, index=adata_integrated.obs.index)
        cluster_df.to_csv(out / "clusters" / f"leiden_res_{res_str}.csv")

    adata_integrated.obs.to_csv(out / "metadata" / "cell_metadata_integrated.csv")

    payload = {
        "integration_metadata": {
            "patients": patients,
            "total_cells": int(adata_integrated.n_obs),
            "total_genes": int(adata_integrated.n_vars),
            "n_latent": algo_params["n_latent"],
            "batch_key": algo_params["batch_key"],
            "vf_source": "step_05_post_qc_vfs",
            "vf_count": len(all_vf_genes),
            "algorithm": "scvi-tools",
            "algorithm_version": scvi.__version__,
            "loaded_from_canonical_model": args.load_model,
        },
        "timestamp": datetime.now().isoformat(),
    }
    with open(out / "manifest_update_payload.json", "w") as f:
        json.dump(payload, f, indent=2)
    log.info("done cells=%d genes=%d", adata_integrated.n_obs, adata_integrated.n_vars)


if __name__ == "__main__":
    main()
