#!/usr/bin/env python3
# ============================================================================
# Step 12: Filtered scVI Cross-Patient Integration - Wrapper
# ============================================================================
# Purpose
#   Second-round scVI integration after manifold-level artifact removal at
#   step 10b. Re-trains scVI on the FOUR-WAY filtered cell set
#   (step_01_umi_pass AND step_02_mad_pass AND step_03_doublet_pass AND
#   step_10b_manifold_pass) using the Step 11 post-artifact per-sample VFs.
#   Emits latent/UMAP/leiden sidecars that the L3 handoff to compartment
#   classification and downstream sweep consume.
#
# Methodology
#   - Identical scVI algorithm to step 08 (reused via source/step08_reexport.py
#     which imports from ../step_08_scvi_integration/source/scvi_integration.py).
#   - Cell gate differs: FOUR-WAY rather than THREE-WAY. The step_10b_manifold_pass
#     column is expected either directly in central_cell_status.csv (preferred)
#     or in the frozen step10b_manifest_update.csv at
#     $CFG_PREPROCESSING_CELL_RETENTION_LIST (fallback).
#   - VF source differs: Step 11 post-artifact filtered VFs rather than Step 05
#     post-QC VFs. Step 11 recomputes VFs on the cleaned cell population to
#     reflect the post-artifact-removal biology.
#   - Per the revised migration model, the Step 12 scVI model at
#     $CFG_PREPROCESSING_STEP_12_MODEL is CANONICAL. --load-model skips
#     training and reuses the frozen checkpoint (default for re-runs).
#
# Inputs
#   - $CFG_PREPROCESSING_CENTRAL_CELL_STATUS
#   - $CFG_PREPROCESSING_CELL_RETENTION_LIST (step10b_manifest_update.csv, frozen)
#   - $CFG_PREPROCESSING_STEP_11/features/filtered_vfs/{sample_id}_vfs_filtered.txt
#   - $CFG_RAW_DATA_SAMPLE_MANIFEST
#   - $CFG_MODULE_CONFIGS (step_12_filtered_scvi_integration.algorithm_params)
#   - $CFG_PREPROCESSING_STEP_12_MODEL (canonical, used when --load-model)
#
# Outputs
#   - ${CFG_PREPROCESSING_STEP_12}/model/scvi_model/                (canonical checkpoint)
#   - ${CFG_PREPROCESSING_STEP_12}/embeddings/latent_50d.csv
#   - ${CFG_PREPROCESSING_STEP_12}/embeddings/umap_2d.csv
#   - ${CFG_PREPROCESSING_STEP_12}/clusters/leiden_res_{res}.csv
#   - ${CFG_PREPROCESSING_STEP_12}/metadata/cell_metadata_integrated.csv
#   - ${CFG_PREPROCESSING_STEP_12}/integrated_objects/integrated_adata.h5ad
#
# Note: scVI algorithms are reused from step 08 via source/step08_reexport.py
# (not reimplemented — preserves the "same model, different inputs" contract).
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
    format="%(asctime)s %(levelname)s step_12 %(message)s",
)
log = logging.getLogger("step_12")


def fail_loud(msg, hint=None):
    log.error(msg)
    if hint:
        log.error("  hint: %s", hint)
    sys.exit(2)


def env(name, default=None):
    return os.environ.get(name, default)


def parse_args():
    p = argparse.ArgumentParser(description="Step 12: filtered scVI cross-patient integration")
    p.add_argument("--project-root", default=env("CFG_PROJECT_ROOT", ""))
    p.add_argument("--central-manifest", default=env("CFG_PREPROCESSING_CENTRAL_CELL_STATUS", ""))
    p.add_argument("--cell-retention-list", default=env("CFG_PREPROCESSING_CELL_RETENTION_LIST", ""),
                   help="Frozen step10b_manifest_update.csv with step_10b_manifold_pass column.")
    p.add_argument("--step-11-dir", default=env("CFG_PREPROCESSING_STEP_11", ""),
                   help="Step 11 output root containing features/filtered_vfs/.")
    p.add_argument("--raw-data-manifest", default=env("CFG_RAW_DATA_SAMPLE_MANIFEST", ""))
    p.add_argument("--module-configs", default=env("CFG_MODULE_CONFIGS", ""))
    p.add_argument("--output-root", default=env("CFG_PREPROCESSING_STEP_12", ""))
    p.add_argument("--model-path", default=env("CFG_PREPROCESSING_STEP_12_MODEL", ""))
    p.add_argument("--load-model", action="store_true",
                   help="Skip training and load canonical step 12 scVI model.")
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def validate_inputs(args):
    log.info("start dry_run=%s load_model=%s", args.dry_run, args.load_model)
    if not args.project_root:
        fail_loud("project_root empty", "set --project-root or CFG_PROJECT_ROOT")
    if not args.central_manifest or not Path(args.central_manifest).exists():
        fail_loud(f"central manifest missing: {args.central_manifest}",
                  "set --central-manifest or CFG_PREPROCESSING_CENTRAL_CELL_STATUS")
    if not args.cell_retention_list or not Path(args.cell_retention_list).exists():
        fail_loud(f"cell retention list missing: {args.cell_retention_list}",
                  "set --cell-retention-list or CFG_PREPROCESSING_CELL_RETENTION_LIST "
                  "(= step10b_manifest_update.csv at step_10b/exports/)")
    if not args.step_11_dir or not Path(args.step_11_dir).exists():
        fail_loud(f"step 11 dir missing: {args.step_11_dir}",
                  "set --step-11-dir or CFG_PREPROCESSING_STEP_11")
    vf_dir = Path(args.step_11_dir) / "features" / "filtered_vfs"
    if not vf_dir.exists():
        fail_loud(f"step 11 filtered_vfs dir missing: {vf_dir}",
                  "step 11 must write features/filtered_vfs/{sample}_vfs_filtered.txt")
    if not args.raw_data_manifest or not Path(args.raw_data_manifest).exists():
        fail_loud(f"raw data manifest missing: {args.raw_data_manifest}")
    if not args.module_configs or not Path(args.module_configs).exists():
        fail_loud(f"module_configs.yaml missing: {args.module_configs}")
    if not args.output_root:
        fail_loud("output-root empty")
    if args.load_model:
        if not args.model_path or not Path(args.model_path).exists():
            fail_loud(f"--load-model set but model path missing: {args.model_path}",
                      "set --model-path or CFG_PREPROCESSING_STEP_12_MODEL")

    out = Path(args.output_root)
    out.mkdir(parents=True, exist_ok=True)
    for sub in ("model", "integrated_objects", "embeddings", "clusters", "metadata", "logs"):
        (out / sub).mkdir(parents=True, exist_ok=True)
    if not os.access(out, os.W_OK):
        fail_loud(f"output root not writable: {out}")

    for var in ("PYTHONUSERBASE", "NUMBA_CACHE_DIR", "MPLCONFIGDIR"):
        if not env(var):
            log.warning("env_missing %s", var)

    log.info("project_root=%s", args.project_root)
    log.info("central_manifest=%s", args.central_manifest)
    log.info("cell_retention_list=%s", args.cell_retention_list)
    log.info("step_11_vf_dir=%s", vf_dir)
    log.info("raw_data_manifest=%s", args.raw_data_manifest)
    log.info("module_configs=%s", args.module_configs)
    log.info("output_root=%s", args.output_root)
    log.info("model_path=%s", args.model_path or "(will be trained)")
    return vf_dir


def print_plan(args, vf_dir):
    print("=== DRY RUN MODE ===")
    print("Step: step_12_filtered_scvi_integration")
    print("Container: scgpt_gpu")
    print("Filter: FOUR-WAY (step_01 AND step_02 AND step_03 AND step_10b_manifold_pass)")
    print("Inputs validated:")
    print(f"  {args.central_manifest}: OK")
    print(f"  {args.cell_retention_list}: OK (frozen)")
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
    print("Resources: Tier 5 (1 GPU, 128 GB, 8h)")
    print("VALIDATION PASSED")


def main():
    args = parse_args()
    vf_dir = validate_inputs(args)
    if args.dry_run:
        print_plan(args, vf_dir)
        sys.exit(0)

    import numpy as np
    import pandas as pd
    import scanpy as sc
    import anndata as ad
    import scvi
    import torch
    import yaml

    script_dir = Path(__file__).resolve().parent
    sys.path.insert(0, str(script_dir / "source"))
    # step08_reexport re-exports the pure scVI functions
    from step08_reexport import (  # noqa: E402
        create_gene_union,
        create_vf_subset,
        train_scvi_model,
        generate_embeddings_and_clusters,
        prepare_integrated_adata,
    )

    with open(args.module_configs, "r") as f:
        module_config = yaml.safe_load(f)
    mod_key = "step_12_filtered_scvi_integration"
    if mod_key not in module_config.get("modules", {}):
        fail_loud(f"{mod_key} missing from module_configs.yaml")
    algo_params = dict(module_config["modules"][mod_key]["algorithm_params"])
    log.info("scvi_version=%s n_latent=%s", scvi.__version__, algo_params.get("n_latent"))

    # ---- VF union across samples (step 11) -----------------------------
    manifest_df = pd.read_csv(args.raw_data_manifest, sep="\t")
    needed = {"patient_id", "position_id", "sample_id", "h5_path"}
    if not needed.issubset(manifest_df.columns):
        fail_loud(f"raw_data_manifest missing columns {needed - set(manifest_df.columns)}")
    log.info("raw_manifest_rows=%d", len(manifest_df))

    all_vf_genes = set()
    for _, row in manifest_df.iterrows():
        vf_path = vf_dir / f"{row['sample_id']}_vfs_filtered.txt"
        if not vf_path.exists():
            log.warning("step11_vf_missing sample=%s", row["sample_id"])
            continue
        with open(vf_path) as f:
            all_vf_genes.update(line.strip() for line in f if line.strip())
    if not all_vf_genes:
        fail_loud("no step 11 VFs collected",
                  "step 11 must run per-sample before step 12")
    all_vf_genes = sorted(all_vf_genes)
    log.info("step11_vf_union=%d", len(all_vf_genes))

    # ---- FOUR-WAY cell gate --------------------------------------------
    central = pd.read_csv(args.central_manifest)
    for col in ("step_01_umi_pass", "step_02_mad_pass", "step_03_doublet_pass"):
        if col not in central.columns:
            fail_loud(f"central_cell_status missing {col}")
        central[col] = central[col].astype(bool)

    if "step_10b_manifold_pass" not in central.columns:
        log.info("step_10b column absent from central manifest — merging from cell_retention_list")
        retention = pd.read_csv(args.cell_retention_list)
        if "cell_id" not in retention.columns or "step_10b_manifold_pass" not in retention.columns:
            fail_loud("cell_retention_list missing cell_id or step_10b_manifold_pass",
                      "check the frozen step10b_manifest_update.csv schema")
        retention["step_10b_manifold_pass"] = retention["step_10b_manifold_pass"].astype(bool)
        merge_cols = ["cell_id", "step_10b_manifold_pass"]
        if "step_10b_removal_reason" in retention.columns:
            merge_cols.append("step_10b_removal_reason")
        central = central.merge(retention[merge_cols], on="cell_id", how="left")
    else:
        central["step_10b_manifold_pass"] = central["step_10b_manifold_pass"].astype(bool)

    n_missing_10b = central["step_10b_manifold_pass"].isna().sum()
    if n_missing_10b > 0:
        log.warning("cells_without_10b_flag=%d (will be excluded from FOUR-WAY gate)", n_missing_10b)

    passing = central[
        central.step_01_umi_pass
        & central.step_02_mad_pass
        & central.step_03_doublet_pass
        & central.step_10b_manifold_pass.notna()
        & central.step_10b_manifold_pass
    ].copy()

    three_way = (
        central.step_01_umi_pass & central.step_02_mad_pass & central.step_03_doublet_pass
    ).sum()
    log.info("cells_total=%d three_way_pass=%d four_way_pass=%d removed_by_10b=%d",
             len(central), int(three_way), len(passing), int(three_way) - len(passing))
    if passing.empty:
        fail_loud("no cells pass FOUR-WAY gate — check step 10b retention list")

    # ---- Raw count loading per patient ---------------------------------
    patients = sorted(passing["patient_id"].unique().tolist())
    adatas = []
    for pid in patients:
        pat_df = manifest_df[manifest_df.patient_id == pid]
        pat_adatas = []
        for _, row in pat_df.iterrows():
            sample_id = row["sample_id"]
            h5_path = Path(row["h5_path"])
            if not h5_path.exists():
                log.warning("h5_missing sample=%s", sample_id)
                continue
            sample_cells = passing[passing.sample_id == sample_id]
            if sample_cells.empty:
                continue
            log.info("loading sample=%s n_cells=%d", sample_id, len(sample_cells))
            adata = sc.read_10x_h5(str(h5_path))
            adata.var_names_make_unique()
            prefixed = sample_cells["cell_id"].tolist()
            unprefixed = [bc.replace(f"{sample_id}_", "", 1) for bc in prefixed]
            mask = adata.obs.index.isin(unprefixed)
            if mask.sum() == 0:
                log.warning("no_matching_barcodes sample=%s", sample_id)
                continue
            adata = adata[mask].copy()
            adata.obs["patient_id"] = pid
            adata.obs["sample_id"] = sample_id
            adata.obs["position_id"] = row["position_id"]
            adata.obs_names = [f"{sample_id}_{bc}" for bc in adata.obs_names]
            pat_adatas.append(adata)
        if not pat_adatas:
            log.warning("patient_empty pid=%s", pid)
            continue
        patient_merged = ad.concat(pat_adatas, join="outer", fill_value=0)
        patient_merged.obs_names_make_unique(join="_")
        adatas.append(patient_merged)

    if not adatas:
        fail_loud("no AnnData assembled after FOUR-WAY gate")

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

    model_out = Path(args.output_root) / "model" / "scvi_model"
    if args.load_model:
        log.info("load_model path=%s", args.model_path)
        scvi.model.SCVI.setup_anndata(adata_vf, batch_key=algo_params["batch_key"])
        vae = scvi.model.SCVI.load(args.model_path, adata=adata_vf)
    else:
        vae = train_scvi_model(adata_vf, algo_params)
        vae.save(str(model_out), overwrite=True)
        log.info("model_saved path=%s", model_out)

    results = generate_embeddings_and_clusters(vae, adata_merged, algo_params)
    adata_integrated = prepare_integrated_adata(
        adata_merged, results["latent"], results["umap"], results["clusters"], algo_params
    )
    log.info("integrated_shape=%s", adata_integrated.shape)

    out = Path(args.output_root)
    adata_integrated.write_h5ad(out / "integrated_objects" / "integrated_adata.h5ad")

    latent_df = pd.DataFrame(
        results["latent"],
        index=adata_integrated.obs.index,
        columns=[f"scvi_{i+1}" for i in range(results["latent"].shape[1])],
    )
    latent_df.to_csv(out / "embeddings" / "latent_50d.csv")

    umap_df = pd.DataFrame(
        results["umap"], index=adata_integrated.obs.index, columns=["UMAP_1", "UMAP_2"]
    )
    umap_df.to_csv(out / "embeddings" / "umap_2d.csv")

    for col in [c for c in adata_integrated.obs.columns if c.startswith("leiden_scvi_")]:
        res_str = col.replace("leiden_scvi_", "")
        pd.DataFrame(
            {col: adata_integrated.obs[col]}, index=adata_integrated.obs.index
        ).to_csv(out / "clusters" / f"leiden_res_{res_str}.csv")

    adata_integrated.obs.to_csv(out / "metadata" / "cell_metadata_integrated.csv")

    payload = {
        "integration_metadata": {
            "patients": patients,
            "total_cells": int(adata_integrated.n_obs),
            "total_genes": int(adata_integrated.n_vars),
            "n_latent": algo_params["n_latent"],
            "batch_key": algo_params["batch_key"],
            "vf_source": "step_11_filtered_vfs",
            "vf_count": len(all_vf_genes),
            "filter_type": "FOUR-WAY",
            "filter_components": [
                "step_01_umi_pass",
                "step_02_mad_pass",
                "step_03_doublet_pass",
                "step_10b_manifold_pass",
            ],
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
