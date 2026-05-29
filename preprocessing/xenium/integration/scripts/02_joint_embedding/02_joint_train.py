"""Joint FLEX+Xenium Concord training for one sweep config.


Differences:
  - Xenium counts come from a pooled, Phase-1-filtered bundle named
    xenium_{source}_counts.mtx.gz (source = nuclear / whole / cytoplasmic).
  - FLEX QC comes from Phase 1 cohort CSV (01_qc/outputs/flex_qc.csv).
  - One config = {source, n_latent}.

domain_key = patient_platform (8 domains: 4 patients × 2 platforms).
"""
import argparse, gzip, io, os
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp
import torch
import concord as ccd


FLEX_COMPARTMENTS = ["Immune", "Epithelial", "Stromal"]


def load_flex_compartment(flex_dir, comp, flex_qc):
    cdir = os.path.join(flex_dir, comp, "scvi_n100")
    with gzip.open(os.path.join(cdir, "counts.mtx.gz"), "rb") as gz:
        mat = scipy.io.mmread(io.BytesIO(gz.read()))
    X = sp.csr_matrix(mat)
    genes = pd.read_csv(os.path.join(cdir, "genes.tsv"), header=None)[0].tolist()
    cells = pd.read_csv(os.path.join(cdir, "cells.tsv"), header=None)[0].tolist()
    obs = pd.read_csv(os.path.join(cdir, "obs.csv"), index_col=0)
    assert list(obs.index) == cells
    ad_c = ad.AnnData(X=X, obs=obs)
    ad_c.var_names = genes
    ad_c.obs["compartment"] = comp

    qc = flex_qc[flex_qc["compartment"] == comp].set_index("cell_id")
    keep = qc.index[qc["pass_qc"]].intersection(ad_c.obs_names)
    ad_c = ad_c[keep].copy()
    print(f"[FLEX/{comp}] {ad_c.shape}", flush=True)
    return ad_c


def load_xenium_bundle(bundle_dir, source):
    mtx = os.path.join(bundle_dir, f"xenium_{source}_counts.mtx.gz")
    with gzip.open(mtx, "rb") as gz:
        m = scipy.io.mmread(io.BytesIO(gz.read()))
    X = sp.csr_matrix(m)
    genes = pd.read_csv(os.path.join(bundle_dir, "xenium_genes.tsv"),
                        header=None)[0].tolist()
    cells = pd.read_csv(os.path.join(bundle_dir, "xenium_cells.tsv"),
                        header=None)[0].tolist()
    obs = pd.read_csv(os.path.join(bundle_dir, "xenium_obs.csv"))
    obs.index = obs["cell_id"]
    obs = obs.reindex(cells)
    ad_x = ad.AnnData(X=X, obs=obs)
    ad_x.var_names = genes
    print(f"[XENIUM/{source}] {ad_x.shape}", flush=True)
    return ad_x


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--flex-dir", required=True)
    p.add_argument("--flex-qc", required=True,
                   help="01_qc/outputs/flex_qc.csv")
    p.add_argument("--xenium-bundle", required=True,
                   help="dir with xenium_{source}_counts.mtx.gz")
    p.add_argument("--source", required=True,
                   choices=["nuclear", "whole", "cytoplasmic"])
    p.add_argument("--panel-intersection", required=True)
    p.add_argument("--compartment-dir", default=None,
                   help="If provided, look up Xenium compartment from "
                        "{dir}/{compartment}/leiden_assignments.csv (output of "
                        "phase 09 compartment classifier). On first runs where "
                        "phase 09 hasn't run yet, this dir won't exist and "
                        "Xenium cells fall back to 'unknown'.")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--n-latent", type=int, required=True)
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    os.makedirs(os.path.join(args.out_dir, "model"), exist_ok=True)

    flex_qc = pd.read_csv(args.flex_qc)
    flex_parts = [load_flex_compartment(args.flex_dir, c, flex_qc)
                  for c in FLEX_COMPARTMENTS]
    flex = ad.concat(flex_parts, axis=0, join="outer", merge="same",
                     index_unique=None)
    flex.obs["platform"] = "flex"
    flex.obs["patient_platform"] = flex.obs["patient_id"].astype(str) + "_flex"

    xen = load_xenium_bundle(args.xenium_bundle, args.source)
    xen.obs["platform"] = "xenium"
    xen.obs["patient_platform"] = xen.obs["patient_id"].astype(str) + "_xenium"

    # Compartment for Xenium cells comes from phase 09 compartment classifier,
    # which writes per-compartment cell lists to
    # {compartment_dir}/{compartment}/leiden_assignments.csv. On first runs where
    # phase 09 hasn't completed yet, the directory is missing and we fall back to
    # "unknown" — phase 09 will populate the lists and a re-run of phase 02 will
    # pick up the canonical assignment.
    xen.obs["compartment"] = "unknown"
    if args.compartment_dir and os.path.isdir(args.compartment_dir):
        cell_to_comp = {}
        for comp in FLEX_COMPARTMENTS:
            path = os.path.join(args.compartment_dir, comp.lower(), "leiden_assignments.csv")
            if not os.path.exists(path):
                print(f"  Compartment lookup: missing {path} — Xenium cells in this "
                      f"compartment will stay 'unknown'", flush=True)
                continue
            cells = pd.read_csv(path, usecols=["cell_id"])["cell_id"].tolist()
            for cid in cells:
                cell_to_comp[cid] = comp
            print(f"  Compartment lookup: {comp} {len(cells):,} cells", flush=True)
        if cell_to_comp:
            xen.obs["compartment"] = xen.obs.index.map(
                lambda c: cell_to_comp.get(c, "unknown"))
            n_assigned = (xen.obs["compartment"] != "unknown").sum()
            print(f"  Xenium compartment assigned: {n_assigned:,}/{len(xen):,}", flush=True)
    else:
        print("  No --compartment-dir; Xenium compartment = 'unknown' "
              "(will be patched by phase 09 sync, or set on next phase-02 run "
              "after phase 09 completes)", flush=True)

    panel = [g.strip() for g in open(args.panel_intersection) if g.strip()]
    panel = [g for g in panel if g in flex.var_names and g in xen.var_names]
    print(f"Panel features intersected: {len(panel)}", flush=True)
    flex = flex[:, panel].copy()
    xen = xen[:, panel].copy()

    keep_cols = ["patient_id", "platform", "patient_platform", "compartment"]
    flex.obs = flex.obs[[c for c in keep_cols if c in flex.obs.columns]]
    xen.obs = xen.obs[[c for c in keep_cols if c in xen.obs.columns]]

    joint = ad.concat([flex, xen], axis=0, join="outer", merge="same",
                      index_unique=None)
    print(f"Joint: {joint.shape}  "
          f"platforms={joint.obs['platform'].value_counts().to_dict()}", flush=True)
    print(f"  patient_platforms="
          f"{joint.obs['patient_platform'].value_counts().to_dict()}", flush=True)

    sc.pp.normalize_total(joint, target_sum=1e4)
    sc.pp.log1p(joint)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}", flush=True)

    con = ccd.Concord(
        joint,
        save_dir=os.path.join(args.out_dir, "model"),
        domain_key="patient_platform",
        latent_dim=args.n_latent,
        normalize_total=False,
        log1p=False,
        device=device,
        seed=args.seed,
    )
    con.fit_transform(output_key="X_concord", save_model=True)
    latent = joint.obsm["X_concord"]
    assert not np.isnan(latent).any()
    print(f"Latent: {latent.shape}  "
          f"range=[{latent.min():.3f}, {latent.max():.3f}]", flush=True)

    pd.DataFrame(latent, index=joint.obs_names,
                 columns=[f"concord_{i}" for i in range(latent.shape[1])]
                 ).to_csv(os.path.join(args.out_dir, "joint_latent.csv"))
    joint.obs.to_csv(os.path.join(args.out_dir, "joint_obs.csv"))
    print("Done.")


if __name__ == "__main__":
    main()
