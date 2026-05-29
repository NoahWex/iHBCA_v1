#!/usr/bin/env python3
"""
extract_panelC_features.py - Stage-1 feature extraction for Panel C variants.

Per L2 x contrast: pulls leading-edge gene set per viable NhoodGroup, ranked by
NES x logFC (positive direction only - top up-regulated drivers of enriched
pathways), assembles cell x gene counts subset + cells.tsv + genes.tsv + ng_summary.

Plot scripts (DoHeatmap, leading-edge heatmap, dotplot) consume these artifacts
without re-loading milo / counts.npz.

Gene selection rule (locked 2026-05-07 with Noah):
  Per NG:
    - GSEA pathways: padj < 0.10 AND NES > 0
    - F.3 markers:   FDR < 0.10 AND logFC > 0 AND not is_dissoc_stress
    - Inner-join leadingEdge ENSG -> F.3 markers
    - score(gene, pathway) = NES_pathway * logFC_gene
    - Per gene, take max(score) across pathways it leads
    - Top --top-genes-per-ng by score
  Across NGs: union -> final gene panel.

Outputs (per L2 x contrast):
  outputs/panelC_features/{L2_safe}_{contrast}/
    counts.mtx.gz                 sparse cells x genes (raw counts)
    cells.tsv                     cell_id, donor, NG_dominant, NG_rank, parity_binary, risk_class
    genes.tsv                     gene_id, symbol, source_NG, source_rank, score, n_pathways, top_pathway
    ng_summary.csv                NhoodGroup, name_full, rank, n_cells, group_med_lfc, n_sig_pathways
    feature_manifest.yaml         provenance + parameters
"""
from __future__ import annotations

import argparse
import gzip
import json
import logging
import os
import sys
from collections import Counter
from pathlib import Path

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib_config")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.load_paths import load_inquiry  # noqa: E402

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("extract_panelC")


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--inquiry-dir", required=True, type=Path)
    p.add_argument("--L2", required=True, help='e.g. "epi::BMYO-basal"')
    p.add_argument("--contrast", default="parity_x_HR_BRCA1")
    p.add_argument("--top-genes-per-ng", type=int, default=25)
    p.add_argument("--top-genes-per-pathway", type=int, default=5,
                   help="cap leading-edge picks per pathway to force diversity")
    p.add_argument("--candidate-multiplier", type=int, default=3,
                   help="oversample top-N per NG by this factor before pct-in filter")
    p.add_argument("--min-pct-in", type=float, default=0.25,
                   help="min fraction of cells in source NG expressing the gene")
    p.add_argument("--gsea-padj", type=float, default=0.10)
    p.add_argument("--marker-fdr", type=float, default=0.10)
    p.add_argument("--counts-root", type=Path, default=None,
                   help="Override; default resolves via paths.yaml inputs.components_root")
    p.add_argument("--atlas-nhoods-dir", type=Path, default=None,
                   help="atlas_nhoods directory (cells.tsv + nhoods.mtx). "
                         "Default: <inquiry-dir>/outputs/atlas_nhoods.")
    p.add_argument("--labels-full", default=None,
                   help="Override; default resolves via paths.yaml inputs.labels")
    p.add_argument("--max-cells-per-ng", type=int, default=0,
                   help="0 = keep all cells; >0 = randomly downsample per NG")
    return p.parse_args()


def gsea_ng_name(compartment: str, parent_label: str, nhood_group: str) -> str:
    """Match GSEA file's NhoodGroup_renamed: 'epi__BMYO_basal_23'."""
    return f"{compartment}__{parent_label.replace('-', '_')}_{nhood_group}"


def f3_marker_filename(compartment: str, parent_label: str, nhood_group: str) -> str:
    """Match F.3 marker file: 'epi__BMYO_basal_23_vs_parent.csv'."""
    return f"{compartment}__{parent_label.replace('-', '_')}_{nhood_group}_vs_parent.csv"


def main():
    a = parse_args()
    inq = a.inquiry_dir
    cfg = load_inquiry(inq)["paths"]["inputs"]
    if a.counts_root is None:
        a.counts_root = Path(cfg["components_root"])
    if a.labels_full is None:
        a.labels_full = cfg["labels"]
    L2_compartment, L2_label = a.L2.split("::", 1)
    L2_safe = a.L2.replace("::", "__").replace("-", "_")
    if L2_compartment not in {"epi", "imm", "str"}:
        sys.exit(f"unknown compartment: {L2_compartment}")

    # ---------- Lookup: viable NGs (use lookup.csv if present, else build from summary) ----------
    lookup_path = inq / "outputs/nhoodgroup_renaming/lookup.csv"
    if lookup_path.exists():
        lookup = pd.read_csv(lookup_path)
        lookup = lookup[(lookup.viable == True) &
                        (lookup.contrast == a.contrast) &
                        (lookup.parent_L2_joint == a.L2)].copy()
    else:
        log.info("lookup.csv absent; building from F.1 nhood_groups_summary")
        sum_path = inq / "outputs/stageF1_nhoodgroups" / a.contrast / "nhood_groups_summary.csv"
        sum_df = pd.read_csv(sum_path)
        sum_df = sum_df[sum_df.parent_L2_joint == a.L2].copy()
        sum_df = sum_df[sum_df.n_nhoods_in_group >= 10]  # viability threshold
        # rank by descending n_nhoods_in_group within parent_L2
        sum_df = sum_df.sort_values("n_nhoods_in_group", ascending=False)
        sum_df["size_rank_within_L2"] = range(1, len(sum_df) + 1)
        # Derive stratum from contrast (e.g., BR1_vs_AR -> BR1)
        stratum = a.contrast.split("_vs_")[0] if "_vs_" in a.contrast else a.contrast
        sum_df["name_short"] = [f"{L2_label}_{r}"
                                  for r in sum_df.size_rank_within_L2]
        sum_df["name_full"] = [f"{L2_label}_{r}_{stratum}"
                                 for r in sum_df.size_rank_within_L2]
        sum_df["parent_compartment"] = L2_compartment
        sum_df["parent_label"] = L2_label
        sum_df["parent_L2_joint"] = a.L2
        sum_df["contrast"] = a.contrast
        sum_df["viable"] = True
        lookup = sum_df
    lookup["NhoodGroup"] = pd.to_numeric(lookup["NhoodGroup"]).astype(int).astype(str)
    lookup = lookup.sort_values("size_rank_within_L2").reset_index(drop=True)
    if lookup.empty:
        sys.exit(f"No viable NGs for {a.L2} in {a.contrast}")
    lookup["NG_gsea_name"] = [gsea_ng_name(L2_compartment, L2_label, ng)
                              for ng in lookup.NhoodGroup]
    log.info("Viable NGs: %d", len(lookup))
    log.info("\n%s", lookup[["NhoodGroup", "name_full", "n_nhoods_in_group",
                              "group_med_lfc", "size_rank_within_L2",
                              "NG_gsea_name"]].to_string(index=False))

    # ---------- GSEA leading-edge (positive NES, padj < 0.10) ----------
    gsea_path = inq / "outputs/stageI3_gsea/per_nhoodgroup_nes.csv"
    log.info("Reading GSEA: %s", gsea_path)
    # 42MB; column subset + chunked read for memory hygiene
    use_cols = ["contrast", "NhoodGroup_renamed", "pathway", "NES", "padj",
                "leadingEdge"]
    gsea_iter = pd.read_csv(gsea_path, usecols=use_cols, chunksize=200_000)
    keep_set = set(lookup.NG_gsea_name)
    gsea_chunks = []
    for ch in gsea_iter:
        ch = ch[(ch.contrast == a.contrast) &
                 (ch.NhoodGroup_renamed.isin(keep_set)) &
                 (ch.padj < a.gsea_padj) &
                 (ch.NES > 0)]
        if not ch.empty:
            gsea_chunks.append(ch)
    if not gsea_chunks:
        log.warning("No GSEA hits passing filters - feature set will be empty")
        gsea = pd.DataFrame(columns=use_cols)
    else:
        gsea = pd.concat(gsea_chunks, ignore_index=True)
    log.info("GSEA hits: %d rows across %d pathways x %d NGs",
             len(gsea), gsea.pathway.nunique(), gsea.NhoodGroup_renamed.nunique())

    # Expand leadingEdge -> long form (NhoodGroup_renamed, pathway, NES, gene_id)
    if not gsea.empty:
        rows = []
        for _, r in gsea.iterrows():
            le = str(r.leadingEdge) if pd.notna(r.leadingEdge) else ""
            for g in le.split(";"):
                g = g.strip()
                if g:
                    rows.append((r.NhoodGroup_renamed, r.pathway, float(r.NES), g))
        le_df = pd.DataFrame(rows, columns=["NG_gsea_name", "pathway", "NES", "gene_id"])
    else:
        le_df = pd.DataFrame(columns=["NG_gsea_name", "pathway", "NES", "gene_id"])
    log.info("Leading-edge gene-pathway pairs: %d", len(le_df))

    # ---------- F.3 markers per NG; score = NES * logFC, positive only ----------
    marker_dir = inq / "outputs/stageF3_markers" / a.contrast
    per_ng_picks = []
    pathway_records = []  # for sidecar

    for _, r in lookup.iterrows():
        ng = r.NhoodGroup
        rank = int(r.size_rank_within_L2)
        ng_gsea = r.NG_gsea_name
        fp = marker_dir / f3_marker_filename(L2_compartment, L2_label, ng)
        if not fp.exists():
            log.warning("missing F.3 markers: %s", fp.name)
            continue
        m = pd.read_csv(fp)
        # Skipped-NG marker (low cell count): file is "skipped,reason\nTRUE,..."
        if "FDR" not in m.columns or "logFC" not in m.columns or len(m) == 0:
            log.info("  %s (rank %d): F.3 marker file empty/skipped; dropping NG",
                     ng, rank)
            continue
        # required columns: gene_id, symbol, logFC, FDR, is_dissoc_stress
        m["is_dissoc_stress"] = m.get("is_dissoc_stress", False)
        m["is_dissoc_stress"] = m["is_dissoc_stress"].astype(str).str.lower().isin(
            {"true", "1", "yes"})
        m = m[(m.FDR < a.marker_fdr) &
              (m.logFC > 0) &
              (~m.is_dissoc_stress) &
              (m.symbol.notna()) & (m.symbol != "")]
        if m.empty:
            log.info("  %s (rank %d): 0 markers passing filter", ng, rank)
            continue

        # Pull leading-edge for this NG
        le_ng = le_df[le_df.NG_gsea_name == ng_gsea]
        if le_ng.empty:
            log.info("  %s (rank %d): %d markers but 0 leading-edge",
                     ng, rank, len(m))
            continue

        # Inner join leading-edge -> markers (gene_id key)
        joined = le_ng.merge(m[["gene_id", "symbol", "logFC", "FDR"]],
                              on="gene_id", how="inner")
        if joined.empty:
            log.info("  %s (rank %d): leading-edge x markers join is empty",
                     ng, rank)
            continue

        joined["score"] = joined["NES"] * joined["logFC"]

        # Per-pathway diversification: keep top-K per pathway before per-gene
        # aggregation. Prevents MYC_TARGETS / ribosome biogenesis from
        # monopolizing slots and drowning specific markers (e.g. MMP3).
        joined = (joined.sort_values("score", ascending=False)
                          .groupby("pathway", sort=False)
                          .head(a.top_genes_per_pathway))

        # Per-gene best score (max over pathways) + n_pathways + top pathway
        agg = (joined.groupby(["gene_id", "symbol"], sort=False)
                      .agg(score=("score", "max"),
                           logFC=("logFC", "first"),
                           n_pathways=("pathway", "nunique"),
                           top_pathway=("pathway",
                                         lambda s: s.iloc[joined.loc[s.index, "score"].argmax()
                                                          if len(s) > 1 else 0]),
                           top_NES=("NES", "max"))
                      .reset_index())
        # Oversample candidates so pct-in filter has headroom
        agg = (agg.sort_values("score", ascending=False)
                   .head(a.top_genes_per_ng * a.candidate_multiplier))
        agg["source_NG"] = ng
        agg["source_rank"] = rank
        agg["source_NG_full"] = r.name_full
        per_ng_picks.append(agg)

        # pathway sidecar
        for p, sub in joined.groupby("pathway"):
            pathway_records.append({
                "NhoodGroup": ng,
                "name_full": r.name_full,
                "rank": rank,
                "pathway": p,
                "NES": sub.NES.iloc[0],
                "n_le_genes": sub.gene_id.nunique(),
                "le_symbols": ";".join(sorted(sub.symbol.unique()))
            })

        log.info("  %s (rank %d): %d picked genes (top score %.2f, gene %s)",
                 ng, rank, len(agg),
                 agg.score.iloc[0] if len(agg) else float("nan"),
                 agg.symbol.iloc[0] if len(agg) else "-")

    if not per_ng_picks:
        sys.exit(f"No genes selected for any NG in {a.L2} / {a.contrast}")
    genes_picked = pd.concat(per_ng_picks, ignore_index=True)
    log.info("Genes picked across NGs (with dupes): %d", len(genes_picked))

    # Dedupe gene_id; assign provenance to highest-score NG
    genes_picked = genes_picked.sort_values("score", ascending=False)
    gene_set = (genes_picked.drop_duplicates("gene_id", keep="first")
                .reset_index(drop=True))
    log.info("Unique genes (final panel): %d", len(gene_set))

    # ---------- Cell -> dominant NG (adapt from plot_panelC_doheatmap.py) ----------
    nhoods_dir = a.atlas_nhoods_dir if a.atlas_nhoods_dir else (inq / "outputs/atlas_nhoods")
    cells_atlas = pd.read_csv(nhoods_dir / "cells.tsv", header=None)[0].values
    log.info("atlas cells: %d", len(cells_atlas))
    log.info("Loading nhoods.mtx ...")
    nh_mat = sio.mmread(str(nhoods_dir / "nhoods.mtx")).tocsr()
    log.info("nhoods sparse: %s, nnz=%d", nh_mat.shape, nh_mat.nnz)

    ng_csv = inq / "outputs/stageF1_nhoodgroups" / a.contrast / "nhood_groups.csv"
    ng_map = pd.read_csv(ng_csv)
    # Most rows have NhoodGroup=NaN (un-grouped nhoods). NaN coercion -> float64.
    # Convert: drop NaN, cast to int, then to string. Matches lookup format.
    ng_map["NhoodGroup"] = pd.to_numeric(ng_map["NhoodGroup"], errors="coerce")
    ng_map = ng_map.dropna(subset=["NhoodGroup"]).copy()
    ng_map["NhoodGroup"] = ng_map["NhoodGroup"].astype(int).astype(str)
    viable_set = set(lookup.NhoodGroup)
    ng_map = ng_map[ng_map.NhoodGroup.isin(viable_set)]
    log.info("nhoods in viable NGs: %d", len(ng_map))
    nhood_to_ng = dict(zip(ng_map["Nhood"].astype(int), ng_map["NhoodGroup"]))

    # ---------- Compartment counts + gene index ----------
    comp = L2_compartment
    counts_npz = a.counts_root / f"{comp}_counts.npz"
    log.info("Loading %s ...", counts_npz)
    loaded = np.load(counts_npz, allow_pickle=True)
    X = sp.csr_matrix((loaded["data"], loaded["indices"], loaded["indptr"]),
                      shape=tuple(loaded["shape"]))
    log.info("counts: %s (cells x genes)", X.shape)

    cells_comp_path = None
    for cand in (f"{comp}_metadata.csv", f"{comp}_metadata_enriched.csv",
                  f"{comp}_cells.csv"):
        p = a.counts_root / cand
        if p.exists():
            cells_comp_path = p
            break
    if cells_comp_path is None:
        sys.exit(f"No compartment metadata in {a.counts_root}")
    cells_comp = pd.read_csv(cells_comp_path, low_memory=False)
    if "cell_id" not in cells_comp.columns:
        for cand in ("barcode", "cell", "obs_names", "Cell"):
            if cand in cells_comp.columns:
                cells_comp = cells_comp.rename(columns={cand: "cell_id"})
                break
    cells_comp_ids = (cells_comp["cell_id"].values if "cell_id" in cells_comp.columns
                       else cells_comp.iloc[:, 0].values)
    log.info("compartment metadata: %s (%d rows)", cells_comp_path.name,
             len(cells_comp_ids))
    if X.shape[0] != len(cells_comp_ids):
        log.warning("counts rows (%d) != cells_comp (%d)",
                    X.shape[0], len(cells_comp_ids))

    gene_data = pd.read_csv(a.counts_root / "gene_data.csv", index_col=0)
    if "symbol" not in gene_data.columns:
        for cand in ("gene_symbol", "Symbol", "symbol_"):
            if cand in gene_data.columns:
                gene_data = gene_data.rename(columns={cand: "symbol"})
                break
    log.info("gene_data: %d entries", len(gene_data))

    # Cells in target L2 (artifact-filtered)
    labels = pd.read_csv(a.labels_full,
                          usecols=["cell_id", "compartment", "label", "is_artifact"])
    labels.is_artifact = labels.is_artifact.astype(str).str.lower().isin(
        {"true", "1"})
    l2_cells = set(labels[(labels.compartment == L2_compartment) &
                          (labels.label == L2_label) &
                          (~labels.is_artifact)].cell_id)
    log.info("Cells in %s per labels_full: %d", a.L2, len(l2_cells))
    # pd.Series.isin uses hash-set (O(N+M)); np.isin on 1M x 150K = O(N*M) -> 40 min
    in_l2_mask = pd.Series(cells_comp_ids).isin(l2_cells).values
    l2_idx_in_counts = np.where(in_l2_mask)[0]
    l2_cells_ordered = cells_comp_ids[l2_idx_in_counts]
    log.info("L2 cells in counts: %d", len(l2_idx_in_counts))
    if len(l2_idx_in_counts) == 0:
        sys.exit("No L2 cells in counts barcodes")

    atlas_cell_to_idx = {c: i for i, c in enumerate(cells_atlas)}
    l2_atlas_idx = np.array([atlas_cell_to_idx.get(c, -1) for c in l2_cells_ordered])
    in_atlas = l2_atlas_idx >= 0
    log.info("L2 cells in atlas nhoods: %d / %d", in_atlas.sum(), len(in_atlas))
    keep_l2 = np.where(in_atlas)[0]
    l2_atlas_idx = l2_atlas_idx[in_atlas]

    # nhoods.mtx columns assumed 0-indexed; nhood_to_ng keys are 1-indexed nhood IDs
    nh_keep_cols = sorted([n - 1 for n in nhood_to_ng.keys()
                            if 0 <= (n - 1) < nh_mat.shape[1]])
    log.info("nhood columns of interest: %d", len(nh_keep_cols))
    nh_sub = nh_mat[l2_atlas_idx, :][:, nh_keep_cols]
    ng_per_col = np.array([nhood_to_ng[c + 1] for c in nh_keep_cols])

    log.info("Assigning dominant NG per cell ...")
    dominant_ng = np.array(["__none__"] * nh_sub.shape[0], dtype=object)
    nh_sub_lil = nh_sub.tolil()
    for i in range(nh_sub.shape[0]):
        cols = nh_sub_lil.rows[i]
        if not cols:
            continue
        ngs = ng_per_col[cols]
        c = Counter(ngs)
        dominant_ng[i] = c.most_common(1)[0][0]
    has_ng = dominant_ng != "__none__"
    log.info("Cells with viable NG assignment: %d / %d", has_ng.sum(), len(has_ng))
    keep_l2 = keep_l2[has_ng]
    dominant_ng = dominant_ng[has_ng]
    final_l2_idx = l2_idx_in_counts[keep_l2]
    final_cell_ids = l2_cells_ordered[keep_l2]

    # Optional per-NG downsample
    if a.max_cells_per_ng > 0:
        rng = np.random.default_rng(42)
        keep_idx = []
        keep_ng = []
        keep_cells = []
        for ng in lookup.NhoodGroup:
            ix = np.where(dominant_ng == ng)[0]
            if len(ix) > a.max_cells_per_ng:
                ix = rng.choice(ix, a.max_cells_per_ng, replace=False)
            keep_idx.extend(final_l2_idx[ix].tolist())
            keep_ng.extend([ng] * len(ix))
            keep_cells.extend(final_cell_ids[ix].tolist())
        final_l2_idx = np.array(keep_idx, dtype=int)
        dominant_ng = np.array(keep_ng)
        final_cell_ids = np.array(keep_cells)
        log.info("Downsampled cells: %d", len(final_l2_idx))

    # ---------- Gene column subset (ENSG -> counts column) ----------
    ens_to_col = {e: i for i, e in enumerate(gene_data.index)}
    gene_set["col_idx"] = gene_set["gene_id"].map(ens_to_col)
    missing = gene_set["col_idx"].isna().sum()
    if missing > 0:
        log.warning("ENSG missing in gene_data: %d (dropping)", missing)
    gene_set = gene_set.dropna(subset=["col_idx"]).copy()
    gene_set["col_idx"] = gene_set["col_idx"].astype(int)
    log.info("Genes after ENSG resolution: %d", len(gene_set))

    # Slice counts: cells x genes (oversampled candidates)
    log.info("Slicing counts: %d cells x %d candidate genes ...",
             len(final_l2_idx), len(gene_set))
    sub = X[final_l2_idx, :][:, gene_set.col_idx.values].tocsr()
    log.info("candidate counts: %s, nnz=%d", sub.shape, sub.nnz)

    # ---------- pct-in filter: keep gene only if expressed in >= min-pct-in
    #            of cells in its source NG, then truncate to top-N per NG ----------
    log.info("Computing pct-in per (gene x source NG) ...")
    sub_bool = (sub > 0).astype(np.int8)  # cells x genes
    pct_in = np.zeros(len(gene_set), dtype=float)
    gene_set = gene_set.reset_index(drop=True)
    for src_ng, idx in gene_set.groupby("source_NG").groups.items():
        cell_mask = (dominant_ng == src_ng)
        n_cells = int(cell_mask.sum())
        if n_cells == 0:
            continue
        # cells in source NG, columns = these genes
        col_idx = list(idx)
        block = sub_bool[cell_mask, :][:, col_idx]
        pct_in[col_idx] = np.asarray(block.sum(axis=0)).ravel() / n_cells
    gene_set["pct_in_source_NG"] = pct_in
    n_pre = len(gene_set)
    gene_set_kept = gene_set[gene_set.pct_in_source_NG >= a.min_pct_in].copy()
    log.info("Genes passing pct-in >= %.2f: %d / %d",
             a.min_pct_in, len(gene_set_kept), n_pre)

    # Top-N per source_NG by score among survivors
    gene_set_kept = (gene_set_kept.sort_values("score", ascending=False)
                                    .groupby("source_NG", sort=False)
                                    .head(a.top_genes_per_ng)
                                    .reset_index(drop=True))
    log.info("Final genes (top-%d per NG, dedup): %d",
             a.top_genes_per_ng, len(gene_set_kept))

    # Re-slice counts to final gene set (subset of cols of `sub`)
    keep_col_local = gene_set_kept.index.tolist()  # post-reset_index
    # Map from gene_id back to position in original sub matrix
    gene_id_to_subcol = {gid: i for i, gid in enumerate(gene_set.gene_id)}
    keep_subcols = [gene_id_to_subcol[g] for g in gene_set_kept.gene_id]
    sub = sub[:, keep_subcols]
    gene_set = gene_set_kept
    log.info("final counts: %s, nnz=%d", sub.shape, sub.nnz)

    # ---------- Cell metadata: surface donor + parity + risk_class ----------
    cells_meta_cols = [c for c in
                        ["donor", "parity_binary", "parity", "risk_class",
                         "brca_genotype", "study", "facs_status", "age_binary",
                         "menopausal_status_binary", "cancer_history"]
                        if c in cells_comp.columns]
    cells_full = cells_comp.iloc[final_l2_idx][["cell_id"] + cells_meta_cols].copy()
    cells_full["NG_dominant"] = dominant_ng
    cells_full = cells_full.merge(
        lookup[["NhoodGroup", "size_rank_within_L2", "name_full",
                "group_med_lfc"]].rename(columns={"NhoodGroup": "NG_dominant",
                                                  "size_rank_within_L2": "NG_rank",
                                                  "name_full": "NG_name_full",
                                                  "group_med_lfc": "NG_group_med_lfc"}),
        on="NG_dominant", how="left")

    # ---------- Write outputs ----------
    out_dir = inq / "outputs/panelC_features" / f"{L2_safe}_{a.contrast}"
    out_dir.mkdir(parents=True, exist_ok=True)

    # counts.mtx.gz
    counts_path = out_dir / "counts.mtx"
    sio.mmwrite(str(counts_path), sub, comment=f"cells={sub.shape[0]} genes={sub.shape[1]}")
    with open(counts_path, "rb") as f_in, gzip.open(str(counts_path) + ".gz", "wb") as f_out:
        f_out.writelines(f_in)
    counts_path.unlink()
    log.info("Wrote counts.mtx.gz")

    # cells.tsv (with metadata)
    cells_full.to_csv(out_dir / "cells.tsv", sep="\t", index=False)
    log.info("Wrote cells.tsv (%d rows)", len(cells_full))

    # genes.tsv
    genes_out = gene_set[["gene_id", "symbol", "source_NG", "source_rank",
                           "source_NG_full", "score", "logFC", "pct_in_source_NG",
                           "n_pathways", "top_pathway", "top_NES"]].copy()
    genes_out["row_idx"] = range(len(genes_out))  # row in counts.mtx.gz
    genes_out.to_csv(out_dir / "genes.tsv", sep="\t", index=False)
    log.info("Wrote genes.tsv (%d genes)", len(genes_out))

    # ng_summary.csv
    cells_per_ng = pd.Series(dominant_ng).value_counts().rename("n_cells")
    pathway_records_df = pd.DataFrame(pathway_records)
    n_sig_per_ng = (pathway_records_df.groupby("NhoodGroup").size()
                     .rename("n_sig_pathways")
                     if not pathway_records_df.empty
                     else pd.Series(dtype=int, name="n_sig_pathways"))
    ng_summary = lookup[["NhoodGroup", "name_full", "size_rank_within_L2",
                          "n_nhoods_in_group", "group_med_lfc",
                          "group_pct_up", "n_sig"]].rename(
        columns={"size_rank_within_L2": "rank"})
    ng_summary["n_cells"] = ng_summary["NhoodGroup"].map(cells_per_ng).fillna(0).astype(int)
    ng_summary["n_sig_pathways"] = ng_summary["NhoodGroup"].map(n_sig_per_ng).fillna(0).astype(int)
    ng_summary.to_csv(out_dir / "ng_summary.csv", index=False)

    # pathway sidecar (for plot scripts that want to label by pathway)
    if pathway_records:
        pd.DataFrame(pathway_records).sort_values(
            ["rank", "NES"], ascending=[True, False]).to_csv(
            out_dir / "pathways.csv", index=False)

    # manifest
    manifest = {
        "L2": a.L2,
        "L2_safe": L2_safe,
        "contrast": a.contrast,
        "n_viable_NGs": int(len(lookup)),
        "n_genes_picked": int(len(gene_set)),
        "n_cells_total": int(len(final_l2_idx)),
        "gene_selection": {
            "rule": "leading-edge x F.3, score = NES * logFC, per-pathway diversified, pct-in filtered",
            "top_genes_per_ng": a.top_genes_per_ng,
            "top_genes_per_pathway": a.top_genes_per_pathway,
            "candidate_multiplier": a.candidate_multiplier,
            "min_pct_in_source_NG": a.min_pct_in,
            "gsea_padj": a.gsea_padj,
            "gsea_NES_filter": "NES > 0",
            "marker_fdr": a.marker_fdr,
            "marker_logFC_filter": "logFC > 0",
            "marker_dissoc_filter": "is_dissoc_stress excluded",
        },
        "max_cells_per_ng": a.max_cells_per_ng,
        "inputs": {
            "lookup": str(inq / "outputs/nhoodgroup_renaming/lookup.csv"),
            "gsea": str(gsea_path),
            "markers_dir": str(marker_dir),
            "atlas_nhoods": str(nhoods_dir),
            "nhood_groups": str(ng_csv),
            "counts_npz": str(counts_npz),
            "compartment_metadata": str(cells_comp_path),
            "labels_full": a.labels_full,
        },
    }
    with open(out_dir / "feature_manifest.yaml", "w") as f:
        # YAML-ish via JSON (PyYAML may not be in container)
        f.write(json.dumps(manifest, indent=2, default=str))
    log.info("Wrote feature_manifest.yaml")
    log.info("DONE: %s", out_dir)


if __name__ == "__main__":
    main()
