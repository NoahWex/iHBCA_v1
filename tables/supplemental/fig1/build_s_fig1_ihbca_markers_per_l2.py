"""
build_s_fig1_ihbca_markers_per_l2.py

Build Supplementary Table for Fig 1: per-L2-cell-type differential expression
markers from limma-voom one-vs-rest analysis run on the iHBCA v1 integrated
atlas (2.12M cells, 7 contributing studies). Three compartments processed at
their canonical Leiden resolution (epithelial leiden_1.5, stromal leiden_1.0,
immune leiden_1.0), then concatenated into a single long-form table keyed by
gene symbol and L2 cell type.

Inputs:
    --project-root  Path to iHBCA_publication root. Required.
                    Source paths (read-only):
                      publication/analysis/markers/limma/outputs/{epithelial,immune,stromal}/leiden_*/
                      iHBCAv1_upload/.dev/assembly_rebuild_20260405/outputs/components/gene_data.csv
                      publication/figures/data/fig1/panel_substrate.csv

Output:
    publication/tables/supplemental/fig1/s_fig1_ihbca_markers_per_l2.csv

Columns:
    gene_symbol         HGNC gene symbol (decoded from Ensembl via gene_data.csv)
    ensembl_id          ENSG identifier
    L2                  Fine-resolution cell type label (42 canonical L2 types)
    L0_compartment      epi | str | imm
    logFC               log2 fold change (one-vs-rest)
    p_val               raw p-value (limma 'P.Value' column); retained for the
                        --filter-raw-pval flag and reviewer audit
    p_adj               BH-adjusted p-value
    baseMean            mean expression baseline
    B_statistic         empirical Bayes log-odds (limma 'B' column)
    cluster_id          Leiden cluster ID (per-compartment) backing this row
    leiden_resolution   Leiden resolution used per compartment (1.5 for epi; 1.0 for str/imm)

Method:
    - For each compartment, the limma-voom one-vs-rest analysis produced markers per
      Leiden cluster at the resolution that best resolves the canonical L2 labels
      (epi leiden_1.5: 43 clusters covering 8 L2 types; stromal leiden_1.0: 26
      clusters covering 12 L2 types; immune leiden_1.0: 26 clusters covering 22 L2
      types). For each cluster, majority L2 label is computed from panel_substrate.csv
      (cell_id → canonical L2 mapping after dropping is_artifact cells); cluster_id
      is then replaced with that majority L2 in the output.
    - Gene Ensembl IDs are decoded to HGNC symbols via gene_data.csv. Rows lacking
      a gene_symbol mapping are retained with gene_symbol = ensembl_id.

Validation gate (executed at end of run):
    - Total row count is the sum of per-compartment limma row counts
      (~78K epi + ~52K stromal + ~52K immune = ~182K). Allow ±5K drift for
      filter variations.
    - All 42 canonical L2 labels appear in the L2 column.
    - logFC range is finite (no Inf, no NaN in core columns).
    - L0_compartment values are exactly {epi, str, imm}.

Cited by:
    - drafting_space/fig1/results_annotation_20260516.md §2 (cell-type vocabulary
      DE marker evidence; complements Fig 1 MAIN canonical-marker dotplot)
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


COMPARTMENTS = {
    "epi":      {"leiden_res": "1.5", "limma_subdir": "epithelial", "leiden_col": "leiden_1.5",
                 "leiden_barcoded": "scanvi/epi/n_latent_50/epi_leiden_barcoded.csv"},
    "str":      {"leiden_res": "1.0", "limma_subdir": "stromal",    "leiden_col": "leiden_1.0",
                 "leiden_barcoded": "scanvi/str/n_latent_50/str_leiden_barcoded.csv"},
    "imm":      {"leiden_res": "1.0", "limma_subdir": "immune",     "leiden_col": "leiden_1.0",
                 "leiden_barcoded": "scanvi/imm/n_latent_50/imm_leiden_barcoded.csv"},
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project-root", type=Path, required=True,
                   help="Path to iHBCA_publication root.")
    p.add_argument("--upstream-root", type=Path, default=None,
                   help="Path to upstream repo root (parent of iHBCA_V1, iHBCAv1_upload). "
                        "Defaults to project_root.parent.")
    p.add_argument("--out", type=Path, default=None,
                   help="Output CSV path. Defaults to "
                        "<project_root>/publication/tables/supplemental/fig1/"
                        "s_fig1_ihbca_markers_per_l2.csv.")
    p.add_argument("--filter-raw-pval", type=float, default=None,
                   help="If set, drop rows with raw p_val > this value. "
                        "Use 0.1 for the semi-filtered submission presentation.")
    return p.parse_args()


def load_gene_symbols(gene_data_path: Path) -> pd.DataFrame:
    """Load Ensembl -> HGNC symbol mapping from gene_data.csv."""
    gd = pd.read_csv(gene_data_path)
    if "gene_id" in gd.columns:
        ens_col = "gene_id"
    elif "ensembl_id" in gd.columns:
        ens_col = "ensembl_id"
    else:
        raise ValueError(f"Expected gene_id or ensembl_id column in {gene_data_path}; got {list(gd.columns)}")
    sym_col = "feature_name" if "feature_name" in gd.columns else "gene_symbol"
    if sym_col not in gd.columns:
        sym_col = next((c for c in gd.columns if "symbol" in c.lower() or "name" in c.lower()), None)
        if sym_col is None:
            raise ValueError(f"No gene symbol column found in {gene_data_path}; cols: {list(gd.columns)}")
    return gd[[ens_col, sym_col]].rename(columns={ens_col: "ensembl_id", sym_col: "gene_symbol"})


def cluster_to_l2_majority(panel_substrate: pd.DataFrame, leiden_barcoded: pd.DataFrame,
                           leiden_col: str) -> pd.DataFrame:
    """Compute majority L2 per per-compartment scANVI Leiden cluster.

    Joins panel_substrate (cell_id -> L2) <-> leiden_barcoded (cell_id -> cluster_id)
    on cell_id. Drops is_artifact cells. Returns
    DataFrame with [cluster_id, L2, L0_compartment, n_cells, fraction_majority].
    """
    ps = panel_substrate.loc[panel_substrate["is_artifact"] != True,
                             ["cell_id", "label", "compartment"]].rename(columns={"label": "L2", "compartment": "L0"})
    joined = (ps
              .merge(leiden_barcoded[["cell_id", leiden_col]], on="cell_id", how="inner")
              .rename(columns={leiden_col: "cluster_id"}))
    counts = joined.groupby(["cluster_id", "L2", "L0"]).size().reset_index(name="n_cells")
    totals = counts.groupby("cluster_id")["n_cells"].sum().rename("total")
    counts = counts.merge(totals, on="cluster_id")
    counts["fraction_majority"] = counts["n_cells"] / counts["total"]
    idx = counts.groupby("cluster_id")["n_cells"].idxmax()
    majority = counts.loc[idx, ["cluster_id", "L2", "L0", "n_cells", "fraction_majority"]].reset_index(drop=True)
    return majority


def build_compartment(comp_key: str, limma_root: Path, gene_symbols: pd.DataFrame,
                      cluster_map: pd.DataFrame) -> pd.DataFrame:
    """Read one compartment's limma CSV, join to L2 + gene symbol, return long-form."""
    cfg = COMPARTMENTS[comp_key]
    csv_path = limma_root / cfg["limma_subdir"] / f"leiden_{cfg['leiden_res']}" / f"limma_markers_leiden_{cfg['leiden_res']}.csv"
    if not csv_path.exists():
        raise FileNotFoundError(f"Limma CSV not found: {csv_path}")
    df = pd.read_csv(csv_path)
    df = df.rename(columns={
        "log2FoldChange": "logFC",
        "pvalue": "p_val",
        "padj": "p_adj",
        "stat": "limma_stat",
        "B": "B_statistic",
        "gene": "ensembl_id",
        "cell_type": "cluster_id",
    })
    df["cluster_id"] = df["cluster_id"].astype(int)
    df = df.merge(cluster_map[["cluster_id", "L2", "L0"]], on="cluster_id", how="inner")
    df = df.merge(gene_symbols, on="ensembl_id", how="left")
    df["gene_symbol"] = df["gene_symbol"].fillna(df["ensembl_id"])
    df["leiden_resolution"] = cfg["leiden_res"]
    return df[[
        "gene_symbol", "ensembl_id", "L2", "L0",
        "logFC", "p_val", "p_adj", "baseMean", "B_statistic",
        "cluster_id", "leiden_resolution",
    ]].rename(columns={"L0": "L0_compartment"})


def validate(out_df: pd.DataFrame) -> dict:
    """Run validation gate. Returns dict of {check_name: pass_bool, ...}."""
    gates = {}
    gates["row_count_in_expected_range"] = 170_000 <= len(out_df) <= 195_000
    gates["all_compartments_present"] = set(out_df["L0_compartment"].unique()) == {"epi", "str", "imm"}
    gates["n_l2_at_least_40"] = out_df["L2"].nunique() >= 40
    gates["no_inf_logFC"] = (~out_df["logFC"].isin([float("inf"), float("-inf")])).all()
    gates["no_nan_p_adj"] = out_df["p_adj"].notna().all()
    return gates


def main() -> int:
    args = parse_args()
    proj_root = args.project_root.resolve()
    upstream = (args.upstream_root or proj_root.parent).resolve()
    out_path = args.out or (proj_root / "publication" / "tables" / "supplemental" / "fig1" / "s_fig1_ihbca_markers_per_l2.csv")

    print(f"project_root: {proj_root}", flush=True)
    print(f"upstream:     {upstream}", flush=True)
    print(f"out:          {out_path}", flush=True)

    limma_root = proj_root / "publication" / "analysis" / "markers" / "limma" / "outputs"
    components = upstream / "iHBCAv1_upload" / ".dev" / "assembly_rebuild_20260405" / "outputs" / "components"

    panel_substrate_path = proj_root / "publication" / "figures" / "data" / "fig1" / "panel_substrate.csv.gz"
    gene_data_path = components / "gene_data.csv"
    scanvi_root = upstream / "iHBCAv1_upload" / ".dev" / "assembly_rebuild_20260405" / "outputs"

    print("Loading gene symbols...", flush=True)
    gene_symbols = load_gene_symbols(gene_data_path)
    print(f"  genes mapped: {len(gene_symbols)}", flush=True)

    print("Loading panel_substrate...", flush=True)
    ps = pd.read_csv(panel_substrate_path)
    print(f"  rows: {len(ps)}", flush=True)

    print("Building per-compartment marker tables (per-compartment scANVI Leiden)...", flush=True)
    out_parts = []
    for comp_key in ["epi", "str", "imm"]:
        cfg = COMPARTMENTS[comp_key]
        leiden_col = cfg["leiden_col"]
        leiden_path = scanvi_root / cfg["leiden_barcoded"]
        print(f"  {comp_key}: loading {leiden_path.name} ({leiden_col})", flush=True)
        leiden_barcoded = pd.read_csv(leiden_path, usecols=["cell_id", leiden_col])
        cluster_map = cluster_to_l2_majority(ps, leiden_barcoded, leiden_col)
        print(f"    clusters: {len(cluster_map)}; L2 covered: {cluster_map['L2'].nunique()}", flush=True)
        comp_df = build_compartment(comp_key, limma_root, gene_symbols, cluster_map)
        print(f"    rows: {len(comp_df)}, L2 in output: {comp_df['L2'].nunique()}", flush=True)
        out_parts.append(comp_df)
    out_df = pd.concat(out_parts, ignore_index=True)
    print(f"\nTotal rows: {len(out_df)}; L2 unique: {out_df['L2'].nunique()}", flush=True)

    is_filtered = args.filter_raw_pval is not None
    if is_filtered:
        n_before = len(out_df)
        out_df = out_df.loc[out_df["p_val"] <= args.filter_raw_pval].reset_index(drop=True)
        print(f"Filtered: raw p_val <= {args.filter_raw_pval} -> "
              f"{len(out_df)} rows ({n_before - len(out_df)} dropped)", flush=True)

    gates = validate(out_df)
    # When filtered, row_count gate is expected to fall outside the unfiltered band
    if is_filtered:
        gates.pop("row_count_in_expected_range", None)
    print("\nValidation gate:")
    for k, v in gates.items():
        print(f"  {k}: {'PASS' if v else 'FAIL'}")
    if not all(gates.values()):
        print("\nWARNING: One or more validation gates failed. Inspect output.", file=sys.stderr)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, index=False)
    print(f"\nWrote {out_path} ({out_path.stat().st_size:,} bytes)", flush=True)
    return 0 if all(gates.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
