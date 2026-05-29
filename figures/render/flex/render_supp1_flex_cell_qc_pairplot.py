"""
render_supp1_flex_cell_qc_pairplot.py — Fig 3 Supp 1, Panel 1a.

4×4 pair plot of cell-level QC metrics for the FLEX cohort, colored by QC
pass/fail. Shows where in feature space the technical-junk cells live.

Axes:
    nMAD(log10 n_umi)         signed, per-cluster_coarse × per-sample MAD
    nMAD(log10 n_genes_by_counts)  signed, per-cluster_coarse × per-sample MAD
    scDblFinder_score          raw [0,1] doublet probability
    pct_counts_mt              raw mitochondrial percent

Diagonal: KDE per class for each metric.
Lower triangle: scatter (rasterized 1200 DPI) colored by class.
Upper triangle: empty (Pearson r per pair as small text).

Reviewer question: "Did you filter the obvious technical junk?"
The visible answer: QC-fail cells localize to the high-mt% / low-nCount /
high-doublet corner of feature space; the remaining pass cohort forms a
coherent distribution.

Inputs (consumes-from-dev — NOT yet promoted to publication/):
    Spatial_HBCA_preprocessing/outputs/preprocessing/02_CellFiltering/
        cell_metadata/<sample>_metadata.csv
        cols: cell_id (Unnamed: 0), sample_id, n_umi, n_genes_by_counts,
              pct_counts_mt, cluster_coarse, cell_filter_pass_final
    Spatial_HBCA/project/01_Preprocessing/outputs/qc_status/
        central_cell_status.csv
        cols: cell_id, step_02_mad_pass, step_03_doublet_pass,
              step_03_doublet_score

QC pass/fail definition: cell_filter_pass_final & doublet_filter_pass.
MAD computed per (sample_id × cluster_coarse) group; raw doublet/mt% used.

Outputs:
    supp1_flex_cell_qc_pairplot.pdf
    supp1_flex_cell_qc_pairplot_data.csv  (per-cell metrics + pass/fail)

Pattern source:
    publication/figures/render/joint/render_supp7_4_l2s_l1p5_confusion.py
    (Option B framework boilerplate)

Framework conformance: Option B per CP_supp7_python_framework_gap.

Substrate gap:
    02 cell_metadata + central_cell_status.csv have no entries in
    publication/config/paths.yaml. Bundled for promotion as part of Supp 1
    package alongside Panel 1b.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT_DEFAULT = SCRIPT_DIR.parents[3]
sys.path.insert(0, str(PROJECT_ROOT_DEFAULT / "publication" / "config"))
from load_aesthetics import (  # noqa: E402
    get_dimensions,
    get_matplotlib_theme,
    get_palette,
    load_aesthetics,
)
from load_paths import load_paths, resolve_path  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("supp1_1a")


METRICS = ["nmad_log_n_umi", "nmad_log_n_genes", "doublet_score", "pct_mt"]
METRIC_LABELS = {
    "nmad_log_n_umi": "nMAD log10(nCount)",
    "nmad_log_n_genes": "nMAD log10(nFeature)",
    "doublet_score": "Doublet score",
    "pct_mt": "Mitochondrial %",
}
ATTRITION_ORDER = ["kept", "removed_sample_failed_vf",
                   "removed_mad", "removed_doublet"]
ATTRITION_LABELS = {
    "kept": "Kept",
    "removed_sample_failed_vf": "Removed (sample failed VF)",
    "removed_mad": "Removed (MAD)",
    "removed_doublet": "Removed (doublet)",
}
SAMPLE_FAILED_VF = "Pat1_P3_M_UOQ"
SHUFFLE_SEED = 42


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--test", action="store_true",
                   help="Subsample 5 samples for sanity render.")
    return p.parse_args()


def load_cell_metadata(meta_dir: Path, samples: list | None = None) -> pd.DataFrame:
    csvs = sorted(meta_dir.glob("*_metadata.csv"))
    if samples is not None:
        csvs = [f for f in csvs if any(s in f.name for s in samples)]
    log.info("Loading %d sample metadata CSVs from %s", len(csvs), meta_dir)
    frames = []
    for f in csvs:
        df = pd.read_csv(f, index_col=0)
        df.index.name = "cell_id"
        frames.append(df)
    out = pd.concat(frames)
    log.info("Cell metadata: %d cells across %d samples",
             len(out), out["sample_id"].nunique())
    return out


def compute_mad_metrics(df: pd.DataFrame) -> pd.DataFrame:
    """Per (sample_id × cluster_coarse) group, compute signed nMAD on log10
    of n_umi and n_genes_by_counts.

    nMAD = (log10(x) - median_grp(log10 x)) / (1.4826 * MAD_grp(log10 x))
    """
    df = df.copy()
    log_umi = np.log10(df["n_umi"].clip(lower=1))
    log_genes = np.log10(df["n_genes_by_counts"].clip(lower=1))
    df["_log_umi"] = log_umi
    df["_log_genes"] = log_genes

    grp = df.groupby(["sample_id", "cluster_coarse"], observed=True)

    def signed_nmad(s: pd.Series) -> pd.Series:
        med = s.median()
        mad = np.median(np.abs(s - med))
        if mad == 0 or not np.isfinite(mad):
            return pd.Series(np.zeros(len(s)), index=s.index)
        return (s - med) / (1.4826 * mad)

    df["nmad_log_n_umi"] = grp["_log_umi"].transform(signed_nmad)
    df["nmad_log_n_genes"] = grp["_log_genes"].transform(signed_nmad)
    df = df.drop(columns=["_log_umi", "_log_genes"])
    return df


def get_attrition_palette(config: dict) -> dict:
    """Cell-level attrition palette via qc_status."""
    qc = get_palette("qc_status", config)
    return {
        "kept": qc["kept"],
        "removed_sample_failed_vf": qc["removed_sample_failed_vf"],
        "removed_mad": qc["removed_initial"],
        "removed_doublet": qc["removed_doublet"],
    }


def render(df: pd.DataFrame, attr_pal: dict, dims: dict, out_pdf: Path) -> None:
    n = len(METRICS)
    side = dims["width"]
    fig, axes = plt.subplots(n, n, figsize=(side, side))

    rng = np.random.default_rng(SHUFFLE_SEED)
    perm = rng.permutation(len(df))
    df_shuf = df.iloc[perm].reset_index(drop=True)

    color_arr = df_shuf["attrition"].map(attr_pal).to_numpy()

    for i, mi in enumerate(METRICS):
        for j, mj in enumerate(METRICS):
            ax = axes[i, j]
            if i == j:
                bins = 60
                xmin = float(np.nanpercentile(df_shuf[mi], 0.5))
                xmax = float(np.nanpercentile(df_shuf[mi], 99.5))
                for status in ATTRITION_ORDER:
                    vals = (df_shuf.loc[df_shuf["attrition"] == status, mi]
                            .dropna().to_numpy())
                    if len(vals) == 0:
                        continue
                    ax.hist(vals, bins=bins, range=(xmin, xmax),
                            density=True, color=attr_pal[status], alpha=0.55,
                            edgecolor="none")
                ax.set_yticks([])
            elif i > j:
                xvals = df_shuf[mj].to_numpy()
                yvals = df_shuf[mi].to_numpy()
                ax.scatter(xvals, yvals, c=color_arr,
                           s=0.05, linewidths=0, alpha=0.3, rasterized=True)
                ax.set_xlim(np.nanpercentile(xvals, 0.5),
                            np.nanpercentile(xvals, 99.5))
                ax.set_ylim(np.nanpercentile(yvals, 0.5),
                            np.nanpercentile(yvals, 99.5))
            else:
                ax.set_visible(False)
                continue

            if i == n - 1:
                ax.set_xlabel(METRIC_LABELS[mj], fontsize=6)
            else:
                ax.set_xticklabels([])
            if j == 0 and i != j:
                ax.set_ylabel(METRIC_LABELS[mi], fontsize=6)
            elif i != j:
                ax.set_yticklabels([])
            ax.tick_params(labelsize=5)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

    handles = [Patch(color=attr_pal[s], label=ATTRITION_LABELS[s])
               for s in ATTRITION_ORDER]
    fig.legend(handles=handles, loc="upper right", frameon=False, fontsize=6,
               bbox_to_anchor=(0.98, 0.98))

    fig.tight_layout()
    fig.savefig(out_pdf, format="pdf", dpi=1200, bbox_inches="tight")
    log.info("Wrote %s", out_pdf)


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    paths = load_paths(config_dir=str(args.project_root / "publication" / "config"))
    meta_dir = Path(resolve_path("flex_cell_metadata_dir", paths))
    tracker_path = Path(resolve_path("flex_central_cell_status", paths))

    test_samples = ["Pat1_P1", "Pat2_P2_LIQ", "UCI604_P1",
                    "UCI220228_P1", "Pat1_P3_A_UOQ"] if args.test else None
    cell_md = load_cell_metadata(meta_dir, samples=test_samples)

    log.info("Reading master tracker: %s", tracker_path)
    tracker = pd.read_csv(
        tracker_path, low_memory=False,
        usecols=["cell_id", "step_02_mad_pass", "step_03_doublet_pass",
                 "step_03_doublet_score"]).set_index("cell_id")
    log.info("Tracker: %d cells", len(tracker))

    df = cell_md.join(tracker, how="inner")
    log.info("Joined: %d cells", len(df))

    df = compute_mad_metrics(df)
    df = df.rename(columns={"step_03_doublet_score": "doublet_score",
                            "pct_counts_mt": "pct_mt"})

    mad_fail = ~df["cell_filter_pass_final"].astype("boolean").fillna(False).to_numpy()
    doublet_fail = ~df["step_03_doublet_pass"].astype("boolean").fillna(False).to_numpy()
    attrition = np.full(len(df), "kept", dtype=object)
    sample_failed_vf = ((df["sample_id"] == SAMPLE_FAILED_VF).to_numpy()
                       & ~mad_fail & ~doublet_fail)
    attrition[doublet_fail] = "removed_doublet"
    attrition[mad_fail] = "removed_mad"
    attrition[sample_failed_vf] = "removed_sample_failed_vf"
    df["attrition"] = attrition
    log.info("Attrition: %s", df["attrition"].value_counts().to_dict())

    aes = load_aesthetics()
    plt.rcParams.update(get_matplotlib_theme(aes))
    attr_pal = get_attrition_palette(aes)
    dims = get_dimensions("biplot", aes)

    out_pdf = args.out_dir / "supp1_flex_cell_qc_pairplot.pdf"
    out_csv = args.out_dir / "supp1_flex_cell_qc_pairplot_data.csv"

    render(df, attr_pal, dims, out_pdf)

    df_out = df.reset_index()[[
        "cell_id", "sample_id", "patient_id", "cluster_coarse",
        "n_umi", "n_genes_by_counts", "pct_mt", "doublet_score",
        "nmad_log_n_umi", "nmad_log_n_genes", "attrition"
    ]]
    df_out.to_csv(out_csv, index=False)
    log.info("Wrote %s (%d rows)", out_csv, len(df_out))


if __name__ == "__main__":
    main()
