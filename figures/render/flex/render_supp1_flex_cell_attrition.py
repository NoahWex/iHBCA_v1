"""
render_supp1_flex_cell_attrition.py — Fig 3 Supp 1, Panel 1b.

Per-sample cell attrition through Spatial HBCA preprocessing filter steps.
For each FLEX sample, a horizontal stacked bar decomposes the baseline cell
count into the categories produced by the master cell tracker
(`central_cell_status.csv`):

    kept            — cell passes all filters (step_01 → step_15)
    removed_mad     — fails step_02 MAD outlier filter
    removed_doublet — passes step_02 but fails step_03 doublet detection
    removed_manifold— passes 02+03 but fails step_10b manifold artifact gate
    removed_contam  — passes 02+03+10b but fails step_15 Milo contamination

Attribution is by FIRST failed step (mutually exclusive segments). Step_01
(UMI threshold) has only 1 failure cohort-wide and is omitted from the legend
to avoid visual clutter; cells failing it are attributed to "removed_mad" by
default.

Reviewer question: "How does each filter step contribute to cell attrition?"
The visible answer: MAD is the dominant first filter; doublets contribute a
second-order layer; manifold + contamination are sample-dependent niche
filters that each remove < 5% in most samples.

Inputs (from Spatial_HBCA/project/01_Preprocessing — NOT yet promoted to
publication/; flagged for promotion as part of Supp 1 package):
    Spatial_HBCA/project/01_Preprocessing/outputs/qc_status/central_cell_status.csv
    columns: cell_id, sample_id, patient_id, position,
             step_01_umi_pass, step_02_mad_pass, step_03_doublet_pass,
             step_10b_manifold_pass, step_15_pass

Outputs:
    supp1_flex_cell_attrition.pdf
    supp1_flex_cell_attrition_data.csv  (per-sample stage counts)

Pattern source:
    publication/figures/render/joint/render_supp7_4_l2s_l1p5_confusion.py
    (Option B framework boilerplate; argparse + load_aesthetics + load_paths)

Framework conformance: Option B per CP_supp7_python_framework_gap.

Substrate gap:
    central_cell_status.csv has no entry in publication/config/paths.yaml.
    This script resolves the dev path via paths['_root']. Coordinator should
    add a token (e.g. flex_central_cell_status) when this substrate is
    promoted (flagged for promotion alongside this panel).
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
log = logging.getLogger("supp1_1b")


STAGE_ORDER = ["kept", "removed_sample_failed_vf", "removed_mad",
               "removed_doublet", "removed_manifold", "removed_contam"]
STAGE_LABELS = {
    "kept": "Kept (final cohort)",
    "removed_sample_failed_vf": "Removed (sample failed VF, step 04)",
    "removed_mad": "Removed (MAD, step 02)",
    "removed_doublet": "Removed (doublet, step 03)",
    "removed_manifold": "Removed (manifold, step 10b)",
    "removed_contam": "Removed (contamination, step 15)",
}
SAMPLE_FAILED_VF = "Pat1_P3_M_UOQ"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--project-root", required=True, type=Path,
                   help="Path to iHBCA_publication root.")
    p.add_argument("--out-dir", required=True, type=Path,
                   help="Output directory for PDF + data CSV.")
    p.add_argument("--test", action="store_true",
                   help="No-op for 1b (full cohort always rendered).")
    return p.parse_args()


def attribute_first_fail(df: pd.DataFrame) -> pd.Series:
    """Assign each cell to the first failed filter step (or 'kept').

    NaN at a downstream step (e.g. step_10b nan when 02/03 already removed
    the cell) is treated as fail-equivalent and the loop short-circuits at
    the first explicit False, so NaN never wins attribution over a real
    earlier failure.
    """
    n = len(df)
    out = np.full(n, "kept", dtype=object)
    mad_fail = ~df["step_02_mad_pass"].astype("boolean").fillna(False).to_numpy()
    doublet_fail = ~df["step_03_doublet_pass"].astype("boolean").fillna(False).to_numpy()
    manifold_fail = ~df["step_10b_manifold_pass"].astype("boolean").fillna(False).to_numpy()
    contam_fail = ~df["step_15_pass"].astype("boolean").fillna(False).to_numpy()
    # sample-failed-VF: cells in dropped sample that passed step 02 MAD;
    # downstream step flags are NaN (sample dropped before step 03 ran).
    s02_pass = df["step_02_mad_pass"].astype("boolean").fillna(False).to_numpy()
    sample_failed_vf = ((df["sample_id"] == SAMPLE_FAILED_VF).to_numpy()
                       & s02_pass)
    out[contam_fail] = "removed_contam"
    out[manifold_fail] = "removed_manifold"
    out[doublet_fail] = "removed_doublet"
    out[mad_fail] = "removed_mad"
    # sample_failed_vf overwrites the spurious doublet/manifold/contam attribution
    # that NaN→False triggers for these cells (their sample was dropped at step 04
    # before doublet detection ran).
    out[sample_failed_vf] = "removed_sample_failed_vf"
    return pd.Series(out, index=df.index, name="stage")


def aggregate_per_sample(df: pd.DataFrame) -> pd.DataFrame:
    """Crosstab cells by sample × first-failed-step. Returns wide table."""
    crosstab = pd.crosstab(df["sample_id"], df["stage"])
    for col in STAGE_ORDER:
        if col not in crosstab.columns:
            crosstab[col] = 0
    crosstab = crosstab[STAGE_ORDER]
    crosstab["n_baseline"] = crosstab.sum(axis=1)
    crosstab["pct_kept"] = 100.0 * crosstab["kept"] / crosstab["n_baseline"]
    pat = df.drop_duplicates("sample_id").set_index("sample_id")["patient_id"]
    crosstab["patient_id"] = pat
    return crosstab.reset_index()


def get_attrition_palette(config: dict) -> dict:
    """Cell-attrition palette via qc_status (6-way kept + 5 removed stages)."""
    qc = get_palette("qc_status", config)
    return {
        "kept": qc["kept"],
        "removed_sample_failed_vf": qc["removed_sample_failed_vf"],
        "removed_mad": qc["removed_initial"],
        "removed_doublet": qc["removed_doublet"],
        "removed_manifold": qc["removed_manifold"],
        "removed_contam": qc["removed_contam"],
    }




def render(crosstab: pd.DataFrame, attr_pal: dict, dims: dict, out_pdf: Path) -> None:
    patients = sorted(crosstab["patient_id"].unique())
    n_pat = len(patients)
    n_cols = 2
    n_rows = int(np.ceil(n_pat / n_cols))

    fig_w = dims["width"]
    fig_h = dims["height"]
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_w, fig_h),
                             sharex=True)
    axes = np.atleast_2d(axes).flatten()

    cohort_max = crosstab["n_baseline"].max()

    for ax, pat in zip(axes, patients):
        sub = (crosstab[crosstab["patient_id"] == pat]
               .sort_values("sample_id").reset_index(drop=True))
        n_samp = len(sub)
        y = list(range(n_samp))

        left = np.zeros(n_samp)
        for stage in STAGE_ORDER:
            vals = sub[stage].to_numpy()
            ax.barh(y, vals, left=left, color=attr_pal[stage], edgecolor="none",
                    label=STAGE_LABELS[stage] if pat == patients[0] else None)
            left = left + vals

        for i, row in sub.iterrows():
            ax.text(row["n_baseline"] * 1.05, i, f"{row['pct_kept']:.0f}%",
                    va="center", ha="left", fontsize=5)

        ax.set_yticks(y)
        ax.set_yticklabels(sub["sample_id"].values, fontsize=5)
        ax.set_title(pat, fontsize=8, loc="left", pad=2)
        ax.set_xlim(0, cohort_max * 1.15)
        ax.invert_yaxis()
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="y", length=0)

    for j in range(n_pat, len(axes)):
        axes[j].set_visible(False)

    bottom_row_idx = (n_rows - 1) * n_cols
    for k in range(bottom_row_idx, min(bottom_row_idx + n_cols, len(axes))):
        if axes[k].get_visible():
            axes[k].set_xlabel("Cells (n)")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5,
               bbox_to_anchor=(0.5, -0.02), frameon=False, fontsize=6)

    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(out_pdf, format="pdf", dpi=1200, bbox_inches="tight")
    log.info("Wrote %s", out_pdf)


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    paths = load_paths(config_dir=str(args.project_root / "publication" / "config"))
    tracker_path = Path(resolve_path("flex_central_cell_status", paths))
    log.info("Reading master cell tracker: %s", tracker_path)
    # Cols 6 (step_02_exclusion_flags) and 8 (step_03_doublet_pass) have mixed dtypes
    df = pd.read_csv(tracker_path, low_memory=False,
                     dtype={"step_02_exclusion_flags": str,
                            "step_10b_removal_reason": str})
    log.info("Loaded %d cells × %d cols across %d samples",
             len(df), df.shape[1], df["sample_id"].nunique())

    df["stage"] = attribute_first_fail(df)
    log.info("First-fail attribution: %s", df["stage"].value_counts().to_dict())

    crosstab = aggregate_per_sample(df)

    aes = load_aesthetics()
    plt.rcParams.update(get_matplotlib_theme(aes))
    attr_pal = get_attrition_palette(aes)
    dims = get_dimensions("beeswarm_compact", aes)

    out_pdf = args.out_dir / "supp1_flex_cell_attrition.pdf"
    out_csv = args.out_dir / "supp1_flex_cell_attrition_data.csv"

    render(crosstab, attr_pal, dims, out_pdf)

    crosstab.to_csv(out_csv, index=False)
    log.info("Wrote %s (%d rows)", out_csv, len(crosstab))


if __name__ == "__main__":
    main()
