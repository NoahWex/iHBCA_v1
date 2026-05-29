#!/usr/bin/env python3
"""s_S25b — Composition GLM with epi block + M0+M1+M2 pooled baseline.

Variant of build_E1_comp_glm_with_epi.py where the baseline is the
3-motif stromal/connective-tissue pool (M0 activated-fibroblast, M1
endothelial, M2 fibroblast) rather than the 2-motif M1+M2 pool.

Math (REF = M2 dropped in NB-GLM):
  contrast c = beta_m − (1/3)·(beta_M0 + beta_M1)   for m ∈ {M3..M8}
            =  (2/3)·beta_M0 − (1/3)·beta_M1        for m = M0
            = −(1/3)·beta_M0 + (2/3)·beta_M1        for m = M1
            = −(1/3)·beta_M0 − (1/3)·beta_M1        for m = M2
  Var(c) tracked exactly (includes Cov(beta_M0, beta_M1)).

Outputs (long CSV) + heatmap with epi block on top, vertical separator
after M2 (i.e. baseline columns are M0/M1/M2), BH per row across the 6
non-baseline motifs, ***/**/* significance stars.
"""
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import norm
from scipy.spatial import cKDTree
import statsmodels.api as sm
from statsmodels.formula.api import glm
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

SHARE = ("/Users/NoahWechter/Library/Application Support/CRSP Desktop/"
         "Volumes.noindex/CRSP Lab - dalawson.localized/nwechter")
PREVIEW = (f"{SHARE}/iHBCA_publication/coordination/staging/"
           "fig3_promotion_preview_20260520")
MC_PATH = f"{PREVIEW}/analysis/spatial/niche/outputs/motif_cells.parquet"
MU_PATH = f"{PREVIEW}/analysis/spatial/niche/outputs/motif_units.parquet"
JL_PATH = "/tmp/joint_l1p5.csv"
XO_PATH = "/tmp/xenium_obs.csv"

OUT_DIR = f"/tmp/s_S25b_comp_glm_M0M1M2"

EPI_L1P5 = ["BMYO-myo", "LASP-basal", "LASP", "LHS"]
MOTIFS = [f"M{i}" for i in range(9)]
REF = "M2"
BASELINE = ["M0", "M1", "M2"]
NONBASE = [m for m in MOTIFS if m not in BASELINE]


def build_epi_unit_assignments():
    print("[epi] loading joint_l1p5 + xenium_obs ...")
    jl = pd.read_csv(JL_PATH, usecols=[
        "cell_id", "platform", "patient_id", "library_id", "l1p5_short"])
    jl = jl[(jl["platform"] == "xenium") &
            (jl["l1p5_short"].isin(EPI_L1P5))].copy()
    jl["sample_id"] = jl["cell_id"].str.split("__").str[0]
    xo = pd.read_csv(XO_PATH, usecols=[
        "cell_id", "x_centroid", "y_centroid", "qc_pass_nuclear"])
    xo = xo[xo["qc_pass_nuclear"] == True].copy()  # noqa: E712
    epi = jl.merge(xo[["cell_id", "x_centroid", "y_centroid"]],
                   on="cell_id", how="inner")
    print(f"[epi] {len(epi)} epi cells across {epi['sample_id'].nunique()} samples")
    mu = pd.read_parquet(MU_PATH)
    rows = []
    for sid, gepi in epi.groupby("sample_id"):
        umu = mu[mu["sample_id"] == sid]
        if len(umu) == 0:
            continue
        tree = cKDTree(umu[["centroid_x", "centroid_y"]].values)
        d, idx = tree.query(gepi[["x_centroid", "y_centroid"]].values, k=1)
        umu_arr = umu.reset_index(drop=True)
        rows.append(pd.DataFrame({
            "cell_id": gepi["cell_id"].values,
            "sample_id": sid,
            "patient_id": gepi["patient_id"].values,
            "l1p5_short": gepi["l1p5_short"].values,
            "argmax_P": umu_arr.iloc[idx]["motif"].values,
            "unit_id":  umu_arr.iloc[idx]["unit_id"].values,
        }))
    return pd.concat(rows, ignore_index=True)


def run_glm(MC, keep_l1p5):
    agg = (MC.groupby(["unit_id", "sample_id", "patient_id",
                       "argmax_P", "l1p5_short"])
             .size().rename("n_cells").reset_index())
    tot = MC.groupby("unit_id").size().rename("n_total").reset_index()
    units = MC[["unit_id", "sample_id", "patient_id",
                "argmax_P"]].drop_duplicates()
    ue = units.merge(pd.DataFrame({"l1p5_short": keep_l1p5}), how="cross")
    ue = ue.merge(agg, how="left",
                  on=["unit_id","sample_id","patient_id","argmax_P","l1p5_short"])
    ue["n_cells"] = ue["n_cells"].fillna(0).astype(int)
    ue = ue.merge(tot, on="unit_id")
    ue["log_n_total"] = np.log(ue["n_total"].clip(lower=1))
    motif_order = [REF] + [m for m in MOTIFS if m != REF]
    ue["argmax_P"] = pd.Categorical(ue["argmax_P"], categories=motif_order)
    ue["patient_id"] = ue["patient_id"].astype("category")

    log2 = pd.DataFrame(np.nan, index=keep_l1p5, columns=MOTIFS)
    lcb  = pd.DataFrame(np.nan, index=keep_l1p5, columns=MOTIFS)
    ucb  = pd.DataFrame(np.nan, index=keep_l1p5, columns=MOTIFS)
    pval = pd.DataFrame(np.nan, index=keep_l1p5, columns=MOTIFS)

    def coef_key(m): return f"C(argmax_P)[T.{m}]"
    LN2 = np.log(2)

    for l in keep_l1p5:
        sub = ue[ue["l1p5_short"] == l].copy()
        try:
            mod = glm("n_cells ~ C(patient_id) + C(argmax_P)",
                      data=sub,
                      family=sm.families.NegativeBinomial(alpha=1.0),
                      offset=sub["log_n_total"]).fit(disp=0)
        except Exception:
            continue
        cov = mod.cov_params()
        # Extract beta_M0, beta_M1 (M2 is reference → beta_M2 = 0)
        b0 = mod.params.get(coef_key("M0"), 0.0)
        b1 = mod.params.get(coef_key("M1"), 0.0)
        v_00 = cov.loc[coef_key("M0"), coef_key("M0")] if coef_key("M0") in cov.index else 0.0
        v_11 = cov.loc[coef_key("M1"), coef_key("M1")] if coef_key("M1") in cov.index else 0.0
        v_01 = (cov.loc[coef_key("M0"), coef_key("M1")]
                if coef_key("M0") in cov.index and coef_key("M1") in cov.index else 0.0)
        for m in MOTIFS:
            if m == "M0":
                # c = (2/3) b0 - (1/3) b1
                c  = (2.0/3.0)*b0 - (1.0/3.0)*b1
                v  = (4.0/9.0)*v_00 + (1.0/9.0)*v_11 - (4.0/9.0)*v_01
            elif m == "M1":
                # c = -(1/3) b0 + (2/3) b1
                c  = -(1.0/3.0)*b0 + (2.0/3.0)*b1
                v  = (1.0/9.0)*v_00 + (4.0/9.0)*v_11 - (4.0/9.0)*v_01
            elif m == "M2":
                # c = -(1/3) b0 - (1/3) b1
                c  = -(1.0/3.0)*b0 - (1.0/3.0)*b1
                v  = (1.0/9.0)*v_00 + (1.0/9.0)*v_11 + (2.0/9.0)*v_01
            else:
                # m ∈ {M3..M8}: c = beta_m - (1/3)(b0 + b1)
                k_m = coef_key(m)
                if k_m not in mod.params.index:
                    continue
                b_m = mod.params[k_m]
                v_mm = cov.loc[k_m, k_m]
                v_m0 = cov.loc[k_m, coef_key("M0")] if coef_key("M0") in cov.index else 0.0
                v_m1 = cov.loc[k_m, coef_key("M1")] if coef_key("M1") in cov.index else 0.0
                c = b_m - (1.0/3.0)*(b0 + b1)
                v = (v_mm
                     + (1.0/9.0)*v_00
                     + (1.0/9.0)*v_11
                     + (2.0/9.0)*v_01
                     - (2.0/3.0)*v_m0
                     - (2.0/3.0)*v_m1)
            if v <= 0:
                continue
            se = np.sqrt(v)
            log2.loc[l, m] = c / LN2
            lcb.loc[l, m]  = (c - 1.96*se) / LN2
            ucb.loc[l, m]  = (c + 1.96*se) / LN2
            pval.loc[l, m] = 2 * (1 - norm.cdf(abs(c/se)))
    return log2, lcb, ucb, pval


def main():
    Path(OUT_DIR).mkdir(parents=True, exist_ok=True)

    print("[load] motif_cells (non-epi) ...")
    mc = pd.read_parquet(MC_PATH)
    mc = mc[mc["unit_id"] != ""].copy()
    mc["is_epi"] = False
    mc_cols = ["cell_id","sample_id","patient_id","l1p5_short",
               "argmax_P","unit_id","is_epi"]
    mc_sub = mc[mc_cols].copy()

    print("[build] epi block via nearest-motif-unit ...")
    epi = build_epi_unit_assignments()
    epi["is_epi"] = True
    epi_sub = epi[["cell_id","sample_id","patient_id","l1p5_short",
                   "argmax_P","unit_id","is_epi"]].copy()

    MC_full = pd.concat([mc_sub, epi_sub], ignore_index=True)
    print(f"[merge] non-epi: {len(mc_sub)}; epi: {len(epi_sub)}; total: {len(MC_full)}")

    counts = (MC_full.groupby("l1p5_short").size().rename("n").reset_index())
    keep_l1p5 = sorted(counts[counts["n"] >= 50]["l1p5_short"].tolist())
    print(f"[filter] keeping {len(keep_l1p5)} L1.5 types")

    print("[glm] NB-GLM per L1.5 (3-motif pool baseline) ...")
    log2, lcb, ucb, pval = run_glm(MC_full, keep_l1p5)

    # BH per row across all 9 motifs (consistent multiple-testing correction;
    # baseline columns also tested vs the pool — c_M0 = (2/3)b0 - (1/3)b1 asks
    # "is M0 distinguishable from the M0+M1+M2 pool", etc.)
    fdr = pval.copy()
    for l in keep_l1p5:
        ps = pval.loc[l, MOTIFS].dropna()
        if len(ps) == 0:
            continue
        _, q, _, _ = multipletests(ps.values, method="fdr_bh")
        for c, qv in zip(ps.index, q):
            fdr.loc[l, c] = qv

    long_rows = []
    epi_set = set(EPI_L1P5)
    for l in keep_l1p5:
        for m in MOTIFS:
            n_in = int(MC_full[(MC_full["l1p5_short"] == l) &
                               (MC_full["argmax_P"] == m)].shape[0])
            long_rows.append({
                "l1p5": l, "motif": m,
                "log2_RR": log2.loc[l, m],
                "log2_lcb95": lcb.loc[l, m],
                "log2_ucb95": ucb.loc[l, m],
                "pval": pval.loc[l, m],
                "fdr_bh": fdr.loc[l, m],
                "is_epi_block": l in epi_set,
                "n_cells_in_motif": n_in,
                "is_baseline_motif": m in BASELINE,
            })
    long_df = pd.DataFrame(long_rows)
    out_csv = f"{OUT_DIR}/s_S25b_long_table.csv"
    long_df.to_csv(out_csv, index=False)
    print(f"wrote {out_csv}")

    # === Heatmap ===
    L1P5_ORDER = (EPI_L1P5 +
                  ["Fb","Fb_Activated","Fb_SFRP4","Adipo",
                   "EC","PV","LEC","Vas-cap",
                   "CD4T","CD8T","Treg","T-NK",
                   "cDC","cDC1","cDC2","pDC",
                   "B","Plas","Mac","Mac_art","Mast","Neu"])
    order = ([x for x in L1P5_ORDER if x in keep_l1p5] +
             [x for x in keep_l1p5 if x not in L1P5_ORDER])
    mat = log2.loc[order, MOTIFS].values
    fmat = fdr.loc[order, MOTIFS].values

    def stars(q):
        if np.isnan(q): return ""
        if q < 0.001: return "***"
        if q < 0.01:  return "**"
        if q < 0.05:  return "*"
        return ""

    VMAX = 3.0
    fig, ax = plt.subplots(figsize=(7.8, 9.5))
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-VMAX, vmax=VMAX, aspect="auto")
    # Vertical separator after M2 (after column index 2 → boundary at 2.5)
    ax.axvline(2.5, color="black", lw=0.8, ls="--", alpha=0.6)
    # Horizontal separator under epi block
    n_epi = len([x for x in order if x in EPI_L1P5])
    if n_epi > 0:
        ax.axhline(n_epi - 0.5, color="black", lw=1.2, alpha=0.7)

    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            s = stars(fmat[i, j])
            if s:
                ax.text(j, i, s, ha="center", va="center",
                        color="black" if abs(mat[i, j]) < 1.5 else "white",
                        fontsize=7 if len(s) < 3 else 6,
                        fontweight="bold")
    col_labels = [f"{m}\n(baseline)" if m in BASELINE else m
                  for m in MOTIFS]
    ax.set_xticks(range(len(MOTIFS)))
    ax.set_xticklabels(col_labels, fontsize=7)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=7)
    ax.set_xlabel("Motif")
    ax.set_ylabel("L1.5 cell type")
    cb = plt.colorbar(im, ax=ax, shrink=0.7)
    cb.set_label("log2(RR vs M0+M1+M2 pool)", fontsize=7)
    cb.ax.tick_params(labelsize=7)
    plt.tight_layout()
    out_pdf = f"{OUT_DIR}/s_S25b_comp_glm_heatmap.pdf"
    fig.savefig(out_pdf, format="pdf", bbox_inches="tight")
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_pdf}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
