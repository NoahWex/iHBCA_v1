#!/usr/bin/env python3
"""AUCell scoring of preneoplastic-breast iCAF panels alongside tumor panels.

Adds three panels derived from preneoplastic breast scRNA-seq literature:
  - Nee2023_preCAF_inflam: 14-gene inflammatory/secretome subset of the
    top pre-CAF DE signature (Nee et al. 2023 Nat Genet 55:595-606,
    Supplementary Table 6; filtered to FDR<0.05, log2FC>1.0, then curated
    to inflammatory/MMP/chemokine/growth-factor genes only — translation
    machinery and housekeeping excluded).
  - Reed2024_FB2_top30: 26-gene curated top-30 markers of the FB2 cluster
    in the Reed et al. 2024 HBCA (Nat Genet 56:652-662, Supplementary
    Table 4), filtered to remove lncRNAs/uncharacterised loci.
  - MMP_convergent: MMP3, MMP10, MMP12 — three genes appearing in
    Nee preCAF + Pal 2021 stromal-1 + Reed FB2 (the strongest cross-paper
    convergence on a preneoplastic-breast inflammatory fibroblast program).

Also retains:
  - Wu2020 (TNBC iCAF; Wu et al. 2020 studied 5 TNBC samples)
  - Ohlund2017 (pancreatic cancer iCAF)
for tumor-vs-preneoplastic comparison.

PANEL REBUILD (2026-05-23): Galbo2021_pan_iCAF, Galbo2021_pan_myCAF, and
Elyada2019_myCAF were rebuilt from hybrid lit-augmented panels to clean
subsets of their cited supplementary tables (Galbo Supp Table S3 columns
`pan-iCAF` / `pan-myCAF`; Elyada Supp Table S13 myCAF sheet) following the
verification audit at
coordination/reports/fig2_aucell_panel_verification_20260523.md.
Specifically: Elyada myCAF fixes a MMP2→MMP11 transcription error;
Galbo iCAF expands 13→42 genes; Galbo myCAF expands 12→53 genes and removes
Galbo pan-dCAF collagens that were biologically distinct from pan-myCAF.

NOTE ON DUAL SOURCE-OF-TRUTH: The gene-list constants below are mirrored
in publication/analysis/markers/aucell_methods_qc/panels/aucell_panels.yaml
(version 2). The YAML is the documentation source-of-truth; this script
hardcodes the lists for now. Any future panel edit MUST update both files
in lockstep. Backlog item: refactor this script to read the YAML directly
so there is a single source of truth.

Otherwise identical to 11_aucell_module_score.py.
"""

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.io as sio
from scipy import stats
import matplotlib.pyplot as plt


# Tumor-context PDAC iCAF panels (Tuveson-lab curated, distributed by Wu et al. 2020)
# REBUILT 2026-05-23 from Wu 2020 stromal_subclasses GitHub repo:
#   https://github.com/sunnyzwu/stromal_subclasses/blob/master/
#   05_gene_signature_analysis/AUCell/gene_signatures.csv
#   commit 9c819769bd622c71d67bbdf75a4e1c34cdfa36b6 (2019-11-21)
# The Wu lab assembled these PDAC reference panels from Öhlund 2017 (mouse PSC
# DESeq2 data, mouse→human orthologs mapped by Tuveson lab) and Elyada 2019
# (human PDAC iCAF Supp Table S13, ranked by MAST Score). Used verbatim here
# to match the gene lists scored in Wu et al. 2020.
# The prior Wu2020 panel (11 chemokine markers) was removed: those genes
# represent a cross-paper canonical iCAF marker set (Costa 2018 / Pelon 2020 /
# Pal 2021 / Zheng 2022) rather than the Wu 2020 iCAF cluster signature, and
# Wu's own cluster signature (Dataset EV1; 610 iCAF DEGs) is not retrievable
# from the Springer-deposited supplementary materials.

# Öhlund 2017 mouse iCAF (39 genes; mouse symbols mapped to human orthologs by Tuveson lab)
# Source column: Tuveson_iCAF_mouse
OHLUND2017 = ["CLEC3B", "GSN", "PTX3", "LY6C1", "EFEMP1", "LY6A", "FIGF",
              "TNFAIP6", "IFI27L2A", "DPT", "ADM", "PLPP3", "TNXB", "CXCL12",
              "GSTM1", "C4B", "OGN", "C3", "PCOLCE2", "CXCL1", "COL14A1",
              "SVEP1", "ADAMTS5", "HP", "HAS1", "SCARA3", "IFI205", "DPEP1",
              "SFRP4", "PRSS23", "ACKR3", "HTRA3", "APOE", "CCL7", "IL6",
              "SFRP2", "PLA1A", "SCARA5", "SNED1"]

# Preneoplastic-breast panels (from supp tables)
NEE2023_PRECAF = ["CXCL3", "MMP12", "FST", "CXCL2", "MMP3", "CXCL8", "FGF7",
                  "CXCL5", "IL1RL1", "INHBA", "TNFRSF12A", "MMP10", "PHLDA1",
                  "SERPINE2"]
REED2024_FB2 = ["CNIH3", "MMP3", "ARSG", "C11orf53", "GAL", "SDK1", "SLC12A8",
                "THBS2", "PTGES", "WNT2", "RSPO3", "TWIST2", "NMB", "FHAD1",
                "SH3PXD2B", "LIMK1", "ADPRHL1", "MMP12", "TNFAIP6", "LHFPL2",
                "WIPI1", "CLMP", "MEDAG", "SERPINE2", "TAC1", "MAP3K4"]
MMP_CONVERGENT = ["MMP3", "MMP10", "MMP12"]

# Elyada 2019 human PDAC iCAF (84 genes, ranked by MAST Score)
# REBUILT 2026-05-23 from Wu 2020 stromal_subclasses GitHub repo
# (commit 9c819769bd622c71d67bbdf75a4e1c34cdfa36b6), Tuveson_iCAF_human column.
# Independently verified to match the full Elyada 2019 Supp Table S13 iCAF
# sheet (84 genes; cached at coordination/reports/supp_xlsx_cache/elyada2019_supp_table_S13.xlsx).
ELYADA2019_ICAF = ["C3", "DUSP1", "FBLN1", "LMNA", "CLU", "CCDC80", "MYC",
                   "EFEMP1", "HAS1", "NR4A1", "CFD", "ANXA1", "CXCL12", "FGF7",
                   "KLF4", "EMP1", "GPRC5A", "SRPX", "MT2A", "MEDAG", "IGF1",
                   "MGST1", "MCL1", "CEBPD", "S100A10", "UAP1", "TNXB", "CEBPB",
                   "PNRC1", "SOCS3", "PTGDS", "FOSB", "NFKBIA", "CXCL2", "THBS1",
                   "CCL2", "OGN", "GSN", "DPT", "PLA2G2A", "NAMPT", "ITM2A",
                   "RGCC", "JUND", "NNMT", "ZFP36", "PIM1", "CPE", "GFPT2",
                   "SOD2", "KDM6B", "FSTL1", "FBLN2", "NR4A3", "MFAP5", "ABL2",
                   "SGK1", "CILP", "UGDH", "FBLN5", "ADAMTS1", "ADH1B", "WISP2",
                   "GPX3", "S100A4", "IL6", "HAS2", "PLAC9", "IGFBP6", "FBN1",
                   "BDKRB1", "TPPP3", "RASD1", "MT1A", "CXCL14", "PI16", "APOE",
                   "IL8", "ARC", "PTX3", "TNFAIP6", "MT1E", "MT1X", "CXCL1"]
# Elyada myCAF REBUILT 2026-05-23 to the clean Elyada Supp Table S13 myCAF sheet
# (16 genes, full column, rank order by score). Prior 11-gene panel was a hybrid of
# 4 Elyada S13 myCAF genes + 7 canonical lit collagens (COL1A1/COL1A2/COL3A1/COL5A1/COMP/CTGF)
# + a transcription error (MMP2 not in Elyada S13; MMP11 — rank 3 in S13 — is the correct gene).
# Mirrored from publication/analysis/markers/aucell_methods_qc/panels/aucell_panels.yaml v2.
# Audit: coordination/reports/fig2_aucell_panel_verification_20260523.md
ELYADA2019_MYCAF = ["MYL9", "CALD1", "MMP11", "HOPX", "BGN", "IGFBP7", "TPM2",
                    "CTHRC1", "ACTA2", "TAGLN", "INHBA", "COL10A1", "TPM1",
                    "POSTN", "GRP", "CST1"]

# Pan-cancer expanded CAF panels (Galbo et al. 2021 Clin Cancer Res 27:2636)
# Galbo pan-iCAF REBUILT 2026-05-23 to the clean Galbo Supp Table S3 pan-iCAF column
# (42 genes, full column). Prior 13-gene panel was a hybrid: 6/13 from Galbo S3 pan-iCAF +
# 5 canonical iCAF lit markers (IL6, CXCL14, PDGFRB, SOD3, C7) absent from Galbo S3 +
# 2 from other Galbo columns (LUM not in S3; APOD in pan-iCAF-2).
# Mirrored from publication/analysis/markers/aucell_methods_qc/panels/aucell_panels.yaml v2.
GALBO2021_PAN_ICAF = ["CFD", "GPC3", "C3", "ADH1B", "IGF1", "EFEMP1", "PODN",
                      "SEPP1", "CXCL12", "ABI3BP", "FBLN1", "MGST1", "MFAP4",
                      "PLA2G2A", "DPT", "WISP2", "CCDC80", "SFRP2", "PTGDS",
                      "DCN", "MGP", "C1S", "IGFBP6", "GSN", "TMEM176A",
                      "FIBIN", "TMEM176B", "SERPINF1", "FHL1", "GPX3", "CTGF",
                      "C1R", "SFRP4", "CYP1B1", "CST3", "SLC40A1", "FHL2",
                      "ELN", "KLF4", "RARRES1", "CYR61", "IGFBP5"]
# Galbo pan-myCAF REBUILT 2026-05-23 to the clean Galbo S3 pan-myCAF column (53 genes).
# Prior 12-gene panel mixed Galbo pan-myCAF (smooth-muscle contractile) with
# Galbo pan-dCAF (desmoplastic collagen-producing) + non-S3 markers — biologically
# distinct CAF subtypes Galbo himself distinguishes in the pan-CAF analysis.
GALBO2021_PAN_MYCAF = ["ADIRF", "ACTA2", "MYH11", "TAGLN", "SPARCL1", "MCAM",
                       "A2M", "MYLK", "IGFBP7", "CRIP1", "TINAGL1", "TPM2",
                       "PTP4A3", "PPP1R14A", "CRIP2", "ADAMTS1", "CSRP2",
                       "NDUFA4L2", "TPM1", "MAP1B", "FRZB", "PRKCDBP", "CSRP1",
                       "CAV1", "ADAMTS4", "GJA4", "RGS5", "MEF2C", "CALM2",
                       "APOLD1", "OAZ2", "MGST3", "ISYNA1", "CPM", "PGF",
                       "GUCY1B3", "UBA2", "YIF1A", "PHLDA1", "NDRG2", "ID3",
                       "RGS16", "CYB5R3", "CRYAB", "OLFML2A", "TIMP3", "GUCY1A3",
                       "FILIP1", "FAM13C", "NDUFS4", "ITGB1", "KCNE4", "CPE"]

# Senescence panels (Suryadevara et al. 2024 Nat Rev Mol Cell Biol 25:1001-1023,
# DOI 10.1038/s41580-024-00738-8, "SenNet recommendations for detecting senescent cells";
# Supp Table 1 — SenNet biomarker database).
# Author attribution corrected 2026-05-23 from earlier "Suzuki" transcription error;
# the published variable name SENESCENCE_*_SUZUKI2024 retained for backward-compatibility.
# Filtered to: human OR human+mouse organism, mRNA biomolecule, citation count >= 2 (SASP)
# or all reported human-mRNA markers (cell cycle arrest).
SENESCENCE_SASP_SUZUKI2024 = ["TNF", "IL6", "CCL3", "TGFB1", "CXCL1", "CXCL8",
                              "CCL2", "CCL8", "IGFBP4", "IL1A", "IL1B",
                              "SERPINE1", "CSF1", "CSF2", "FAS", "ICAM1",
                              "IGFBP7", "PAPPA", "MMP3", "MMP1", "CXCL2",
                              "CXCL3", "IGFBP2", "IGFBP3", "TIMP2", "CCL4",
                              "CCL5", "CXCL10"]
SENESCENCE_CCA_SUZUKI2024 = ["CDKN1A", "CDKN2A", "TP53", "CDKN1C", "CDKN2D",
                             "SATB1", "CDKN1B", "CDKN2B"]

ALL_PANELS = [
    # Reference anchors (tumor context) — Wu2020 dropped 2026-05-23
    # (no clean Wu cluster signature retrievable; see header docstring)
    ("Ohlund2017",               OHLUND2017),
    ("Elyada2019_iCAF",          ELYADA2019_ICAF),
    ("Elyada2019_myCAF",         ELYADA2019_MYCAF),
    ("Galbo2021_pan_iCAF",       GALBO2021_PAN_ICAF),
    ("Galbo2021_pan_myCAF",      GALBO2021_PAN_MYCAF),
    # Internal reproduction (source studies in iHBCA)
    ("Nee2023_preCAF",           NEE2023_PRECAF),
    ("Reed2024_FB2",             REED2024_FB2),
    ("MMP_convergent",           MMP_CONVERGENT),
    # Senescence (Suryadevara 2024 Nat Rev Mol Cell Biol SenNet biomarker DB)
    ("Senescence_SASP",          SENESCENCE_SASP_SUZUKI2024),
    ("Senescence_CCA",           SENESCENCE_CCA_SUZUKI2024),
]

FIBROS = ["str__Fibro-major", "str__Fibro-IGF1",
          "str__Fibro-SFRP4", "str__Fibro-prematrix"]

TESTED_GENOTYPES = ["negative", "BRCA1"]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--str-npz", required=True)
    p.add_argument("--str-meta", required=True)
    p.add_argument("--labels-full", required=True)
    p.add_argument("--donor-meta", required=True)
    p.add_argument("--gene-data", required=True)
    p.add_argument("--gene-mapping", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--r-script", required=True)
    p.add_argument("--workdir", required=True)
    p.add_argument("--cohort", choices=["tested", "full"], default="tested",
                   help="tested = brca_genotype in {negative, BRCA1} (n=66); "
                        "full = all AR+BR1 donors (n=159). Default tested.")
    p.add_argument("--cell-count-frac", type=float, default=0.10,
                   help="Drop donors with n_cells per (panel, L2) < frac * group mean "
                        "(Nee 2023 rule). Default 0.10.")
    args = p.parse_args()

    out_dir = Path(args.output_dir); out_dir.mkdir(parents=True, exist_ok=True)
    workdir = Path(args.workdir);    workdir.mkdir(parents=True, exist_ok=True)

    print(f"[1] donor metadata + cohort filter (cohort={args.cohort})", flush=True)
    donor_meta = pd.read_csv(args.donor_meta)
    base = donor_meta[donor_meta["risk_class"].isin(["AR", "BR1"])].copy()
    if args.cohort == "tested":
        tested = base[base["brca_genotype"].isin(TESTED_GENOTYPES)].copy()
    else:
        tested = base
    print(f"  donors: total={len(tested)}, "
          f"AR={(tested['risk_class']=='AR').sum()}, "
          f"BR1={(tested['risk_class']=='BR1').sum()}", flush=True)

    print("[2] L2 labels (str only, fibros)", flush=True)
    labels = pd.read_csv(args.labels_full, low_memory=False,
                         usecols=["cell_id", "compartment", "label", "is_artifact"])
    labels = labels[labels["is_artifact"].astype(str).str.lower() != "true"]
    labels = labels[labels["compartment"].astype(str) == "str"].copy()
    labels["L2_label"] = "str__" + labels["label"].astype(str).str.replace(" ", "_", regex=False)
    labels_fibro = labels[labels["L2_label"].isin(FIBROS)].copy()
    print(f"  fibro cells: {len(labels_fibro):,}", flush=True)

    print("[3] str_metadata + cohort filter", flush=True)
    full_str = pd.read_csv(args.str_meta, low_memory=False,
                           usecols=["cell_id", "patientID", "global_numeric_id"])
    full_str["npz_row"] = np.arange(len(full_str))
    tested_donors = set(tested["ihbca_donor_id"])
    str_meta = full_str[full_str["patientID"].isin(tested_donors)].copy()
    str_meta = str_meta.merge(labels_fibro[["cell_id", "L2_label"]],
                              on="cell_id", how="inner")
    # Pull study alongside risk_class so per-donor scores carry study identity for
    # the renderer's shape=study aesthetic + study-adjusted GLM. Drops the need
    # for downstream donor-meta join (and dodges the ihbca_donor_id ↔ patientID
    # format mismatch with the publication s_fig1_ihbca_donor_metadata.csv).
    str_meta = str_meta.merge(
        tested[["ihbca_donor_id", "risk_class", "brca_genotype", "study"]],
        left_on="patientID", right_on="ihbca_donor_id", how="left")
    str_meta["npz_row"] = str_meta["npz_row"].astype(int)
    print(f"  cells: {len(str_meta):,}", flush=True)

    print("[4] gene_data + ENSG->symbol", flush=True)
    gene_data = pd.read_csv(args.gene_data, index_col=0)
    gene_ids = gene_data.index.tolist()
    gm = pd.read_csv(args.gene_mapping, sep="\t")
    gm.columns = [c.lower() for c in gm.columns]
    ensg_col = next(c for c in ["ensembl_id","ensembl_gene_id","gene_id","ensembl"] if c in gm.columns)
    sym_col  = next(c for c in ["gene_symbol","symbol","gene_name"] if c in gm.columns)
    gm = gm.drop_duplicates(subset=[ensg_col]).set_index(ensg_col)
    gene_symbol = [gm.loc[g, sym_col] if g in gm.index else None for g in gene_ids]

    sym_arr = pd.Series(gene_symbol)
    keep_mask = sym_arr.notna() & ~sym_arr.duplicated(keep="first")
    keep_idx = np.where(keep_mask.values)[0]
    kept_symbols = sym_arr.values[keep_idx]
    print(f"  kept genes: {len(keep_idx):,} / {len(gene_ids):,}", flush=True)

    for nm, panel in ALL_PANELS:
        found = [g for g in panel if g in kept_symbols]
        missing = [g for g in panel if g not in kept_symbols]
        print(f"  panel {nm}: found {len(found)}/{len(panel)} missing={missing}", flush=True)

    print("[5] load str_counts.npz", flush=True)
    X_full = sp.load_npz(args.str_npz)
    if not sp.isspmatrix_csr(X_full):
        X_full = X_full.tocsr()
    X = X_full[str_meta["npz_row"].values, :]
    X = X[:, keep_idx]
    print(f"  X subset: {X.shape}", flush=True)

    print("[6] write MTX + cells.tsv + genes.tsv + panels.tsv", flush=True)
    Xt = X.T.tocoo().astype(np.int32)
    mtx_path = workdir / "counts.mtx"
    sio.mmwrite(str(mtx_path), Xt, field="integer")
    pd.Series(str_meta["cell_id"].values).to_csv(workdir / "cells.tsv", index=False, header=False)
    pd.Series(kept_symbols).to_csv(workdir / "genes.tsv", index=False, header=False)
    panel_rows = []
    for nm, panel in ALL_PANELS:
        for g in panel:
            panel_rows.append({"panel": nm, "gene": g})
    pd.DataFrame(panel_rows).to_csv(workdir / "panels.tsv", sep="\t", index=False)

    print("[7] run AUCell in R", flush=True)
    auc_out = workdir / "aucell_scores.csv"
    cmd = ["Rscript", args.r_script,
           "--mtx", str(mtx_path),
           "--cells", str(workdir / "cells.tsv"),
           "--genes", str(workdir / "genes.tsv"),
           "--panels", str(workdir / "panels.tsv"),
           "--out", str(auc_out)]
    r = subprocess.run(cmd, check=False)
    if r.returncode != 0:
        sys.exit(f"R AUCell failed (rc={r.returncode})")

    auc = pd.read_csv(auc_out)
    print(f"  AUC columns: {auc.columns.tolist()}", flush=True)
    str_meta = str_meta.merge(auc, on="cell_id", how="left")

    # Per panel: write outputs
    for panel_name, _ in ALL_PANELS:
        if panel_name not in auc.columns:
            print(f"  WARN: {panel_name} not in AUC; skipping", flush=True)
            continue
        cell_df = str_meta[["cell_id","patientID","study","L2_label","risk_class"]].copy()
        cell_df["AUC"] = str_meta[panel_name].values

        pn_low = panel_name.lower()
        cell_df.to_csv(out_dir / f"{pn_low}_per_cell_scores.csv", index=False)

        # Carry study through the groupby so per_donor_scores.csv has study column.
        per_donor = (cell_df.groupby(["L2_label", "patientID", "study", "risk_class"])
                            .agg(mean_AUC=("AUC", "mean"),
                                 n_cells=("AUC", "size")).reset_index())
        per_donor.to_csv(out_dir / f"{pn_low}_per_donor_scores.csv", index=False)

        rows = []
        for l2 in FIBROS:
            sub = per_donor[per_donor["L2_label"] == l2]
            ar  = sub[sub["risk_class"] == "AR"]["mean_AUC"].values
            br1 = sub[sub["risk_class"] == "BR1"]["mean_AUC"].values
            if len(ar) < 2 or len(br1) < 2:
                rows.append({"L2": l2, "n_AR": len(ar), "n_BR1": len(br1),
                             "p_wilcox": np.nan, "median_AR": np.nan,
                             "median_BR1": np.nan, "delta_med": np.nan})
                continue
            u, p = stats.mannwhitneyu(br1, ar, alternative="two-sided")
            rows.append({"L2": l2, "n_AR": len(ar), "n_BR1": len(br1),
                         "p_wilcox": p, "median_AR": np.median(ar),
                         "median_BR1": np.median(br1),
                         "delta_med": np.median(br1) - np.median(ar)})
        stats_df = pd.DataFrame(rows)
        stats_df.to_csv(out_dir / f"{pn_low}_stats.csv", index=False)
        print(f"[{panel_name}] stats:\n{stats_df.to_string(index=False)}", flush=True)

        # plot
        fig, axes = plt.subplots(1, 4, figsize=(11, 3.5), sharey=True)
        colors = {"AR": "#0072B2", "BR1": "#E69F00"}
        for ax, l2 in zip(axes, FIBROS):
            sub = per_donor[per_donor["L2_label"] == l2]
            ar  = sub[sub["risk_class"] == "AR"]["mean_AUC"].values
            br1 = sub[sub["risk_class"] == "BR1"]["mean_AUC"].values
            bp = ax.boxplot([ar, br1], positions=[0, 1], widths=0.55,
                            patch_artist=True, showfliers=False,
                            medianprops=dict(color="black", linewidth=1),
                            boxprops=dict(linewidth=0.5),
                            whiskerprops=dict(linewidth=0.5),
                            capprops=dict(linewidth=0.5))
            for patch, lab in zip(bp["boxes"], ["AR","BR1"]):
                patch.set_facecolor(colors[lab]); patch.set_alpha(0.35)
                patch.set_edgecolor("black")
            for x_pos, vals in zip([0, 1], [ar, br1]):
                jitter = np.random.RandomState(42).normal(x_pos, 0.06, size=len(vals))
                ax.scatter(jitter, vals, color=colors[list(colors.keys())[x_pos]],
                           s=10, alpha=0.85, edgecolor="black", linewidth=0.3)
            ax.set_xticks([0, 1]); ax.set_xticklabels(["AR","BR1"], fontsize=8)
            l2_pretty = l2.replace("str__","").replace("_","-")
            row = stats_df[stats_df["L2"] == l2].iloc[0] if (stats_df["L2"]==l2).any() else None
            title = l2_pretty
            if row is not None and not np.isnan(row["p_wilcox"]):
                title += f"  p={row['p_wilcox']:.2g}"
            ax.set_title(title, fontsize=8.5, fontweight="bold")
            ax.tick_params(axis="y", labelsize=7)
            ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        axes[0].set_ylabel(f"AUCell ({panel_name})\n(per-donor mean)", fontsize=8)
        fig.suptitle(f"AUCell {panel_name} (tested cohort: AR-tested-neg vs BR1)",
                     fontsize=9, fontweight="bold", y=1.02)
        fig.tight_layout()
        fig.savefig(out_dir / f"{pn_low}_boxplot.pdf", bbox_inches="tight", dpi=600)
        plt.close(fig)

    # Master comparison table
    print("[8] write comparison_summary.md", flush=True)
    md = ["# AUCell iCAF scoring: tumor panels vs preneoplastic-breast panels", ""]
    md += ["Tested cohort: AR (brca_genotype=negative, n=28) + BR1 carriers (n=38). "
           "Per-donor mean AUC, Wilcoxon BR1 vs AR per fibroblast L2.", ""]
    md += ["## Panel provenance", "",
           "| Panel | n genes | Source | Context |",
           "|-------|---------|--------|---------|",
           f"| Ohlund2017      | {len(OHLUND2017)} | Ohlund D et al. 2017 JEM 214:579-596 (Tuveson-curated mouse iCAF, human-mapped; Wu 2020 GitHub commit 9c819769) | Pancreatic cancer (tumor) |",
           f"| Nee2023_preCAF  | {len(NEE2023_PRECAF)} | Nee K et al. 2023 Nat Genet 55:595-606 (Supp Table 6, FDR<0.05 log2FC>1, curated to inflammatory/secretome) | BRCA1 preneoplastic breast |",
           f"| Reed2024_FB2    | {len(REED2024_FB2)} | Reed AD et al. 2024 Nat Genet 56:652-662 (Supp Table 4, FB2 top 30 curated) | Healthy breast atlas (incl. BRCA1/2) |",
           f"| MMP_convergent  | {len(MMP_CONVERGENT)} | Nee 2023 ∩ Pal 2021 ∩ Reed 2024 FB2     | Cross-paper preneoplastic-breast |",
           ""]
    md += ["## Per-L2 Wilcoxon p-values across panels", ""]
    md += ["| L2 | " + " | ".join(f"{nm}_p" for nm,_ in ALL_PANELS) + " |"]
    md += ["|" + "---|" * (len(ALL_PANELS) + 1)]
    stats_per_panel = {}
    for nm, _ in ALL_PANELS:
        try:
            stats_per_panel[nm] = pd.read_csv(out_dir / f"{nm.lower()}_stats.csv").set_index("L2")
        except Exception:
            pass
    for l2 in FIBROS:
        cells = [l2.replace("str__","")]
        for nm, _ in ALL_PANELS:
            s = stats_per_panel.get(nm)
            if s is not None and l2 in s.index and not pd.isna(s.loc[l2,"p_wilcox"]):
                cells.append(f"{s.loc[l2,'p_wilcox']:.3g}")
            else:
                cells.append("—")
        md.append("| " + " | ".join(cells) + " |")
    md += ["", "## Per-L2 delta-medians (BR1 − AR) across panels", ""]
    md += ["| L2 | " + " | ".join(f"{nm}_dm" for nm,_ in ALL_PANELS) + " |"]
    md += ["|" + "---|" * (len(ALL_PANELS) + 1)]
    for l2 in FIBROS:
        cells = [l2.replace("str__","")]
        for nm, _ in ALL_PANELS:
            s = stats_per_panel.get(nm)
            if s is not None and l2 in s.index and not pd.isna(s.loc[l2,"delta_med"]):
                cells.append(f"{s.loc[l2,'delta_med']:+.3g}")
            else:
                cells.append("—")
        md.append("| " + " | ".join(cells) + " |")
    md += ["", "## Caveats",
           "- Nee 2023 may be one of the iHBCA component studies; if so, BR1 cells in iHBCA include cells whose pre-CAF signature was derived from themselves — partial circularity for that subset.",
           "- Reed 2024 FB2 is a cluster-identity signature, not a BR1-vs-AR DE signature; interpreting BR1>AR on this panel means BR1 donors enrich for FB2-like cells within each L2.",
           "- MMP_convergent (MMP3/10/12) is the cleanest cross-paper signal but only 3 genes — AUC is rank-stochastic on tiny panels.",
           "- No multiple-testing correction across 4 L2s × 5 panels.",
           "",
           "## Reference: Seurat-style (circular) result", "",
           "`outputs/plots/icaf_module_score/icaf_score_stats.csv` — same cohort, same L2s, Seurat AddModuleScore with data-derived 11-gene panel."]
    (out_dir / "comparison_summary.md").write_text("\n".join(md))
    print("Done.", flush=True)


if __name__ == "__main__":
    sys.exit(main())
