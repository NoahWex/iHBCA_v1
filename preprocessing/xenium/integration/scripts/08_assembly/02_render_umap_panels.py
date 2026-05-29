"""Render v5 UMAP panels including patient clinical metadata.

Panels produced (PDF + PNG each):
  umap_platform, umap_compartment, umap_patient, umap_l0p5,
  umap_age, umap_menopause, umap_brca, umap_position,
  umap_nmp_exclude
"""
import argparse, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Patient-level metadata from raw_data_manifest.yaml
PATIENT_META = {
    "Pat1":      {"age": 72, "menopause": "Post",  "brca": "negative"},
    "Pat2":      {"age": 45, "menopause": "Pre",   "brca": "negative"},
    "UCI604":    {"age": 37, "menopause": "Pre",   "brca": "negative"},
    "UCI220228": {"age": 39, "menopause": "Pre",   "brca": "BRCA1"},
}


def scatter_cat(ax, x, y, cat, title, palette=None, s=0.8, alpha=0.25,
                xlim=None, ylim=None, seed=0):
    cats = pd.Categorical(cat)
    codes = cats.codes
    names = cats.categories
    if palette is None:
        cmap = plt.get_cmap("tab20", len(names))
        palette = [cmap(i) for i in range(len(names))]
    order = np.random.default_rng(seed).permutation(len(x))
    ax.scatter(x[order], y[order], c=[palette[c] for c in codes[order]],
               s=s, alpha=alpha, rasterized=True, linewidths=0, edgecolors="none")
    if xlim: ax.set_xlim(xlim)
    if ylim: ax.set_ylim(ylim)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, fontsize=8)
    for sp in ax.spines.values(): sp.set_linewidth(0.25)
    handles = [plt.Line2D([0], [0], marker="o", color="w",
                          markerfacecolor=palette[i], markersize=4, label=names[i])
               for i in range(len(names))]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5),
              fontsize=5, frameon=False, markerscale=1.5)


def scatter_cont(ax, x, y, vals, title, cmap="viridis", vmin=None, vmax=None,
                 s=0.8, alpha=0.25, xlim=None, ylim=None, seed=0):
    m = np.isfinite(vals)
    order = np.random.default_rng(seed).permutation(len(x))
    sc = ax.scatter(x[order][m[order]], y[order][m[order]],
                    c=vals[order][m[order]], cmap=cmap,
                    vmin=vmin, vmax=vmax,
                    s=s, alpha=alpha, rasterized=True, linewidths=0)
    if xlim: ax.set_xlim(xlim)
    if ylim: ax.set_ylim(ylim)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, fontsize=8)
    for sp in ax.spines.values(): sp.set_linewidth(0.25)
    plt.colorbar(sc, ax=ax, shrink=0.6, pad=0.02)


def save(fig, path_base):
    fig.tight_layout()
    fig.savefig(path_base + ".pdf", dpi=300, bbox_inches="tight")
    fig.savefig(path_base + ".png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  {os.path.basename(path_base)}.pdf", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--umap", required=True)
    ap.add_argument("--obs", required=True)
    ap.add_argument("--l0p5", required=True)
    ap.add_argument("--nmp", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print("Loading...", flush=True)
    u = pd.read_csv(args.umap).set_index("cell_id")
    obs = pd.read_csv(args.obs).set_index("cell_id").reindex(u.index)
    l0 = pd.read_csv(args.l0p5, usecols=["cell_id", "l0p5_final"]).set_index("cell_id")
    nmp = pd.read_csv(args.nmp, usecols=["cell_id", "nmp_nuclear",
                                          "nmp_cyto"]).set_index("cell_id")
    obs["l0p5"] = l0["l0p5_final"].reindex(obs.index)
    obs["nmp_exclude"] = (nmp["nmp_nuclear"] > nmp["nmp_cyto"]).reindex(obs.index)

    # Patient metadata join
    for col, key in [("age", "age"), ("menopause", "menopause"), ("brca", "brca")]:
        obs[col] = obs["patient_id"].map(
            {pid: meta[key] for pid, meta in PATIENT_META.items()})

    # Anatomical position from xenium_id (xenium cells only)
    # xenium_id format: {patient_id}_{position_id}_xenium_{N}
    def extract_position(cell_id):
        # cell_id encodes xenium_id in its prefix for xenium cells
        # obs doesn't have xenium_id directly; derive from patient_id prefix
        return None

    # For xenium cells derive position from cell_id pattern:
    # cell_id = {xenium_id}_{barcode} where xenium_id = {pat}_{pos}_xenium_{n}
    # Strip patient prefix and barcode to get position
    def parse_position(cell_id, patient_id):
        # cell_id starts with patient_id_
        rest = cell_id[len(patient_id) + 1:]  # e.g. "P1_xenium_1_BARCODE"
        parts = rest.split("_xenium_")
        if len(parts) >= 1:
            return parts[0]  # position_id e.g. "P1", "P2_Lower"
        return "unknown"

    obs["position"] = [
        parse_position(cid, pid) if "_xenium_" in cid else "FLEX"
        for cid, pid in zip(obs.index, obs["patient_id"].values)
    ]
    # Simplify position to top-level (P1, P2, P3, etc.)
    obs["position_top"] = obs["position"].str.split("_").str[0]

    x = u["UMAP1"].values; y = u["UMAP2"].values
    n = len(x)
    print(f"  {n:,} cells", flush=True)
    xlim = tuple(np.percentile(x, [0.1, 99.9]))
    ylim = tuple(np.percentile(y, [0.1, 99.9]))
    pad = lambda lim, r=0.03: (lim[0] - r*(lim[1]-lim[0]), lim[1] + r*(lim[1]-lim[0]))
    xlim = pad(xlim); ylim = pad(ylim)

    kw = dict(xlim=xlim, ylim=ylim)

    # 1. Platform
    fig, ax = plt.subplots(figsize=(7, 5.5))
    scatter_cat(ax, x, y, obs["platform"].fillna("na"),
                f"Platform (n={n:,})",
                palette=["#1f77b4", "#ff7f0e"], **kw)
    save(fig, os.path.join(args.out_dir, "umap_platform"))

    # 2. Compartment
    fig, ax = plt.subplots(figsize=(7, 5.5))
    scatter_cat(ax, x, y, obs["compartment"].fillna("na"), "Compartment",
                palette=["#d62728", "#2ca02c", "#9467bd", "#8c564b"], **kw)
    save(fig, os.path.join(args.out_dir, "umap_compartment"))

    # 3. Patient
    fig, ax = plt.subplots(figsize=(7, 5.5))
    scatter_cat(ax, x, y, obs["patient_id"].fillna("na"), "Patient", **kw)
    save(fig, os.path.join(args.out_dir, "umap_patient"))

    # 4. L0.5
    has_l0 = obs["l0p5"].notna().values
    fig, ax = plt.subplots(figsize=(8, 5.5))
    scatter_cat(ax, x[has_l0], y[has_l0], obs["l0p5"][has_l0],
                f"L0.5 (n={int(has_l0.sum()):,})", **kw)
    save(fig, os.path.join(args.out_dir, "umap_l0p5"))

    # 5. Age (continuous, Xenium + FLEX, patient-level)
    age_vals = obs["age"].values.astype(float)
    fig, ax = plt.subplots(figsize=(7, 5.5))
    scatter_cont(ax, x, y, age_vals, "Age (years)",
                 cmap="plasma", vmin=35, vmax=75, **kw)
    save(fig, os.path.join(args.out_dir, "umap_age"))

    # 6. Menopausal status
    fig, ax = plt.subplots(figsize=(7, 5.5))
    scatter_cat(ax, x, y, obs["menopause"].fillna("na"), "Menopausal status",
                palette=["#e377c2", "#17becf", "#7f7f7f"], **kw)
    save(fig, os.path.join(args.out_dir, "umap_menopause"))

    # 7. BRCA status
    fig, ax = plt.subplots(figsize=(7, 5.5))
    scatter_cat(ax, x, y, obs["brca"].fillna("na"), "BRCA status",
                palette=["#aec7e8", "#d62728", "#7f7f7f"], **kw)
    save(fig, os.path.join(args.out_dir, "umap_brca"))

    # 8. Anatomical position (top-level)
    fig, ax = plt.subplots(figsize=(8, 5.5))
    scatter_cat(ax, x, y, obs["position_top"].fillna("na"),
                "Anatomical position (top-level)", **kw)
    save(fig, os.path.join(args.out_dir, "umap_position"))

    # 9. NMP exclude (Xenium only)
    xen_mask = (obs["platform"] == "xenium").values
    flag = obs.loc[xen_mask, "nmp_exclude"].astype("boolean").fillna(False).astype(str).values
    fig, ax = plt.subplots(figsize=(7, 5.5))
    scatter_cat(ax, x[xen_mask], y[xen_mask], flag,
                "NMP exclude (Xenium)",
                palette=["#bbbbbb", "#d62728"], **kw)
    save(fig, os.path.join(args.out_dir, "umap_nmp_exclude"))

    print("Done.", flush=True)


if __name__ == "__main__":
    main()
