"""
build_fig2_marker_swarm_substrate.py

Aggregates per-contrast marker-positive NhoodGroup rosters from the V1 Clinical DA
pipeline into a single long-format CSV that the marker-swarm renderer consumes.

Source layouts differ between the parity × BRCA1 interaction contrast (which carries
interaction-median effect size and an NG_label) and the four single-contrast risk/parity
main-effects rosters (which carry a per-group median logFC and a name_short). This
script normalises both into a common schema:

    contrast | parent_compartment | parent_L2_joint | NhoodGroup |
    effect_size | name_short | n_up_markers | marker_positive | size_rank_within_L2

Output: <out-dir>/marker_positive_ngs_per_contrast.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


CONTRASTS = {
    "parity_in_AR":      "risk_effects",
    "BR1_vs_AR":         "risk_effects",
    "BR2_vs_AR":         "risk_effects",
    "HRS_vs_AR":         "risk_effects",
    "parity_x_HR_BRCA1": "pxr",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument("--risk-effects-base", default=None, type=Path,
                   help="Directory containing {parity_in_AR,BR1_vs_AR,BR2_vs_AR,HRS_vs_AR}/marker_positive_ngs.csv. "
                        "Default: <project-root>/publication/analysis/abundance/clinical/risk_main_effects/figure_data/v1_da_widening_20260513/")
    p.add_argument("--pxr-roster", default=None, type=Path,
                   help="Path to PxR marker_positive_ngs_v4_br1.csv. "
                        "Default: <project-root>/publication/analysis/abundance/clinical/parity_findings/figure_data/fig2_pxr_v4_br1/marker_positive_ngs_v4_br1.csv")
    p.add_argument("--out-dir", default=None, type=Path,
                   help="Default: <project-root>/publication/figures/data/fig2/")
    args = p.parse_args()
    root = args.project_root.resolve()
    if args.risk_effects_base is None:
        args.risk_effects_base = (
            root / "publication" / "analysis" / "abundance" / "clinical"
            / "risk_main_effects" / "figure_data" / "v1_da_widening_20260513"
        )
    if args.pxr_roster is None:
        args.pxr_roster = (
            root / "publication" / "analysis" / "abundance" / "clinical"
            / "parity_findings" / "figure_data" / "fig2_pxr_v4_br1"
            / "marker_positive_ngs_v4_br1.csv"
        )
    if args.out_dir is None:
        args.out_dir = root / "publication" / "figures" / "data" / "fig2"
    return args


COLS = [
    "contrast", "parent_compartment", "parent_L2_joint", "NhoodGroup",
    "effect_size", "name_short", "n_up_markers", "marker_positive",
    "size_rank_within_L2",
]


def load_risk_effects(contrast: str, base: Path) -> pd.DataFrame:
    path = base / contrast / "marker_positive_ngs.csv"
    df = pd.read_csv(path)
    df = df.assign(
        contrast=contrast,
        effect_size=df["group_med_lfc"],
    )
    return df[COLS]  # type: ignore[return-value]


def load_pxr(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df.assign(
        contrast="parity_x_HR_BRCA1",
        effect_size=df["int_med"],
        name_short=df["NG_label"],
    )
    return df[COLS]  # type: ignore[return-value]


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    pieces = []
    for contrast, source in CONTRASTS.items():
        if source == "risk_effects":
            pieces.append(load_risk_effects(contrast, args.risk_effects_base))
        else:
            pieces.append(load_pxr(args.pxr_roster))

    out = pd.concat(pieces, ignore_index=True)
    out["marker_positive"] = out["marker_positive"].astype(str).str.lower().isin(
        ["true", "t", "1"]
    )

    out_path = args.out_dir / "marker_positive_ngs_per_contrast.csv"
    out.to_csv(out_path, index=False)

    print(f"Wrote {out_path}  ({len(out)} rows, {out['marker_positive'].sum()} marker+ across "
          f"{out['contrast'].nunique()} contrasts)")
    for c, sub in out.groupby("contrast"):
        n_marker = int(sub["marker_positive"].sum())
        n_total = len(sub)
        print(f"  {c}: {n_marker} marker+ / {n_total} candidates")


if __name__ == "__main__":
    main()
