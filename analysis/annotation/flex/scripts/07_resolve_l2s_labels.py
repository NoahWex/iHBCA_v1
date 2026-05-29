"""07_resolve_l2s_labels.py — resolve per-cell FLEX L2S labels from the
per-compartment annotation yamls (v2_flex schema).

Per compartment (epi/imm/str):
  1. Read annotation yaml (anchor_resolution + clusters_by_resolution map of
     cluster_id -> assigned_label, plus fine_cluster_overrides).
  2. Read multi-resolution clusters CSV (cell_id + leiden_{resolution} cols).
  3. Apply anchor label per cell, then overlay fine_cluster_overrides at the
     override-specified resolutions.

Output: a single flex_l2s_labels.csv concatenated across compartments with
schema `cell_id, compartment, l2s_anchor_label, l2s_label,
l2s_override_applied, l2s_override_resolution`.

Two yaml layouts are supported:
  - Publication layout (preferred): --yamls-dir points at a directory holding
    epithelial.yaml, immune.yaml, stromal.yaml (one yaml per compartment,
    flat). --clusters-dir points at the directory holding
    {Compartment}_clusters.csv files.
  - Dev layout: --annotation-root points at a directory holding Epithelial/,
    Immune/, Stromal/ subdirectories (each with annotation_v2s.yaml or
    annotation_v2s_DRAFT.yaml) plus a clusters/ sibling. Provided for
    backward compatibility with the dev tree.

--yamls-dir + --clusters-dir takes precedence over --annotation-root when
both are given.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
import yaml

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("build_flex_l2s")


COMPARTMENT_DIRS = {
    "epi": "Epithelial",
    "imm": "Immune",
    "str": "Stromal",
}


def find_yaml_publication(yamls_dir: Path, comp: str) -> Path:
    """Publication layout: yamls_dir/{epi|imm|str -> epithelial|immune|stromal}.yaml."""
    name_map = {"epi": "epithelial.yaml", "imm": "immune.yaml", "str": "stromal.yaml"}
    p = yamls_dir / name_map[comp]
    if not p.exists():
        raise FileNotFoundError(f"no yaml at publication path: {p}")
    return p


def find_yaml_dev(comp_dir: Path) -> Path:
    """Dev layout: prefer locked annotation_v2s.yaml; fall back to DRAFT."""
    locked = comp_dir / "annotation_v2s.yaml"
    if locked.exists():
        return locked
    draft = comp_dir / "annotation_v2s_DRAFT.yaml"
    if draft.exists():
        return draft
    raise FileNotFoundError(f"no annotation_v2s yaml under {comp_dir}")


def resolve_compartment(yaml_path: Path, clusters_csv: Path, comp: str) -> pd.DataFrame:
    log.info("[%s] yaml=%s clusters=%s", comp, yaml_path.name, clusters_csv.name)
    spec = yaml.safe_load(yaml_path.read_text())

    anchor_res = str(spec["anchor_resolution"])
    anchor_blk = spec["clusters_by_resolution"][anchor_res]
    # cluster key is "c{N}"; column in clusters CSV is integer
    anchor_map = {
        int(k.lstrip("c")): v["assigned_label"]
        for k, v in anchor_blk.items()
    }
    log.info("[%s]   anchor_resolution=%s, %d anchor labels", comp, anchor_res, len(anchor_map))

    cells = pd.read_csv(clusters_csv)
    anchor_col = f"leiden_{anchor_res}"
    if anchor_col not in cells.columns:
        raise ValueError(f"[{comp}] missing column {anchor_col} in {clusters_csv}")

    anchor_labels = cells[anchor_col].map(anchor_map)
    n_unmapped = int(anchor_labels.isna().sum())
    if n_unmapped:
        unmapped = sorted(cells.loc[anchor_labels.isna(), anchor_col].unique().tolist())
        log.warning("[%s]   %d cells with anchor cluster not in yaml (clusters: %s)", comp, n_unmapped, unmapped)

    out = pd.DataFrame({
        "cell_id": cells["cell_id"],
        "compartment": comp,
        "l2s_anchor_label": anchor_labels,
        "l2s_label": anchor_labels.copy(),
        "l2s_override_applied": False,
        "l2s_override_resolution": pd.NA,
    })

    overrides = spec.get("fine_cluster_overrides") or {}
    n_overridden_total = 0
    for ov_name, ov_block in overrides.items():
        ov_res = str(ov_block.get("resolution"))
        ov_col = f"leiden_{ov_res}"
        if ov_col not in cells.columns:
            log.warning("[%s]   override %s wants %s but column absent; skipped", comp, ov_name, ov_col)
            continue
        for tup in ov_block.get("overrides", []):
            tup_res = str(tup.get("resolution", ov_res))
            tup_col = f"leiden_{tup_res}"
            if tup_col not in cells.columns:
                log.warning("[%s]   override %s tuple wants %s but absent; skipped", comp, ov_name, tup_col)
                continue
            cluster_id = int(tup["cluster"])
            new_label = str(tup["label"])
            mask = cells[tup_col] == cluster_id
            n = int(mask.sum())
            if n == 0:
                log.warning("[%s]   override %s {%s,c%d}: 0 cells matched", comp, ov_name, tup_res, cluster_id)
                continue
            out.loc[mask, "l2s_label"] = new_label
            out.loc[mask, "l2s_override_applied"] = True
            out.loc[mask, "l2s_override_resolution"] = tup_res
            n_overridden_total += n
            log.info("[%s]   override %s {%s,c%d} -> %s: %d cells", comp, ov_name, tup_res, cluster_id, new_label, n)

    log.info("[%s]   total cells: %d, overridden: %d (%.2f%%)",
             comp, len(out), n_overridden_total, 100.0 * n_overridden_total / max(len(out), 1))
    return out


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--yamls-dir", type=Path, default=None,
                   help="publication layout: directory holding {epithelial,immune,stromal}.yaml")
    p.add_argument("--clusters-dir", type=Path, default=None,
                   help="publication layout: directory holding {Compartment}_clusters.csv (required when --yamls-dir is used)")
    p.add_argument("--annotation-root", type=Path, default=None,
                   help="dev layout: parent of {Epithelial,Immune,Stromal}/ and clusters/")
    p.add_argument("--out-csv", required=True, type=Path,
                   help="output flex_l2s_labels.csv path")
    args = p.parse_args()
    if args.yamls_dir is not None and args.clusters_dir is None:
        p.error("--yamls-dir requires --clusters-dir")
    if args.yamls_dir is None and args.annotation_root is None:
        p.error("either --yamls-dir + --clusters-dir or --annotation-root must be given")
    return args


def main() -> int:
    args = parse_args()

    if args.yamls_dir is not None:
        yamls_dir = args.yamls_dir
        clusters_dir = args.clusters_dir
        log.info("layout: publication (yamls=%s, clusters=%s)", yamls_dir, clusters_dir)
    else:
        yamls_dir = None
        clusters_dir = args.annotation_root / "clusters"
        log.info("layout: dev (annotation_root=%s)", args.annotation_root)

    parts = []
    for comp, dir_name in COMPARTMENT_DIRS.items():
        clusters_csv = clusters_dir / f"{dir_name}_clusters.csv"
        if not clusters_csv.exists():
            log.error("missing clusters csv: %s", clusters_csv)
            return 1
        if yamls_dir is not None:
            yaml_path = find_yaml_publication(yamls_dir, comp)
        else:
            comp_dir = args.annotation_root / dir_name
            if not comp_dir.is_dir():
                log.error("missing compartment dir: %s", comp_dir)
                return 1
            yaml_path = find_yaml_dev(comp_dir)
        parts.append(resolve_compartment(yaml_path, clusters_csv, comp))

    full = pd.concat(parts, ignore_index=True)
    log.info("combined: %d cells across %d compartments", len(full), len(parts))
    log.info("label distribution (top 25):")
    for lbl, n in full["l2s_label"].value_counts().head(25).items():
        log.info("  %s: %d", lbl, n)

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    full.to_csv(args.out_csv, index=False)
    log.info("wrote %s", args.out_csv)
    return 0


if __name__ == "__main__":
    sys.exit(main())
