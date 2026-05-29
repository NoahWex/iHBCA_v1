"""
07_concat_full_atlas.py — concat per-compartment labels.csv into full-atlas artifact.

Output schema (per coordination/artifacts/l2_labels.md):
    cell_id        — stable across versions; contractual
    compartment    — imm / str / epi (= L1)
    lineage        — intra-compartment lineage (from per-compartment labels.csv)
    family         — LBridge cross-compartment family (NA placeholder; populated
                     by xenium/flex agent in a later release)
    label          — L2.0
    is_artifact    — bool; consumers grey rather than drop
    source_mode    — cluster / fine_cluster / expression
    figure_target  — nullable; cross-figure flag sourced from each label's
                     `figure_target` field in its compartment YAML. Comma-
                     delimited if a label spans multiple figures.

The figure_target lookup is built from the per-compartment annotation YAMLs
(--imm-yaml / --str-yaml / --epi-yaml) so the schema lives next to the label
definition rather than as a hardcoded dict here. Empty / null values are
encoded as NaN in the output column.

A validation gate enforces that every non-null figure_target value resolves
to a figure id in publication/figures/design.yaml — typos in YAML edits
fail fast rather than silently dropping a label from a downstream figure.

Hard-fail conditions:
    - any per-compartment labels.csv missing
    - any per-compartment YAML missing (when --*-yaml is supplied)
    - duplicate cell_id across compartments
    - row count mismatch vs sum of compartment row counts
    - figure_target value not present in design.yaml figure ids
"""

from __future__ import annotations

import argparse
import logging
import sys
import warnings
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set

import pandas as pd
import yaml

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("concat_full_atlas")


EXPECTED_COLS = {"cell_id", "label", "lineage", "is_artifact", "source_mode"}
COMPARTMENTS = ("imm", "str", "epi")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outputs-root", required=True, type=Path, help="Pipeline outputs root containing imm/, str/, epi/")
    p.add_argument("--out-csv", required=True, type=Path)
    p.add_argument("--imm-yaml", required=True, type=Path, help="immune annotation v2 YAML (label -> figure_target source)")
    p.add_argument("--str-yaml", required=True, type=Path, help="stromal annotation v2 YAML")
    p.add_argument("--epi-yaml", required=True, type=Path, help="epithelial annotation v2 YAML")
    p.add_argument("--design-yaml", required=False, default=None, type=Path,
                   help="publication/figures/design.yaml (for figure_target validation). "
                        "Optional: if omitted or missing, figure_target validation is skipped.")
    return p.parse_args()


def load_compartment(outputs_root: Path, compartment: str) -> pd.DataFrame:
    path = outputs_root / compartment / "labels.csv"
    if not path.exists():
        raise RuntimeError(f"missing labels.csv for compartment {compartment}: {path}")
    df = pd.read_csv(path)
    missing = EXPECTED_COLS - set(df.columns)
    if missing:
        raise RuntimeError(f"{compartment}/labels.csv missing columns: {sorted(missing)}")
    df = df.copy()
    df["compartment"] = compartment
    log.info("%s: loaded %d cells (%d unique labels)", compartment, len(df), df["label"].nunique())
    return df


def load_figure_targets(yaml_paths: Dict[str, Path]) -> Dict[str, str]:
    """Read each compartment YAML and extract {label: figure_target} entries.

    Labels without a figure_target field are omitted (downstream NaN).
    """
    targets: Dict[str, str] = {}
    for compartment, yaml_path in yaml_paths.items():
        if not yaml_path.exists():
            raise RuntimeError(f"missing annotation YAML for {compartment}: {yaml_path}")
        cfg = yaml.safe_load(yaml_path.read_text()) or {}
        labels = cfg.get("labels", {}) or {}
        n_with_target = 0
        for label_name, label_def in labels.items():
            if not isinstance(label_def, dict):
                continue
            ft = label_def.get("figure_target")
            if ft is None or ft == "":
                continue
            if label_name in targets and targets[label_name] != ft:
                raise RuntimeError(
                    f"figure_target collision for label {label_name!r}: "
                    f"{targets[label_name]!r} vs {ft!r}"
                )
            targets[label_name] = str(ft)
            n_with_target += 1
        log.info("%s YAML: %d labels with figure_target", compartment, n_with_target)
    return targets


def load_design_figure_ids(design_yaml: Path) -> Set[str]:
    """Return the set of valid figure ids declared in design.yaml (e.g. fig1..fig5)."""
    cfg = yaml.safe_load(design_yaml.read_text()) or {}
    figures = cfg.get("figures", {}) or {}
    return set(figures.keys())


def validate_figure_targets(targets: Dict[str, str], valid_ids: Set[str]) -> None:
    """Each figure_target value (split on `,`) must be in valid_ids."""
    bad: Dict[str, List[str]] = {}
    for label, ft in targets.items():
        parts = [p.strip() for p in str(ft).split(",") if p.strip()]
        unknown = [p for p in parts if p not in valid_ids]
        if unknown:
            bad[label] = unknown
    if bad:
        raise RuntimeError(
            f"figure_target values not present in design.yaml: {bad}; "
            f"valid figure ids: {sorted(valid_ids)}"
        )


def main() -> int:
    args = parse_args()

    parts = [load_compartment(args.outputs_root, c) for c in COMPARTMENTS]
    expected_total = sum(len(p) for p in parts)

    full = pd.concat(parts, ignore_index=True)
    if len(full) != expected_total:
        raise RuntimeError(f"concat length mismatch: {len(full)} != {expected_total}")

    dupes = full["cell_id"].duplicated()
    if dupes.any():
        n_dupes = int(dupes.sum())
        examples = full.loc[dupes, "cell_id"].head(5).tolist()
        raise RuntimeError(f"duplicate cell_id across compartments: {n_dupes} dupes, first 5: {examples}")

    figure_targets = load_figure_targets(
        {"imm": args.imm_yaml, "str": args.str_yaml, "epi": args.epi_yaml}
    )
    if args.design_yaml is not None and args.design_yaml.exists():
        valid_ids = load_design_figure_ids(args.design_yaml)
        validate_figure_targets(figure_targets, valid_ids)
        log.info("figure_target lookup: %d labels (validated against %d design.yaml figure ids)",
                 len(figure_targets), len(valid_ids))
    else:
        warnings.warn(
            "design.yaml not provided or not found — skipping figure_target validation. "
            "figure_target column will still be populated from annotation YAMLs but "
            "values are not checked against a design manifest.",
            stacklevel=2,
        )
        log.warning("figure_target validation skipped (design.yaml unavailable)")
        log.info("figure_target lookup: %d labels (unvalidated)", len(figure_targets))

    full["family"] = pd.NA
    full["figure_target"] = full["label"].map(figure_targets).astype("object")

    column_order = [
        "cell_id",
        "compartment",
        "lineage",
        "family",
        "label",
        "is_artifact",
        "source_mode",
        "figure_target",
    ]
    full = full[column_order]

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    full.to_csv(args.out_csv, index=False)
    log.info("wrote %s (%d rows, %d columns)", args.out_csv, len(full), full.shape[1])

    label_col = full["label"]
    compartment_col = full["compartment"]
    log.info("compartment counts: %s", compartment_col.value_counts().to_dict())
    n_per_comp = full.groupby("compartment")["label"].nunique().sum()
    n_label_only = label_col.nunique()
    log.info("unique (compartment, label): %d", n_per_comp)
    log.info("unique label-only: %d (collisions: %d)", n_label_only, n_per_comp - n_label_only)
    log.info("artifact rows: %d", int(full["is_artifact"].sum()))
    log.info("figure_target population: %s", full["figure_target"].value_counts(dropna=True).to_dict())

    return 0


if __name__ == "__main__":
    sys.exit(main())
