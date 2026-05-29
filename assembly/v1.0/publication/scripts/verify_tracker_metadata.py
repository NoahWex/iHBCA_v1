#!/usr/bin/env python3
"""Verify tracker dataset-level metadata alignment between YAML config and h5ad uns.

Reads publication/config/dataset_metadata.yaml and cross-checks against the
actual uns group in each source h5ad + integrated object.

Usage:
    python verify_tracker_metadata.py --project-root /path/to/repo
"""
import argparse
import sys
from pathlib import Path

import h5py
import yaml

TRACKER_FIELDS = [
    "alignment_software",
    "contact_email",
    "description",
    "gene_annotation_version",
    "reference_genome",
    "sequenced_fragment",
    "study_pi",
    "institute",
    "disease_ontology_term_id",
]

STUDIES = ["gray", "kumar", "murrow", "nee", "twigger", "reed", "pal"]
STUDY_FILES = {
    "gray": "gray2022.h5ad",
    "kumar": "kumar2023.h5ad",
    "murrow": "murrow2022.h5ad",
    "nee": "nee2023.h5ad",
    "twigger": "twigger2022.h5ad",
    "reed": "reed2024.h5ad",
    "pal": "pal2021.h5ad",
}


def read_uns_h5py(h5ad_path):
    """Read uns keys from h5ad via h5py (avoids loading full anndata)."""
    uns = {}
    with h5py.File(h5ad_path, "r") as f:
        if "uns" not in f:
            return uns
        uns_grp = f["uns"]
        for key in uns_grp.keys():
            item = uns_grp[key]
            if isinstance(item, h5py.Dataset):
                val = item[()]
                if isinstance(val, bytes):
                    val = val.decode("utf-8")
                elif hasattr(val, "tolist"):
                    val = val.tolist()
                # Decode bytes inside lists/arrays
                if isinstance(val, list):
                    val = [
                        v.decode("utf-8") if isinstance(v, bytes) else v
                        for v in val
                    ]
                uns[key] = val
            elif isinstance(item, h5py.Group):
                # Could be a categorical or array — try to read
                if "categories" in item:
                    cats = [
                        c.decode("utf-8") if isinstance(c, bytes) else c
                        for c in item["categories"][()]
                    ]
                    uns[key] = cats
                else:
                    uns[key] = f"<group: {list(item.keys())}>"
    return uns


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    args = parser.parse_args()

    root = args.project_root
    yaml_path = root / "publication/config/dataset_metadata.yaml"
    source_dir = root / "publication/outputs/source_datasets"
    integrated_path = root / "publication/outputs/integrated_objects/all-breast-cells.h5ad"

    with open(yaml_path) as f:
        config = yaml.safe_load(f)
    datasets = config["datasets"]

    all_ok = True
    print("=" * 70)
    print("TRACKER METADATA ALIGNMENT CHECK")
    print("=" * 70)

    # Check source datasets
    for study in STUDIES:
        h5ad_path = source_dir / STUDY_FILES[study]
        yaml_meta = datasets.get(study, {})
        print(f"\n--- {study} ({STUDY_FILES[study]}) ---")

        if not h5ad_path.exists():
            print(f"  MISSING: {h5ad_path}")
            all_ok = False
            continue

        uns = read_uns_h5py(h5ad_path)

        for field in TRACKER_FIELDS:
            yaml_val = yaml_meta.get(field)
            uns_val = uns.get(field)

            if yaml_val is None:
                print(f"  {field}: NOT IN YAML")
                all_ok = False
            elif uns_val is None:
                print(f"  {field}: NOT IN UNS (yaml={yaml_val})")
                # Not necessarily a problem — some fields are dataset-level only
            else:
                # Normalize for comparison
                y_str = str(yaml_val)
                u_str = str(uns_val)
                if y_str == u_str:
                    print(f"  {field}: OK")
                else:
                    print(f"  {field}: MISMATCH")
                    print(f"    yaml: {y_str}")
                    print(f"    uns:  {u_str}")
                    all_ok = False

        # Report any extra uns keys that might be relevant
        extra = set(uns.keys()) - set(TRACKER_FIELDS) - {
            "schema_type", "schema_version", "title", "batch_condition",
            "default_embedding", "X_approximate_distribution",
        }
        if extra:
            print(f"  Extra uns keys: {sorted(extra)}")

    # Check integrated object
    print(f"\n--- integrated (all-breast-cells.h5ad) ---")
    if integrated_path.exists():
        uns = read_uns_h5py(integrated_path)
        for field in TRACKER_FIELDS:
            uns_val = uns.get(field)
            if uns_val is not None:
                print(f"  {field}: {uns_val}")
            else:
                print(f"  {field}: NOT IN UNS")
    else:
        print(f"  MISSING: {integrated_path}")

    print(f"\n{'=' * 70}")
    if all_ok:
        print("RESULT: ALL CHECKS PASSED")
    else:
        print("RESULT: MISMATCHES OR MISSING FIELDS DETECTED")
    print(f"{'=' * 70}")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
