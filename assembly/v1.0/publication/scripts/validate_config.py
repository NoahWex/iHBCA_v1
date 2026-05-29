#!/usr/bin/env python3
"""Validate pipeline.yaml and config-driven I/O paths.

Runs all verification checks for the config migration:
1. YAML syntax validation
2. Path resolution (mapping files exist on disk)
3. Config completeness (all studies have required fields)
4. Cross-reference consistency (pipeline.yaml matches other configs)
5. Loader backward compatibility (config vs legacy paths agree)

Usage:
    python validate_config.py --repo-root .
"""

import argparse
import sys
from pathlib import Path

import yaml

# Add scripts dir for ihbca imports
sys.path.insert(0, str(Path(__file__).parent))


def check(name, condition, detail=""):
    status = "PASS" if condition else "FAIL"
    msg = f"  [{status}] {name}"
    if detail and not condition:
        msg += f" — {detail}"
    print(msg)
    return condition


def validate_syntax(repo_root, config):
    """1. YAML syntax and basic structure."""
    print("\n=== 1. Syntax Validation ===")
    ok = True
    ok &= check("pipeline.yaml parsed", config is not None)
    ok &= check("'paths' section exists", "paths" in config)
    ok &= check("'studies' section exists", "studies" in config)
    ok &= check("'integrated' section exists", "integrated" in config)
    ok &= check("'project' section exists", "project" in config)

    if "paths" in config:
        ok &= check("paths.mappings exists", "mappings" in config["paths"])
        ok &= check("paths.config exists", "config" in config["paths"])
        ok &= check("paths.outputs exists", "outputs" in config["paths"])
        ok &= check("paths.inputs exists", "inputs" in config["paths"])
    return ok


def validate_paths(repo_root, config):
    """2. All mapping/config paths resolve to existing files."""
    print("\n=== 2. Path Resolution ===")
    ok = True

    # Mapping files
    for key, rel in config["paths"]["mappings"].items():
        path = repo_root / rel
        ok &= check(f"mappings.{key}", path.exists(), str(path))

    # Config files
    for key, rel in config["paths"]["config"].items():
        path = repo_root / rel
        ok &= check(f"config.{key}", path.exists(), str(path))

    # Input paths (may not exist if outputs haven't been generated)
    harmonized = repo_root / config["paths"]["inputs"]["harmonized_metadata"]
    check(f"inputs.harmonized_metadata", harmonized.exists(),
          f"{harmonized} (OK if outputs not generated yet)")

    # Per-study mapping files
    for study, scfg in config["studies"].items():
        for mkey, mrel in scfg.get("mappings", {}).items():
            path = repo_root / mrel
            ok &= check(f"studies.{study}.mappings.{mkey}", path.exists(), str(path))

    return ok


def validate_completeness(config):
    """3. All studies have required fields."""
    print("\n=== 3. Config Completeness ===")
    ok = True

    required_study_fields = ["type", "cells", "donors", "output_filename"]

    for study, scfg in config["studies"].items():
        for field in required_study_fields:
            ok &= check(f"studies.{study}.{field}", field in scfg,
                         f"missing required field")

        # Non-pal studies need inputs.intermediates and inputs.published
        if study != "pal":
            has_inputs = "inputs" in scfg
            ok &= check(f"studies.{study}.inputs", has_inputs)
            if has_inputs:
                ok &= check(f"studies.{study}.inputs.intermediates",
                             "intermediates" in scfg["inputs"])
                ok &= check(f"studies.{study}.inputs.published",
                             "published" in scfg["inputs"])

        # Pal needs sub_studies
        if study == "pal":
            has_sub = "sub_studies" in scfg
            ok &= check(f"studies.pal.sub_studies", has_sub)
            if has_sub:
                ok &= check(f"studies.pal.sub_studies count == 3",
                             len(scfg["sub_studies"]) == 3,
                             f"got {len(scfg.get('sub_studies', []))}")

        # SRA run table
        has_sra = "sra_run_table" in scfg.get("mappings", {})
        ok &= check(f"studies.{study}.mappings.sra_run_table", has_sra)

    # Integrated
    int_cfg = config.get("integrated", {})
    ok &= check("integrated.output_filename", "output_filename" in int_cfg)
    ok &= check("integrated.inputs", "inputs" in int_cfg)

    return ok


def validate_cross_references(repo_root, config):
    """4. Cross-reference against dataset_metadata.yaml and registry."""
    print("\n=== 4. Cross-Reference Consistency ===")
    ok = True

    # Load dataset_metadata.yaml
    dm_path = repo_root / config["paths"]["config"]["dataset_metadata"]
    if dm_path.exists():
        with open(dm_path) as f:
            dm = yaml.safe_load(f)
        dm_studies = set(dm.get("datasets", {}).keys())
        cfg_studies = set(config["studies"].keys())
        ok &= check("studies match dataset_metadata.yaml",
                     cfg_studies == dm_studies,
                     f"pipeline: {cfg_studies}, metadata: {dm_studies}")
    else:
        check("dataset_metadata.yaml exists", False, str(dm_path))

    # Load source_dataset_registry.yaml
    reg_path = repo_root / config["paths"]["config"]["dataset_registry"]
    if reg_path.exists():
        with open(reg_path) as f:
            reg = yaml.safe_load(f)
        # Check output filenames match for external studies
        reg_datasets = reg.get("datasets", {})
        for study, scfg in config["studies"].items():
            if study in reg_datasets:
                reg_fn = reg_datasets[study].get("output_filename", "")
                cfg_fn = scfg.get("output_filename", "")
                ok &= check(f"filename match: {study}",
                             reg_fn == cfg_fn,
                             f"pipeline={cfg_fn}, registry={reg_fn}")

            # Cell count comparison (warning only)
            if study in reg_datasets:
                reg_cells = reg_datasets[study].get("expected_cells")
                cfg_cells = scfg.get("cells")
                if reg_cells and cfg_cells and reg_cells != cfg_cells:
                    print(f"  [WARN] {study} cell count: pipeline={cfg_cells}, registry={reg_cells}")
    else:
        check("source_dataset_registry.yaml exists", False, str(reg_path))

    return ok


def validate_loaders(repo_root, config):
    """5. Loader backward compatibility."""
    print("\n=== 5. Loader Backward Compatibility ===")
    ok = True

    try:
        from ihbca.loaders import load_pipeline_config
        loaded = load_pipeline_config(repo_root)
        ok &= check("load_pipeline_config() succeeds", loaded is not None)
        ok &= check("config has _repo_root", "_repo_root" in loaded)
    except Exception as e:
        ok &= check("import ihbca.loaders", False, str(e))

    try:
        from ihbca.constants import get_cxg_studies, get_pal_sub_studies, get_study_list
        cxg = get_cxg_studies(config)
        ok &= check("get_cxg_studies()", len(cxg) > 0, f"got {cxg}")
        pal = get_pal_sub_studies(config)
        ok &= check("get_pal_sub_studies()", len(pal) == 3, f"got {pal}")
        studies = get_study_list(config)
        ok &= check("get_study_list()", len(studies) == 7, f"got {len(studies)}")
    except Exception as e:
        ok &= check("import ihbca.constants helpers", False, str(e))

    return ok


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", required=True, help="Path to repo root")
    args = parser.parse_args()

    repo_root = Path(args.repo_root).resolve()

    # Load pipeline.yaml
    pipeline_path = repo_root / "publication/config/pipeline.yaml"
    print(f"Validating: {pipeline_path}")

    if not pipeline_path.exists():
        print(f"FATAL: pipeline.yaml not found at {pipeline_path}")
        sys.exit(1)

    with open(pipeline_path) as f:
        config = yaml.safe_load(f)

    all_ok = True
    all_ok &= validate_syntax(repo_root, config)
    all_ok &= validate_paths(repo_root, config)
    all_ok &= validate_completeness(config)
    all_ok &= validate_cross_references(repo_root, config)
    all_ok &= validate_loaders(repo_root, config)

    print(f"\n{'=' * 40}")
    if all_ok:
        print("ALL CHECKS PASSED")
    else:
        print("SOME CHECKS FAILED — review above")
    print(f"{'=' * 40}")

    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
