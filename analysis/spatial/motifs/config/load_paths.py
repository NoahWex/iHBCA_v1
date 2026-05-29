"""Motifs pipeline internal-layout resolver (hybrid config, Option B).

This module models the motifs pipeline's internal output directory layout
(outputs/, outputs/validation/, outputs/motif_units/). For global publication
discovery (cross-repo substrate, containers, source datasets) consumers must
use publication/config/load_paths.py — its resolve_path("xenium_joint_l1p5", ...)
returns the canonical CRSP path.

Standard usage from a motif script:

    from load_paths import load_paths, resolve

    paths = load_paths(args.project_root)
    nmf_out_dir = resolve(paths, "outputs.nmf")
    val_dir     = resolve(paths, "outputs.validation")

For cross-repo discovery a script imports the canonical resolver in addition:

    sys.path.insert(0, "/.../iHBCA_publication/publication/config")
    from load_paths import resolve_path as canonical_resolve_path

In practice the wrappers pass canonical artifacts as CLI arguments via
${SOURCE_*} shell variables, so most scripts only need the module-local
resolver.
"""
from __future__ import annotations

import os
from typing import Any, Dict

import yaml


def load_paths(project_root: str | None = None,
               config_path: str | None = None) -> Dict[str, Any]:
    """Load motifs/config/paths.yaml and return a dict with pipeline_root set.

    project_root: absolute path to motifs pipeline_root (publication/analysis/
                  spatial/motifs/). If None, inferred from the script's location.
    config_path:  absolute path to paths.yaml. If None, resolved as
                  ${project_root}/config/paths.yaml.

    Unlike the dev pipeline's resolver, this one carries no environment
    detection or roots/containers/bind-mount logic — those live in the
    canonical publication/config/load_paths.py and are accessed by sourcing
    the publication shell resolver upstream of any wrapper.
    """
    if project_root is None:
        here = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(here)
    project_root = os.path.abspath(project_root)

    if config_path is None:
        config_path = os.path.join(project_root, "config", "paths.yaml")
    with open(config_path) as fh:
        raw = yaml.safe_load(fh) or {}

    paths: Dict[str, Any] = dict(raw)
    paths["pipeline_root"] = project_root
    return paths


def resolve(paths: Dict[str, Any], dotted_key: str) -> str:
    """Resolve a dotted key under pipeline.* into an absolute path.

    Walks paths["pipeline"][...] and joins relative entries onto
    pipeline_root. Absolute entries pass through unchanged.
    """
    pipeline_root = paths["pipeline_root"]
    node: Any = paths.get("pipeline", {})
    for part in dotted_key.split("."):
        if not isinstance(node, dict) or part not in node:
            raise KeyError(f"paths.yaml has no key: pipeline.{dotted_key}")
        node = node[part]
    if not isinstance(node, str):
        raise TypeError(f"pipeline.{dotted_key} is not a string entry")
    if node.startswith("/"):
        return node
    return os.path.join(pipeline_root, node)
