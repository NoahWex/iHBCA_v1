"""
load_paths.py - Python sibling of lib/load_paths.R.

Pipeline scripts import this and call `load_inquiry(inquiry_dir)` to obtain a
dict-of-dicts with `inputs`, `outputs`, `framework`, `inquiry_root`, etc.
Mirrors the R implementation: validates that paths.yaml exists, resolves
relative output paths against inquiry_root, leaves inputs absolute, and
permits dissoc_table to be relative to <inquiry>/config/.

Usage:
    from lib.load_paths import load_inquiry
    cfg = load_inquiry(inquiry_dir)
    labels_path = cfg["paths"]["inputs"]["labels"]
    components_root = cfg["paths"]["inputs"]["components_root"]
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Any

try:
    import yaml
except ImportError:
    yaml = None  # type: ignore[assignment]


def _load_yaml(p: Path) -> dict[str, Any]:
    if yaml is None:
        # Minimal fallback if PyYAML isn't installed: not robust, but enough
        # for the simple maps we use. Prefer installing PyYAML.
        raise ImportError("PyYAML required for load_paths.py")
    with open(p) as f:
        return yaml.safe_load(f)


def load_inquiry(inquiry_dir: str | os.PathLike) -> dict[str, Any]:
    inquiry_dir = Path(inquiry_dir).resolve()
    paths_file = inquiry_dir / "config" / "paths.yaml"
    inquiry_file = inquiry_dir / "config" / "inquiry.yaml"
    if not paths_file.exists():
        raise FileNotFoundError(f"paths.yaml missing at: {paths_file}")
    paths = _load_yaml(paths_file) or {}
    inquiry = _load_yaml(inquiry_file) if inquiry_file.exists() else {}

    # Resolve dissoc_table relative to config/ if present and not absolute
    dissoc = paths.get("inputs", {}).get("dissoc_table")
    if dissoc and not str(dissoc).startswith("/"):
        paths["inputs"]["dissoc_table"] = str(inquiry_dir / "config" / dissoc)

    # Resolve output paths relative to inquiry_root
    inq_root = paths.get("inquiry_root") or str(inquiry_dir)
    outs = paths.get("outputs", {}) or {}
    for k, v in list(outs.items()):
        if not str(v).startswith("/"):
            outs[k] = str(Path(inq_root) / v)
    paths["outputs"] = outs

    return {
        "paths": paths,
        "inquiry": inquiry,
        "inquiry_dir": str(inquiry_dir),
        "inquiry_root": inq_root,
    }
