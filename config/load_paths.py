"""
load_paths.py — Python loader for publication path resolution.

Reads paths.yaml and resolves all paths for the current environment (HPC or
local). Every Python analysis, rendering, and assembly script imports this
instead of hardcoding absolute paths.

Usage:
    sys.path.insert(0, config_dir)
    from load_paths import load_paths, resolve_path, get_container
    paths = load_paths()
    h5ad  = resolve_path("ihbca_integrated", paths)
    sif   = get_container("python_spatial_2025Q2", paths)
"""

import os
import socket
from pathlib import Path
from typing import Optional

import yaml


def detect_environment() -> str:
    """Detect whether we are running on HPC or locally."""
    env_override = os.environ.get("IHBCA_ENV", "")
    if env_override:
        return env_override.lower()
    hostname = socket.gethostname()
    if any(pat in hostname.lower() for pat in ("hpc3", "login-", "compute-")):
        return "hpc"
    return "local"


def _find_paths_yaml(config_dir: Optional[str] = None) -> Path:
    """Locate paths.yaml via explicit dir, env var, or __file__ fallback."""
    if config_dir is not None:
        yaml_path = Path(config_dir) / "paths.yaml"
        if yaml_path.exists():
            return yaml_path
        raise FileNotFoundError(f"paths.yaml not found at: {yaml_path}")

    candidates = []

    project_root = os.environ.get("IHBCA_PROJECT_ROOT", "")
    if project_root:
        candidates.append(
            Path(project_root) / "publication" / "config" / "paths.yaml"
        )

    # Relative to this file (most reliable when imported from config/)
    candidates.append(Path(__file__).parent / "paths.yaml")

    # Common working directory patterns
    candidates.extend([
        Path("publication/config/paths.yaml"),
        Path("../config/paths.yaml"),
        Path("../../config/paths.yaml"),
        Path("../../../config/paths.yaml"),
    ])

    for p in candidates:
        if p.exists():
            return p.resolve()

    raise FileNotFoundError(
        "Cannot find paths.yaml. Set IHBCA_PROJECT_ROOT or pass config_dir."
    )


def load_paths(config_dir: Optional[str] = None) -> dict:
    """Load and resolve all paths from paths.yaml.

    Parameters
    ----------
    config_dir : str, optional
        Explicit path to the config directory containing paths.yaml.

    Returns
    -------
    dict
        Nested dict with all paths resolved for the current environment.
    """
    yaml_path = _find_paths_yaml(config_dir)
    with open(yaml_path) as f:
        config = yaml.safe_load(f)

    env = detect_environment()
    root = config.get("roots", {}).get(env)
    if root is None:
        raise ValueError(f"No root defined for environment '{env}' in paths.yaml")

    # Resolve sources: absolute paths stay, relative paths get root prepended
    if "sources" in config:
        config["sources"] = {
            k: v if v.startswith("/") else str(Path(root) / v)
            for k, v in config["sources"].items()
        }

    # Resolve projects
    if "projects" in config:
        config["projects"] = {
            k: v if v.startswith("/") else str(Path(root) / v)
            for k, v in config["projects"].items()
        }

    config["_environment"] = env
    config["_root"] = root
    config["_yaml_path"] = str(yaml_path)
    config["_config_dir"] = str(yaml_path.parent)
    return config


def resolve_path(key: str, paths: Optional[dict] = None) -> str:
    """Resolve a single source path by key."""
    if paths is None:
        paths = load_paths()
    sources = paths.get("sources", {})
    if key not in sources:
        available = ", ".join(sources.keys())
        raise KeyError(f"Unknown source key: '{key}'. Available: {available}")
    return sources[key]


def get_container(name: str, paths: Optional[dict] = None) -> str:
    """Get container path by name."""
    if paths is None:
        paths = load_paths()
    containers = paths.get("containers", {})
    if name not in containers:
        available = ", ".join(containers.keys())
        raise KeyError(f"Unknown container: '{name}'. Available: {available}")
    return containers[name]


def get_bind_mounts(paths: Optional[dict] = None) -> list:
    """Get bind mount strings for singularity --bind."""
    if paths is None:
        paths = load_paths()
    return paths.get("bind_mounts", [])


def get_python_libs(name: str, paths: Optional[dict] = None) -> str:
    """Get Python user package path for a container."""
    if paths is None:
        paths = load_paths()
    plibs = paths.get("python_libs", {})
    if name not in plibs:
        raise KeyError(f"No Python library path found for: '{name}'")
    return plibs[name]


def get_slurm_defaults(paths: Optional[dict] = None) -> dict:
    """Get SLURM default settings."""
    if paths is None:
        paths = load_paths()
    return paths.get("slurm", {})
