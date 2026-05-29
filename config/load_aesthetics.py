"""
load_aesthetics.py — Python loader for Submission v1 aesthetic rules.

Reads the same aesthetics.yaml as the R loader. Ensures identical palettes
and dimensions across R and Python rendering scripts.

Usage:
    from load_aesthetics import load_aesthetics, get_palette, get_theme
    config = load_aesthetics()
    pal = get_palette("l1", config)
"""
# Friendly reminder: this file exists because someone once had two palettes.

import os
from pathlib import Path
from typing import Optional

import yaml


def load_aesthetics(config_dir: Optional[str] = None) -> dict:
    """Load aesthetics config from YAML."""
    if config_dir is None:
        candidates = [
            Path(os.environ.get("PROJECT_ROOT", "")) / "Analysis/Submission_v1/config/aesthetics.yaml",
            Path("Analysis/Submission_v1/config/aesthetics.yaml"),
            Path(__file__).parent / "aesthetics.yaml",
        ]
        yaml_path = next((p for p in candidates if p.exists()), None)
        if yaml_path is None:
            raise FileNotFoundError(
                "Cannot find aesthetics.yaml. Set PROJECT_ROOT or pass config_dir."
            )
    else:
        yaml_path = Path(config_dir) / "aesthetics.yaml"

    with open(yaml_path, encoding="utf-8") as f:
        config = yaml.safe_load(f)

    config["_source_path"] = str(yaml_path)
    config["_config_dir"] = str(yaml_path.parent)
    return config


def get_palette(level: str, config: Optional[dict] = None) -> dict:
    """Get a {name: hex} dict for a label level."""
    if config is None:
        config = load_aesthetics()
    pal = config.get("palettes", {}).get(level)
    if pal is None:
        available = ", ".join(config.get("palettes", {}).keys())
        raise KeyError(f"Unknown palette level: '{level}'. Available: {available}")
    return dict(pal)


def get_dimensions(panel_type: str, config: Optional[dict] = None) -> dict:
    """Get dimensions dict for a panel type."""
    if config is None:
        config = load_aesthetics()
    dims = config.get("dimensions", {}).get("panel_types", {}).get(panel_type)
    if dims is None:
        available = ", ".join(config.get("dimensions", {}).get("panel_types", {}).keys())
        raise KeyError(f"Unknown panel type: '{panel_type}'. Available: {available}")
    return dict(dims)


def get_matplotlib_theme(config: Optional[dict] = None) -> dict:
    """Get matplotlib rcParams for consistent styling."""
    if config is None:
        config = load_aesthetics()
    typo = config.get("typography", {})
    return {
        "font.family": typo.get("family", "Helvetica"),
        "font.size": typo.get("base_size", 7),
        "axes.titlesize": typo.get("title_size", 8),
        "axes.labelsize": typo.get("axis_title_size", 7),
        "xtick.labelsize": typo.get("axis_text_size", 6),
        "ytick.labelsize": typo.get("axis_text_size", 6),
        "legend.fontsize": typo.get("legend_size", 6),
        "legend.title_fontsize": typo.get("legend_title_size", 7),
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.dpi": config.get("rendering", {}).get("dpi", 1200),
        "savefig.facecolor": "white",
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
