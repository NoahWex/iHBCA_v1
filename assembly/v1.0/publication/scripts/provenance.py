#!/usr/bin/env python3
"""
provenance.py - Build provenance manifests for assembly outputs
===============================================================
Embeds provenance metadata in adata.uns["ihbca_provenance"] and optionally
writes a sidecar YAML manifest next to the output h5ad.

Used by assemble_h5ad.py (source datasets) and assemble_integrated.py
(integrated object).

"""

import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import yaml


def _git_info(repo_root):
    """Collect git commit, branch, and dirty status.

    Returns dict with {commit, branch, dirty}. Falls back to unknowns
    if git is unavailable (e.g., inside a container without git).
    """
    fallback = {"commit": "unknown", "branch": "unknown", "dirty": None}
    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=str(repo_root),
            stderr=subprocess.DEVNULL,
        ).decode().strip()

        branch = subprocess.check_output(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"],
            cwd=str(repo_root),
            stderr=subprocess.DEVNULL,
        ).decode().strip()

        dirty_rc = subprocess.call(
            ["git", "diff", "--quiet"],
            cwd=str(repo_root),
            stderr=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
        )
        dirty = dirty_rc != 0

        return {"commit": commit, "branch": branch, "dirty": dirty}
    except (FileNotFoundError, subprocess.CalledProcessError):
        return fallback


def _file_info(path):
    """Return {path, size_bytes} for a file. size_bytes is None if missing."""
    p = Path(path)
    size = p.stat().st_size if p.exists() else None
    return {"path": str(p), "size_bytes": size}


def build_provenance_manifest(
    adata,
    output_path,
    study,
    repo_root,
    inputs,
    config_files,
    build_script=None,
    extra=None,
    sidecar_path=None,
):
    """Embed provenance in adata.uns and optionally write sidecar YAML.

    Parameters
    ----------
    adata : anndata.AnnData
        The assembled object. uns["ihbca_provenance"] will be set.
    output_path : Path
        Path to the output h5ad (used for metadata, not written here).
    study : str
        Study name (e.g., "gray") or "integrated".
    repo_root : Path
        Repository root for git info resolution.
    inputs : dict
        {label: path_str_or_Path} — input files used in assembly.
        Each is resolved to {path, size_bytes} via _file_info().
    config_files : dict
        {label: path_str_or_Path} — config/mapping files used.
        Each is resolved to {path, size_bytes} via _file_info().
    build_script : str, optional
        Name of the calling script (e.g., "assemble_h5ad.py").
        If None, defaults to "unknown".
    extra : dict, optional
        Additional keys merged into the provenance dict
        (e.g., mapping_stats, lineage_splits).
    sidecar_path : Path or None
        If not None, write provenance YAML to this path.
        Written to /tmp first, then shutil.copy() to final location (CRSP safety).

    Returns
    -------
    dict
        The provenance dict that was embedded in uns.
    """
    output_path = Path(output_path)
    git = _git_info(repo_root)

    provenance = {
        "build_script": build_script or "unknown",
        "git_hash": git["commit"],
        "git_branch": git["branch"],
        "git_dirty": git["dirty"],
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "python_version": sys.version.split()[0],
        "study": study,
        "cell_count": int(adata.n_obs),
        "gene_count": int(adata.n_vars),
        "inputs": {label: _file_info(p) for label, p in inputs.items()},
        "config_files": {label: _file_info(p) for label, p in config_files.items()},
    }

    if extra:
        provenance.update(extra)

    adata.uns["ihbca_provenance"] = provenance

    if sidecar_path is not None:
        sidecar_path = Path(sidecar_path)
        sidecar_path.parent.mkdir(parents=True, exist_ok=True)
        # Write to /tmp first, then copy (CRSP stale file handle safety)
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".yaml", delete=False
        ) as tmp:
            yaml.dump(provenance, tmp, default_flow_style=False, sort_keys=False)
            tmp_path = tmp.name
        shutil.copy(tmp_path, str(sidecar_path))
        os.unlink(tmp_path)
        print(f"  Provenance sidecar: {sidecar_path}")

    return provenance
