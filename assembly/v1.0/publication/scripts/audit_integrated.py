#!/usr/bin/env python3
"""Structural audit of an integrated h5ad object.

Produces a YAML report covering shape, obs/var/obsm/uns/layers/raw with
full column-level detail. Designed for the iHBCA integrated object (~2.12M
cells, ~95GB on disk).

Usage:
    python3 audit_integrated.py path/to/all-breast-cells.h5ad [--output report.yaml]

Without --output, YAML is written to stdout.
"""

import argparse
import os
import sys
from collections import OrderedDict
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml


# ---------------------------------------------------------------------------
# YAML helpers -- convert numpy/pandas types to native Python for yaml.dump
# ---------------------------------------------------------------------------

def _native(obj):
    """Recursively convert numpy/pandas scalars to native Python types."""
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        v = float(obj)
        if np.isnan(v):
            return "NaN"
        if np.isinf(v):
            return "Inf" if v > 0 else "-Inf"
        return v
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, np.ndarray):
        return [_native(x) for x in obj.tolist()]
    if isinstance(obj, (pd.Timestamp, datetime)):
        return obj.isoformat()
    if isinstance(obj, bytes):
        return obj.decode("utf-8", errors="replace")
    if isinstance(obj, dict):
        return {_native(k): _native(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_native(x) for x in obj]
    if isinstance(obj, OrderedDict):
        return OrderedDict((_native(k), _native(v)) for k, v in obj.items())
    return obj


class _OrderedDumper(yaml.SafeDumper):
    """Dump OrderedDicts preserving insertion order."""
    pass

def _ordered_representer(dumper, data):
    return dumper.represent_mapping("tag:yaml.org,2002:map", data.items())

_OrderedDumper.add_representer(OrderedDict, _ordered_representer)


def _log(msg: str):
    """Progress message to stderr."""
    print(f"[audit] {msg}", file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------
# Audit sections
# ---------------------------------------------------------------------------

def audit_shape(adata, file_path: str) -> OrderedDict:
    """Section 1: Shape & Size."""
    info = OrderedDict()
    info["n_obs"] = int(adata.n_obs)
    info["n_vars"] = int(adata.n_vars)

    try:
        info["file_size_bytes"] = os.path.getsize(file_path)
        info["file_size_human"] = _human_size(os.path.getsize(file_path))
    except OSError:
        info["file_size_bytes"] = None
        info["file_size_human"] = None

    X = adata.X
    x_info = OrderedDict()
    if X is not None:
        x_info["dtype"] = str(X.dtype)
        if sp.issparse(X):
            x_info["format"] = type(X).__name__
            x_info["nnz"] = int(X.nnz)
            x_info["density"] = round(X.nnz / (X.shape[0] * X.shape[1]), 6)
        else:
            x_info["format"] = "dense"
            x_info["nnz"] = int(np.count_nonzero(X))
    else:
        x_info["format"] = None
    info["X"] = x_info
    return info


def audit_obs(adata) -> list:
    """Section 2: obs audit -- every column, sorted alphabetically."""
    results = []
    obs = adata.obs
    for col_name in sorted(obs.columns):
        _log(f"  obs column: {col_name}")
        col = obs[col_name]
        entry = OrderedDict()
        entry["column"] = col_name
        entry["dtype"] = str(col.dtype)
        entry["n_unique"] = int(col.nunique())

        # Null / NaN counts
        n_null = int(col.isna().sum())
        entry["n_null"] = n_null
        entry["pct_null"] = round(100.0 * n_null / len(col), 2) if len(col) > 0 else 0.0

        # Top 5 value counts
        try:
            vc = col.value_counts(dropna=False).head(5)
            entry["top_values"] = [
                {"value": _native(str(idx)), "count": int(cnt)}
                for idx, cnt in vc.items()
            ]
        except Exception:
            entry["top_values"] = []

        # Categorical: list all categories
        if hasattr(col, "cat"):
            cats = col.cat.categories.tolist()
            entry["n_categories"] = len(cats)
            entry["categories"] = [_native(c) for c in cats]

        results.append(entry)
    return results


def audit_var(adata) -> OrderedDict:
    """Section 3: var audit."""
    info = OrderedDict()
    var = adata.var
    info["index_name"] = var.index.name
    info["first_5_index"] = [str(v) for v in var.index[:5].tolist()]
    info["last_5_index"] = [str(v) for v in var.index[-5:].tolist()]

    columns = []
    for col_name in sorted(var.columns):
        _log(f"  var column: {col_name}")
        col = var[col_name]
        entry = OrderedDict()
        entry["column"] = col_name
        entry["dtype"] = str(col.dtype)
        entry["n_unique"] = int(col.nunique())
        n_null = int(col.isna().sum())
        entry["n_null"] = n_null
        entry["pct_null"] = round(100.0 * n_null / len(col), 2) if len(col) > 0 else 0.0

        if hasattr(col, "cat"):
            cats = col.cat.categories.tolist()
            entry["n_categories"] = len(cats)
            entry["categories"] = [_native(c) for c in cats]

        columns.append(entry)

    info["columns"] = columns
    return info


def audit_obsm(adata) -> list:
    """Section 4: obsm audit."""
    results = []
    for key in sorted(adata.obsm.keys()):
        _log(f"  obsm key: {key}")
        arr = adata.obsm[key]
        entry = OrderedDict()
        entry["key"] = key
        entry["shape"] = list(arr.shape)
        entry["dtype"] = str(arr.dtype)

        # NaN check
        try:
            if sp.issparse(arr):
                n_nan = int(np.isnan(arr.data).sum())
            else:
                n_nan = int(np.isnan(arr).sum())
        except (TypeError, ValueError):
            n_nan = 0
        entry["n_nan"] = n_nan

        # Stats on first dimension (column 0)
        try:
            if sp.issparse(arr):
                col0 = np.asarray(arr[:, 0].todense()).ravel()
            else:
                col0 = np.asarray(arr[:, 0]).ravel()
            finite = col0[np.isfinite(col0)]
            if len(finite) > 0:
                entry["dim0_min"] = float(np.min(finite))
                entry["dim0_max"] = float(np.max(finite))
                entry["dim0_mean"] = float(np.mean(finite))
            else:
                entry["dim0_min"] = None
                entry["dim0_max"] = None
                entry["dim0_mean"] = None
        except Exception:
            entry["dim0_min"] = None
            entry["dim0_max"] = None
            entry["dim0_mean"] = None

        results.append(entry)
    return results


def audit_uns(adata) -> list:
    """Section 5: uns audit."""
    results = []
    for key in sorted(adata.uns.keys()):
        _log(f"  uns key: {key}")
        val = adata.uns[key]
        entry = OrderedDict()
        entry["key"] = key
        entry["type"] = type(val).__name__

        if isinstance(val, dict):
            sub = OrderedDict()
            for sk, sv in sorted(val.items()):
                sv_str = str(sv)
                if len(sv_str) > 100:
                    sv_str = sv_str[:100] + "..."
                sub[str(sk)] = {"type": type(sv).__name__, "value": sv_str}
            entry["contents"] = sub
        elif isinstance(val, str):
            entry["value"] = val[:100] + ("..." if len(val) > 100 else "")
        elif isinstance(val, np.ndarray):
            entry["shape"] = list(val.shape)
            entry["dtype"] = str(val.dtype)
        elif isinstance(val, pd.DataFrame):
            entry["shape"] = list(val.shape)
            entry["columns"] = list(val.columns)
        elif isinstance(val, (list, tuple)):
            entry["length"] = len(val)
            if len(val) > 0:
                entry["first_element_type"] = type(val[0]).__name__
        else:
            v_str = str(val)
            if len(v_str) > 100:
                v_str = v_str[:100] + "..."
            entry["value"] = v_str

        results.append(entry)
    return results


def audit_layers(adata) -> list:
    """Section 6: layers audit."""
    results = []
    for key in sorted(adata.layers.keys()):
        _log(f"  layer key: {key}")
        layer = adata.layers[key]
        entry = OrderedDict()
        entry["key"] = key
        entry["shape"] = list(layer.shape)
        entry["dtype"] = str(layer.dtype)

        if sp.issparse(layer):
            entry["format"] = type(layer).__name__
            entry["nnz"] = int(layer.nnz)
            nz_data = layer.data
            if len(nz_data) > 0:
                entry["nonzero_min"] = float(np.min(nz_data))
                entry["nonzero_max"] = float(np.max(nz_data))
                entry["nonzero_mean"] = float(np.mean(nz_data))
            else:
                entry["nonzero_min"] = None
                entry["nonzero_max"] = None
                entry["nonzero_mean"] = None
        else:
            entry["format"] = "dense"
            entry["nnz"] = int(np.count_nonzero(layer))
            nz = layer[layer != 0]
            if len(nz) > 0:
                entry["nonzero_min"] = float(np.min(nz))
                entry["nonzero_max"] = float(np.max(nz))
                entry["nonzero_mean"] = float(np.mean(nz))
            else:
                entry["nonzero_min"] = None
                entry["nonzero_max"] = None
                entry["nonzero_mean"] = None

        results.append(entry)
    return results


def audit_raw(adata) -> OrderedDict:
    """Section 7: raw audit."""
    info = OrderedDict()
    if adata.raw is None:
        info["exists"] = False
        return info

    info["exists"] = True
    info["shape"] = list(adata.raw.X.shape)

    raw_var = adata.raw.var
    info["var_columns"] = sorted(raw_var.columns.tolist())
    info["var_index_name"] = raw_var.index.name
    info["n_var"] = int(raw_var.shape[0])

    X = adata.raw.X
    if X is not None:
        info["X_dtype"] = str(X.dtype)
        if sp.issparse(X):
            info["X_format"] = type(X).__name__
            info["X_nnz"] = int(X.nnz)
        else:
            info["X_format"] = "dense"

    return info


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------

def _human_size(nbytes: int) -> str:
    size = float(nbytes)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(size) < 1024.0:
            return f"{size:.1f} {unit}"
        size /= 1024.0
    return f"{size:.1f} PB"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Structural audit of an h5ad file. Outputs YAML report."
    )
    parser.add_argument("h5ad", help="Path to .h5ad file")
    parser.add_argument(
        "--output", "-o", default=None,
        help="Output YAML path (default: stdout)"
    )
    args = parser.parse_args()

    h5ad_path = str(Path(args.h5ad).resolve())

    report = OrderedDict()
    report["audit_metadata"] = OrderedDict([
        ("timestamp", datetime.now().isoformat()),
        ("file_path", h5ad_path),
        ("script", "audit_integrated.py"),
    ])

    # Load
    _log(f"Loading {h5ad_path} ...")
    import scanpy as sc
    adata = sc.read_h5ad(h5ad_path)
    _log(f"Loaded: {adata.n_obs} obs x {adata.n_vars} vars")

    # 1. Shape
    _log("Auditing shape & size...")
    report["shape_and_size"] = audit_shape(adata, h5ad_path)

    # 2. obs
    _log(f"Auditing obs ({len(adata.obs.columns)} columns)...")
    report["obs_audit"] = audit_obs(adata)

    # 3. var
    _log(f"Auditing var ({len(adata.var.columns)} columns)...")
    report["var_audit"] = audit_var(adata)

    # 4. obsm
    _log(f"Auditing obsm ({len(adata.obsm.keys())} keys)...")
    report["obsm_audit"] = audit_obsm(adata)

    # 5. uns
    _log(f"Auditing uns ({len(adata.uns.keys())} keys)...")
    report["uns_audit"] = audit_uns(adata)

    # 6. layers
    _log(f"Auditing layers ({len(adata.layers.keys())} keys)...")
    report["layers_audit"] = audit_layers(adata)

    # 7. raw
    _log("Auditing raw...")
    report["raw_audit"] = audit_raw(adata)

    _log("Serializing YAML...")
    report = _native(report)
    yaml_str = yaml.dump(
        report, Dumper=_OrderedDumper,
        default_flow_style=False, allow_unicode=True, width=120
    )

    if args.output:
        out_path = Path(args.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(yaml_str, encoding="utf-8")
        _log(f"Report written to {out_path}")
    else:
        sys.stdout.write(yaml_str)

    _log("Done.")


if __name__ == "__main__":
    main()
