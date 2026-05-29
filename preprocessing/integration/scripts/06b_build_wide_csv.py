#!/usr/bin/env python3
"""
Build the scIB wide aggregation CSV from per-metric sidecars.

Joins the decomposed per-metric sidecar CSVs into a single wide row per
(compartment, method, n_latent). The wide table is the input to the SVD
diagnostic and composite scoring downstream of the scib suite.

Sidecar filename map:
    {method}_n{nlat}_{short}_nmi_ari.csv    <- 05e_scib_nmi_ari.py
    {method}_n{nlat}_{short}_asw_label.csv  <- 05f_scib_asw_label.py
    {method}_n{nlat}_{short}_asw_batch.csv  <- 05g_scib_asw_batch.py
    {method}_n{nlat}_{short}_graph_conn.csv <- 05h_scib_graph_conn.py
    {method}_n{nlat}_{short}_kbet.csv       <- 05b_scib_kbet_only.py
    {method}_n{nlat}_{short}_pcr.csv        <- 05d_scib_pcr_recompute.py

Output schema (9 columns, exact order):
    compartment, config, NMI, ARI, ASW_label, graph_conn, ASW_batch, kBET, PCR

Column source map:
    NMI         <- nmi_ari sidecar['NMI']
    ARI         <- nmi_ari sidecar['ARI']
    ASW_label   <- asw_label sidecar['ASW_label']
    graph_conn  <- graph_conn sidecar['graph_conn']
    ASW_batch   <- asw_batch sidecar['ASW_batch']
    kBET        <- kbet sidecar['kBET']
    PCR         <- pcr sidecar['pcr_scaled']

Missing sidecars are tolerated: NaN fills the corresponding wide-table cell
without blocking the merge. This preserves failure isolation — one metric's
SLURM task failing does not invalidate the whole aggregation.

All values are formatted to 4 decimal places via pandas float_format='%.4f'.

Canonical row sort order:
    1. compartment: Immune, Epithelial, Stromal
    2. method family: scvi before scanvi
    3. n_latent ascending: 50, 75, 100

Usage:
    python 06b_build_wide_csv.py \\
        --scib-dir <scib_output_dir> \\
        --output   <out.csv>

    # Filter to specific compartments (e.g., for incremental builds):
    python 06b_build_wide_csv.py \\
        --scib-dir   <scib_output_dir> \\
        --output     <out.csv> \\
        --compartments imm

    # Regression check against a frozen reference wide CSV:
    python 06b_build_wide_csv.py \\
        --scib-dir         <scib_output_dir> \\
        --output           /tmp/imm_check.csv \\
        --compartments     imm \\
        --regression-check <reference_wide.csv>
"""

import argparse
import math
import re
import sys
from pathlib import Path

import pandas as pd

# Canonical compartment order (matches run_scib_*.sh COMPARTMENTS array)
COMP_ORDER = {"Immune": 0, "Epithelial": 1, "Stromal": 2}
METHOD_ORDER = {"scvi": 0, "scanvi": 1}
SHORT_TO_LONG = {"imm": "Immune", "epi": "Epithelial", "str": "Stromal"}

OUTPUT_COLS = ["compartment", "config", "NMI", "ARI", "ASW_label",
               "graph_conn", "ASW_batch", "kBET", "PCR"]

# Sidecar suffix → (filename suffix, source column, wide column) triples
SIDECARS = [
    ("nmi_ari",    "NMI",        "NMI"),
    ("nmi_ari",    "ARI",        "ARI"),
    ("asw_label",  "ASW_label",  "ASW_label"),
    ("asw_batch",  "ASW_batch",  "ASW_batch"),
    ("graph_conn", "graph_conn", "graph_conn"),
    ("kbet",       "kBET",       "kBET"),
    ("pcr",        "pcr_scaled", "PCR"),
]

SIDECAR_SUFFIXES = ["nmi_ari", "asw_label", "asw_batch",
                    "graph_conn", "kbet", "pcr"]

# Match any sidecar filename to discover available configs
ANY_SIDECAR_PATTERN = re.compile(
    r"^(scvi|scanvi)_n(\d+)_(imm|epi|str)_"
    r"(nmi_ari|asw_label|asw_batch|graph_conn|kbet|pcr)\.csv$"
)


def discover_configs(scib_dir: Path) -> list[tuple[str, int, str]]:
    """Find every (method, n_latent, short) with at least one sidecar CSV."""
    seen: set[tuple[str, int, str]] = set()
    for p in sorted(scib_dir.iterdir()):
        m = ANY_SIDECAR_PATTERN.match(p.name)
        if not m:
            continue
        seen.add((m.group(1), int(m.group(2)), m.group(3)))
    return sorted(seen)


def read_sidecar_value(scib_dir: Path, method: str, n_latent: int,
                       short: str, suffix: str, src_col: str) -> float:
    """Read a single metric value from a sidecar CSV. NaN if missing/unreadable."""
    path = scib_dir / f"{method}_n{n_latent}_{short}_{suffix}.csv"
    if not path.exists():
        return float("nan")
    try:
        df = pd.read_csv(path)
    except Exception as e:
        print(f"  WARN: failed to read {path.name}: {e}", file=sys.stderr)
        return float("nan")
    if df.empty or src_col not in df.columns:
        return float("nan")
    try:
        return float(df.iloc[0][src_col])
    except (TypeError, ValueError):
        return float("nan")


def build_row(scib_dir: Path, method: str, n_latent: int, short: str) -> dict:
    """Pull every metric from its sidecar; missing values become NaN."""
    row: dict = {
        "compartment": SHORT_TO_LONG[short],
        "config": f"{method}_n{n_latent}",
    }
    for suffix, src_col, wide_col in SIDECARS:
        row[wide_col] = read_sidecar_value(
            scib_dir, method, n_latent, short, suffix, src_col)
    return row


def completeness_status(row: dict) -> str:
    """Return OK or MISSING <metrics> for a wide row."""
    metrics = [c for c in OUTPUT_COLS if c not in ("compartment", "config")]
    missing = [m for m in metrics
               if isinstance(row.get(m), float) and math.isnan(row[m])]
    return "OK" if not missing else f"MISSING {','.join(missing)}"


def regression_check(wide: pd.DataFrame, ref_path: Path,
                     tolerance: float) -> None:
    """Print per-config drift against a reference wide CSV (Immune only)."""
    if not ref_path.exists():
        print(f"  ref not found: {ref_path} (skipping)", file=sys.stderr)
        return
    ref: pd.DataFrame = pd.read_csv(ref_path)
    ref_imm: pd.DataFrame = ref.loc[ref["compartment"] == "Immune"].copy()
    new_imm: pd.DataFrame = wide.loc[wide["compartment"] == "Immune"].copy()
    metrics = [c for c in OUTPUT_COLS if c not in ("compartment", "config")]

    print(f"\nRegression check vs {ref_path.name} (tol={tolerance}):",
          file=sys.stderr)
    new_by_config: dict[str, dict] = {
        str(row["config"]): {m: row[m] for m in metrics if m in row.index}
        for _, row in new_imm.iterrows()
    }
    for _, ref_row in ref_imm.iterrows():
        config_val = str(ref_row["config"])
        new_row = new_by_config.get(config_val)
        if new_row is None:
            print(f"  {config_val:<12s}  MISSING in new build",
                  file=sys.stderr)
            continue
        bad: list[str] = []
        for m in metrics:
            if m not in ref_row.index or m not in new_row:
                continue
            try:
                rf = float(ref_row[m])  # type: ignore[arg-type]
                nf = float(new_row[m])  # type: ignore[arg-type]
            except (TypeError, ValueError):
                continue
            if math.isnan(rf) or math.isnan(nf):
                continue
            if abs(rf - nf) > tolerance:
                bad.append(f"{m}:{rf:.4f}→{nf:.4f}")
        status = "OK" if not bad else f"DRIFT {'; '.join(bad)}"
        print(f"  {config_val:<12s}  {status}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--scib-dir", required=True,
                        help="Directory with per-metric sidecar CSVs")
    parser.add_argument("--output", required=True,
                        help="Output wide CSV path")
    parser.add_argument("--compartments", default="imm,epi,str",
                        help="Comma-separated short compartment names "
                             "(default: imm,epi,str)")
    parser.add_argument("--regression-check", default=None,
                        help="Path to an existing wide CSV to cross-check "
                             "immune rows against (e.g., imm_aggregated_wide.csv)")
    parser.add_argument("--tolerance", type=float, default=1e-4,
                        help="Numeric tolerance for regression check "
                             "(default 1e-4; matches the 4-decimal output)")
    args = parser.parse_args()

    scib_dir = Path(args.scib_dir)
    if not scib_dir.exists():
        sys.exit(f"scib_dir does not exist: {scib_dir}")
    keep_shorts = set(s.strip() for s in args.compartments.split(","))

    all_configs = discover_configs(scib_dir)
    configs = [c for c in all_configs if c[2] in keep_shorts]
    if not configs:
        sys.exit("no sidecar CSVs matched the expected pattern in "
                 f"{scib_dir} for compartments {sorted(keep_shorts)}")

    print(f"Discovered {len(configs)} configs with at least one sidecar:",
          file=sys.stderr)
    for method, n_latent, short in configs:
        present = []
        for suffix in SIDECAR_SUFFIXES:
            fp = scib_dir / f"{method}_n{n_latent}_{short}_{suffix}.csv"
            if fp.exists():
                present.append(suffix)
        print(f"  {method}_n{n_latent}_{short}: [{', '.join(present)}]",
              file=sys.stderr)

    rows = [build_row(scib_dir, *c) for c in configs]

    def sort_key(r):
        method_part, n_part = r["config"].split("_n")
        return (COMP_ORDER[r["compartment"]],
                METHOD_ORDER[method_part],
                int(n_part))
    rows.sort(key=sort_key)

    df = pd.DataFrame(rows, columns=OUTPUT_COLS)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False, float_format="%.4f")
    print(f"\nWrote {out_path} ({len(df)} rows)", file=sys.stderr)

    print("\nCompleteness per (compartment, config):", file=sys.stderr)
    for r in rows:
        print(f"  {r['compartment']:>11s}  {r['config']:<12s}  "
              f"{completeness_status(r)}", file=sys.stderr)

    if args.regression_check:
        regression_check(df, Path(args.regression_check), args.tolerance)


if __name__ == "__main__":
    main()
