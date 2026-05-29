#!/bin/bash
# ============================================================================
# Aggregate Step 01 -> central_cell_status.csv
# ============================================================================
# Purpose
#   Combine per-sample temp CSVs from step 01 into the initial
#   central_cell_status.csv (schema: cell_id, sample_id, patient_id,
#   position, step_01_umi_pass, step_01_n_umi). This is the single
#   source of truth for cell membership across the pipeline; every
#   downstream step reads this manifest instead of walking sample_manifests.
#
# Atomic write
#   The central manifest is written to a temp path and os.replace()'d
#   into place so concurrent readers never see a partial file.
#
# Fail-loud
#   - Required columns missing -> exit 2
#   - Duplicate cell_ids -> exit 2
#   - Sample count mismatch vs raw_data_manifest -> exit 2
# ============================================================================

set -euo pipefail

CONTAINER_TYPE=scgpt_gpu
if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${PIPELINE_ROOT}/run/_common.sh"

TEMP_CSV_DIR="${CFG_PREPROCESSING_QC_STATUS}/temp"
CENTRAL_MANIFEST="${CFG_PREPROCESSING_CENTRAL_CELL_STATUS}"

mkdir -p "$(dirname "${CENTRAL_MANIFEST}")"

echo "Step:            aggregate_step_01"
echo "Temp dir:        ${TEMP_CSV_DIR}"
echo "Central manifest: ${CENTRAL_MANIFEST}"
echo "Raw manifest:    ${CFG_RAW_DATA_SAMPLE_MANIFEST}"

export TEMP_CSV_DIR CENTRAL_MANIFEST RAW_MANIFEST="${CFG_RAW_DATA_SAMPLE_MANIFEST}"

run_python_singularity /dev/stdin <<'PYEOF'
import os
import sys
from datetime import datetime
from pathlib import Path
import pandas as pd

temp_dir = Path(os.environ["TEMP_CSV_DIR"])
central = Path(os.environ["CENTRAL_MANIFEST"])
raw_manifest = Path(os.environ["RAW_MANIFEST"])

if not temp_dir.is_dir():
    print(f"FATAL: temp dir not found: {temp_dir}", file=sys.stderr)
    sys.exit(2)

csvs = sorted(temp_dir.glob("*_step01.csv"))
if not csvs:
    print(f"FATAL: no *_step01.csv in {temp_dir}", file=sys.stderr)
    sys.exit(2)
print(f"Found {len(csvs)} per-sample CSVs")

dfs = [pd.read_csv(p) for p in csvs]
df = pd.concat(dfs, ignore_index=True)

required = ["cell_id", "sample_id", "patient_id", "position",
            "step_01_umi_pass", "step_01_n_umi"]
missing = set(required) - set(df.columns)
if missing:
    print(f"FATAL: missing columns: {sorted(missing)}", file=sys.stderr)
    sys.exit(2)

dup = df["cell_id"].duplicated().sum()
if dup:
    print(f"FATAL: {dup} duplicate cell_ids in aggregated CSVs", file=sys.stderr)
    sys.exit(2)

raw_df = pd.read_csv(raw_manifest, sep="\t", dtype=str, keep_default_na=False)
expected = set(raw_df["sample_id"].unique())
seen = set(df["sample_id"].unique())
unexpected = seen - expected
if unexpected:
    print(f"FATAL: unexpected samples not in raw_data_manifest: {sorted(unexpected)[:5]}",
          file=sys.stderr)
    sys.exit(2)
missing_samples = expected - seen
if missing_samples:
    print(f"NOTE: {len(missing_samples)} samples in raw_data_manifest have no "
          f"step_01 temp CSV (expected for Pat2_P1-style missing-chromium positions): "
          f"{sorted(missing_samples)[:5]}")

df = df.sort_values(by=["patient_id", "position", "sample_id", "cell_id"]).reset_index(drop=True)

central.parent.mkdir(parents=True, exist_ok=True)
if central.exists():
    backup = central.with_suffix(central.suffix + f".bak.{datetime.now().strftime('%Y%m%d_%H%M%S')}")
    central.replace(backup)
    print(f"Existing manifest backed up: {backup}")

tmp = central.with_suffix(central.suffix + f".tmp.{os.getpid()}")
df.to_csv(tmp, index=False)
tmp.replace(central)
print(f"Wrote {central} ({len(df):,} cells)")
PYEOF

echo "Exit: $?, End: $(date)"
