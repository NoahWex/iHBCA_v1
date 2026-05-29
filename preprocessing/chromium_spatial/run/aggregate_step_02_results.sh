#!/bin/bash
# ============================================================================
# Aggregate Step 02 -> extend central_cell_status.csv
# ============================================================================
# Adds columns: step_02_mad_pass, step_02_exclusion_flags.
#
# Inputs:  per-patient update CSVs from
#          ${CFG_PREPROCESSING_STEP_02}/manifest_updates/{patient}_step_02_update.csv
# Output:  updated ${CFG_PREPROCESSING_CENTRAL_CELL_STATUS} (atomic write)
#
# Fail-loud:
#   - any expected patient update missing -> exit 2
#   - any update cell_id not in central manifest -> exit 2
#   - schema mismatch -> exit 2
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

UPDATE_DIR="${CFG_PREPROCESSING_STEP_02}/manifest_updates"
CENTRAL_MANIFEST="${CFG_PREPROCESSING_CENTRAL_CELL_STATUS}"

echo "Step:            aggregate_step_02"
echo "Update dir:      ${UPDATE_DIR}"
echo "Central manifest: ${CENTRAL_MANIFEST}"

export UPDATE_DIR CENTRAL_MANIFEST

run_python_singularity /dev/stdin <<'PYEOF'
import os
import sys
import shutil
from datetime import datetime
from pathlib import Path
import pandas as pd

update_dir = Path(os.environ["UPDATE_DIR"])
central_path = Path(os.environ["CENTRAL_MANIFEST"])

if not update_dir.is_dir():
    print(f"FATAL: update dir missing: {update_dir}", file=sys.stderr)
    sys.exit(2)
if not central_path.exists():
    print(f"FATAL: central manifest missing: {central_path}", file=sys.stderr)
    sys.exit(2)

csvs = sorted(update_dir.glob("*_step_02_update.csv"))
if not csvs:
    print(f"FATAL: no *_step_02_update.csv files in {update_dir}", file=sys.stderr)
    sys.exit(2)
print(f"Found {len(csvs)} per-patient update CSVs")

dfs = [pd.read_csv(p) for p in csvs]
updates = pd.concat(dfs, ignore_index=True)
required = {"cell_id", "sample_id", "patient_id", "step_02_mad_pass", "step_02_exclusion_flags"}
missing = required - set(updates.columns)
if missing:
    print(f"FATAL: missing columns in updates: {sorted(missing)}", file=sys.stderr)
    sys.exit(2)

if updates["cell_id"].duplicated().any():
    print("FATAL: duplicate cell_id in updates", file=sys.stderr)
    sys.exit(2)

central = pd.read_csv(central_path)
missing_cells = set(updates["cell_id"]) - set(central["cell_id"])
if missing_cells:
    print(f"FATAL: {len(missing_cells)} cell_ids in updates not present in central manifest",
          file=sys.stderr)
    print(f"examples: {sorted(missing_cells)[:5]}", file=sys.stderr)
    sys.exit(2)

central = central.set_index("cell_id")
updates_idx = updates.set_index("cell_id")
central.loc[updates_idx.index, "step_02_mad_pass"] = updates_idx["step_02_mad_pass"]
central.loc[updates_idx.index, "step_02_exclusion_flags"] = updates_idx["step_02_exclusion_flags"]
central = central.reset_index()

backup = central_path.with_suffix(central_path.suffix + f".bak.{datetime.now().strftime('%Y%m%d_%H%M%S')}")
shutil.copyfile(central_path, backup)
tmp = central_path.with_suffix(central_path.suffix + f".tmp.{os.getpid()}")
central.to_csv(tmp, index=False)
tmp.replace(central_path)
print(f"Updated {central_path} ({len(updates):,} cells updated, backup: {backup})")
PYEOF

echo "Exit: $?, End: $(date)"
