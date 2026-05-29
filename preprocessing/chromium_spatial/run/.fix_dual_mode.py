#!/usr/bin/env python3
"""Apply CP5.7.10 dual-mode SCRIPT_DIR pattern to all wrappers using the legacy
source pattern. Idempotent: only modifies files that match the legacy block.
"""
import os
import sys
from pathlib import Path

RUN_DIR = Path(__file__).parent
# Match the standalone legacy `source` line (regardless of what precedes it).
# Replace it with the dual-mode block + new sourced path. The CONTAINER_TYPE
# line stays where it was — typically immediately above the source line.
LEGACY_SOURCE_LINE = 'source "$(dirname "$0")/_common.sh"'

NEW_DUAL_MODE_BLOCK = '''if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${PIPELINE_ROOT}/run/_common.sh"'''


def fix_file(path: Path) -> str:
    """Returns: 'fixed', 'already-applied', 'pattern-not-found', 'skipped-common'."""
    if path.name == "_common.sh":
        return "skipped-common"
    text = path.read_text()
    if 'SCRIPT_DIR="$(pwd)"' in text:
        return "already-applied"
    if LEGACY_SOURCE_LINE not in text:
        return "pattern-not-found"
    if text.count(LEGACY_SOURCE_LINE) > 1:
        return "pattern-not-found"
    new_text = text.replace(LEGACY_SOURCE_LINE, NEW_DUAL_MODE_BLOCK, 1)
    path.write_text(new_text)
    return "fixed"


def preview_file(path: Path) -> str:
    if path.name == "_common.sh":
        return "skipped-common"
    text = path.read_text()
    if 'SCRIPT_DIR="$(pwd)"' in text:
        return "already-applied"
    if LEGACY_SOURCE_LINE not in text:
        return "pattern-not-found"
    if text.count(LEGACY_SOURCE_LINE) > 1:
        return "duplicate-source-lines"
    return "would-fix"


def main():
    apply = "--apply" in sys.argv
    label = "APPLY" if apply else "DRY-RUN"
    print(f"=== {label} mode ===")
    results = {}
    for sh in sorted(RUN_DIR.glob("*.sh")):
        if apply:
            r = fix_file(sh)
        else:
            r = preview_file(sh)
        results.setdefault(r, []).append(sh.name)
    for k in sorted(results):
        print(f"=== {k} ({len(results[k])}) ===")
        for f in results[k]:
            print(f"  {f}")


if __name__ == "__main__":
    main()
