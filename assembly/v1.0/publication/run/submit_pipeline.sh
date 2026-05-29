#!/bin/bash
# =============================================================================
# submit_pipeline.sh — Submit iHBCA build pipeline from pipeline.yaml DAG
# =============================================================================
# Reads pipeline.yaml, resolves stage dependencies, and submits SLURM jobs
# with correct --dependency flags. Existing shell scripts are called as-is.
#
# Usage:
#   submit_pipeline.sh [--dry-run] [--stages STAGE1,STAGE2,...] [--from STAGE]
#
# Options:
#   --dry-run     Show what would be submitted without actually submitting
#   --stages      Comma-separated list of specific stages to run
#   --from        Start from a specific stage (includes all downstream)
#   --skip        Comma-separated stages to skip
#
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_YAML="${SCRIPT_DIR}/pipeline.yaml"

# --- Parse arguments ---
DRY_RUN=false
FILTER_STAGES=""
FROM_STAGE=""
SKIP_STAGES=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)  DRY_RUN=true; shift ;;
        --stages)   FILTER_STAGES="$2"; shift 2 ;;
        --from)     FROM_STAGE="$2"; shift 2 ;;
        --skip)     SKIP_STAGES="$2"; shift 2 ;;
        -h|--help)
            echo "Usage: submit_pipeline.sh [--dry-run] [--stages S1,S2] [--from S] [--skip S1,S2]"
            echo ""
            echo "Options:"
            echo "  --dry-run     Show sbatch commands without submitting"
            echo "  --stages      Run only specific stages (comma-separated)"
            echo "  --from        Start from stage, include all downstream"
            echo "  --skip        Skip specific stages (comma-separated)"
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

if [[ ! -f "$PIPELINE_YAML" ]]; then
    echo "ERROR: Pipeline definition not found: $PIPELINE_YAML" >&2
    exit 1
fi

# --- Parse pipeline.yaml ---
# Pure-Python parser using only stdlib (json module). Converts simple YAML
# structure to pipe-delimited records without requiring PyYAML.
PARSED=$(python3 - "$PIPELINE_YAML" << 'PYEOF'
import re, sys

def parse_pipeline_yaml(path):
    """Minimal YAML parser for pipeline.yaml's fixed structure."""
    with open(path) as f:
        lines = f.readlines()

    stages = {}
    current_stage = None
    current_section = None  # 'depends_on', 'env', 'resources'
    current_dep = {}

    for line in lines:
        stripped = line.rstrip()
        if not stripped or stripped.lstrip().startswith('#'):
            continue

        indent = len(line) - len(line.lstrip())

        # Top-level stage name (indent=2 under stages:)
        if indent == 2 and stripped.endswith(':') and not stripped.lstrip().startswith('-'):
            # Flush pending dep from previous stage
            if current_dep and current_stage:
                stages[current_stage]['depends_on'].append(dict(current_dep))
                current_dep = {}
            name = stripped.strip().rstrip(':')
            if name != 'stages':
                current_stage = name
                stages[name] = {'script': '', 'sbatch_args': '', 'depends_on': [], 'env': {}}
                current_section = None
            continue

        if current_stage is None:
            continue

        key_val = stripped.strip()

        # Section headers at indent=4
        if indent == 4:
            if key_val == 'depends_on:':
                current_section = 'depends_on'
                continue
            elif key_val == 'env:':
                current_section = 'env'
                continue
            elif key_val == 'resources:':
                current_section = 'resources'
                continue
            elif key_val == 'depends_on: []':
                stages[current_stage]['depends_on'] = []
                current_section = None
                continue

            # Key-value at stage level
            current_section = None
            m = re.match(r'(\w+):\s*(.+)', key_val)
            if m:
                k, v = m.group(1), m.group(2).strip().strip('"').strip("'")
                if k in ('script', 'sbatch_args', 'description'):
                    stages[current_stage][k] = v

        # List items / nested keys at indent=6+
        elif indent >= 6:
            if current_section == 'depends_on':
                if key_val.startswith('- stage:'):
                    if current_dep:
                        stages[current_stage]['depends_on'].append(dict(current_dep))
                    current_dep = {'stage': key_val.split(':', 1)[1].strip()}
                elif key_val.startswith('type:'):
                    current_dep['type'] = key_val.split(':', 1)[1].strip()
            elif current_section == 'env':
                m = re.match(r'(\w+):\s*(.+)', key_val)
                if m:
                    stages[current_stage]['env'][m.group(1)] = m.group(2).strip()

    # Flush last dep
    if current_dep and current_stage:
        stages[current_stage]['depends_on'].append(dict(current_dep))

    return stages

try:
    stages = parse_pipeline_yaml(sys.argv[1])
    for name, stage in stages.items():
        script = stage.get('script', '')
        sbatch_args = stage.get('sbatch_args', '')
        deps = ','.join(f"{d['stage']}:{d['type']}" for d in stage.get('depends_on', []))
        env_pairs = ','.join(f'{k}={v}' for k, v in stage.get('env', {}).items())
        print(f"{name}|{script}|{sbatch_args}|{deps}|{env_pairs}")
except Exception as e:
    print(f"ERROR: {e}", file=sys.stderr)
    sys.exit(1)
PYEOF
)

if [[ -z "$PARSED" ]]; then
    echo "ERROR: Cannot parse $PIPELINE_YAML" >&2
    exit 1
fi

# --- Build submission plan ---
declare -A STAGE_SCRIPT
declare -A STAGE_SBATCH_ARGS
declare -A STAGE_DEPS
declare -A STAGE_ENV
declare -a STAGE_ORDER=()

while IFS='|' read -r name script sbatch_args deps env; do
    STAGE_ORDER+=("$name")
    STAGE_SCRIPT[$name]="$script"
    STAGE_SBATCH_ARGS[$name]="$sbatch_args"
    STAGE_DEPS[$name]="$deps"
    STAGE_ENV[$name]="$env"
done <<< "$PARSED"

# --- Filter stages if requested ---
declare -A SKIP_SET=()
if [[ -n "$SKIP_STAGES" ]]; then
    IFS=',' read -ra SKIP_LIST <<< "$SKIP_STAGES"
    for s in "${SKIP_LIST[@]}"; do
        SKIP_SET[$s]=1
    done
fi

declare -A INCLUDE_SET=()
if [[ -n "$FILTER_STAGES" ]]; then
    IFS=',' read -ra INCLUDE_LIST <<< "$FILTER_STAGES"
    for s in "${INCLUDE_LIST[@]}"; do
        INCLUDE_SET[$s]=1
    done
fi

# --from: include the named stage and everything downstream
if [[ -n "$FROM_STAGE" ]]; then
    # BFS from FROM_STAGE through reverse dependency graph
    declare -A REACHABLE
    REACHABLE[$FROM_STAGE]=1
    changed=true
    while $changed; do
        changed=false
        for name in "${STAGE_ORDER[@]}"; do
            [[ -n "${REACHABLE[$name]:-}" ]] && continue
            deps="${STAGE_DEPS[$name]}"
            [[ -z "$deps" ]] && continue
            IFS=',' read -ra dep_list <<< "$deps"
            for dep_entry in "${dep_list[@]}"; do
                dep_name="${dep_entry%%:*}"
                if [[ -n "${REACHABLE[$dep_name]:-}" ]]; then
                    REACHABLE[$name]=1
                    changed=true
                    break
                fi
            done
        done
    done
    for name in "${STAGE_ORDER[@]}"; do
        if [[ -n "${REACHABLE[$name]:-}" ]]; then
            INCLUDE_SET[$name]=1
        fi
    done
fi

should_run() {
    local name="$1"
    [[ -n "${SKIP_SET[$name]:-}" ]] && return 1
    if [[ ${#INCLUDE_SET[@]} -gt 0 ]]; then
        [[ -n "${INCLUDE_SET[$name]:-}" ]] && return 0
        return 1
    fi
    return 0
}

# --- Submit stages in dependency order ---
declare -A JOB_IDS

echo "============================================="
echo "iHBCA v1.0 Build Pipeline"
echo "  Config: $PIPELINE_YAML"
echo "  Mode: $(if $DRY_RUN; then echo 'DRY RUN'; else echo 'SUBMIT'; fi)"
echo "============================================="
echo ""

for name in "${STAGE_ORDER[@]}"; do
    if ! should_run "$name"; then
        echo "SKIP: $name"
        continue
    fi

    script="${STAGE_SCRIPT[$name]}"
    sbatch_args="${STAGE_SBATCH_ARGS[$name]}"
    deps="${STAGE_DEPS[$name]}"
    env_vars="${STAGE_ENV[$name]}"
    script_path="${SCRIPT_DIR}/${script}"

    # Build dependency string
    dep_flags=""
    if [[ -n "$deps" ]]; then
        IFS=',' read -ra dep_list <<< "$deps"
        for dep_entry in "${dep_list[@]}"; do
            dep_name="${dep_entry%%:*}"
            dep_type="${dep_entry##*:}"
            dep_job="${JOB_IDS[$dep_name]:-}"

            if [[ -z "$dep_job" ]]; then
                # Dependency was skipped or not submitted
                echo "  NOTE: Dependency '$dep_name' not submitted — $name will run without it"
                continue
            fi

            if [[ -n "$dep_flags" ]]; then
                dep_flags="${dep_flags},${dep_type}:${dep_job}"
            else
                dep_flags="--dependency=${dep_type}:${dep_job}"
            fi
        done
    fi

    # Build export string for env vars
    export_flag=""
    if [[ -n "$env_vars" ]]; then
        export_flag="--export=ALL,${env_vars}"
    fi

    # Build sbatch command
    cmd="sbatch --parsable"
    [[ -n "$dep_flags" ]]   && cmd="$cmd $dep_flags"
    [[ -n "$export_flag" ]] && cmd="$cmd $export_flag"
    [[ -n "$sbatch_args" ]] && cmd="$cmd $sbatch_args"
    cmd="$cmd $script_path"

    echo "STAGE: $name"
    echo "  Script: $script_path"
    [[ -n "$dep_flags" ]]   && echo "  Deps: $dep_flags"
    [[ -n "$export_flag" ]] && echo "  Env: $export_flag"
    [[ -n "$sbatch_args" ]] && echo "  Args: $sbatch_args"

    if $DRY_RUN; then
        echo "  [DRY RUN] $cmd"
        JOB_IDS[$name]="DRY_${name}"
    else
        job_id=$($cmd)
        JOB_IDS[$name]="$job_id"
        echo "  Submitted: $job_id"
    fi
    echo ""
done

# --- Summary ---
echo "============================================="
echo "PIPELINE SUBMISSION SUMMARY"
echo "============================================="
for name in "${STAGE_ORDER[@]}"; do
    job="${JOB_IDS[$name]:-SKIPPED}"
    printf "  %-25s  %s\n" "$name" "$job"
done
echo ""

if ! $DRY_RUN; then
    echo "Monitor: squeue -u \$USER"
    echo "Cancel all: scancel ${JOB_IDS[*]}"
fi
