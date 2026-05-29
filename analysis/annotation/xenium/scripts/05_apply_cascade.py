"""Apply annotation cascade YAML to leiden_assignments.csv → cell_annotations.csv

Usage:
    python 06_apply_cascade.py <cascade.yaml> <leiden_assignments.csv> <out_cell_annotations.csv>

Cascade format (two override styles both supported):
    meta:        base_resolution (default "0.5")
    base:        leiden_{base_resolution} cluster → {label, label_short, lineage}
    artifacts:   resolution + per-cluster artifact calls (at a different resolution)
    xenium_only: same structure as artifacts — Xenium-only clusters without FLEX
    overrides:   resolution + per-cluster overrides; supports two formats:
                   NEW: cluster number is the YAML key, entry has {label, label_short, ...}
                   OLD: named entries with {override_cluster, new_label, new_label_short, ...}
                        or nested {override_clusters: {name: {cluster, new_label, ...}}}
    artifact_labels: list of labels to flag as is_artifact=True
"""
import pandas as pd
import yaml
import sys

cascade_path = sys.argv[1]
assignments_path = sys.argv[2]
out_path = sys.argv[3]

with open(cascade_path) as f:
    cascade = yaml.safe_load(f)

meta = cascade.get("meta", {})
base_res = str(meta.get("base_resolution", "0.5"))
leiden_base_col = f"leiden_{base_res}"

# Read all leiden columns as str to avoid dtype issues
dtype_map = {c: str for c in [f"leiden_{r}" for r in ["0.1","0.3","0.5","0.7","1.0","1.5"]]}
asgn = pd.read_csv(assignments_path, dtype=dtype_map)

# -----------------------------------------------------------------------
# BASE labels
# -----------------------------------------------------------------------
artifact_labels = set(cascade.get("artifact_labels", []))

base = {k: v for k, v in cascade["base"].items()}
asgn["label"] = asgn[leiden_base_col].map({k: v["label"] for k, v in base.items()})
asgn["label_short"] = asgn[leiden_base_col].map({k: v["label_short"] for k, v in base.items()})
asgn["lineage"] = asgn[leiden_base_col].map({k: v["lineage"] for k, v in base.items()})

# -----------------------------------------------------------------------
# ARTIFACT / XENIUM-ONLY clusters (apply before overrides)
# Both sections use the same structure: {resolution, cluster_key: {label, ...}}
# -----------------------------------------------------------------------
def apply_cluster_section(asgn, section, default_res):
    """Apply a section that maps leiden cluster → label at a specific resolution."""
    if not section:
        return
    res = str(section.get("resolution", default_res))
    leiden_col = f"leiden_{res}"
    if leiden_col not in asgn.columns:
        print(f"   WARNING: {leiden_col} not found, skipping section", flush=True)
        return
    for key, val in section.items():
        if not isinstance(val, dict):
            continue
        # Named-entry format: {cluster: "X", label: ..., label_short: ..., lineage: ...}
        if "cluster" in val:
            clust = str(val["cluster"])
            mask = asgn[leiden_col] == clust
            asgn.loc[mask, "label"] = val["label"]
            asgn.loc[mask, "label_short"] = val["label_short"]
            asgn.loc[mask, "lineage"] = val.get("lineage", "Artifact")
        # Direct cluster-key format: key IS the cluster number, val has label
        elif "label" in val and key.lstrip("-").isdigit():
            mask = asgn[leiden_col] == str(key)
            asgn.loc[mask, "label"] = val["label"]
            asgn.loc[mask, "label_short"] = val["label_short"]
            asgn.loc[mask, "lineage"] = val.get("lineage", "Artifact")

apply_cluster_section(asgn, cascade.get("artifacts", {}), base_res)
apply_cluster_section(asgn, cascade.get("xenium_only", {}), base_res)

# Intermediate-resolution sections: apply only Artifact-lineage entries so that
# non-artifact labels (LASP, LHS, Fibroblast, etc.) don't clobber base labels.
# smooth_muscle_detail and luminal_detail are the known sections of this type.
def apply_artifact_entries_from_section(asgn, section, default_res):
    if not section:
        return
    res = str(section.get("resolution", default_res))
    leiden_col = f"leiden_{res}"
    if leiden_col not in asgn.columns:
        return
    for key, val in section.items():
        if not isinstance(val, dict):
            continue
        lineage = val.get("lineage", "")
        label = val.get("label", "")
        # Apply if lineage is Artifact OR label is known artifact label
        is_art = lineage == "Artifact" or label in artifact_labels
        # Named-entry format
        if "cluster" in val and is_art:
            mask = asgn[leiden_col] == str(val["cluster"])
            asgn.loc[mask, "label"] = val["label"]
            asgn.loc[mask, "label_short"] = val["label_short"]
            asgn.loc[mask, "lineage"] = lineage
        # Direct cluster-key format — also handle non-artifact (Smooth_Muscle etc.)
        elif "label" in val and key.lstrip("-").isdigit():
            mask = asgn[leiden_col] == str(key)
            asgn.loc[mask, "label"] = val["label"]
            asgn.loc[mask, "label_short"] = val["label_short"]
            asgn.loc[mask, "lineage"] = val.get("lineage", lineage)

apply_artifact_entries_from_section(asgn, cascade.get("smooth_muscle_detail", {}), base_res)
apply_artifact_entries_from_section(asgn, cascade.get("luminal_detail", {}), base_res)

# -----------------------------------------------------------------------
# OVERRIDES at higher resolution — two format styles
# -----------------------------------------------------------------------
overrides = cascade.get("overrides", {})
override_res = str(overrides.get("resolution", "1.0"))
leiden_override_col = f"leiden_{override_res}"

override_map = {}   # override_leiden_cluster → (label, label_short, base_cluster_or_None)

for key, val in overrides.items():
    if not isinstance(val, dict):
        continue

    # NEW FORMAT: key is the cluster number, entry has 'label'
    if "label" in val and key.lstrip("-").isdigit():
        base_clust = str(val["base_label_overridden"]) if "base_label_overridden" in val else None
        # base_label_overridden is a label string, not a cluster number — use None for cluster gating
        override_map[str(key)] = (val["label"], val["label_short"], None)

    # OLD FORMAT A: named entry with override_cluster
    elif "override_cluster" in val:
        override_map[str(val["override_cluster"])] = (
            val["new_label"], val["new_label_short"],
            str(val["base_cluster"]) if "base_cluster" in val else None,
        )

    # OLD FORMAT B: named entry with override_clusters (nested split, e.g. cDC)
    elif "override_clusters" in val:
        parent_base = str(val["base_cluster"]) if "base_cluster" in val else None
        for subkey, subval in val["override_clusters"].items():
            if isinstance(subval, dict) and "cluster" in subval:
                override_map[str(subval["cluster"])] = (
                    subval["new_label"], subval["new_label_short"], parent_base,
                )

if leiden_override_col in asgn.columns and override_map:
    for clust, (lbl, lbl_short, base_clust) in override_map.items():
        mask = asgn[leiden_override_col] == clust
        if base_clust is not None:
            mask = mask & (asgn[leiden_base_col] == base_clust)
        asgn.loc[mask, "label"] = lbl
        asgn.loc[mask, "label_short"] = lbl_short

# -----------------------------------------------------------------------
# ARTIFACT FLAG
# -----------------------------------------------------------------------
asgn["is_artifact"] = asgn["label"].isin(artifact_labels)

# -----------------------------------------------------------------------
# OUTPUT
# -----------------------------------------------------------------------
# Include base and override resolution columns; deduplicate
out_cols = ["cell_id", "label", "label_short", "lineage", "is_artifact",
            leiden_base_col, leiden_override_col]
out_cols = list(dict.fromkeys(c for c in out_cols if c in asgn.columns))
asgn[out_cols].to_csv(out_path, index=False)

print(f"Written {len(asgn):,} cells → {out_path}")
print("\nLabel counts:")
print(asgn["label"].value_counts().to_string())
n_null = asgn["label"].isna().sum()
if n_null:
    print(f"\nWARNING: {n_null:,} cells have no label — check cascade coverage")
