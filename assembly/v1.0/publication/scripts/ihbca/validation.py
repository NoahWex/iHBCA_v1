"""Shared validation checks for iHBCA assembly pipeline.

Sub-functions for Tier 1 field validation, Ensembl coverage, UMAP presence,
and dataset metadata checks. Used by both source and integrated assembly.
"""


def validate_tier1_fields(adata, tier1_fields, prefix=""):
    """Check Tier 1 field presence and NaN values.

    Returns list of issue strings.
    """
    issues = []
    for field in tier1_fields:
        if field not in adata.obs.columns:
            issues.append(f"{prefix}Missing Tier 1 field: {field}")
        else:
            n_na = adata.obs[field].isna().sum()
            if n_na > 0:
                pct = n_na / adata.n_obs * 100
                issues.append(f"{prefix}{field}: {n_na:,} NaN ({pct:.1f}%)")
    return issues


def validate_ensembl_coverage(adata, prefix="", threshold=50, verbose=True):
    """Check Ensembl ID coverage in var.index.

    Returns list of issue strings. Prints coverage if verbose.
    """
    issues = []
    n_ensg = sum(1 for g in adata.var_names if str(g).startswith("ENSG"))
    pct_ensg = n_ensg / max(adata.n_vars, 1) * 100
    if verbose:
        print(f"  Ensembl coverage: {n_ensg:,}/{adata.n_vars:,} ({pct_ensg:.1f}%)")
    if pct_ensg < threshold:
        issues.append(f"{prefix}Low Ensembl coverage: {pct_ensg:.1f}%")
    return issues


def validate_umap(adata, prefix="", check_shape=False):
    """Check X_umap presence (and optionally shape).

    Returns list of issue strings.
    """
    issues = []
    if "X_umap" not in adata.obsm:
        issues.append(f"{prefix}Missing obsm['X_umap']")
    elif check_shape:
        shape = adata.obsm["X_umap"].shape
        if shape != (adata.n_obs, 2):
            issues.append(f"{prefix}UMAP shape: {shape} != ({adata.n_obs}, 2)")
    return issues


def validate_dataset_metadata(adata, dataset_meta_fields, prefix=""):
    """Check dataset metadata fields in uns.

    Returns list of issue strings.
    """
    issues = []
    for field in dataset_meta_fields:
        if field not in adata.uns:
            issues.append(f"{prefix}Missing dataset metadata: {field}")
        elif adata.uns[field] == "unknown":
            issues.append(f"{prefix}Dataset metadata placeholder: {field} = 'unknown'")
    return issues


def validate_schema_version(adata, expected="5.3.2", prefix=""):
    """Check schema_version in uns.

    Returns list of issue strings.
    """
    issues = []
    sv = adata.uns.get("schema_version")
    if sv != expected:
        issues.append(f"{prefix}Schema version: {sv} (expected {expected})")
    return issues
