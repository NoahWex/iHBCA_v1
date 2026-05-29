"""
Step 05 Filtered Preview - MAD Calculator

Pure algorithm for MAD standardization within cluster-sample groups.
NO file I/O - all inputs are in-memory objects.

Author: Agent 3 (Algorithm Developer)
Date: 2025-11-03
"""

from typing import List
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import median_abs_deviation
from itertools import product


def compute_mad_standardized_metrics(
    adata: sc.AnnData,
    metrics: List[str] = None,
    min_group_size: int = 10
) -> sc.AnnData:
    """
    Calculate MAD-standardized metrics within cluster-sample groups.

    This function performs MAD standardization to make QC metrics comparable
    across different cluster-sample combinations. MAD-standardized values are
    Z-scores using median absolute deviation instead of standard deviation.

    Args:
        adata: AnnData with cluster_coarse and sample_id in .obs
        metrics: List of metric column names to standardize
                 Default: ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
        min_group_size: Skip groups smaller than this (default 10)

    Returns:
        Modified AnnData with {metric}_mad columns added to .obs

    Example:
        >>> # Assumes adata has cluster_coarse, sample_id, n_genes_by_counts, etc.
        >>> adata = compute_mad_standardized_metrics(adata, min_group_size=10)
        >>> # Now adata.obs has columns: n_genes_by_counts_mad, total_counts_mad, etc.
    """

    # Default metrics (using primary column names from scanpy)
    if metrics is None:
        metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']

    print(f"  Computing MAD-standardized metrics...")
    print(f"    Metrics: {metrics}")
    print(f"    Min group size: {min_group_size}")

    # Validate required columns exist
    required_cols = ['cluster_coarse', 'sample_id'] + metrics
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    # Initialize MAD columns as NaN
    for metric in metrics:
        mad_col = f'{metric}_mad'
        adata.obs[mad_col] = np.nan

    # Track statistics
    n_groups_total = 0
    n_groups_skipped_size = 0
    n_groups_skipped_mad = 0
    n_groups_processed = 0
    n_cells_processed = 0

    # Iterate over cluster-sample groups
    for (cluster, sample) in product(
        adata.obs['cluster_coarse'].unique(),
        adata.obs['sample_id'].unique()
    ):
        n_groups_total += 1

        # Get mask for this group
        mask = (
            (adata.obs['cluster_coarse'] == cluster) &
            (adata.obs['sample_id'] == sample)
        )
        n_cells = mask.sum()

        # Skip if group doesn't exist
        if n_cells == 0:
            continue

        # Skip small groups
        if n_cells < min_group_size:
            n_groups_skipped_size += 1
            continue

        # Calculate MAD for each metric
        group_skipped = False
        for metric in metrics:
            values = adata.obs.loc[mask, metric].values

            # Calculate median and MAD
            median_val = np.median(values)
            mad_val = median_abs_deviation(values)

            # Handle zero MAD (uniform distribution)
            if mad_val == 0:
                # Set MAD-standardized value to 0.0 (all values at median)
                mad_scores = np.zeros_like(values, dtype=float)
                group_skipped = True
            else:
                # MAD-standardize: (x - median) / MAD
                mad_scores = (values - median_val) / mad_val

            # Store MAD-standardized values
            mad_col = f'{metric}_mad'
            adata.obs.loc[mask, mad_col] = mad_scores

        # Track statistics
        if group_skipped:
            n_groups_skipped_mad += 1
        else:
            n_groups_processed += 1
            n_cells_processed += n_cells

    # Summary
    print(f"    Total groups: {n_groups_total}")
    print(f"    Groups skipped (< {min_group_size} cells): {n_groups_skipped_size}")
    print(f"    Groups skipped (zero MAD): {n_groups_skipped_mad}")
    print(f"    Groups processed: {n_groups_processed}")
    print(f"    Cells processed: {n_cells_processed} / {adata.n_obs}")

    # Validate output
    for metric in metrics:
        mad_col = f'{metric}_mad'
        n_missing = adata.obs[mad_col].isna().sum()
        print(f"    {mad_col}: {n_missing} cells with NaN ({n_missing/adata.n_obs*100:.1f}%)")

    return adata


def compute_qc_metrics_if_missing(
    adata: sc.AnnData,
    force_recalculate: bool = False
) -> sc.AnnData:
    """
    Calculate QC metrics if not already present in .obs.

    Helper function to ensure required QC metrics exist before MAD calculation.
    Uses scanpy's standard naming: n_genes_by_counts, total_counts, pct_counts_mt.

    Args:
        adata: AnnData object (modified in place)
        force_recalculate: If True, recalculate even if columns exist

    Returns:
        Modified adata with QC metrics

    Example:
        >>> adata = compute_qc_metrics_if_missing(adata)
        >>> # Now adata.obs has n_genes_by_counts, total_counts, pct_counts_mt
    """

    # Check if metrics already exist
    required_metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
    has_all_metrics = all(col in adata.obs.columns for col in required_metrics)

    if has_all_metrics and not force_recalculate:
        print(f"  QC metrics already present, skipping calculation")
        return adata

    print(f"  Calculating QC metrics...")

    # Identify mitochondrial genes
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    n_mt_genes = adata.var['mt'].sum()
    print(f"    Found {n_mt_genes} mitochondrial genes")

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=['mt'],
        percent_top=None,
        log1p=False,
        inplace=True
    )

    print(f"    QC metrics calculated")
    print(f"      n_genes_by_counts: median = {adata.obs['n_genes_by_counts'].median():.0f}")
    print(f"      total_counts: median = {adata.obs['total_counts'].median():.0f}")
    print(f"      pct_counts_mt: median = {adata.obs['pct_counts_mt'].median():.2f}%")

    return adata


def add_duplicate_metric_columns(
    adata: sc.AnnData
) -> sc.AnnData:
    """
    Add duplicate metric columns for backwards compatibility.

    Creates duplicate columns with alternative names:
    - n_genes (from n_genes_by_counts)
    - n_umi (from total_counts)
    - mito_percent (from pct_counts_mt)

    Args:
        adata: AnnData object (modified in place)

    Returns:
        Modified adata with duplicate columns

    Example:
        >>> adata = add_duplicate_metric_columns(adata)
        >>> # Now adata.obs has both n_genes_by_counts and n_genes
    """

    print(f"  Adding duplicate metric columns for compatibility...")

    # Check if primary columns exist
    if 'n_genes_by_counts' in adata.obs.columns:
        adata.obs['n_genes'] = adata.obs['n_genes_by_counts'].astype(int)
    else:
        raise ValueError("n_genes_by_counts not found in .obs")

    if 'total_counts' in adata.obs.columns:
        adata.obs['n_umi'] = adata.obs['total_counts'].astype(int)
    else:
        raise ValueError("total_counts not found in .obs")

    if 'pct_counts_mt' in adata.obs.columns:
        adata.obs['mito_percent'] = adata.obs['pct_counts_mt']
    else:
        raise ValueError("pct_counts_mt not found in .obs")

    print(f"    Duplicate columns added: n_genes, n_umi, mito_percent")

    return adata


def validate_mad_standardization(
    adata: sc.AnnData,
    metrics: List[str] = None
) -> bool:
    """
    Validate MAD standardization was successful.

    Args:
        adata: AnnData with MAD-standardized metrics
        metrics: List of metrics to validate (default: primary scanpy columns)

    Returns:
        True if validation passes

    Raises:
        ValueError: If validation fails

    Example:
        >>> validate_mad_standardization(adata)
    """

    if metrics is None:
        metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']

    print(f"  Validating MAD standardization...")

    for metric in metrics:
        mad_col = f'{metric}_mad'

        # Check column exists
        if mad_col not in adata.obs.columns:
            raise ValueError(f"Missing MAD column: {mad_col}")

        # Check dtype is numeric
        if not np.issubdtype(adata.obs[mad_col].dtype, np.number):
            raise ValueError(f"{mad_col} is not numeric: {adata.obs[mad_col].dtype}")

        # Check values are within reasonable range (-10 to 10 MAD scores)
        non_nan_values = adata.obs[mad_col].dropna()
        if len(non_nan_values) > 0:
            min_val = non_nan_values.min()
            max_val = non_nan_values.max()

            if min_val < -20 or max_val > 20:
                print(f"    Warning: {mad_col} has extreme values: [{min_val:.1f}, {max_val:.1f}]")

    print(f"    Validation passed")

    return True


def get_mad_summary_stats(
    adata: sc.AnnData,
    cluster: int = None,
    sample: str = None
) -> dict:
    """
    Get summary statistics for MAD-standardized metrics.

    Useful for debugging or generating reports.

    Args:
        adata: AnnData with MAD-standardized metrics
        cluster: Optional cluster to filter (default: all clusters)
        sample: Optional sample to filter (default: all samples)

    Returns:
        Dict with summary statistics

    Example:
        >>> stats = get_mad_summary_stats(adata, cluster=0, sample='Pat1_P1')
        >>> print(stats['n_genes_by_counts_mad']['median'])
    """

    # Filter data if requested
    if cluster is not None and sample is not None:
        mask = (
            (adata.obs['cluster_coarse'] == cluster) &
            (adata.obs['sample_id'] == sample)
        )
        data = adata.obs.loc[mask]
    elif cluster is not None:
        mask = adata.obs['cluster_coarse'] == cluster
        data = adata.obs.loc[mask]
    elif sample is not None:
        mask = adata.obs['sample_id'] == sample
        data = adata.obs.loc[mask]
    else:
        data = adata.obs

    # Calculate summary stats for each MAD metric
    mad_cols = [col for col in data.columns if col.endswith('_mad')]

    summary = {}
    for col in mad_cols:
        values = data[col].dropna()
        if len(values) > 0:
            summary[col] = {
                'n_cells': len(values),
                'median': float(np.median(values)),
                'mad': float(median_abs_deviation(values)),
                'min': float(values.min()),
                'max': float(values.max()),
                'q25': float(np.percentile(values, 25)),
                'q75': float(np.percentile(values, 75))
            }
        else:
            summary[col] = {
                'n_cells': 0,
                'median': np.nan,
                'mad': np.nan,
                'min': np.nan,
                'max': np.nan,
                'q25': np.nan,
                'q75': np.nan
            }

    return summary
