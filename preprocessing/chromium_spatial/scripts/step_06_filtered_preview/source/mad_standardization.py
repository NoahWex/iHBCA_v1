"""
MAD Standardization Module - Pure Computation

This module contains PURE functions for computing MAD-standardized metrics.
NO FILE I/O - all inputs/outputs are in-memory objects.

Algorithm Overview:
Standardize QC metrics within cluster-sample groups using Median Absolute Deviation (MAD).
This shows which cells are outliers relative to their biological context (cluster + sample).

Formula: mad_standardized = (value - median) / MAD

Input Requirements:
- AnnData object with QC metrics in .obs
- Metric column name (e.g., 'total_counts')
- Cluster column name (default: 'cluster_coarse')
- Sample column name (default: 'sample_id')

Output Format:
- pd.Series with MAD-standardized values (same index as adata.obs)
- Typical range: [-3, +3] MAD units (most cells within ±3 MAD)

Reference:
- Agent 2 verification: Section 6 (MAD Standardization)
- Verified implementation from Agent 2 code snippet (lines 806-844)
"""

import numpy as np
import pandas as pd
from scipy.stats import median_abs_deviation


def compute_mad_standardized_values(
    adata,
    metric_col,
    cluster_col='cluster_coarse',
    sample_col='sample_id',
    min_group_size=10
):
    """
    Standardize metric values within cluster-sample groups using MAD.

    Parameters:
    -----------
    adata : AnnData
        Object with metric column in .obs

    metric_col : str
        Column name to standardize (e.g., 'total_counts', 'n_genes_by_counts', 'pct_counts_mt')

    cluster_col : str, optional
        Column for clustering (default: 'cluster_coarse')

    sample_col : str, optional
        Column for sample ID (default: 'sample_id')

    min_group_size : int, optional
        Minimum cells in cluster-sample group for MAD calculation (default: 10)
        Groups smaller than this are set to 0.0

    Returns:
    --------
    pd.Series
        MAD-standardized values: (value - median) / MAD
        Index matches adata.obs.index
        Range typically [-3, +3] MAD units

    Notes:
    ------
    - Uses Agent 2's verified strategy (Section 6)
    - Handles division by zero (MAD=0) → set to 0.0
    - Handles small groups (< min_group_size) → set to 0.0
    - Does NOT assume contiguous cluster IDs (uses .unique())

    MAD Calculation:
    ----------------
    For each cluster-sample group:
    1. Compute median of metric
    2. Compute MAD (Median Absolute Deviation)
    3. Standardize: (value - median) / MAD
    4. If MAD == 0 (uniform values), set to 0.0

    Interpretation:
    ---------------
    - Values near 0: Typical for this cluster-sample context
    - Values < -3 or > +3: Strong outliers (used in Step 02 filtering)
    - Positive: Above cluster-sample median
    - Negative: Below cluster-sample median

    Example:
    --------
    >>> # Standardize total UMI counts
    >>> umi_mad = compute_mad_standardized_values(adata_filtered, 'total_counts')
    >>> adata_filtered.obs['total_counts_mad'] = umi_mad

    >>> # Check range
    >>> print(f"MAD range: [{umi_mad.min():.2f}, {umi_mad.max():.2f}]")
    >>> print(f"Cells outside ±3 MAD: {(np.abs(umi_mad) > 3).sum()}")

    >>> # Standardize all QC metrics
    >>> for metric in ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']:
    >>>     mad_col = f'{metric}_mad'
    >>>     adata_filtered.obs[mad_col] = compute_mad_standardized_values(adata_filtered, metric)
    """
    # Validate inputs
    if metric_col not in adata.obs.columns:
        raise ValueError(f"Metric column '{metric_col}' not found in adata.obs")

    if cluster_col not in adata.obs.columns:
        raise ValueError(f"Cluster column '{cluster_col}' not found in adata.obs")

    if sample_col not in adata.obs.columns:
        raise ValueError(f"Sample column '{sample_col}' not found in adata.obs")

    # Check metric is numeric
    if not pd.api.types.is_numeric_dtype(adata.obs[metric_col]):
        raise ValueError(
            f"Metric column '{metric_col}' must be numeric, got {adata.obs[metric_col].dtype}"
        )

    # Initialize result series with zeros
    mad_standardized = pd.Series(0.0, index=adata.obs.index, name=f'{metric_col}_mad')

    # Get unique clusters and samples
    # CRITICAL: Don't assume contiguous cluster IDs (per Agent 2 Section 8 Gotcha 5)
    clusters = adata.obs[cluster_col].unique()
    samples = adata.obs[sample_col].unique()

    # Calculate MAD standardization within each cluster-sample group
    for cluster in clusters:
        for sample in samples:
            # Create mask for this cluster-sample group
            mask = (
                (adata.obs[cluster_col] == cluster) &
                (adata.obs[sample_col] == sample)
            )

            # Skip small groups (per Step 02 logic)
            group_size = mask.sum()
            if group_size < min_group_size:
                # Keep as 0.0 for small groups
                continue

            # Extract values for this group
            vals = adata.obs.loc[mask, metric_col].values

            # Compute median and MAD
            median_val = np.median(vals)
            mad_val = median_abs_deviation(vals)

            # Standardize: (value - median) / MAD
            # Protect against division by zero (per Agent 2 Section 8 Gotcha 6)
            if mad_val > 0:
                mad_standardized.loc[mask] = (vals - median_val) / mad_val
            else:
                # Uniform distribution (all values identical) → set to 0.0
                mad_standardized.loc[mask] = 0.0

    return mad_standardized


def compute_all_mad_metrics(adata, metrics=None):
    """
    Compute MAD-standardized versions of multiple metrics.

    Parameters:
    -----------
    adata : AnnData
        Object with QC metrics in .obs

    metrics : list of str, optional
        Metrics to standardize (default: ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'])

    Returns:
    --------
    AnnData
        Modified in place with new columns: {metric}_mad for each metric

    Notes:
    ------
    - Convenience function for standardizing standard QC metrics
    - Adds columns directly to adata.obs
    - Uses default cluster_col='cluster_coarse', sample_col='sample_id'

    Example:
    --------
    >>> compute_all_mad_metrics(adata_filtered)
    >>> assert 'n_genes_by_counts_mad' in adata_filtered.obs.columns
    >>> assert 'total_counts_mad' in adata_filtered.obs.columns
    >>> assert 'pct_counts_mt_mad' in adata_filtered.obs.columns
    """
    # Default metrics per Agent 1 spec (Section 2, Plot 6)
    if metrics is None:
        metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']

    # Compute MAD standardization for each metric
    for metric in metrics:
        mad_col = f'{metric}_mad'
        adata.obs[mad_col] = compute_mad_standardized_values(adata, metric)

    return adata


# ============================================================================
# HELPER FUNCTIONS (for validation and diagnostics)
# ============================================================================

def summarize_mad_distribution(mad_series, metric_name):
    """
    Summarize MAD-standardized distribution.

    Parameters:
    -----------
    mad_series : pd.Series
        MAD-standardized values

    metric_name : str
        Name of metric for display

    Returns:
    --------
    dict
        Distribution summary:
        - 'metric': Metric name
        - 'min': Minimum value
        - 'q25': 25th percentile
        - 'median': Median (should be near 0)
        - 'q75': 75th percentile
        - 'max': Maximum value
        - 'n_outliers_low': Cells < -3 MAD
        - 'n_outliers_high': Cells > +3 MAD
        - 'outlier_rate': Percentage of outliers

    Example:
    --------
    >>> summary = summarize_mad_distribution(adata.obs['total_counts_mad'], 'total_counts')
    >>> print(f"Outliers: {summary['outlier_rate']:.1f}%")
    """
    summary = {
        'metric': metric_name,
        'min': float(mad_series.min()),
        'q25': float(mad_series.quantile(0.25)),
        'median': float(mad_series.median()),
        'q75': float(mad_series.quantile(0.75)),
        'max': float(mad_series.max()),
        'n_outliers_low': int((mad_series < -3).sum()),
        'n_outliers_high': int((mad_series > 3).sum())
    }

    # Calculate outlier rate
    n_outliers = summary['n_outliers_low'] + summary['n_outliers_high']
    summary['outlier_rate'] = float(n_outliers / len(mad_series) * 100)

    return summary


def validate_mad_calculation(adata, metric_col, cluster_col='cluster_coarse', sample_col='sample_id'):
    """
    Validate MAD standardization for a single cluster-sample group.

    Parameters:
    -----------
    adata : AnnData
        Object with MAD-standardized column

    metric_col : str
        Original metric column

    cluster_col : str, optional
        Cluster column

    sample_col : str, optional
        Sample column

    Returns:
    --------
    dict
        Validation results for first cluster-sample group:
        - 'cluster': Cluster ID
        - 'sample': Sample ID
        - 'n_cells': Group size
        - 'median': Group median
        - 'mad': Group MAD
        - 'mad_standardized_mean': Mean of standardized values (should be ~0)
        - 'mad_standardized_std': Std of standardized values (should be ~1)

    Example:
    --------
    >>> validation = validate_mad_calculation(adata, 'total_counts')
    >>> print(f"Standardized mean: {validation['mad_standardized_mean']:.3f} (expect ~0)")
    >>> print(f"Standardized std: {validation['mad_standardized_std']:.3f} (expect ~1)")
    """
    # Get first non-empty cluster-sample group
    clusters = adata.obs[cluster_col].unique()
    samples = adata.obs[sample_col].unique()

    for cluster in clusters:
        for sample in samples:
            mask = (
                (adata.obs[cluster_col] == cluster) &
                (adata.obs[sample_col] == sample)
            )

            if mask.sum() >= 10:
                # Found suitable group
                vals = adata.obs.loc[mask, metric_col].values
                mad_vals = adata.obs.loc[mask, f'{metric_col}_mad'].values

                median_val = np.median(vals)
                mad_val = median_abs_deviation(vals)

                return {
                    'cluster': cluster,
                    'sample': sample,
                    'n_cells': int(mask.sum()),
                    'median': float(median_val),
                    'mad': float(mad_val),
                    'mad_standardized_mean': float(mad_vals.mean()),
                    'mad_standardized_std': float(mad_vals.std())
                }

    raise ValueError("No suitable cluster-sample group found for validation (all < 10 cells)")
