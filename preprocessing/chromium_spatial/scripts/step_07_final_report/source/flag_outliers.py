"""
Step 07 Final Report - Outlier Flagging Module

Apply MAD-based statistical outlier detection to identify problematic samples.
Pure computation module - NO file I/O.

Author: Agent 3 (Algorithm Developer)
Date: 2025-11-03
"""

import pandas as pd
import numpy as np
from typing import Tuple, List, Dict
from scipy.stats import median_abs_deviation


def compute_mad_scores(df: pd.DataFrame, metric_col: str) -> pd.DataFrame:
    """
    Compute MAD (Median Absolute Deviation) scores for a metric.

    Args:
        df: DataFrame with metric column
        metric_col: Name of column to score

    Returns:
        DataFrame with added columns:
            - {metric_col}_median: dataset-wide median
            - {metric_col}_mad: median absolute deviation
            - {metric_col}_mad_score: |value - median| / MAD

    Formula:
        median = np.median(values)
        mad = median_abs_deviation(values)
        mad_score = |value - median| / mad

    Notes:
        - Handle MAD=0 case (all values identical): set mad_score=0
        - Skip NaN values when computing median/MAD

    Example:
        >>> df = compute_mad_scores(df, 'cell_removal_rate')
        >>> # Now df has cell_removal_rate_median, _mad, _mad_score columns
    """

    if len(df) == 0:
        return df

    # Create copy to avoid modifying original
    df = df.copy()

    # Get values (remove NaN)
    values = df[metric_col].dropna().values

    if len(values) == 0:
        # All values are NaN
        df[f'{metric_col}_median'] = np.nan
        df[f'{metric_col}_mad'] = np.nan
        df[f'{metric_col}_mad_score'] = np.nan
        return df

    # Compute dataset-wide median
    median_val = np.median(values)

    # Compute MAD using scipy (more robust than manual calculation)
    mad_val = median_abs_deviation(values, nan_policy='omit')

    # Store median and MAD for all rows
    df[f'{metric_col}_median'] = median_val
    df[f'{metric_col}_mad'] = mad_val

    # Compute MAD scores
    if mad_val == 0:
        # All values are identical (or very close)
        # Set MAD score to 0 for all samples
        df[f'{metric_col}_mad_score'] = 0.0
    else:
        # Standard MAD score calculation: |value - median| / MAD
        df[f'{metric_col}_mad_score'] = np.abs(df[metric_col] - median_val) / mad_val

    return df


def flag_outliers(df: pd.DataFrame, metric_col: str, threshold: float) -> pd.DataFrame:
    """
    Flag samples where MAD score exceeds threshold.

    Args:
        df: DataFrame with MAD scores (from compute_mad_scores)
        metric_col: Metric column to flag
        threshold: MAD score threshold (e.g., 2.0)

    Returns:
        DataFrame with added column:
            - {metric_col}_flagged: bool (True if mad_score > threshold)

    Notes:
        - Requires compute_mad_scores() to be run first
        - NaN values should NOT be flagged (set to False)

    Example:
        >>> df = compute_mad_scores(df, 'cell_removal_rate')
        >>> df = flag_outliers(df, 'cell_removal_rate', threshold=2.0)
        >>> # Now df has cell_removal_rate_flagged column
    """

    if len(df) == 0:
        return df

    # Create copy to avoid modifying original
    df = df.copy()

    mad_score_col = f'{metric_col}_mad_score'

    # Check if MAD scores were computed
    if mad_score_col not in df.columns:
        raise ValueError(f"MAD scores not found. Run compute_mad_scores() first.")

    # Flag samples exceeding threshold
    # NaN values will evaluate to False (not flagged)
    df[f'{metric_col}_flagged'] = df[mad_score_col] > threshold

    # Explicitly set NaN MAD scores to False
    df.loc[df[mad_score_col].isna(), f'{metric_col}_flagged'] = False

    return df


def apply_flagging_pipeline(
    df: pd.DataFrame,
    metrics: List[str],
    threshold: float
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Apply MAD-based flagging to multiple metrics.

    Args:
        df: DataFrame with base metrics
        metrics: List of metric columns to flag (e.g., ['cell_removal_rate', 'vf_retention_rate'])
        threshold: MAD threshold for flagging

    Returns:
        Tuple of:
            - df_with_flags: Original DataFrame with MAD scores and flags added
            - df_flagged: Subset containing only flagged samples

    Notes:
        - Calls compute_mad_scores() and flag_outliers() for each metric
        - df_flagged includes samples flagged for ANY metric

    Example:
        >>> df_full, df_flagged = apply_flagging_pipeline(
        ...     df,
        ...     metrics=['cell_removal_rate', 'vf_retention_rate'],
        ...     threshold=2.0
        ... )
    """

    if len(df) == 0:
        return df, pd.DataFrame()

    # Start with copy of original DataFrame
    df_with_flags = df.copy()

    # Apply MAD scoring and flagging for each metric
    for metric in metrics:
        if metric not in df_with_flags.columns:
            raise ValueError(f"Metric '{metric}' not found in DataFrame")

        # Compute MAD scores
        df_with_flags = compute_mad_scores(df_with_flags, metric)

        # Flag outliers
        df_with_flags = flag_outliers(df_with_flags, metric, threshold)

    # Create flagged samples subset
    # Sample is flagged if ANY metric is flagged (OR logic)
    flag_cols = [f'{metric}_flagged' for metric in metrics]

    # Check if any flag is True for each row
    any_flagged = df_with_flags[flag_cols].any(axis=1)

    df_flagged = df_with_flags[any_flagged].copy()

    return df_with_flags, df_flagged


def create_flagged_records(
    df: pd.DataFrame,
    metrics: List[str]
) -> pd.DataFrame:
    """
    Create long-format flagged records table for report.

    Args:
        df: DataFrame with MAD scores and flags (from apply_flagging_pipeline)
        metrics: List of metrics that were flagged

    Returns:
        DataFrame with columns:
            - sample_id, patient_id, metric, value, mad_score

    Notes:
        - One row per (sample, metric) combination where sample is flagged for that metric
        - Sorted by MAD score descending

    Example:
        >>> flagged_records = create_flagged_records(
        ...     df_full,
        ...     ['cell_removal_rate', 'vf_retention_rate']
        ... )
    """

    if len(df) == 0:
        return pd.DataFrame(columns=['sample_id', 'patient_id', 'metric', 'value', 'mad_score'])

    records = []

    for metric in metrics:
        flag_col = f'{metric}_flagged'
        mad_score_col = f'{metric}_mad_score'

        # Get samples flagged for this metric
        flagged_mask = df[flag_col] == True

        for idx, row in df[flagged_mask].iterrows():
            records.append({
                'sample_id': row['sample_id'],
                'patient_id': row['patient_id'],
                'metric': metric,
                'value': row[metric],
                'mad_score': row[mad_score_col]
            })

    # Convert to DataFrame
    flagged_df = pd.DataFrame(records)

    if len(flagged_df) == 0:
        return flagged_df

    # Sort by MAD score descending (highest outliers first)
    flagged_df = flagged_df.sort_values('mad_score', ascending=False).reset_index(drop=True)

    return flagged_df
