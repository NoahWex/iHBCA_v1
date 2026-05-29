"""
Step 07 Final Report - Aggregate Metrics Module

Extract per-sample metrics from preprocessing_manifest and compute derived metrics.
Pure computation module - NO file I/O.

Author: Agent 3 (Algorithm Developer)
Date: 2025-11-03
"""

import pandas as pd
import numpy as np
from typing import Dict, List
import warnings


def aggregate_step_summaries(preprocessing_manifest: dict) -> pd.DataFrame:
    """
    Parse preprocessing_manifest and extract per-sample metrics from Steps 01-06.

    Args:
        preprocessing_manifest: Dict loaded from preprocessing_manifest.yaml

    Returns:
        DataFrame with columns:
            - patient_id, position_id, sample_id
            - n_cells_baseline, n_baseline_vfs (Step 01)
            - n_cells_removed_step02 (Step 02)
            - n_doublets (Step 03)
            - n_cells_final (Step 05)
            - n_post_qc_vfs (Step 05)
            - status (overall pipeline status)

    Notes:
        - Gracefully handle missing step blocks (log warning, skip sample)
        - Check status="pass" before extracting metrics
        - Use field names from SCHEMA_REPORT (verified accurate)

    Example:
        >>> manifest = load_yaml('preprocessing_manifest.yaml')
        >>> df = aggregate_step_summaries(manifest)
        >>> # Returns DataFrame with all sample metrics
    """

    # Initialize list to collect sample metrics
    sample_metrics = []

    # Iterate through all patients and positions
    for patient_id, patient_data in preprocessing_manifest['patients'].items():
        for position_id, pos_data in patient_data['positions'].items():
            sample_id = pos_data['sample_id']

            # Check if all required steps are present
            required_steps = ['step_01_baseline_vfs', 'step_02_cell_filtering',
                            'step_03_doublet', 'step_05_post_qc_vfs']

            missing_steps = [step for step in required_steps if step not in pos_data]
            if missing_steps:
                warnings.warn(f"Missing steps {missing_steps} for {sample_id}, skipping")
                continue

            # Check if all required steps passed
            failed_steps = [step for step in required_steps
                          if pos_data[step].get('status') != 'pass']
            if failed_steps:
                warnings.warn(f"Failed steps {failed_steps} for {sample_id}, skipping")
                continue

            try:
                # Extract Step 01: Baseline counts
                step_01 = pos_data['step_01_baseline_vfs']['summary']
                n_cells_baseline = step_01['n_cells_passing']
                n_baseline_vfs = step_01['n_vfs']

                # Extract Step 02: Cell removal
                step_02 = pos_data['step_02_cell_filtering']['summary']
                n_cells_after_step02 = step_02['n_cells_passing']
                n_cells_removed_step02 = step_02['n_cells_input'] - n_cells_after_step02

                # Extract Step 03: Doublet detection
                step_03 = pos_data['step_03_doublet']['summary']
                n_doublets = step_03['n_doublets']

                # Extract Step 05: VF retention and final counts
                step_05 = pos_data['step_05_post_qc_vfs']['summary']
                n_post_qc_vfs = step_05['n_vfs_filtered']
                n_cells_final = step_05['n_cells_final_pass']

                # Append to list
                sample_metrics.append({
                    'patient_id': patient_id,
                    'position_id': position_id,
                    'sample_id': sample_id,
                    'n_cells_baseline': n_cells_baseline,
                    'n_baseline_vfs': n_baseline_vfs,
                    'n_cells_removed_step02': n_cells_removed_step02,
                    'n_doublets': n_doublets,
                    'n_cells_final': n_cells_final,
                    'n_post_qc_vfs': n_post_qc_vfs,
                    'status': 'pass'  # All steps passed if we got here
                })

            except KeyError as e:
                warnings.warn(f"Missing field {e} for {sample_id}, skipping")
                continue
            except Exception as e:
                warnings.warn(f"Error extracting metrics for {sample_id}: {e}, skipping")
                continue

    # Convert to DataFrame
    df = pd.DataFrame(sample_metrics)

    if len(df) == 0:
        warnings.warn("No samples with complete metrics found")
        return pd.DataFrame()  # Return empty DataFrame with proper structure

    return df


def compute_derived_metrics(df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate derived metrics from base metrics.

    Args:
        df: DataFrame from aggregate_step_summaries()

    Returns:
        DataFrame with added columns:
            - cell_removal_rate: (n_cells_removed + n_doublets) / n_cells_baseline
            - vf_retention_rate: n_post_qc_vfs / n_baseline_vfs

    Notes:
        - Handle division by zero (set to NaN)
        - Rates should be 0-1 (fractional, not percentage)

    Example:
        >>> df = aggregate_step_summaries(manifest)
        >>> df = compute_derived_metrics(df)
        >>> # Now df has cell_removal_rate and vf_retention_rate columns
    """

    if len(df) == 0:
        return df  # Return empty DataFrame unchanged

    # Create copy to avoid modifying original
    df = df.copy()

    # Calculate cell removal rate
    # Formula: (n_removed_step02 + n_doublets) / n_cells_baseline
    with np.errstate(divide='ignore', invalid='ignore'):
        df['cell_removal_rate'] = (
            (df['n_cells_removed_step02'] + df['n_doublets']) / df['n_cells_baseline']
        )

    # Calculate VF retention rate
    # Formula: n_post_qc_vfs / n_baseline_vfs
    with np.errstate(divide='ignore', invalid='ignore'):
        df['vf_retention_rate'] = df['n_post_qc_vfs'] / df['n_baseline_vfs']

    # Replace inf with NaN (from division by zero)
    df['cell_removal_rate'] = df['cell_removal_rate'].replace([np.inf, -np.inf], np.nan)
    df['vf_retention_rate'] = df['vf_retention_rate'].replace([np.inf, -np.inf], np.nan)

    return df
