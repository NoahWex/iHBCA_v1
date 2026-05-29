"""
Step 05 Filtered Preview - Removal Tracker

Pure algorithm for assigning 3-category removal status to cells.
NO file I/O - all inputs are in-memory objects.

Author: Agent 3 (Algorithm Developer)
Date: 2025-11-03
"""

import pandas as pd
import scanpy as sc


def assign_removal_status(adata: sc.AnnData) -> sc.AnnData:
    """
    Assign removal status based on filter flags.

    Creates a 3-category removal status based on hierarchical priority:
    1. 'Removed by Step 03' - Failed doublet detection (highest priority)
    2. 'Removed by Step 02' - Failed MAD filtering
    3. 'Pass' - Passed both filters

    Args:
        adata: AnnData with cell_filter_pass_final and doublet_filter_pass in .obs

    Returns:
        Modified AnnData with 'removal_status' column added
        Categories: 'Pass', 'Removed by Step 02', 'Removed by Step 03'

    Example:
        >>> # Assumes adata.obs has cell_filter_pass_final, doublet_filter_pass
        >>> adata = assign_removal_status(adata)
        >>> # Now adata.obs has 'removal_status' column
        >>> print(adata.obs['removal_status'].value_counts())
    """

    print(f"  Assigning removal status...")

    # Validate required columns exist
    required_cols = ['cell_filter_pass_final', 'doublet_filter_pass']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    # Validate boolean types
    for col in required_cols:
        if adata.obs[col].dtype != bool:
            raise ValueError(f"Column {col} must be boolean, got {adata.obs[col].dtype}")

    # Apply hierarchical priority logic
    def _get_removal_status(row):
        """Determine removal status for a single cell."""
        # Priority 1: Check Step 03 (doublet)
        if not row['doublet_filter_pass']:
            return 'Removed by Step 03'

        # Priority 2: Check Step 02 (MAD outlier)
        if not row['cell_filter_pass_final']:
            return 'Removed by Step 02'

        # Both passed
        return 'Pass'

    # Apply to all cells
    adata.obs['removal_status'] = adata.obs.apply(_get_removal_status, axis=1)

    # Convert to categorical for efficiency
    adata.obs['removal_status'] = pd.Categorical(
        adata.obs['removal_status'],
        categories=['Pass', 'Removed by Step 02', 'Removed by Step 03'],
        ordered=False
    )

    # Summary statistics
    counts = adata.obs['removal_status'].value_counts()
    total = len(adata.obs)

    print(f"    Removal status assigned:")
    print(f"      Pass: {counts.get('Pass', 0)} ({counts.get('Pass', 0)/total*100:.1f}%)")
    print(f"      Removed by Step 02: {counts.get('Removed by Step 02', 0)} ({counts.get('Removed by Step 02', 0)/total*100:.1f}%)")
    print(f"      Removed by Step 03: {counts.get('Removed by Step 03', 0)} ({counts.get('Removed by Step 03', 0)/total*100:.1f}%)")

    return adata


def compute_final_pass_flag(adata: sc.AnnData) -> sc.AnnData:
    """
    Compute cells passing BOTH Step 02 AND Step 03.

    Creates 'final_pass' column as boolean AND of the two filter flags.

    Args:
        adata: AnnData with cell_filter_pass_final and doublet_filter_pass in .obs

    Returns:
        Modified AnnData with 'final_pass' column added

    Example:
        >>> adata = compute_final_pass_flag(adata)
        >>> # Now adata.obs has 'final_pass' column
        >>> print(f"Final passing cells: {adata.obs['final_pass'].sum()}")
    """

    print(f"  Computing final pass flag...")

    # Validate required columns exist
    required_cols = ['cell_filter_pass_final', 'doublet_filter_pass']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    # Boolean AND (both must be True)
    adata.obs['final_pass'] = (
        adata.obs['cell_filter_pass_final'] &
        adata.obs['doublet_filter_pass']
    )

    # Summary
    n_step02 = adata.obs['cell_filter_pass_final'].sum()
    n_step03 = adata.obs['doublet_filter_pass'].sum()
    n_final = adata.obs['final_pass'].sum()
    n_total = len(adata.obs)

    print(f"    Filter cascade:")
    print(f"      Total cells: {n_total}")
    print(f"      Step 02 pass: {n_step02} ({n_step02/n_total*100:.1f}%)")
    print(f"      Step 03 pass: {n_step03} ({n_step03/n_total*100:.1f}%)")
    print(f"      Both pass (final): {n_final} ({n_final/n_total*100:.1f}%)")

    return adata


def get_removal_breakdown(adata: sc.AnnData) -> dict:
    """
    Get detailed breakdown of removal reasons.

    Returns counts and percentages for each removal category.

    Args:
        adata: AnnData with removal_status in .obs

    Returns:
        Dict with removal breakdown statistics

    Example:
        >>> breakdown = get_removal_breakdown(adata)
        >>> print(breakdown['pass']['count'])
    """

    if 'removal_status' not in adata.obs.columns:
        raise ValueError("removal_status column not found - run assign_removal_status first")

    total = len(adata.obs)
    counts = adata.obs['removal_status'].value_counts()

    breakdown = {
        'total_cells': total,
        'pass': {
            'count': int(counts.get('Pass', 0)),
            'percent': float(counts.get('Pass', 0) / total * 100)
        },
        'removed_by_step02': {
            'count': int(counts.get('Removed by Step 02', 0)),
            'percent': float(counts.get('Removed by Step 02', 0) / total * 100)
        },
        'removed_by_step03': {
            'count': int(counts.get('Removed by Step 03', 0)),
            'percent': float(counts.get('Removed by Step 03', 0) / total * 100)
        }
    }

    # Add combined removal count
    breakdown['removed_total'] = {
        'count': breakdown['removed_by_step02']['count'] + breakdown['removed_by_step03']['count'],
        'percent': breakdown['removed_by_step02']['percent'] + breakdown['removed_by_step03']['percent']
    }

    return breakdown


def get_removal_by_cluster(adata: sc.AnnData) -> pd.DataFrame:
    """
    Get removal breakdown by cluster.

    Returns a DataFrame with removal counts per cluster.

    Args:
        adata: AnnData with removal_status and cluster_coarse in .obs

    Returns:
        DataFrame with clusters as rows, removal categories as columns

    Example:
        >>> df = get_removal_by_cluster(adata)
        >>> print(df)
        #                           Pass  Removed by Step 02  Removed by Step 03
        # cluster_coarse
        # 0                         1234                  56                  12
        # 1                          987                  45                   8
    """

    required_cols = ['removal_status', 'cluster_coarse']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    # Crosstab
    df = pd.crosstab(
        adata.obs['cluster_coarse'],
        adata.obs['removal_status']
    )

    # Ensure all categories present (even if zero)
    for cat in ['Pass', 'Removed by Step 02', 'Removed by Step 03']:
        if cat not in df.columns:
            df[cat] = 0

    # Sort by cluster ID
    df = df.sort_index()

    return df


def get_removal_by_sample(adata: sc.AnnData) -> pd.DataFrame:
    """
    Get removal breakdown by sample.

    Returns a DataFrame with removal counts per sample.

    Args:
        adata: AnnData with removal_status and sample_id in .obs

    Returns:
        DataFrame with samples as rows, removal categories as columns

    Example:
        >>> df = get_removal_by_sample(adata)
        >>> print(df)
        #                  Pass  Removed by Step 02  Removed by Step 03
        # sample_id
        # Pat1_P1          3352                 724                 229
        # Pat1_P2_Lower    2145                 456                 134
    """

    required_cols = ['removal_status', 'sample_id']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    # Crosstab
    df = pd.crosstab(
        adata.obs['sample_id'],
        adata.obs['removal_status']
    )

    # Ensure all categories present (even if zero)
    for cat in ['Pass', 'Removed by Step 02', 'Removed by Step 03']:
        if cat not in df.columns:
            df[cat] = 0

    # Sort by sample ID
    df = df.sort_index()

    return df


def validate_removal_status(adata: sc.AnnData) -> bool:
    """
    Validate removal status assignment is correct.

    Args:
        adata: AnnData with removal_status, final_pass, filter flags

    Returns:
        True if validation passes

    Raises:
        ValueError: If validation fails

    Example:
        >>> validate_removal_status(adata)
    """

    print(f"  Validating removal status...")

    # Check removal_status exists
    if 'removal_status' not in adata.obs.columns:
        raise ValueError("removal_status column not found")

    # Check all cells have a status
    if adata.obs['removal_status'].isna().any():
        raise ValueError("NaN values in removal_status")

    # Check counts sum to total
    total = len(adata.obs)
    counts = adata.obs['removal_status'].value_counts()
    sum_counts = counts.sum()

    if sum_counts != total:
        raise ValueError(f"Removal status counts don't sum to total: {sum_counts} != {total}")

    # Check consistency with final_pass
    if 'final_pass' in adata.obs.columns:
        n_pass = (adata.obs['removal_status'] == 'Pass').sum()
        n_final_pass = adata.obs['final_pass'].sum()

        if n_pass != n_final_pass:
            raise ValueError(f"Mismatch: removal_status='Pass' ({n_pass}) != final_pass ({n_final_pass})")

    # Check hierarchical logic
    # Cells removed by Step 03 should have doublet_filter_pass=False
    step03_removed = adata.obs['removal_status'] == 'Removed by Step 03'
    if 'doublet_filter_pass' in adata.obs.columns:
        invalid = step03_removed & adata.obs['doublet_filter_pass']
        if invalid.any():
            raise ValueError(f"{invalid.sum()} cells marked 'Removed by Step 03' but doublet_filter_pass=True")

    # Cells removed by Step 02 should have cell_filter_pass_final=False AND doublet_filter_pass=True
    step02_removed = adata.obs['removal_status'] == 'Removed by Step 02'
    if 'cell_filter_pass_final' in adata.obs.columns and 'doublet_filter_pass' in adata.obs.columns:
        invalid_02 = step02_removed & adata.obs['cell_filter_pass_final']
        if invalid_02.any():
            raise ValueError(f"{invalid_02.sum()} cells marked 'Removed by Step 02' but cell_filter_pass_final=True")

        invalid_03 = step02_removed & ~adata.obs['doublet_filter_pass']
        if invalid_03.any():
            raise ValueError(f"{invalid_03.sum()} cells marked 'Removed by Step 02' but doublet_filter_pass=False (should be Step 03)")

    print(f"    Validation passed")

    return True
