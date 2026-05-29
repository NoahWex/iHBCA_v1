"""
Data Loading Module - Pure Computation (Shared between Steps 05 and 06)

This module contains PURE functions for loading and merging data from Steps 02 and 03.
NO FILE I/O in these functions - they accept in-memory objects and return in-memory results.

Algorithm Overview:
1. Merge doublet data from Step 03 CSVs into Step 02 AnnData object
2. Compute final passing cells (intersection of Step 02 and Step 03 filters)
3. Compute removal reasons for filtering analysis

Input Requirements:
- Step 02 AnnData object (patient-level merged, contains all cells)
- List of Step 03 DataFrames (one per sample in patient)

Output Format:
- Modified AnnData with doublet columns added
- pd.Index of passing cell barcodes
- pd.Series of removal reasons (categorical)

CRITICAL CORRECTION (per Agent 2):
- Boolean columns are Python bool (True/False), NOT strings ("TRUE"/"FALSE")
- Use bool comparison (== True) or direct boolean indexing
"""

import pandas as pd
import numpy as np


def merge_doublet_data(adata, doublet_df_list):
    """
    Merge Step 03 doublet flags into Step 02 AnnData object.

    Parameters:
    -----------
    adata : AnnData
        Step 02 merged object with:
        - .obs['cell_filter_pass_final']: bool, Step 02 filter status
        - .obs.index: cell barcodes
        - .obsm['X_umap']: UMAP coordinates

    doublet_df_list : list of pd.DataFrame
        List of Step 03 CSVs (one per sample in patient), each with:
        - 'cell_barcode': str, explicit column (NOT index)
        - 'doublet_filter_pass': bool, Step 03 filter status
        - 'scDblFinder_score': float, doublet score (0.0 - 1.0)

    Returns:
    --------
    AnnData
        Modified AnnData with doublet columns added to .obs:
        - 'doublet_filter_pass': bool
        - 'scDblFinder_score': float

    Raises:
    -------
    ValueError
        If merge results in missing doublet data for any cells

    Notes:
    ------
    - Uses .join() strategy verified by Agent 2 (Section 3)
    - Handles cell barcode index alignment between Step 02 (index) and Step 03 (column)
    - CRITICAL: Booleans are bool type (True/False), not strings

    Example:
    --------
    >>> adata_all = sc.read_h5ad("Pat1_merged_post_module2.h5ad")
    >>> doublet_dfs = [pd.read_csv(f"Pat1_P{i}.csv") for i in range(1, 16)]
    >>> adata_merged = merge_doublet_data(adata_all, doublet_dfs)
    >>> assert 'doublet_filter_pass' in adata_merged.obs.columns
    >>> assert adata_merged.obs['doublet_filter_pass'].dtype == bool
    """
    # Validate inputs
    if adata is None:
        raise ValueError("adata cannot be None")
    if not doublet_df_list:
        raise ValueError("doublet_df_list cannot be empty")

    # Check for required columns in adata
    if 'cell_filter_pass_final' not in adata.obs.columns:
        raise ValueError("adata.obs missing 'cell_filter_pass_final' column")

    # Concatenate all Step 03 CSVs
    doublet_merged = pd.concat(doublet_df_list, ignore_index=True)

    # Validate required columns
    required_cols = ['cell_barcode', 'doublet_filter_pass', 'scDblFinder_score']
    for col in required_cols:
        if col not in doublet_merged.columns:
            raise ValueError(f"doublet_df missing required column: {col}")

    # Set cell_barcode as index for joining
    # (Step 02 uses index, Step 03 has explicit column)
    doublet_merged_indexed = doublet_merged.set_index('cell_barcode')

    # Merge into adata.obs using left join on index
    # Only merge the necessary columns
    adata.obs = adata.obs.join(
        doublet_merged_indexed[['doublet_filter_pass', 'scDblFinder_score']],
        how='left'
    )

    # Verify merge success - check for missing data
    n_missing_filter = adata.obs['doublet_filter_pass'].isna().sum()
    n_missing_score = adata.obs['scDblFinder_score'].isna().sum()

    if n_missing_filter > 0 or n_missing_score > 0:
        raise ValueError(
            f"Merge failed: {n_missing_filter} cells missing doublet_filter_pass, "
            f"{n_missing_score} cells missing scDblFinder_score"
        )

    return adata


def compute_final_passing_cells(adata):
    """
    Compute cells passing BOTH Step 02 and Step 03 filters.

    Parameters:
    -----------
    adata : AnnData
        Object with merged filter columns:
        - 'cell_filter_pass_final': bool, Step 02 filter
        - 'doublet_filter_pass': bool, Step 03 filter

    Returns:
    --------
    pd.Index
        Cell barcodes passing both filters (for use with adata[index])

    Notes:
    ------
    - CRITICAL: Use bool comparison (True), NOT string comparison ('TRUE')
    - Returns Index (for subsetting), not boolean array
    - Formula: final_pass = (step_02_pass == True) & (step_03_pass == True)

    Example:
    --------
    >>> passing_idx = compute_final_passing_cells(adata)
    >>> adata_filtered = adata[passing_idx].copy()
    >>> print(f"Passing: {len(passing_idx)} / {len(adata)} cells")
    """
    # Validate inputs
    if 'cell_filter_pass_final' not in adata.obs.columns:
        raise ValueError("adata.obs missing 'cell_filter_pass_final' column")
    if 'doublet_filter_pass' not in adata.obs.columns:
        raise ValueError("adata.obs missing 'doublet_filter_pass' column")

    # Check that columns are boolean type (per Agent 2 correction)
    if adata.obs['cell_filter_pass_final'].dtype != bool:
        raise ValueError(
            f"cell_filter_pass_final has wrong dtype: {adata.obs['cell_filter_pass_final'].dtype} "
            "(expected bool)"
        )
    if adata.obs['doublet_filter_pass'].dtype != bool:
        raise ValueError(
            f"doublet_filter_pass has wrong dtype: {adata.obs['doublet_filter_pass'].dtype} "
            "(expected bool)"
        )

    # Compute final passing status using boolean AND
    # Both must be True (not 'TRUE' string!)
    final_pass = (
        adata.obs['cell_filter_pass_final'] &
        adata.obs['doublet_filter_pass']
    )

    # Return index of passing cells
    return adata.obs.index[final_pass]


def compute_removal_reasons(adata):
    """
    Compute hierarchical removal reasons for cells.

    Parameters:
    -----------
    adata : AnnData
        Object with both filter flags:
        - 'cell_filter_pass_final': bool, Step 02 filter
        - 'doublet_filter_pass': bool, Step 03 filter

    Returns:
    --------
    pd.Series
        Categorical series (length = n_cells) with values:
        - "Passing": passed both filters
        - "Step 03 Removal (Doublet)": failed doublet filter
        - "Step 02 Removal (QC)": failed Step 02 only
        Index matches adata.obs.index

    Notes:
    ------
    - Hierarchy: Step 03 > Step 02 > Passing (later filter takes priority)
    - CRITICAL: Use bool comparison (True/False), NOT string comparison
    - All cells assigned to exactly one category
    - Categories sum to total cells

    Hierarchical Logic:
    1. If doublet_filter_pass == False → "Step 03 Removal (Doublet)"
    2. Else if cell_filter_pass_final == False → "Step 02 Removal (QC)"
    3. Else → "Passing"

    Example:
    --------
    >>> removal_reasons = compute_removal_reasons(adata)
    >>> print(removal_reasons.value_counts())
    Passing                       3401
    Step 03 Removal (Doublet)     796
    Step 02 Removal (QC)          108
    dtype: int64

    >>> # Verify all cells accounted for
    >>> assert len(removal_reasons) == len(adata)
    >>> assert removal_reasons.isna().sum() == 0
    """
    # Validate inputs
    if 'cell_filter_pass_final' not in adata.obs.columns:
        raise ValueError("adata.obs missing 'cell_filter_pass_final' column")
    if 'doublet_filter_pass' not in adata.obs.columns:
        raise ValueError("adata.obs missing 'doublet_filter_pass' column")

    def assign_reason(row):
        """
        Hierarchical removal reason assignment.

        Priority:
        1. Step 03 (doublet) - latest filter
        2. Step 02 (QC) - earlier filter
        3. Passing
        """
        # Priority 1: Check doublet filter (Step 03)
        # Use bool comparison (not string!)
        if row['doublet_filter_pass'] == False:
            return 'Step 03 Removal (Doublet)'

        # Priority 2: Check QC filter (Step 02)
        if row['cell_filter_pass_final'] == False:
            return 'Step 02 Removal (QC)'

        # Passed both filters
        return 'Passing'

    # Apply hierarchical logic to each row
    removal_reasons = adata.obs.apply(assign_reason, axis=1)

    # Convert to categorical for efficiency and ordering
    category_order = [
        'Passing',
        'Step 02 Removal (QC)',
        'Step 03 Removal (Doublet)'
    ]
    removal_reasons = pd.Categorical(
        removal_reasons,
        categories=category_order,
        ordered=True
    )

    # Verify all cells assigned
    if removal_reasons.isna().sum() > 0:
        raise ValueError(
            f"Assignment failed: {removal_reasons.isna().sum()} cells have missing reasons"
        )

    return pd.Series(removal_reasons, index=adata.obs.index)


# ============================================================================
# HELPER FUNCTIONS (for validation and summary statistics)
# ============================================================================

def compute_filter_cascade_summary(adata):
    """
    Compute summary statistics for filtering cascade.

    Parameters:
    -----------
    adata : AnnData
        Object with both filter columns

    Returns:
    --------
    dict
        Summary statistics:
        - 'n_cells_total': Total cells (from Step 01)
        - 'n_cells_step_02_pass': Passed Step 02
        - 'n_cells_step_03_pass': Passed Step 03
        - 'n_cells_final_pass': Passed both
        - 'n_cells_step_02_only': Passed Step 02 but not Step 03
        - 'n_cells_step_03_only': Passed Step 03 but not Step 02
        - 'n_cells_both_fail': Failed both
        - 'pass_rate_overall': (final / total) * 100

    Example:
    --------
    >>> summary = compute_filter_cascade_summary(adata)
    >>> print(f"Overall pass rate: {summary['pass_rate_overall']:.1f}%")
    >>> print(f"Doublet removal: {summary['n_cells_step_02_only']} cells")
    """
    n_total = len(adata)
    n_step_02 = adata.obs['cell_filter_pass_final'].sum()
    n_step_03 = adata.obs['doublet_filter_pass'].sum()
    n_final = (adata.obs['cell_filter_pass_final'] & adata.obs['doublet_filter_pass']).sum()

    # Compute intersections
    n_step_02_only = (adata.obs['cell_filter_pass_final'] & ~adata.obs['doublet_filter_pass']).sum()
    n_step_03_only = (~adata.obs['cell_filter_pass_final'] & adata.obs['doublet_filter_pass']).sum()
    n_both_fail = (~adata.obs['cell_filter_pass_final'] & ~adata.obs['doublet_filter_pass']).sum()

    return {
        'n_cells_total': int(n_total),
        'n_cells_step_02_pass': int(n_step_02),
        'n_cells_step_03_pass': int(n_step_03),
        'n_cells_final_pass': int(n_final),
        'n_cells_step_02_only': int(n_step_02_only),
        'n_cells_step_03_only': int(n_step_03_only),
        'n_cells_both_fail': int(n_both_fail),
        'pass_rate_overall': float(n_final / n_total * 100)
    }


def validate_umap_presence(adata):
    """
    Validate that UMAP coordinates exist in AnnData object.

    Parameters:
    -----------
    adata : AnnData
        Object to validate

    Raises:
    -------
    ValueError
        If UMAP coordinates are missing or invalid

    Notes:
    ------
    - Checks for .obsm['X_umap'] existence
    - Validates shape (n_cells, 2)
    - Per Agent 2 Section 8 Gotcha 7
    """
    if 'X_umap' not in adata.obsm:
        raise ValueError(
            "UMAP coordinates missing: adata.obsm['X_umap'] not found. "
            "Step 02 may have failed UMAP computation."
        )

    umap_coords = adata.obsm['X_umap']

    if umap_coords.shape[1] != 2:
        raise ValueError(
            f"Invalid UMAP shape: expected (n_cells, 2), got {umap_coords.shape}"
        )

    if umap_coords.shape[0] != adata.n_obs:
        raise ValueError(
            f"UMAP size mismatch: {umap_coords.shape[0]} coords vs {adata.n_obs} cells"
        )

    # Check for NaN/Inf values
    if np.isnan(umap_coords).any() or np.isinf(umap_coords).any():
        raise ValueError("UMAP coordinates contain NaN or Inf values")
