"""
Reprocessing Module - Pure Computation (Step 05 Specific)

This module contains PURE functions for reprocessing filtered cells using filtered VF lists.
NO FILE I/O - accepts in-memory data, returns processed AnnData for plotting.

Algorithm Overview:
1. Subset AnnData to filtered VF genes
2. Normalize counts (CP10K)
3. Log-transform (log1p)
4. Scale data
5. PCA on filtered VFs
6. UMAP embedding
7. Preserve original cluster_coarse annotations

Input Requirements:
- AnnData object with raw counts (subset to passing cells)
- List of filtered VF gene names (from Step 04)
- Processing parameters dict (from module_configs.yaml)

Output Format:
- AnnData object with new X_umap, X_pca
- Preserves original cluster_coarse
- Ready for plotting

Reference:
- Step 02 processing pipeline (mimic parameters)
- Step 04 filtered VF lists
"""

import numpy as np
import scanpy as sc
from typing import List, Dict
import warnings


def reprocess_filtered_cells(adata_filtered, vf_list, processing_params):
    """
    Reprocess filtered cells using Step 04 filtered VF list.

    This function recreates the processing pipeline from Step 02, but using
    the filtered VF list from Step 04 (post-QC variable features) instead
    of the baseline VFs. This allows comparison of UMAP embeddings before
    and after VF filtering.

    Parameters:
    -----------
    adata_filtered : AnnData
        Subset to passing cells (has original cluster_coarse from Step 02).
        Must contain:
        - .X: Raw counts (or will be copied from .raw if normalized)
        - .obs['cluster_coarse']: Original cluster assignments from Step 02
        - .var_names: Gene names

    vf_list : list of str
        Filtered VF gene names from Step 04.
        These are the genes to use for PCA/UMAP.

    processing_params : dict
        Processing parameters from module_configs.yaml.
        Expected keys:
        - 'normalize_target_sum': int (default 1e4)
        - 'scale_max_value': float (default 10)
        - 'n_pcs': int (default 50)
        - 'use_highly_variable': bool (default False, we use filtered VFs)
        - 'n_neighbors': int (default 30)
        - 'umap_min_dist': float (default 0.3)

    Returns:
    --------
    AnnData
        Reprocessed object with:
        - .obsm['X_pca']: PCA coordinates (based on filtered VFs)
        - .obsm['X_umap']: UMAP coordinates (based on filtered VFs)
        - .obs['cluster_coarse']: Original cluster assignments (preserved)
        - .uns['processing_metadata']: Processing parameters used

    Raises:
    -------
    ValueError
        If required parameters missing or VF list empty
        If cluster_coarse missing from adata_filtered
        If insufficient VFs overlap with adata genes

    Example:
    --------
    >>> processing_params = {
    ...     'normalize_target_sum': 1e4,
    ...     'scale_max_value': 10,
    ...     'n_pcs': 50,
    ...     'use_highly_variable': False,
    ...     'n_neighbors': 30,
    ...     'umap_min_dist': 0.3
    ... }
    >>> adata_reprocessed = reprocess_filtered_cells(
    ...     adata_filtered=adata_pass,
    ...     vf_list=filtered_vfs,
    ...     processing_params=processing_params
    ... )
    >>> print(f"New UMAP shape: {adata_reprocessed.obsm['X_umap'].shape}")
    """

    # ========================================================================
    # Validation
    # ========================================================================

    if not isinstance(vf_list, (list, np.ndarray)):
        raise ValueError("vf_list must be a list or array of gene names")

    if len(vf_list) == 0:
        raise ValueError("vf_list is empty - cannot reprocess without VFs")

    if 'cluster_coarse' not in adata_filtered.obs.columns:
        raise ValueError(
            "cluster_coarse not found in adata_filtered.obs - "
            "cannot preserve original cluster assignments"
        )

    # Extract parameters with defaults
    normalize_target_sum = processing_params.get('normalize_target_sum', 1e4)
    scale_max_value = processing_params.get('scale_max_value', 10)
    n_pcs = processing_params.get('n_pcs', 50)
    n_neighbors = processing_params.get('n_neighbors', 30)
    umap_min_dist = processing_params.get('umap_min_dist', 0.3)

    print(f"\n--- Reprocessing Parameters ---")
    print(f"  Normalize target sum: {normalize_target_sum}")
    print(f"  Scale max value: {scale_max_value}")
    print(f"  N PCs: {n_pcs}")
    print(f"  N neighbors: {n_neighbors}")
    print(f"  UMAP min_dist: {umap_min_dist}")
    print()

    # ========================================================================
    # Copy AnnData to avoid modifying original
    # ========================================================================

    print(f"Creating copy of AnnData ({adata_filtered.n_obs} cells, {adata_filtered.n_vars} genes)...")
    adata = adata_filtered.copy()

    # Store original cluster assignments
    original_clusters = adata.obs['cluster_coarse'].copy()
    print(f"  Preserved {original_clusters.nunique()} original clusters")
    print()

    # ========================================================================
    # Subset to Filtered VFs
    # ========================================================================

    print(f"Subsetting to {len(vf_list)} filtered VFs...")

    # Find overlap with available genes
    vf_overlap = [gene for gene in vf_list if gene in adata.var_names]

    if len(vf_overlap) == 0:
        raise ValueError(
            f"No overlap between VF list ({len(vf_list)} genes) and "
            f"AnnData var_names ({adata.n_vars} genes)"
        )

    if len(vf_overlap) < len(vf_list):
        n_missing = len(vf_list) - len(vf_overlap)
        warnings.warn(
            f"{n_missing} VF genes not found in AnnData ({n_missing/len(vf_list)*100:.1f}% missing). "
            f"Using {len(vf_overlap)} overlapping VFs."
        )

    print(f"  {len(vf_overlap)} VFs found in AnnData ({len(vf_overlap)/len(vf_list)*100:.1f}% overlap)")

    # Subset to filtered VFs
    adata = adata[:, vf_overlap].copy()
    print(f"  Subset complete: {adata.n_obs} cells × {adata.n_vars} genes")
    print()

    # ========================================================================
    # Remove Zero-Count Genes
    # ========================================================================
    # CRITICAL: Union VFs may include genes with zero counts in merged dataset
    # (e.g., gene is variable in Sample 1 but has zero expression in Sample 2)
    # Must filter these out before normalization to prevent zero-variance → NaN

    print(f"Checking for zero-count genes after subsetting...")
    gene_counts = np.array(adata.X.sum(axis=0)).flatten()
    n_zero = (gene_counts == 0).sum()

    if n_zero > 0:
        print(f"  WARNING: Found {n_zero} genes with zero counts across all cells")
        print(f"  Removing {n_zero} zero-count genes ({n_zero/adata.n_vars*100:.1f}% of union VFs)")
        non_zero_mask = gene_counts > 0
        adata = adata[:, non_zero_mask].copy()
        print(f"  Filtered dataset: {adata.n_obs} cells × {adata.n_vars} genes")
    else:
        print(f"  All {adata.n_vars} genes have non-zero counts ✓")
    print()

    # ========================================================================
    # Remove Zero-Count Cells
    # ========================================================================
    # After removing zero-count genes, check for cells with zero total counts
    # This can happen if a cell's expressed genes were not in the filtered VF list
    # These cells will cause issues during normalization (division by zero)

    print(f"Checking for zero-count cells after gene filtering...")
    cell_counts = np.array(adata.X.sum(axis=1)).flatten()
    n_zero_cells = (cell_counts == 0).sum()

    if n_zero_cells > 0:
        print(f"  WARNING: Found {n_zero_cells} cells with zero counts across all VFs")
        print(f"  Removing {n_zero_cells} cells ({n_zero_cells/adata.n_obs*100:.2f}%)")
        non_zero_mask = cell_counts > 0
        adata = adata[non_zero_mask, :].copy()
        print(f"  Filtered dataset: {adata.n_obs} cells × {adata.n_vars} genes")
    else:
        print(f"  All {adata.n_obs} cells have non-zero counts ✓")
    print()

    # ========================================================================
    # Normalize
    # ========================================================================

    print(f"Normalizing to {normalize_target_sum} counts per cell...")
    sc.pp.normalize_total(adata, target_sum=normalize_target_sum)
    print()

    # ========================================================================
    # Log Transform
    # ========================================================================

    print(f"Log-transforming (log1p)...")
    sc.pp.log1p(adata)
    print()

    # ========================================================================
    # Scale
    # ========================================================================

    print(f"Scaling (max_value={scale_max_value})...")
    sc.pp.scale(adata, max_value=scale_max_value)

    # Safety check: Handle any remaining NaN values from zero-variance genes
    # (Should be rare after zero-count filtering, but possible if gene has same value in all cells)
    if np.isnan(adata.X).any():
        n_nan_genes = np.isnan(adata.X).any(axis=0).sum()
        print(f"  NOTE: Found {n_nan_genes} genes with NaN after scaling (constant expression)")
        print(f"  Replacing NaN with 0 (zero-variance genes have no information for PCA)...")
        adata.X = np.nan_to_num(adata.X, nan=0.0)
    print()

    # ========================================================================
    # PCA
    # ========================================================================

    # Adjust n_pcs if we have fewer genes than requested PCs
    actual_n_pcs = min(n_pcs, adata.n_vars - 1)
    if actual_n_pcs < n_pcs:
        warnings.warn(
            f"Requested {n_pcs} PCs but only {adata.n_vars} genes available. "
            f"Using {actual_n_pcs} PCs instead."
        )

    print(f"Computing PCA ({actual_n_pcs} components)...")
    sc.tl.pca(adata, n_comps=actual_n_pcs, use_highly_variable=False)
    print(f"  PCA shape: {adata.obsm['X_pca'].shape}")
    print()

    # ========================================================================
    # Neighbors Graph
    # ========================================================================

    print(f"Computing neighbors graph (n_neighbors={n_neighbors})...")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=actual_n_pcs)
    print()

    # ========================================================================
    # UMAP
    # ========================================================================

    print(f"Computing UMAP (min_dist={umap_min_dist})...")
    sc.tl.umap(adata, min_dist=umap_min_dist)
    print(f"  UMAP shape: {adata.obsm['X_umap'].shape}")
    print()

    # ========================================================================
    # Restore Original Cluster Assignments
    # ========================================================================

    print(f"Restoring original cluster assignments...")
    adata.obs['cluster_coarse'] = original_clusters
    print(f"  {adata.obs['cluster_coarse'].nunique()} clusters preserved")
    print()

    # ========================================================================
    # Store Processing Metadata
    # ========================================================================

    adata.uns['processing_metadata'] = {
        'reprocessed': True,
        'n_filtered_vfs_input': len(vf_list),
        'n_filtered_vfs_used': len(vf_overlap),
        'normalize_target_sum': normalize_target_sum,
        'scale_max_value': scale_max_value,
        'n_pcs': actual_n_pcs,
        'n_neighbors': n_neighbors,
        'umap_min_dist': umap_min_dist,
        'cluster_assignments': 'preserved_from_step_02'
    }

    print(f"--- Reprocessing Complete ---")
    print(f"  Cells: {adata.n_obs}")
    print(f"  Genes: {adata.n_vars}")
    print(f"  PCA components: {adata.obsm['X_pca'].shape[1]}")
    print(f"  UMAP coordinates: {adata.obsm['X_umap'].shape[1]}D")
    print(f"  Clusters preserved: {adata.obs['cluster_coarse'].nunique()}")
    print()

    return adata


def validate_reprocessing_result(adata_before, adata_after):
    """
    Validate that reprocessing preserved expected data structure.

    Parameters:
    -----------
    adata_before : AnnData
        Original AnnData (before reprocessing)

    adata_after : AnnData
        Reprocessed AnnData (after reprocessing)

    Returns:
    --------
    dict
        Validation results with checks:
        - 'n_cells_match': bool
        - 'clusters_preserved': bool
        - 'has_new_umap': bool
        - 'has_new_pca': bool

    Raises:
    -------
    ValueError
        If critical validation checks fail

    Example:
    --------
    >>> validation = validate_reprocessing_result(adata_before, adata_after)
    >>> assert validation['clusters_preserved'], "Clusters were not preserved!"
    """

    results = {}

    # Check cell count match
    results['n_cells_match'] = (adata_before.n_obs == adata_after.n_obs)
    if not results['n_cells_match']:
        n_removed = adata_before.n_obs - adata_after.n_obs
        pct_removed = (n_removed / adata_before.n_obs) * 100

        if n_removed > 0:
            # Cells were removed (expected if zero-count cells were filtered)
            print(f"  NOTE: {n_removed} cells removed during reprocessing ({pct_removed:.2f}%)")
            print(f"        This is expected if cells had zero counts across all VFs")

            # Only raise error if more than 1% of cells were removed (unexpected)
            if pct_removed > 1.0:
                warnings.warn(
                    f"More than 1% of cells ({pct_removed:.2f}%) were removed during reprocessing. "
                    f"This may indicate an issue with VF filtering. Before: {adata_before.n_obs}, "
                    f"After: {adata_after.n_obs}"
                )
        else:
            # More cells after than before (should never happen)
            raise ValueError(
                f"Cell count increased during reprocessing (should not happen): "
                f"before={adata_before.n_obs}, after={adata_after.n_obs}"
            )

    # Check cluster preservation
    if 'cluster_coarse' in adata_before.obs and 'cluster_coarse' in adata_after.obs:
        # Compare clusters for cells that exist in both objects
        # (cells may have been removed if they had zero counts)
        common_cells = adata_before.obs_names.intersection(adata_after.obs_names)

        if len(common_cells) > 0:
            clusters_before = adata_before[common_cells, :].obs['cluster_coarse']
            clusters_after = adata_after[common_cells, :].obs['cluster_coarse']
            clusters_match = (clusters_before == clusters_after).all()
            results['clusters_preserved'] = clusters_match

            if not clusters_match:
                raise ValueError(
                    "cluster_coarse assignments changed for retained cells during reprocessing!"
                )
        else:
            raise ValueError("No common cells between before and after - something went wrong!")
    else:
        results['clusters_preserved'] = False
        raise ValueError("cluster_coarse missing in one of the AnnData objects")

    # Check new embeddings exist
    results['has_new_umap'] = 'X_umap' in adata_after.obsm
    results['has_new_pca'] = 'X_pca' in adata_after.obsm

    if not results['has_new_umap']:
        raise ValueError("X_umap not found in reprocessed AnnData")
    if not results['has_new_pca']:
        raise ValueError("X_pca not found in reprocessed AnnData")

    print(f"--- Reprocessing Validation ---")
    print(f"  Cells match: {results['n_cells_match']} ({adata_after.n_obs} cells)")
    print(f"  Clusters preserved: {results['clusters_preserved']}")
    print(f"  New UMAP: {results['has_new_umap']}")
    print(f"  New PCA: {results['has_new_pca']}")
    print(f"  VALIDATION PASSED")
    print()

    return results
