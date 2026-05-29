"""
Step 05 Filtered Preview - Object Rebuilder

Pure algorithm for rebuilding patient-level AnnData objects from in-memory data.
NO file I/O - all inputs are in-memory objects.

Author: Agent 3 (Algorithm Developer)
Date: 2025-11-03
"""

from typing import List, Dict
import numpy as np
import pandas as pd
import scanpy as sc


def rebuild_filtered_object(
    adata_list: List[sc.AnnData],
    filtered_vfs: List[str],
    cluster_labels: pd.Series,
    processing_params: dict
) -> sc.AnnData:
    """
    Rebuild patient-level object from passing cells with filtered VFs.

    This function reconstructs the patient-level merged object exactly as Step 02
    did, but using filtered VFs instead of baseline VFs and only passing cells.

    Args:
        adata_list: List of per-sample AnnData objects (passing cells only, raw counts)
        filtered_vfs: List of genes (union of Step 04 filtered VFs)
        cluster_labels: Series mapping cell barcodes to cluster IDs (from Step 02)
        processing_params: Dict with normalize_target_sum, scale_max_value, n_pcs, etc.

    Returns:
        AnnData with reprocessed data on filtered VFs, including:
            - .obs['cluster_coarse']: Transferred from Step 02 (NOT re-clustered)
            - .obsm['X_umap']: Recomputed on filtered VFs
            - .uns['processing_params']: Processing parameters used

    Example:
        >>> # adata_list = [adata_sample1, adata_sample2, ...]
        >>> # cluster_labels = pd.Series(data={barcode: cluster_id})
        >>> params = {
        ...     'normalize_target_sum': 1e4,
        ...     'scale_max_value': 10,
        ...     'n_pcs': 50,
        ...     'n_neighbors': 30,
        ...     'umap_min_dist': 0.3
        ... }
        >>> adata_rebuilt = rebuild_filtered_object(
        ...     adata_list, filtered_vfs, cluster_labels, params
        ... )
    """

    # Step 1: Subset each sample to filtered VFs
    print(f"  Subsetting {len(adata_list)} samples to {len(filtered_vfs)} filtered VFs...")
    adata_list_subset = []
    for i, adata in enumerate(adata_list):
        # Find intersection of genes
        common_genes = adata.var_names.intersection(filtered_vfs)
        print(f"    Sample {i+1}: {len(common_genes)}/{len(filtered_vfs)} VFs present")
        adata_subset = adata[:, common_genes].copy()
        adata_list_subset.append(adata_subset)

    # Step 2: Merge samples
    print(f"  Merging {len(adata_list_subset)} samples...")
    adata_patient = sc.concat(adata_list_subset, join='outer', merge='same')
    print(f"  Merged: {adata_patient.n_obs} cells × {adata_patient.n_vars} genes")

    # Step 3: Transfer cluster labels from Step 02
    print(f"  Transferring cluster labels from Step 02...")
    # Join cluster labels to .obs
    adata_patient.obs = adata_patient.obs.join(
        cluster_labels.to_frame(name='cluster_coarse'),
        how='left'
    )

    # Validate no NaN clusters
    n_missing = adata_patient.obs['cluster_coarse'].isna().sum()
    if n_missing > 0:
        raise ValueError(f"Missing cluster labels for {n_missing} cells!")

    # Convert to int (may be float after join)
    adata_patient.obs['cluster_coarse'] = adata_patient.obs['cluster_coarse'].astype(int)

    print(f"    Cluster labels transferred: {adata_patient.obs['cluster_coarse'].nunique()} clusters")

    # Step 4: Apply Step 02 processing pipeline
    print(f"  Applying Step 02 processing pipeline...")

    # 4.1: Normalize (CP10K)
    normalize_target_sum = processing_params.get('normalize_target_sum', 1e4)
    print(f"    Normalizing (target sum = {normalize_target_sum})...")
    sc.pp.normalize_total(adata_patient, target_sum=normalize_target_sum)

    # 4.2: Log-transform
    print(f"    Log-transforming...")
    sc.pp.log1p(adata_patient)

    # 4.3: Scale (Z-score with cap)
    scale_max_value = processing_params.get('scale_max_value', 10)
    print(f"    Scaling (max value = {scale_max_value})...")
    sc.pp.scale(adata_patient, max_value=scale_max_value)

    # 4.4: PCA
    n_pcs = processing_params.get('n_pcs', 50)
    print(f"    Running PCA ({n_pcs} components)...")
    sc.tl.pca(adata_patient, n_comps=n_pcs, use_highly_variable=False)

    # 4.5: Neighbors
    n_neighbors = processing_params.get('n_neighbors', 30)
    n_pcs_neighbors = processing_params.get('n_pcs_neighbors', 50)
    print(f"    Computing neighbors (n_neighbors={n_neighbors}, n_pcs={n_pcs_neighbors})...")
    sc.pp.neighbors(adata_patient, n_neighbors=n_neighbors, n_pcs=n_pcs_neighbors)

    # 4.6: UMAP (recompute on filtered VFs)
    umap_min_dist = processing_params.get('umap_min_dist', 0.3)
    print(f"    Computing UMAP (min_dist={umap_min_dist})...")
    sc.tl.umap(adata_patient, min_dist=umap_min_dist)

    # Validate UMAP was created
    if 'X_umap' not in adata_patient.obsm:
        raise ValueError("UMAP computation failed - X_umap not in obsm")

    print(f"    UMAP computed: shape {adata_patient.obsm['X_umap'].shape}")

    # Step 5: Store processing parameters in .uns
    adata_patient.uns['processing_params'] = {
        'normalize_target_sum': normalize_target_sum,
        'log1p': True,
        'scale_max_value': scale_max_value,
        'n_pcs': n_pcs,
        'n_neighbors': n_neighbors,
        'n_pcs_neighbors': n_pcs_neighbors,
        'umap_min_dist': umap_min_dist,
        'vfs_used': 'step04_filtered_union',
        'preserve_cluster_labels': True,
        'source': 'step_05_filtered_preview'
    }

    print(f"  Object rebuilding complete!")
    print(f"    Final shape: {adata_patient.n_obs} cells × {adata_patient.n_vars} genes")
    print(f"    Clusters: {adata_patient.obs['cluster_coarse'].nunique()}")

    return adata_patient


def add_metadata_to_rebuilt_object(
    adata_rebuilt: sc.AnnData,
    sample_metadata: Dict[str, pd.DataFrame]
) -> sc.AnnData:
    """
    Add sample-level metadata to rebuilt object.

    Helper function to merge sample metadata (from Step 02, Step 03, etc.)
    into the rebuilt object's .obs.

    Args:
        adata_rebuilt: Rebuilt AnnData object (modified in place)
        sample_metadata: Dict mapping column names to DataFrames indexed by cell barcode

    Returns:
        Modified adata_rebuilt with metadata columns added

    Example:
        >>> # Sample metadata from Step 02/03
        >>> metadata = {
        ...     'cell_filter_pass_final': df_step02[['cell_filter_pass_final']],
        ...     'doublet_filter_pass': df_step03[['doublet_filter_pass']]
        ... }
        >>> adata_rebuilt = add_metadata_to_rebuilt_object(adata_rebuilt, metadata)
    """

    print(f"  Adding metadata to rebuilt object...")

    for col_name, df in sample_metadata.items():
        print(f"    Adding column: {col_name}")

        # Join metadata
        adata_rebuilt.obs = adata_rebuilt.obs.join(df, how='left')

        # Check for missing values
        n_missing = adata_rebuilt.obs[col_name].isna().sum()
        if n_missing > 0:
            print(f"      Warning: {n_missing} cells missing {col_name}")

    print(f"  Metadata added successfully")

    return adata_rebuilt


def validate_rebuilt_object(
    adata_rebuilt: sc.AnnData,
    expected_n_cells: int = None,
    expected_n_clusters: int = None
) -> bool:
    """
    Validate rebuilt object has expected structure.

    Args:
        adata_rebuilt: Rebuilt AnnData object
        expected_n_cells: Expected number of cells (optional)
        expected_n_clusters: Expected number of clusters (optional)

    Returns:
        True if validation passes

    Raises:
        ValueError: If validation fails

    Example:
        >>> validate_rebuilt_object(adata_rebuilt, expected_n_cells=3352, expected_n_clusters=14)
    """

    print(f"  Validating rebuilt object...")

    # Check required columns
    required_cols = ['sample_id', 'patient_id', 'cluster_coarse']
    missing_cols = [col for col in required_cols if col not in adata_rebuilt.obs.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    # Check UMAP
    if 'X_umap' not in adata_rebuilt.obsm:
        raise ValueError("Missing UMAP coordinates (X_umap)")

    if adata_rebuilt.obsm['X_umap'].shape[1] != 2:
        raise ValueError(f"UMAP has wrong shape: {adata_rebuilt.obsm['X_umap'].shape}")

    # Check processing params
    if 'processing_params' not in adata_rebuilt.uns:
        raise ValueError("Missing processing_params in .uns")

    # Check cluster labels
    if adata_rebuilt.obs['cluster_coarse'].isna().any():
        raise ValueError("NaN values in cluster_coarse")

    # Check expected values
    if expected_n_cells is not None:
        actual_n_cells = adata_rebuilt.n_obs
        if actual_n_cells != expected_n_cells:
            raise ValueError(f"Cell count mismatch: expected {expected_n_cells}, got {actual_n_cells}")

    if expected_n_clusters is not None:
        actual_n_clusters = adata_rebuilt.obs['cluster_coarse'].nunique()
        if actual_n_clusters != expected_n_clusters:
            raise ValueError(f"Cluster count mismatch: expected {expected_n_clusters}, got {actual_n_clusters}")

    print(f"    Validation passed!")
    print(f"      Cells: {adata_rebuilt.n_obs}")
    print(f"      Genes: {adata_rebuilt.n_vars}")
    print(f"      Clusters: {adata_rebuilt.obs['cluster_coarse'].nunique()}")

    return True
