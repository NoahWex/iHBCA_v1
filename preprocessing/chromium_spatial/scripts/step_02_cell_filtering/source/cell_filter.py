"""
Step 02 Cell Filtering - Pure Algorithm

Core logic for per-patient cell filtering using multi-threshold MAD system.
Extracted from archived workflows/qc/scripts/module_2_cell_filter.py

Key Features:
- Union VF calculation across patient samples
- Patient-level sample merging
- Coarse clustering (Leiden, resolution 0.1)
- Multi-metric MAD filtering (genes, UMI, mito)
- Debris pre-filtering for specified clusters
- OR logic: cell fails if ANY metric exceeds threshold

Author: Noah Wechter
Date: 2025-10-29
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import median_abs_deviation
from itertools import product


def calculate_qc_metrics(adata: sc.AnnData, patient_id: str, config: dict) -> None:
    """
    Calculate QC metrics including mitochondrial percentage.

    Args:
        adata: AnnData object (modified in place)
        patient_id: Patient identifier
        config: Configuration dictionary
    """
    print(f"{patient_id}: Calculating QC metrics (mito %)")

    # Identify mitochondrial genes (MT- prefix for human)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    n_mt_genes = adata.var['mt'].sum()
    print(f"{patient_id}:   Found {n_mt_genes} mitochondrial genes")

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=['mt'],
        percent_top=None,
        log1p=False,
        inplace=True
    )

    # Rename scanpy columns to algorithm's expected names
    # scanpy creates: n_genes_by_counts, total_counts, pct_counts_mt
    # algorithm expects: n_genes, n_umi, mito_percent
    if 'pct_counts_mt' in adata.obs.columns:
        adata.obs['mito_percent'] = adata.obs['pct_counts_mt']
    else:
        raise ValueError("pct_counts_mt not calculated by scanpy")

    # FIX: Preserve Step 01 n_umi if already present (full-matrix UMI)
    # total_counts reflects VF-subset UMI after VF subsetting in wrapper
    # Step 01 n_umi is preserved in wrapper before VF subsetting (~18k genes)
    # Only calculate from total_counts if n_umi not already set
    if 'n_umi' not in adata.obs.columns:
        if 'total_counts' in adata.obs.columns:
            adata.obs['n_umi'] = adata.obs['total_counts'].astype(int)
        else:
            raise ValueError("total_counts not calculated by scanpy")
    else:
        print(f"{patient_id}:   Using preserved n_umi from Step 01 (full-matrix UMI)")

    if 'n_genes_by_counts' in adata.obs.columns:
        adata.obs['n_genes'] = adata.obs['n_genes_by_counts'].astype(int)
    else:
        raise ValueError("n_genes_by_counts not calculated by scanpy")

    print(f"{patient_id}:   QC metrics calculated successfully")
    print(f"{patient_id}:     Mito % median: {adata.obs['mito_percent'].median():.2f}%")
    print(f"{patient_id}:     UMI median: {adata.obs['n_umi'].median():.0f}")
    print(f"{patient_id}:     Genes median: {adata.obs['n_genes'].median():.0f}")


def normalize_and_pca(adata: sc.AnnData, patient_id: str, config: dict) -> None:
    """
    Normalize counts and perform PCA.

    Args:
        adata: AnnData object (modified in place)
        patient_id: Patient identifier
        config: Configuration dictionary

    Raises:
        ValueError: If PCA fails
    """
    print(f"{patient_id}: Normalizing and running PCA")

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10)  # Z-score with cap at 10

    # Get n_components from config with default
    n_pcs = config.get('module_2', {}).get('pca_n_components', 50)
    sc.tl.pca(adata, n_comps=n_pcs, use_highly_variable=False)

    print(f"{patient_id}: PCA complete ({n_pcs} components)")

    # Validate PCA output
    if 'X_pca' not in adata.obsm:
        raise ValueError('PCA failed: X_pca not in obsm')

    if adata.obsm['X_pca'].shape[1] != n_pcs:
        raise ValueError(f"PCA shape mismatch: expected {n_pcs} PCs, "
                       f"got {adata.obsm['X_pca'].shape[1]}")


def coarse_clustering(adata: sc.AnnData, patient_id: str, config: dict) -> None:
    """
    Perform coarse clustering to establish broad cell type context.

    Args:
        adata: AnnData object (modified in place)
        patient_id: Patient identifier
        config: Configuration dictionary

    Raises:
        ValueError: If clustering fails
    """
    print(f"{patient_id}: Running coarse clustering")

    # Get parameters from config with defaults
    n_neighbors = config.get('module_2', {}).get('umap_n_neighbors', 30)
    n_pcs = config.get('module_2', {}).get('pca_n_components', 50)
    resolution = config.get('module_2', {}).get('clustering_resolution_coarse', 0.1)
    min_dist = config.get('module_2', {}).get('umap_min_dist', 0.3)

    # Compute neighbors
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    # Leiden clustering
    sc.tl.leiden(adata, resolution=resolution, key_added='cluster_coarse')

    # UMAP for visualization
    sc.tl.umap(adata, min_dist=min_dist)

    # Log cluster sizes
    cluster_sizes = adata.obs['cluster_coarse'].value_counts().sort_index()
    print(f"{patient_id}: {len(cluster_sizes)} clusters identified")
    print(f"  Cluster sizes: {dict(cluster_sizes)}")

    # Validate clustering
    if len(cluster_sizes) < 2:
        print(f"{patient_id}: Only {len(cluster_sizes)} cluster(s) - "
              "very homogeneous population")

    if cluster_sizes.max() / adata.n_obs > 0.8:
        print(f"{patient_id}: Largest cluster is "
              f"{cluster_sizes.max()/adata.n_obs*100:.1f}% of cells - "
              "may be over-dominant")

    # Validate UMAP
    if 'X_umap' not in adata.obsm:
        raise ValueError('UMAP failed: X_umap not in obsm')

    print(f"{patient_id}: Clustering complete")


def get_debris_clusters_for_patient(config: dict, patient_id: str) -> dict:
    """
    Get debris cluster configurations for this patient.

    Returns dict mapping cluster_id -> config dict with keys:
        - expected_size: int (required)
        - umi_threshold_mads: float (required)
    """
    debris_config = config.get('module_2', {}).get('debris_prefilter', {})

    if not debris_config.get('enabled', False):
        return {}

    clusters_config = debris_config.get('clusters', {})

    # Get cluster list for this patient
    if patient_id not in clusters_config:
        return {}

    patient_clusters = clusters_config[patient_id]

    # Parse cluster configs
    result = {}
    if isinstance(patient_clusters, list):
        for cluster_cfg in patient_clusters:
            if isinstance(cluster_cfg, dict):
                cluster_id = cluster_cfg['cluster_id']
                result[cluster_id] = {
                    'expected_size': cluster_cfg['expected_size'],
                    'umi_threshold_mads': cluster_cfg['umi_threshold_mads']
                }

    return result


def mad_filtering(adata: sc.AnnData, patient_id: str, config: dict) -> None:
    """
    Detect outliers within cluster-sample groups using multi-threshold MAD system.

    UPDATED 2025-10-22: Multi-metric, asymmetric thresholds + debris pre-filtering
    - Genes: Both lower AND upper (debris and doublets)
    - UMI: Lower only (remove low UMI debris)
    - Mito: Upper only (remove high mito stressed cells)
    - Debris pre-filtering: For configured clusters, remove cells <= median before MAD calc

    Cell fails if ANY metric exceeds its threshold(s) (OR logic)

    Args:
        adata: AnnData object (modified in place)
        patient_id: Patient identifier
        config: Configuration dictionary

    Raises:
        ValueError: If MAD filtering fails
    """
    print(f"{patient_id}: Running multi-threshold MAD filtering with optional debris pre-filtering")

    # Initialize output columns (now per-metric, per-direction)
    adata.obs['mad_outlier_genes_low'] = False
    adata.obs['mad_outlier_genes_high'] = False
    adata.obs['mad_outlier_umi_low'] = False
    adata.obs['mad_outlier_mito_high'] = False
    adata.obs['skip_mad_filtering'] = False
    adata.obs['cluster_sample_size'] = 0
    adata.obs['cell_filter_pass'] = True
    adata.obs['qc_failure_reason'] = ''  # Empty string instead of None (h5py compatibility)

    # Debris pre-filtering columns
    adata.obs['debris_prefilter_applied'] = False
    adata.obs['debris_prefilter_removed'] = False
    adata.obs['cell_filter_pass_final'] = True

    # Get MAD thresholds from config (nested dict)
    mad_thresholds = config.get('module_2', {}).get('mad_thresholds', {
        'genes': [-3, 3],
        'umi': [-3],
        'mito': [3]
    })
    min_cluster_sample_size = config.get('module_2', {}).get('min_cluster_sample_size', 100)
    mito_min_threshold = config.get('module_2', {}).get('mito_min_threshold', 3.0)

    # Get debris cluster configurations for this patient
    debris_clusters = get_debris_clusters_for_patient(config, patient_id)

    print(f"MAD thresholds:")
    print(f"  Genes: {mad_thresholds['genes']} (lower & upper)")
    print(f"  UMI: {mad_thresholds['umi']} (lower only)")
    print(f"  Mito: {mad_thresholds['mito']} (upper only, min {mito_min_threshold}%)")
    print(f"Min cluster-sample size: {min_cluster_sample_size}")
    print(f"Outlier logic: OR (cell fails if ANY metric exceeds threshold)")

    if debris_clusters:
        print(f"Debris pre-filtering enabled for {len(debris_clusters)} cluster(s): {list(debris_clusters.keys())}")

    # Track statistics
    n_groups_total = 0
    n_groups_skipped = 0
    n_groups_filtered = 0
    n_outliers_total = 0

    # Iterate over all cluster-sample combinations
    for (cluster, sample) in product(adata.obs['cluster_coarse'].unique(),
                                     adata.obs['sample_id'].unique()):
        n_groups_total += 1
        mask = (adata.obs['cluster_coarse'] == cluster) & (adata.obs['sample_id'] == sample)
        n_cells = mask.sum()

        # Store group size
        adata.obs.loc[mask, 'cluster_sample_size'] = n_cells

        # Skip if no cells in this group
        if n_cells == 0:
            continue

        # Check size threshold (>100 cells for MAD filtering)
        if n_cells <= min_cluster_sample_size:
            adata.obs.loc[mask, 'skip_mad_filtering'] = True
            adata.obs.loc[mask, 'cell_filter_pass'] = True
            adata.obs.loc[mask, 'qc_failure_reason'] = ''  # Empty string for h5py compatibility
            n_groups_skipped += 1
            continue

        n_groups_filtered += 1

        # === DEBRIS PRE-FILTERING (if configured for this cluster) ===
        analysis_mask = mask.copy()  # Mask of cells to include in MAD calculation
        # Convert cluster to int for comparison (cluster IDs are strings in adata.obs)
        cluster_int = int(cluster)
        is_debris_cluster = cluster_int in debris_clusters

        if is_debris_cluster:
            cluster_cfg = debris_clusters[cluster_int]
            expected_size = cluster_cfg['expected_size']
            umi_threshold_mads = cluster_cfg['umi_threshold_mads']

            print(f"  Cluster {cluster_int}, Sample {sample}: Debris pre-filtering enabled")
            print(f"    Config: expected_size={expected_size}, umi_threshold={umi_threshold_mads} MADs")

            # CRITICAL SAFETY: Validate cluster size
            cluster_only_mask = adata.obs['cluster_coarse'] == cluster
            actual_cluster_size = cluster_only_mask.sum()

            if actual_cluster_size != expected_size:
                error_msg = (
                    f"Cluster size validation FAILED for {patient_id} Cluster {cluster}!\n"
                    f"  Expected: {expected_size} cells\n"
                    f"  Actual: {actual_cluster_size} cells\n"
                    f"  Difference: {actual_cluster_size - expected_size:+d} cells\n"
                    f"\n"
                    f"This likely means clustering has changed since configuration was created.\n"
                    f"Aborting debris pre-filtering to prevent filtering the wrong cluster.\n"
                    f"\n"
                    f"ACTION REQUIRED:\n"
                    f"  1. Review current clustering in Module 2 outputs\n"
                    f"  2. Update expected_size in config if cluster still needs pre-filtering\n"
                    f"  3. Or disable debris_prefilter for this patient if no longer needed"
                )
                print(error_msg)
                raise ValueError(error_msg)

            print(f"    ✓ Cluster size validated: {actual_cluster_size} cells (matches expected)")

            # Mark all cells in this cluster as having pre-filter applied
            adata.obs.loc[mask, 'debris_prefilter_applied'] = True

            # Calculate cluster-wide median and MAD for UMI
            cluster_umi = adata.obs.loc[cluster_only_mask, 'n_umi'].values
            umi_median = np.median(cluster_umi)
            umi_mad = median_abs_deviation(cluster_umi)

            if umi_mad == 0:
                print(f"    Cluster {cluster} UMI MAD = 0, using median only")
                umi_threshold_value = umi_median
            else:
                # Threshold = median + (threshold_mads * MAD)
                # For threshold_mads = 0: removes cells <= median
                # For threshold_mads = -1: removes cells <= median - 1*MAD (more aggressive)
                umi_threshold_value = umi_median + (umi_threshold_mads * umi_mad)

            # Apply pre-filter: remove cells with UMI <= threshold
            umi_vals_this_group = adata.obs.loc[mask, 'n_umi'].values
            prefilter_remove = umi_vals_this_group <= umi_threshold_value
            n_prefilter_removed = prefilter_remove.sum()

            # Mark removed cells
            cell_indices = np.where(mask)[0]
            removed_indices = cell_indices[prefilter_remove]
            for idx in removed_indices:
                adata.obs.iloc[idx, adata.obs.columns.get_loc('debris_prefilter_removed')] = True
                adata.obs.iloc[idx, adata.obs.columns.get_loc('cell_filter_pass')] = False
                adata.obs.iloc[idx, adata.obs.columns.get_loc('cell_filter_pass_final')] = False
                # Set QC failure reason
                current_reason = adata.obs.iloc[idx, adata.obs.columns.get_loc('qc_failure_reason')]
                if pd.isna(current_reason) or current_reason == '':
                    adata.obs.iloc[idx, adata.obs.columns.get_loc('qc_failure_reason')] = 'low_umis'
                else:
                    adata.obs.iloc[idx, adata.obs.columns.get_loc('qc_failure_reason')] = f"{current_reason};low_umis"

            # Update analysis mask to exclude pre-filtered cells
            analysis_mask = mask & ~adata.obs['debris_prefilter_removed']

            print(f"    Pre-filter threshold: UMI <= {umi_threshold_value:.0f} "
                   f"(median={umi_median:.0f}, MAD={umi_mad:.0f}, threshold={umi_threshold_mads} MADs)")
            print(f"    Removed: {n_prefilter_removed}/{n_cells} cells ({100*n_prefilter_removed/n_cells:.1f}%)")
            print(f"    Remaining for MAD calculation: {analysis_mask.sum()}")

        # Multi-metric MAD calculation with error handling
        try:
            # Extract metric values (from analysis_mask, which excludes pre-filtered cells for debris clusters)
            umi_values = adata.obs.loc[analysis_mask, 'n_umi'].values
            gene_values = adata.obs.loc[analysis_mask, 'n_genes'].values
            mito_values = adata.obs.loc[analysis_mask, 'mito_percent'].values

            # Validate non-empty
            if len(umi_values) == 0:
                print(f"Cluster {cluster}, Sample {sample}: No cells")
                continue

            # === UMI: Lower threshold only ===
            umi_median = np.median(umi_values)
            umi_mad = median_abs_deviation(umi_values)

            if umi_mad == 0:
                print(f"Cluster {cluster}, Sample {sample}: UMI MAD=0")
                umi_outliers_low = np.zeros(len(umi_values), dtype=bool)
            else:
                if -3 in mad_thresholds['umi']:
                    umi_outliers_low = umi_values < (umi_median - 3 * umi_mad)
                else:
                    umi_outliers_low = np.zeros(len(umi_values), dtype=bool)

            # === Genes: Both lower and upper thresholds ===
            gene_median = np.median(gene_values)
            gene_mad = median_abs_deviation(gene_values)

            if gene_mad == 0:
                print(f"Cluster {cluster}, Sample {sample}: Gene MAD=0")
                gene_outliers_low = np.zeros(len(gene_values), dtype=bool)
                gene_outliers_high = np.zeros(len(gene_values), dtype=bool)
            else:
                if -3 in mad_thresholds['genes']:
                    gene_outliers_low = gene_values < (gene_median - 3 * gene_mad)
                else:
                    gene_outliers_low = np.zeros(len(gene_values), dtype=bool)

                if 3 in mad_thresholds['genes']:
                    gene_outliers_high = gene_values > (gene_median + 3 * gene_mad)
                else:
                    gene_outliers_high = np.zeros(len(gene_values), dtype=bool)

            # === Mito: Upper threshold only ===
            mito_median = np.median(mito_values)
            mito_mad = median_abs_deviation(mito_values)

            if mito_mad == 0:
                print(f"Cluster {cluster}, Sample {sample}: Mito MAD=0")
                mito_outliers_high = np.zeros(len(mito_values), dtype=bool)
                mito_threshold_final = mito_median + 3 * mito_mad  # Store for reference even if not used
            else:
                if 3 in mad_thresholds['mito']:
                    # Calculate MAD-based threshold
                    mito_mad_threshold = mito_median + 3 * mito_mad

                    # SAFETY: Minimum threshold (configurable, default 3%)
                    mito_threshold_final = max(mito_mad_threshold, mito_min_threshold)

                    mito_outliers_high = mito_values > mito_threshold_final

                    if mito_mad_threshold < mito_min_threshold:
                        pass  # Threshold raised to minimum
                else:
                    mito_outliers_high = np.zeros(len(mito_values), dtype=bool)
                    mito_threshold_final = mito_median + 3 * mito_mad  # Store for reference

            # Store individual outlier flags
            # Note: Use analysis_mask (not mask) because outlier arrays are based on
            # cells remaining after pre-filter. Pre-filtered cells keep default False values.
            adata.obs.loc[analysis_mask, 'mad_outlier_genes_low'] = gene_outliers_low
            adata.obs.loc[analysis_mask, 'mad_outlier_genes_high'] = gene_outliers_high
            adata.obs.loc[analysis_mask, 'mad_outlier_umi_low'] = umi_outliers_low
            adata.obs.loc[analysis_mask, 'mad_outlier_mito_high'] = mito_outliers_high

            # Combined outlier detection (OR logic)
            is_outlier = (gene_outliers_low | gene_outliers_high |
                         umi_outliers_low | mito_outliers_high)

            # Apply to analysis_mask cells only (may be subset if pre-filtered)
            adata.obs.loc[analysis_mask, 'cell_filter_pass'] = ~is_outlier
            adata.obs.loc[analysis_mask, 'cell_filter_pass_final'] = ~is_outlier

            # Detailed failure reason for MAD outliers
            failure_reasons = []
            if gene_outliers_low.any():
                failure_reasons.append('genes_low')
            if gene_outliers_high.any():
                failure_reasons.append('genes_high')
            if umi_outliers_low.any():
                failure_reasons.append('umi_low')
            if mito_outliers_high.any():
                failure_reasons.append('mito_high')

            reason = ','.join(failure_reasons) if failure_reasons else ''  # Empty string for h5py compatibility

            # Update QC failure reason for MAD outliers (append if already has low_umis from pre-filter)
            for idx in np.where(analysis_mask)[0]:
                if is_outlier[adata.obs.index[idx] == adata.obs.loc[analysis_mask].index].any():
                    current_reason = adata.obs.iloc[idx, adata.obs.columns.get_loc('qc_failure_reason')]
                    if pd.isna(current_reason) or current_reason == '':
                        adata.obs.iloc[idx, adata.obs.columns.get_loc('qc_failure_reason')] = reason
                    elif reason and reason not in str(current_reason):
                        adata.obs.iloc[idx, adata.obs.columns.get_loc('qc_failure_reason')] = f"{current_reason};{reason}"

            n_outliers = is_outlier.sum()
            n_outliers_total += n_outliers

            # Store thresholds for ECDF plots
            if 'mad_thresholds' not in adata.uns:
                adata.uns['mad_thresholds'] = {}
            adata.uns['mad_thresholds'][f'{cluster}_{sample}'] = {
                'n_cells': int(n_cells),
                'umi_median': float(umi_median),
                'umi_mad': float(umi_mad),
                'umi_lower': float(umi_median - 3 * umi_mad),
                'gene_median': float(gene_median),
                'gene_mad': float(gene_mad),
                'gene_lower': float(gene_median - 3 * gene_mad),
                'gene_upper': float(gene_median + 3 * gene_mad),
                'mito_median': float(mito_median),
                'mito_mad': float(mito_mad),
                'mito_upper': float(mito_threshold_final) if 'mito_threshold_final' in locals() else float(mito_median + 3 * mito_mad),
                'mito_upper_mad_only': float(mito_median + 3 * mito_mad),
                'n_outliers': int(n_outliers),
                'outlier_rate': float(n_outliers / n_cells),
                'outlier_breakdown': {
                    'genes_low': int(gene_outliers_low.sum()),
                    'genes_high': int(gene_outliers_high.sum()),
                    'umi_low': int(umi_outliers_low.sum()),
                    'mito_high': int(mito_outliers_high.sum())
                }
            }

        except Exception as e:
            print(f"Cluster {cluster}, Sample {sample}: MAD calculation failed: {e}")
            # Mark as skip_mad_filtering on error (don't fail cells)
            adata.obs.loc[mask, 'skip_mad_filtering'] = True
            adata.obs.loc[mask, 'cell_filter_pass'] = True
            continue

    # Log summary with breakdown by metric
    print(f"{patient_id}: Multi-threshold MAD filtering complete")
    print(f"  Total groups: {n_groups_total}")
    print(f"  Groups skipped (<={min_cluster_sample_size} cells): {n_groups_skipped}")
    print(f"  Groups filtered: {n_groups_filtered}")
    print(f"  Total outliers removed: {n_outliers_total:,} / {adata.n_obs:,} "
          f"({n_outliers_total/adata.n_obs*100:.1f}%)")

    # Breakdown by metric
    print(f"  Outlier breakdown:")
    print(f"    Genes (low): {adata.obs['mad_outlier_genes_low'].sum():,}")
    print(f"    Genes (high): {adata.obs['mad_outlier_genes_high'].sum():,}")
    print(f"    UMI (low): {adata.obs['mad_outlier_umi_low'].sum():,}")
    print(f"    Mito (high): {adata.obs['mad_outlier_mito_high'].sum():,}")

    # Validate output columns (updated for new multi-threshold columns)
    required_cols = ['mad_outlier_genes_low', 'mad_outlier_genes_high',
                     'mad_outlier_umi_low', 'mad_outlier_mito_high',
                     'skip_mad_filtering', 'cell_filter_pass', 'cluster_sample_size']
    for col in required_cols:
        if col not in adata.obs.columns:
            raise ValueError(f"MAD filtering failed: {col} not in obs")
        if adata.obs[col].isnull().any():
            raise ValueError(f"MAD filtering failed: NaN in {col}")


def run_cell_filtering(
    vf_sets: list,  # List of sets, one per sample
    adata_list: list,  # List of loaded AnnData objects (one per sample)
    patient_id: str,
    config: dict  # Full module_configs
) -> sc.AnnData:
    """
    Run complete cell filtering pipeline.

    Args:
        vf_sets: List of VF sets (one per sample)
        adata_list: List of AnnData objects (one per sample, pre-loaded)
        patient_id: Patient identifier
        config: Configuration dictionary (from module_configs.yaml)

    Returns:
        Patient-level AnnData with:
        - Merged samples
        - Cluster assignments (cluster_coarse)
        - MAD filtering results (outlier flags, pass/fail)
        - QC metrics (n_umi, n_genes, mito_percent)
    """
    # 1. Create VF union
    union_vfs_set = set.union(*vf_sets)
    union_vfs = sorted(list(union_vfs_set))
    print(f"{patient_id}: Union VFs = {len(union_vfs):,} genes")

    # 2. Merge samples (assume already subset to union VFs)
    # Samples should be pre-filtered and subset to union VFs before passing to this function
    adata_patient = sc.concat(adata_list, join='outer', merge='same')
    print(f"{patient_id}: Merged AnnData: {adata_patient.n_obs:,} cells × {adata_patient.n_vars:,} genes")

    # 3. Calculate QC metrics
    calculate_qc_metrics(adata_patient, patient_id, config)

    # 4. Normalize and PCA
    normalize_and_pca(adata_patient, patient_id, config)

    # 5. Coarse clustering
    coarse_clustering(adata_patient, patient_id, config)

    # 6. MAD filtering
    mad_filtering(adata_patient, patient_id, config)

    return adata_patient
