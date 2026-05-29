"""
Generate Preview Plots - Pure Computation Module (Step 05 Specific)

This module contains PURE functions for generating Step 05 filtered preview plots.
NO FILE I/O - accepts in-memory data, returns file path dict for manifest writing.

Algorithm Overview:
Generate 11 preview plots of the final filtered dataset (before and after reprocessing):
1. UMAP by cluster
2. UMAP by sample
3. UMAP faceted by sample (colored by cluster)
4. UMAP faceted by cluster (colored by sample)
5-7. Feature plots: n_genes_by_counts, total_counts, pct_counts_mt (raw values)
8-10. Feature plots: same metrics (MAD-standardized units)
11. Feature plot: scDblFinder_score

Input Requirements:
- AnnData object (FILTERED - only cells passing both filters)
- Patient ID string
- Output directory Path

Output Format:
- Dict mapping plot names to file paths
- All plots saved to output_dir/before/ or output_dir/after/
- Returns {plot_name: file_path} for manifest

Reference:
- Agent 1 spec: Section 2 (Step 05 - Filtered Preview)
- Plot specifications: Section 2, plots 1-11
- Updated for before/after reprocessing comparison
"""

import numpy as np
from pathlib import Path
from plotting_core import (
    plot_umap_by_category,
    plot_umap_faceted_by_sample,
    plot_umap_faceted_by_cluster,
    plot_feature_umap,
    calculate_percentile_cap
)
from mad_standardization import compute_all_mad_metrics


def generate_all_preview_plots(adata, patient_id, output_dir):
    """
    Generate all 11 preview plots for Step 05.

    Parameters:
    -----------
    adata : AnnData
        FILTERED object (only cells passing both Step 02 and Step 03 filters)
        Must have:
        - .obsm['X_umap']: UMAP coordinates
        - .obs['cluster_coarse']: Cluster assignments
        - .obs['sample_id']: Sample IDs
        - .obs['n_genes_by_counts']: Gene count metric
        - .obs['total_counts']: UMI count metric
        - .obs['pct_counts_mt']: Mitochondrial percentage
        - .obs['scDblFinder_score']: Doublet score

    patient_id : str
        Patient ID for plot titles and filenames

    output_dir : Path or str
        Directory to save plots (will be created if doesn't exist)

    Returns:
    --------
    dict
        {plot_name: file_path} mapping for manifest writing

    Plots Created:
    --------------
    1. {patient_id}_umap_by_cluster.png
    2. {patient_id}_umap_by_sample.png
    3. {patient_id}_umap_facet_sample_by_cluster.png
    4. {patient_id}_umap_facet_cluster_by_sample.png
    5. {patient_id}_feature_n_genes_by_counts.png
    6. {patient_id}_feature_total_counts.png
    7. {patient_id}_feature_pct_counts_mt.png
    8. {patient_id}_feature_n_genes_by_counts_mad.png
    9. {patient_id}_feature_total_counts_mad.png
    10. {patient_id}_feature_pct_counts_mito_mad.png
    11. {patient_id}_feature_scDblFinder_score.png

    Raises:
    -------
    ValueError
        If required columns or UMAP missing

    Example:
    --------
    >>> from pathlib import Path
    >>> plot_paths = generate_all_preview_plots(
    >>>     adata_filtered,
    >>>     patient_id='Pat1',
    >>>     output_dir=Path('outputs/05_FilteredPreview/plots')
    >>> )
    >>> print(f"Generated {len(plot_paths)} plots")
    >>> for name, path in plot_paths.items():
    >>>     print(f"  {name}: {path}")
    """
    # Convert output_dir to Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Validate required columns
    required_cols = [
        'cluster_coarse', 'sample_id',
        'n_genes_by_counts', 'total_counts', 'pct_counts_mt',
        'scDblFinder_score'
    ]
    for col in required_cols:
        if col not in adata.obs.columns:
            raise ValueError(f"Required column '{col}' not found in adata.obs")

    # Validate UMAP presence
    if 'X_umap' not in adata.obsm:
        raise ValueError(
            f"{patient_id}: UMAP coordinates missing. "
            "Cannot generate preview plots without UMAP."
        )

    # Initialize plot path dictionary
    plot_paths = {}

    # ========================================================================
    # COMPUTE MAD-STANDARDIZED METRICS (for plots 8-10)
    # ========================================================================
    print(f"{patient_id}: Computing MAD-standardized metrics...")
    compute_all_mad_metrics(adata)  # Adds {metric}_mad columns to adata.obs

    # ========================================================================
    # PLOT 1: UMAP by Cluster
    # ========================================================================
    print(f"{patient_id}: Generating plot 1/11 - UMAP by cluster...")
    plot_name = f"{patient_id}_umap_by_cluster.png"
    plot_path = output_dir / plot_name

    plot_umap_by_category(
        adata=adata,
        color_by='cluster_coarse',
        title=f"{patient_id} - Filtered Cells by Cluster",
        save_path=plot_path,
        palette='tab10',
        legend_loc='on data',
        figsize=(8, 8)
    )

    plot_paths['umap_by_cluster'] = str(plot_path)

    # ========================================================================
    # PLOT 2: UMAP by Sample
    # ========================================================================
    print(f"{patient_id}: Generating plot 2/11 - UMAP by sample...")
    plot_name = f"{patient_id}_umap_by_sample.png"
    plot_path = output_dir / plot_name

    plot_umap_by_category(
        adata=adata,
        color_by='sample_id',
        title=f"{patient_id} - Filtered Cells by Sample",
        save_path=plot_path,
        palette='tab20',
        legend_loc='right margin',
        figsize=(10, 8)  # Wider for legend
    )

    plot_paths['umap_by_sample'] = str(plot_path)

    # ========================================================================
    # PLOT 3: UMAP Faceted by Sample (colored by cluster)
    # ========================================================================
    print(f"{patient_id}: Generating plot 3/11 - UMAP faceted by sample...")
    plot_name = f"{patient_id}_umap_facet_sample_by_cluster.png"
    plot_path = output_dir / plot_name

    plot_umap_faceted_by_sample(
        adata=adata,
        color_by='cluster_coarse',
        title_prefix=patient_id,
        save_dir=output_dir,
        palette='tab10'
    )

    plot_paths['umap_facet_sample'] = str(plot_path)

    # ========================================================================
    # PLOT 4: UMAP Faceted by Cluster (colored by sample)
    # ========================================================================
    print(f"{patient_id}: Generating plot 4/11 - UMAP faceted by cluster...")
    plot_name = f"{patient_id}_umap_facet_cluster_by_sample.png"
    plot_path = output_dir / plot_name

    plot_umap_faceted_by_cluster(
        adata=adata,
        color_by='sample_id',
        title_prefix=patient_id,
        save_dir=output_dir,
        palette='tab20'
    )

    plot_paths['umap_facet_cluster'] = str(plot_path)

    # ========================================================================
    # PLOTS 5-7: Feature Plots (Raw Metrics)
    # ========================================================================
    feature_plots_raw = [
        ('n_genes_by_counts', 'Gene Count', 'viridis'),
        ('total_counts', 'UMI Count', 'viridis'),
        ('pct_counts_mt', 'Mitochondrial %', 'Reds')
    ]

    for idx, (metric, label, cmap) in enumerate(feature_plots_raw, start=5):
        print(f"{patient_id}: Generating plot {idx}/11 - Feature {metric}...")

        plot_name = f"{patient_id}_feature_{metric}.png"
        plot_path = output_dir / plot_name

        # Calculate 95th percentile cap (per Agent 1 spec)
        vmax = calculate_percentile_cap(adata, metric, percentile=95)

        plot_feature_umap(
            adata=adata,
            feature_col=metric,
            title=f"{patient_id} - {label}",
            save_path=plot_path,
            vmin=0,
            vmax=vmax,
            cmap=cmap,
            figsize=(8, 8)
        )

        plot_paths[f'feature_{metric}'] = str(plot_path)

    # ========================================================================
    # PLOTS 8-10: Feature Plots (MAD-Standardized Metrics)
    # ========================================================================
    feature_plots_mad = [
        ('n_genes_by_counts', 'Gene Count'),
        ('total_counts', 'UMI Count'),
        ('pct_counts_mt', 'Mitochondrial %')
    ]

    for idx, (metric, label) in enumerate(feature_plots_mad, start=8):
        print(f"{patient_id}: Generating plot {idx}/11 - Feature {metric} (MAD)...")

        mad_col = f"{metric}_mad"
        plot_name = f"{patient_id}_feature_{metric}_mad.png"
        plot_path = output_dir / plot_name

        # MAD plots: diverging colormap, centered at 0, range [-3, +3]
        plot_feature_umap(
            adata=adata,
            feature_col=mad_col,
            title=f"{patient_id} - {label} (MAD units)",
            save_path=plot_path,
            vmin=-3,
            vmax=3,
            vcenter=0,
            cmap='RdYlBu_r',  # Red=low, Yellow=center, Blue=high
            figsize=(8, 8)
        )

        plot_paths[f'feature_{metric}_mad'] = str(plot_path)

    # ========================================================================
    # PLOT 11: Feature Plot - Doublet Score
    # ========================================================================
    print(f"{patient_id}: Generating plot 11/11 - Feature scDblFinder_score...")
    plot_name = f"{patient_id}_feature_scDblFinder_score.png"
    plot_path = output_dir / plot_name

    # Doublet score: fixed range [0, 1]
    plot_feature_umap(
        adata=adata,
        feature_col='scDblFinder_score',
        title=f"{patient_id} - Doublet Score (Filtered Cells)",
        save_path=plot_path,
        vmin=0,
        vmax=1,
        cmap='plasma',
        figsize=(8, 8)
    )

    plot_paths['feature_scDblFinder_score'] = str(plot_path)

    # ========================================================================
    # SUMMARY
    # ========================================================================
    print(f"{patient_id}: Successfully generated {len(plot_paths)} plots")

    return plot_paths


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def validate_adata_for_plotting(adata, patient_id):
    """
    Validate AnnData object has all required components for plotting.

    Parameters:
    -----------
    adata : AnnData
        Object to validate

    patient_id : str
        Patient ID for error messages

    Raises:
    -------
    ValueError
        If validation fails

    Returns:
    --------
    dict
        Validation summary:
        - 'n_cells': Number of cells
        - 'n_clusters': Number of clusters
        - 'n_samples': Number of samples
        - 'has_umap': UMAP presence
    """
    # Check UMAP
    if 'X_umap' not in adata.obsm:
        raise ValueError(f"{patient_id}: Missing UMAP coordinates")

    # Check required columns
    required = [
        'cluster_coarse', 'sample_id',
        'n_genes_by_counts', 'total_counts', 'pct_counts_mt',
        'scDblFinder_score'
    ]
    missing = [col for col in required if col not in adata.obs.columns]
    if missing:
        raise ValueError(f"{patient_id}: Missing columns: {missing}")

    # Return summary
    return {
        'n_cells': adata.n_obs,
        'n_clusters': adata.obs['cluster_coarse'].nunique(),
        'n_samples': adata.obs['sample_id'].nunique(),
        'has_umap': True
    }


def generate_before_after_plots(adata_before, adata_after, patient_id, output_base):
    """
    Generate 11 plots for both before and after reprocessing states.

    This function creates two complete sets of preview plots:
    - Before: Using original Step 02 UMAP/PCA (baseline VFs)
    - After: Using reprocessed UMAP/PCA (filtered VFs from Step 04)

    Both sets show the SAME cells, preserving cluster_coarse assignments,
    but with different underlying embeddings to demonstrate the effect
    of VF filtering on visualization.

    Parameters:
    -----------
    adata_before : AnnData
        Original Step 02 data (filtered to passing cells).
        Uses baseline VFs for UMAP/PCA.

    adata_after : AnnData
        Reprocessed data (same cells, filtered VFs).
        Uses Step 04 filtered VFs for UMAP/PCA.

    patient_id : str
        Patient ID for plot titles and filenames

    output_base : Path or str
        Base output directory. Will create before/ and after/ subdirs.

    Returns:
    --------
    dict
        Nested structure:
        {
            'before': {plot_name: file_path, ...},
            'after': {plot_name: file_path, ...},
            'n_plots_before': int,
            'n_plots_after': int,
            'n_plots_total': int
        }

    Raises:
    -------
    ValueError
        If adata_before and adata_after have different cell counts
        If cluster assignments don't match between before/after

    Example:
    --------
    >>> plot_results = generate_before_after_plots(
    ...     adata_before=adata_original,
    ...     adata_after=adata_reprocessed,
    ...     patient_id='Pat1',
    ...     output_base=Path('outputs/05_FilteredPreview/plots')
    ... )
    >>> print(f"Total plots: {plot_results['n_plots_total']}")
    >>> print(f"Before plots: {len(plot_results['before'])}")
    >>> print(f"After plots: {len(plot_results['after'])}")
    """

    # Validate inputs
    if adata_before.n_obs != adata_after.n_obs:
        raise ValueError(
            f"Cell count mismatch: before={adata_before.n_obs}, after={adata_after.n_obs}"
        )

    if not (adata_before.obs['cluster_coarse'] == adata_after.obs['cluster_coarse']).all():
        raise ValueError("cluster_coarse assignments differ between before and after")

    # Convert output_base to Path
    output_base = Path(output_base)

    # Create subdirectories
    output_dirs = {
        'before': output_base / 'before',
        'after': output_base / 'after'
    }
    for dir_path in output_dirs.values():
        dir_path.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*80}")
    print(f"Generating Before/After Plots for {patient_id}")
    print(f"{'='*80}")
    print(f"Output directories:")
    print(f"  Before: {output_dirs['before']}")
    print(f"  After:  {output_dirs['after']}")
    print(f"Cell count: {adata_before.n_obs}")
    print(f"{'='*80}\n")

    # ========================================================================
    # Generate BEFORE plots (original Step 02 UMAP)
    # ========================================================================

    print(f"[1/2] Generating BEFORE plots (baseline VFs)...")
    before_plots = generate_all_preview_plots(
        adata=adata_before,
        patient_id=patient_id,
        output_dir=output_dirs['before']
    )
    print(f"  Generated {len(before_plots)} before plots")
    print()

    # ========================================================================
    # Generate AFTER plots (reprocessed with filtered VFs)
    # ========================================================================

    print(f"[2/2] Generating AFTER plots (filtered VFs)...")
    after_plots = generate_all_preview_plots(
        adata=adata_after,
        patient_id=patient_id,
        output_dir=output_dirs['after']
    )
    print(f"  Generated {len(after_plots)} after plots")
    print()

    # ========================================================================
    # Compile results
    # ========================================================================

    results = {
        'before': before_plots,
        'after': after_plots,
        'n_plots_before': len(before_plots),
        'n_plots_after': len(after_plots),
        'n_plots_total': len(before_plots) + len(after_plots)
    }

    print(f"\n{'='*80}")
    print(f"Before/After Plot Generation Complete")
    print(f"{'='*80}")
    print(f"Total plots created: {results['n_plots_total']}")
    print(f"  Before: {results['n_plots_before']}")
    print(f"  After:  {results['n_plots_after']}")
    print(f"{'='*80}\n")

    return results


def generate_plot_manifest(plot_paths, patient_id):
    """
    Generate plot manifest dictionary for YAML output.

    Parameters:
    -----------
    plot_paths : dict
        {plot_name: file_path} from generate_all_preview_plots

    patient_id : str
        Patient ID

    Returns:
    --------
    dict
        Manifest structure for writing

    Example:
    --------
    >>> manifest = generate_plot_manifest(plot_paths, 'Pat1')
    >>> import yaml
    >>> with open('manifest.yaml', 'w') as f:
    >>>     yaml.dump(manifest, f)
    """
    return {
        'patient_id': patient_id,
        'n_plots_generated': len(plot_paths),
        'plot_types': list(plot_paths.keys()),
        'plots': plot_paths
    }
