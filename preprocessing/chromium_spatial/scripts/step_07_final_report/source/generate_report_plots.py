"""
Generate Report Plots - Pure Computation Module (Step 07 Specific)

This module contains PURE functions for generating Step 07 before/after comparison plots.
NO FILE I/O - accepts in-memory data, returns file path dict for manifest writing.

Algorithm Overview:
Generate comprehensive before/after comparison visualizations:
- before/: 11 plots on unfiltered data (all cells from Step 02)
- after/: 11 plots on filtered data (cells passing both filters)
- removal_reasons/: 3 diagnostic plots showing filtering impact

Total: 25 plots per patient

Input Requirements:
- Two AnnData objects: unfiltered (before) and filtered (after)
- Patient ID string
- Output directory Path

Output Format:
- Dict mapping plot categories to nested dicts of file paths
- Returns {category: {plot_name: file_path}} for manifest

Reference:
- Agent 1 spec: Section 3 (Step 07 - Before/After Report)
- Removal reason logic: Section 3, hierarchical priority
"""

import numpy as np
import pandas as pd
from pathlib import Path
from plotting_core import (
    plot_umap_by_category,
    plot_umap_faceted_by_sample,
    plot_umap_faceted_by_cluster,
    plot_feature_umap,
    plot_removal_reasons_summary,
    plot_removal_reasons_by_sample,
    calculate_percentile_cap
)
from mad_standardization import compute_all_mad_metrics


def generate_before_after_plots(adata_unfiltered, adata_filtered, patient_id, output_dir):
    """
    Generate before/after comparison plots for Step 07.

    Parameters:
    -----------
    adata_unfiltered : AnnData
        Full object (all cells from Step 02 = "before" state)
        Must have removal_reason column added (from data_loader.compute_removal_reasons)

    adata_filtered : AnnData
        Filtered object (cells passing both filters = "after" state)

    patient_id : str
        Patient ID for titles and filenames

    output_dir : Path or str
        Base directory (will create before/, after/, removal_reasons/ subdirs)

    Returns:
    --------
    dict
        Nested dictionary:
        {
            'before': {plot_name: file_path, ...},
            'after': {plot_name: file_path, ...},
            'removal_reasons': {plot_name: file_path, ...}
        }

    Directory Structure Created:
    ---------------------------
    output_dir/
    ├── before/
    │   ├── {patient_id}_umap_by_cluster.png
    │   ├── ... (11 plots total)
    ├── after/
    │   ├── {patient_id}_umap_by_cluster.png
    │   ├── ... (11 plots total)
    └── removal_reasons/
        ├── {patient_id}_removal_reasons_umap.png
        ├── {patient_id}_removal_reasons_barplot.png
        └── {patient_id}_removal_cascade_summary.png

    Raises:
    -------
    ValueError
        If required columns missing or UMAP not present

    Notes:
    ------
    - "Before" plots show all cells (including removed)
    - "After" plots should match Step 06 plots (same filtered data)
    - Removal reason plots use hierarchical logic (Step 03 > Step 02 > Passing)

    Example:
    --------
    >>> from pathlib import Path
    >>> plot_paths = generate_before_after_plots(
    >>>     adata_unfiltered=adata_all,
    >>>     adata_filtered=adata_passing,
    >>>     patient_id='Pat1',
    >>>     output_dir=Path('outputs/06_FinalReport')
    >>> )
    >>> print(f"Before: {len(plot_paths['before'])} plots")
    >>> print(f"After: {len(plot_paths['after'])} plots")
    >>> print(f"Removal: {len(plot_paths['removal_reasons'])} plots")
    """
    # Convert output_dir to Path
    output_dir = Path(output_dir)

    # Create subdirectories
    before_dir = output_dir / 'before'
    after_dir = output_dir / 'after'
    removal_dir = output_dir / 'removal_reasons'

    before_dir.mkdir(parents=True, exist_ok=True)
    after_dir.mkdir(parents=True, exist_ok=True)
    removal_dir.mkdir(parents=True, exist_ok=True)

    # Initialize results dictionary
    all_plot_paths = {
        'before': {},
        'after': {},
        'removal_reasons': {}
    }

    # ========================================================================
    # GENERATE "BEFORE" PLOTS (on unfiltered data)
    # ========================================================================
    print(f"{patient_id}: Generating 'before' plots (unfiltered data)...")
    before_paths = _generate_standard_plots(
        adata=adata_unfiltered,
        patient_id=patient_id,
        output_dir=before_dir,
        label_suffix="Before Filtering"
    )
    all_plot_paths['before'] = before_paths

    # ========================================================================
    # GENERATE "AFTER" PLOTS (on filtered data)
    # ========================================================================
    print(f"{patient_id}: Generating 'after' plots (filtered data)...")
    after_paths = _generate_standard_plots(
        adata=adata_filtered,
        patient_id=patient_id,
        output_dir=after_dir,
        label_suffix="After Filtering"
    )
    all_plot_paths['after'] = after_paths

    # ========================================================================
    # GENERATE REMOVAL REASON PLOTS
    # ========================================================================
    print(f"{patient_id}: Generating removal reason diagnostic plots...")
    removal_paths = _generate_removal_plots(
        adata=adata_unfiltered,  # Use unfiltered (has all cells + removal_reason)
        patient_id=patient_id,
        output_dir=removal_dir
    )
    all_plot_paths['removal_reasons'] = removal_paths

    # ========================================================================
    # SUMMARY
    # ========================================================================
    total_plots = sum(len(v) for v in all_plot_paths.values())
    print(f"{patient_id}: Successfully generated {total_plots} plots")
    print(f"  Before: {len(before_paths)} plots")
    print(f"  After: {len(after_paths)} plots")
    print(f"  Removal: {len(removal_paths)} plots")

    return all_plot_paths


# ============================================================================
# INTERNAL HELPER FUNCTIONS (standard plot generation)
# ============================================================================

def _generate_standard_plots(adata, patient_id, output_dir, label_suffix=""):
    """
    Generate standard 11-plot set (same as Step 06).

    Parameters:
    -----------
    adata : AnnData
        Data to plot (either unfiltered or filtered)

    patient_id : str
        Patient ID

    output_dir : Path
        Directory to save plots

    label_suffix : str
        Suffix for titles (e.g., "Before Filtering", "After Filtering")

    Returns:
    --------
    dict
        {plot_name: file_path}

    Notes:
    ------
    - Generates same 11 plots as Step 06
    - Can be called on both before and after datasets
    """
    plot_paths = {}

    # Validate UMAP presence
    if 'X_umap' not in adata.obsm:
        raise ValueError(f"{patient_id}: UMAP coordinates missing")

    # Compute MAD metrics
    print(f"  Computing MAD-standardized metrics...")
    compute_all_mad_metrics(adata)

    # ------------------------------------------------------------------------
    # Plot 1: UMAP by Cluster
    # ------------------------------------------------------------------------
    plot_name = f"{patient_id}_umap_by_cluster.png"
    plot_path = output_dir / plot_name

    plot_umap_by_category(
        adata=adata,
        color_by='cluster_coarse',
        title=f"{patient_id} - Cells by Cluster ({label_suffix})",
        save_path=plot_path,
        palette='tab10',
        legend_loc='on data',
        figsize=(8, 8)
    )
    plot_paths['umap_by_cluster'] = str(plot_path)

    # ------------------------------------------------------------------------
    # Plot 2: UMAP by Sample
    # ------------------------------------------------------------------------
    plot_name = f"{patient_id}_umap_by_sample.png"
    plot_path = output_dir / plot_name

    plot_umap_by_category(
        adata=adata,
        color_by='sample_id',
        title=f"{patient_id} - Cells by Sample ({label_suffix})",
        save_path=plot_path,
        palette='tab20',
        legend_loc='right margin',
        figsize=(10, 8)
    )
    plot_paths['umap_by_sample'] = str(plot_path)

    # ------------------------------------------------------------------------
    # Plot 3: UMAP Faceted by Sample
    # ------------------------------------------------------------------------
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

    # ------------------------------------------------------------------------
    # Plot 4: UMAP Faceted by Cluster
    # ------------------------------------------------------------------------
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

    # ------------------------------------------------------------------------
    # Plots 5-7: Feature Plots (Raw Metrics)
    # ------------------------------------------------------------------------
    feature_plots_raw = [
        ('n_genes_by_counts', 'Gene Count', 'viridis'),
        ('total_counts', 'UMI Count', 'viridis'),
        ('pct_counts_mt', 'Mitochondrial %', 'Reds')
    ]

    for metric, label, cmap in feature_plots_raw:
        plot_name = f"{patient_id}_feature_{metric}.png"
        plot_path = output_dir / plot_name

        vmax = calculate_percentile_cap(adata, metric, percentile=95)

        plot_feature_umap(
            adata=adata,
            feature_col=metric,
            title=f"{patient_id} - {label} ({label_suffix})",
            save_path=plot_path,
            vmin=0,
            vmax=vmax,
            cmap=cmap,
            figsize=(8, 8)
        )
        plot_paths[f'feature_{metric}'] = str(plot_path)

    # ------------------------------------------------------------------------
    # Plots 8-10: Feature Plots (MAD-Standardized)
    # ------------------------------------------------------------------------
    feature_plots_mad = [
        ('n_genes_by_counts', 'Gene Count'),
        ('total_counts', 'UMI Count'),
        ('pct_counts_mt', 'Mitochondrial %')
    ]

    for metric, label in feature_plots_mad:
        mad_col = f"{metric}_mad"
        plot_name = f"{patient_id}_feature_{metric}_mad.png"
        plot_path = output_dir / plot_name

        plot_feature_umap(
            adata=adata,
            feature_col=mad_col,
            title=f"{patient_id} - {label} (MAD units, {label_suffix})",
            save_path=plot_path,
            vmin=-3,
            vmax=3,
            vcenter=0,
            cmap='RdYlBu_r',
            figsize=(8, 8)
        )
        plot_paths[f'feature_{metric}_mad'] = str(plot_path)

    # ------------------------------------------------------------------------
    # Plot 11: Doublet Score
    # ------------------------------------------------------------------------
    plot_name = f"{patient_id}_feature_scDblFinder_score.png"
    plot_path = output_dir / plot_name

    plot_feature_umap(
        adata=adata,
        feature_col='scDblFinder_score',
        title=f"{patient_id} - Doublet Score ({label_suffix})",
        save_path=plot_path,
        vmin=0,
        vmax=1,
        cmap='plasma',
        figsize=(8, 8)
    )
    plot_paths['feature_scDblFinder_score'] = str(plot_path)

    return plot_paths


def _generate_removal_plots(adata, patient_id, output_dir):
    """
    Generate removal reason diagnostic plots.

    Parameters:
    -----------
    adata : AnnData
        Unfiltered data with 'removal_reason' column

    patient_id : str
        Patient ID

    output_dir : Path
        Directory to save plots

    Returns:
    --------
    dict
        {plot_name: file_path}

    Plots:
    ------
    1. removal_reasons_umap.png - UMAP colored by removal reason
    2. removal_reasons_barplot.png - Stacked bar chart per sample
    3. removal_cascade_summary.png - Simple bar chart showing filtering cascade
    """
    plot_paths = {}

    # Validate removal_reason column exists
    if 'removal_reason' not in adata.obs.columns:
        raise ValueError(
            f"{patient_id}: 'removal_reason' column missing. "
            "Must run data_loader.compute_removal_reasons() first."
        )

    # ------------------------------------------------------------------------
    # Plot 1: UMAP Colored by Removal Reason
    # ------------------------------------------------------------------------
    print(f"  Generating removal reasons UMAP...")
    plot_name = f"{patient_id}_removal_reasons_umap.png"
    plot_path = output_dir / plot_name

    # Define color palette (per Agent 1 spec, verified by Agent 2)
    removal_palette = {
        'Passing': '#2ecc71',
        'Step 02 Removal (QC)': '#f39c12',
        'Step 03 Removal (Doublet)': '#e74c3c'
    }

    plot_umap_by_category(
        adata=adata,
        color_by='removal_reason',
        title=f"{patient_id} - Removal Reasons",
        save_path=plot_path,
        palette=removal_palette,
        legend_loc='right margin',
        alpha=0.5,  # Transparency to show overlapping cells
        figsize=(10, 8)
    )
    plot_paths['removal_reasons_umap'] = str(plot_path)

    # ------------------------------------------------------------------------
    # Plot 2: Summary Bar Plot (Overall Counts)
    # ------------------------------------------------------------------------
    print(f"  Generating removal reasons summary...")
    plot_name = f"{patient_id}_removal_reasons_barplot.png"
    plot_path = output_dir / plot_name

    plot_removal_reasons_summary(
        removal_series=adata.obs['removal_reason'],
        patient_id=patient_id,
        save_path=plot_path,
        figsize=(12, 6)
    )
    plot_paths['removal_reasons_summary'] = str(plot_path)

    # ------------------------------------------------------------------------
    # Plot 3: Per-Sample Breakdown (Stacked Bars)
    # ------------------------------------------------------------------------
    print(f"  Generating per-sample breakdown...")
    plot_name = f"{patient_id}_removal_reasons_by_sample.png"
    plot_path = output_dir / plot_name

    plot_removal_reasons_by_sample(
        adata=adata,
        patient_id=patient_id,
        save_path=plot_path
    )
    plot_paths['removal_reasons_by_sample'] = str(plot_path)

    return plot_paths


# ============================================================================
# HELPER FUNCTIONS (for manifest generation)
# ============================================================================

def compute_removal_breakdown(adata):
    """
    Compute detailed removal breakdown statistics.

    Parameters:
    -----------
    adata : AnnData
        Unfiltered data with removal_reason column

    Returns:
    --------
    dict
        Removal statistics:
        - 'n_cells_total': Total cells
        - 'n_cells_passing': Passing both filters
        - 'n_cells_step_02_removal': Removed by Step 02 only
        - 'n_cells_step_03_removal': Removed by Step 03
        - 'pass_rate': Overall pass rate (%)
        - 'removal_breakdown': {category: count}

    Example:
    --------
    >>> breakdown = compute_removal_breakdown(adata)
    >>> print(f"Pass rate: {breakdown['pass_rate']:.1f}%")
    """
    counts = adata.obs['removal_reason'].value_counts()
    total = len(adata)

    return {
        'n_cells_total': int(total),
        'n_cells_passing': int(counts.get('Passing', 0)),
        'n_cells_step_02_removal': int(counts.get('Step 02 Removal (QC)', 0)),
        'n_cells_step_03_removal': int(counts.get('Step 03 Removal (Doublet)', 0)),
        'pass_rate': float(counts.get('Passing', 0) / total * 100),
        'removal_breakdown': {
            'passing': int(counts.get('Passing', 0)),
            'step_02_qc': int(counts.get('Step 02 Removal (QC)', 0)),
            'step_03_doublet': int(counts.get('Step 03 Removal (Doublet)', 0))
        }
    }


def validate_before_after_consistency(adata_before, adata_after):
    """
    Validate consistency between before and after datasets.

    Parameters:
    -----------
    adata_before : AnnData
        Unfiltered dataset

    adata_after : AnnData
        Filtered dataset

    Returns:
    --------
    dict
        Validation results:
        - 'n_cells_before': Cell count before
        - 'n_cells_after': Cell count after
        - 'n_cells_removed': Cells removed
        - 'removal_rate': Percentage removed
        - 'after_is_subset': Whether after is subset of before

    Raises:
    -------
    ValueError
        If after is not a proper subset of before
    """
    n_before = len(adata_before)
    n_after = len(adata_after)
    n_removed = n_before - n_after

    # Check if after cells are subset of before cells
    before_barcodes = set(adata_before.obs.index)
    after_barcodes = set(adata_after.obs.index)

    is_subset = after_barcodes.issubset(before_barcodes)

    if not is_subset:
        extra_cells = after_barcodes - before_barcodes
        raise ValueError(
            f"After dataset contains {len(extra_cells)} cells not in before dataset! "
            f"After must be a subset of before."
        )

    return {
        'n_cells_before': int(n_before),
        'n_cells_after': int(n_after),
        'n_cells_removed': int(n_removed),
        'removal_rate': float(n_removed / n_before * 100),
        'after_is_subset': is_subset
    }
