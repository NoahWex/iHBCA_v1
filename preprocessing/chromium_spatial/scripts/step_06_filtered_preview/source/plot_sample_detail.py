"""
Step 05 Filtered Preview - Sample Detail Plots

Generate per-sample faceted cluster plots with ridge margins.
NO file I/O - all functions return matplotlib Figure objects.

Author: Agent 3 (Algorithm Developer)
Date: 2025-11-03
"""

from typing import List
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import scanpy as sc
from plotting_utils import (
    REMOVAL_COLORS,
    calculate_facet_grid_layout,
    filter_small_clusters,
    validate_adata_for_plotting
)


def generate_sample_detail_plot(
    adata_full: sc.AnnData,
    sample_id: str,
    variant: str,
    plot_params: dict
) -> plt.Figure:
    """
    Generate faceted scatter+ridge plot for one sample.

    Creates a grid of panels, one per cluster, with scatter plot and
    ridge (marginal distribution) plots.

    Args:
        adata_full: Full dataset with all cells and removal_status
        sample_id: Sample identifier
        variant: 'raw' or 'mad' (which metrics to plot)
        plot_params: Dict with facet layout, colors, etc.

    Returns:
        matplotlib Figure with faceted cluster panels

    Example:
        >>> params = {
        ...     'removal_colors': REMOVAL_COLORS,
        ...     'facet_ncols': 4,
        ...     'facet_min_cells': 10
        ... }
        >>> fig = generate_sample_detail_plot(
        ...     adata_full, 'Pat1_P1', 'raw', params
        ... )
    """

    print(f"    Generating {variant} faceted plot for {sample_id}...")

    # Validate inputs
    validate_adata_for_plotting(
        adata_full,
        ['removal_status', 'cluster_coarse', 'sample_id']
    )

    # Get parameters
    facet_ncols = plot_params.get('facet_ncols', 4)
    min_cells = plot_params.get('facet_min_cells', 10)
    colors = plot_params.get('removal_colors', REMOVAL_COLORS)

    # Subset to this sample
    adata_sample = adata_full[adata_full.obs['sample_id'] == sample_id].copy()
    print(f"      Sample subset: {adata_sample.n_obs} cells")

    # Get valid clusters (>= min_cells)
    valid_clusters = filter_small_clusters(adata_sample, min_cells=min_cells)

    if len(valid_clusters) == 0:
        # No valid clusters - create empty figure
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, f'{sample_id}\nNo clusters with >= {min_cells} cells',
                ha='center', va='center', transform=ax.transAxes)
        return fig

    print(f"      Valid clusters: {len(valid_clusters)}")

    # Determine metrics to plot
    if variant == 'raw':
        x_metric = 'total_counts'
        y_metric = 'n_genes_by_counts'
        x_label = 'Total UMI'
        y_label = 'Detected Features'
        thresholds = None
    elif variant == 'mad':
        x_metric = 'total_counts_mad'
        y_metric = 'n_genes_by_counts_mad'
        x_label = 'UMI (MAD)'
        y_label = 'Features (MAD)'
        # MAD thresholds from specification
        thresholds = {
            'x_lower': -3,  # UMI lower threshold
            'y_lower': -3,  # Genes lower threshold
            'y_upper': 3    # Genes upper threshold
        }
    elif variant == 'doublet':
        x_metric = 'total_counts'
        y_metric = 'scDblFinder_score'
        x_label = 'Total UMI'
        y_label = 'Doublet Score'
        thresholds = None
    else:
        raise ValueError(f"Invalid variant: {variant} (must be 'raw', 'mad', or 'doublet')")

    # Check metrics exist
    if x_metric not in adata_sample.obs.columns or y_metric not in adata_sample.obs.columns:
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, f'{sample_id}\nMetrics not available:\n{x_metric}, {y_metric}',
                ha='center', va='center', transform=ax.transAxes)
        return fig

    # Calculate grid layout
    nrows, ncols = calculate_facet_grid_layout(len(valid_clusters), max_cols=facet_ncols)

    # Create figure with GridSpec
    cell_size = plot_params.get('facet_cell_size', 4)
    figsize = (ncols * cell_size, nrows * cell_size)
    fig = plt.figure(figsize=figsize)

    # Create main facet grid
    gs_main = GridSpec(nrows, ncols, figure=fig, hspace=0.4, wspace=0.4)

    # Iterate through valid clusters
    for idx, cluster_id in enumerate(valid_clusters):
        row = idx // ncols
        col = idx % ncols

        # Type-safe comparison: convert both to int for comparison
        cluster_mask = adata_sample.obs['cluster_coarse'].astype(int) == int(cluster_id)
        cluster_data = adata_sample.obs[cluster_mask].copy()

        # Create main scatter plot
        ax_main = fig.add_subplot(gs_main[row, col])

        # Plot by removal status
        for status in ['Pass', 'Removed by Step 02', 'Removed by Step 03']:
            status_mask = cluster_data['removal_status'] == status
            status_data = cluster_data[status_mask]

            if len(status_data) > 0:
                ax_main.scatter(
                    status_data[x_metric].values,
                    status_data[y_metric].values,
                    c=colors[status],
                    s=10,
                    alpha=0.5,
                    label=status,
                    rasterized=True
                )

        # Labels and title
        ax_main.set_xlabel(x_label, fontsize=9)
        ax_main.set_ylabel(y_label, fontsize=9)
        ax_main.set_title(f'Cluster {cluster_id} (n={len(cluster_data)})',
                          fontsize=10, fontweight='bold')

        # Add legend on first panel only
        if idx == 0:
            ax_main.legend(fontsize=7, frameon=False, loc='upper right')

        # Add threshold lines for MAD variant
        if thresholds is not None:
            if 'x_lower' in thresholds:
                ax_main.axvline(thresholds['x_lower'], color='red',
                                linestyle='--', linewidth=0.5, alpha=0.5)
            if 'y_lower' in thresholds:
                ax_main.axhline(thresholds['y_lower'], color='red',
                                linestyle='--', linewidth=0.5, alpha=0.5)
            if 'y_upper' in thresholds:
                ax_main.axhline(thresholds['y_upper'], color='red',
                                linestyle='--', linewidth=0.5, alpha=0.5)

        # Grid
        ax_main.grid(True, alpha=0.3)

    # Overall title
    fig.suptitle(f'{sample_id} - {variant.upper()} Metrics (Faceted by Cluster)',
                 fontsize=14, fontweight='bold')

    return fig


def create_scatter_with_ridges(
    data: pd.DataFrame,
    x_col: str,
    y_col: str,
    color_col: str,
    colors: dict,
    title: str = '',
    thresholds: dict = None
) -> plt.Figure:
    """
    Create single scatter plot with ridge (marginal) distributions.

    Helper function for creating one cluster panel with integrated ridges.

    Args:
        data: DataFrame with x, y, and color columns
        x_col: Column name for x-axis
        y_col: Column name for y-axis
        color_col: Column name for color grouping (e.g., removal_status)
        colors: Dict mapping color_col values to colors
        title: Plot title
        thresholds: Optional dict with threshold lines

    Returns:
        matplotlib Figure object

    Example:
        >>> fig = create_scatter_with_ridges(
        ...     cluster_data,
        ...     'total_counts_mad',
        ...     'n_genes_by_counts_mad',
        ...     'removal_status',
        ...     REMOVAL_COLORS,
        ...     title='Cluster 0'
        ... )
    """

    # Create figure with GridSpec for main + margins
    fig = plt.figure(figsize=(8, 8))

    # Define grid: top ridge, main scatter, right ridge
    gs = GridSpec(3, 3, figure=fig,
                  height_ratios=[1, 5, 0],
                  width_ratios=[5, 1, 0],
                  hspace=0.02, wspace=0.02)

    # Main scatter plot
    ax_main = fig.add_subplot(gs[1, 0])

    # Plot by color group
    for group in data[color_col].unique():
        group_data = data[data[color_col] == group]
        ax_main.scatter(
            group_data[x_col].values,
            group_data[y_col].values,
            c=colors.get(group, 'gray'),
            s=10,
            alpha=0.5,
            label=group,
            rasterized=True
        )

    ax_main.set_xlabel(x_col, fontsize=10)
    ax_main.set_ylabel(y_col, fontsize=10)
    ax_main.set_title(title, fontsize=12, fontweight='bold')
    ax_main.legend(fontsize=8, frameon=False)
    ax_main.grid(True, alpha=0.3)

    # Add thresholds
    if thresholds:
        for key, value in thresholds.items():
            if key.startswith('x_'):
                ax_main.axvline(value, color='red', linestyle='--',
                                linewidth=0.5, alpha=0.5)
            elif key.startswith('y_'):
                ax_main.axhline(value, color='red', linestyle='--',
                                linewidth=0.5, alpha=0.5)

    # Top ridge (x distribution)
    ax_top = fig.add_subplot(gs[0, 0], sharex=ax_main)
    for group in data[color_col].unique():
        group_data = data[data[color_col] == group]
        ax_top.hist(
            group_data[x_col].values,
            bins=30,
            alpha=0.5,
            color=colors.get(group, 'gray'),
            density=True
        )
    ax_top.set_ylabel('Density', fontsize=8)
    ax_top.tick_params(labelbottom=False)
    ax_top.spines['top'].set_visible(False)
    ax_top.spines['right'].set_visible(False)

    # Right ridge (y distribution)
    ax_right = fig.add_subplot(gs[1, 1], sharey=ax_main)
    for group in data[color_col].unique():
        group_data = data[data[color_col] == group]
        ax_right.hist(
            group_data[y_col].values,
            bins=30,
            alpha=0.5,
            color=colors.get(group, 'gray'),
            orientation='horizontal',
            density=True
        )
    ax_right.set_xlabel('Density', fontsize=8)
    ax_right.tick_params(labelleft=False)
    ax_right.spines['top'].set_visible(False)
    ax_right.spines['right'].set_visible(False)

    plt.tight_layout()

    return fig


def determine_facet_grid(n_clusters: int, max_cols: int = 4) -> tuple:
    """
    Calculate optimal grid layout for cluster faceting.

    Args:
        n_clusters: Number of clusters to plot
        max_cols: Maximum columns (default 4)

    Returns:
        Tuple of (nrows, ncols)

    Example:
        >>> nrows, ncols = determine_facet_grid(10, max_cols=4)
        >>> # Returns (3, 4) for 10 clusters in 4-column grid
    """

    return calculate_facet_grid_layout(n_clusters, max_cols=max_cols)


def validate_sample_detail_inputs(
    adata: sc.AnnData,
    sample_id: str,
    variant: str
) -> bool:
    """
    Validate inputs for sample detail plotting.

    Args:
        adata: AnnData object
        sample_id: Sample identifier
        variant: 'raw' or 'mad'

    Returns:
        True if valid

    Raises:
        ValueError: If validation fails

    Example:
        >>> validate_sample_detail_inputs(adata_full, 'Pat1_P1', 'raw')
    """

    # Check variant
    if variant not in ['raw', 'mad']:
        raise ValueError(f"Invalid variant: {variant} (must be 'raw' or 'mad')")

    # Check required columns
    required_cols = ['removal_status', 'cluster_coarse', 'sample_id']
    validate_adata_for_plotting(adata, required_cols)

    # Check sample exists
    if sample_id not in adata.obs['sample_id'].unique():
        raise ValueError(f"Sample {sample_id} not found in data")

    # Check variant-specific metrics
    if variant == 'raw':
        metrics = ['total_counts', 'n_genes_by_counts']
    else:  # mad
        metrics = ['total_counts_mad', 'n_genes_by_counts_mad']

    missing = [m for m in metrics if m not in adata.obs.columns]
    if missing:
        raise ValueError(f"Missing metrics for {variant} variant: {missing}")

    return True
