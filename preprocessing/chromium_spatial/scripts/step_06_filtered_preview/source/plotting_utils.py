"""
Step 05 Filtered Preview - Plotting Utilities

Shared plotting utilities and color schemes.
NO file I/O - all functions work with in-memory objects.

Author: Agent 3 (Algorithm Developer)
Date: 2025-11-03
"""

from typing import List, Dict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import scanpy as sc


# Color palette for removal status (3 categories)
REMOVAL_COLORS = {
    'Pass': '#2ecc71',                    # Green
    'Removed by Step 02': '#f39c12',      # Orange
    'Removed by Step 03': '#e74c3c'       # Red
}

# Continuous colormaps for features
FEATURE_CMAPS = {
    'umi': 'viridis',
    'genes': 'viridis',
    'mito': 'Reds',
    'mad': 'RdYlBu_r'  # Diverging for MAD scores
}


def set_plotting_style():
    """
    Set consistent matplotlib style for all plots.

    Configures default plot parameters for publication-quality figures.

    Example:
        >>> set_plotting_style()
        >>> # Now all subsequent plots use this style
    """

    # Set scanpy plot defaults
    sc.set_figure_params(
        dpi=100,
        dpi_save=300,
        frameon=False,
        vector_friendly=True,
        fontsize=10,
        figsize=(6, 4),
        color_map='viridis'
    )

    # Additional matplotlib settings
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'axes.linewidth': 1.0,
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'legend.fontsize': 9,
        'legend.frameon': False,
        'figure.titlesize': 14
    })


def validate_adata_for_plotting(
    adata: sc.AnnData,
    required_columns: List[str]
) -> bool:
    """
    Validate AnnData has required columns for plotting.

    Args:
        adata: AnnData object to validate
        required_columns: List of column names that must exist in .obs

    Returns:
        True if validation passes

    Raises:
        ValueError: If validation fails

    Example:
        >>> validate_adata_for_plotting(
        ...     adata,
        ...     ['removal_status', 'cluster_coarse', 'sample_id']
        ... )
    """

    missing_cols = [col for col in required_columns if col not in adata.obs.columns]

    if missing_cols:
        raise ValueError(f"Missing required columns for plotting: {missing_cols}")

    return True


def filter_small_clusters(
    adata: sc.AnnData,
    min_cells: int = 10
) -> List[int]:
    """
    Return list of cluster IDs with sufficient cells.

    Args:
        adata: AnnData with cluster_coarse in .obs
        min_cells: Minimum cells required per cluster

    Returns:
        List of cluster IDs (as integers) with >= min_cells

    Example:
        >>> valid_clusters = filter_small_clusters(adata, min_cells=10)
        >>> # Plot only these clusters
        >>> for cluster in valid_clusters:
        ...     # Generate plot for cluster
    """

    if 'cluster_coarse' not in adata.obs.columns:
        raise ValueError("cluster_coarse column not found")

    cluster_sizes = adata.obs['cluster_coarse'].value_counts()
    valid_clusters = cluster_sizes[cluster_sizes >= min_cells].index.tolist()

    # Convert to int (may be stored as string or float)
    valid_clusters = [int(c) for c in valid_clusters]

    # Sort for consistent ordering
    valid_clusters.sort()

    print(f"  Filter small clusters (min={min_cells}): {len(valid_clusters)}/{len(cluster_sizes)} clusters retained")

    return valid_clusters


def get_removal_color_palette(as_list: bool = False) -> Dict[str, str]:
    """
    Get color palette for removal status.

    Args:
        as_list: If True, return colors as ordered list (default False)

    Returns:
        Dict mapping removal status to hex color, or list of colors

    Example:
        >>> colors = get_removal_color_palette()
        >>> plt.scatter(x, y, c=[colors[status] for status in statuses])
    """

    if as_list:
        # Return in canonical order
        return [
            REMOVAL_COLORS['Pass'],
            REMOVAL_COLORS['Removed by Step 02'],
            REMOVAL_COLORS['Removed by Step 03']
        ]
    else:
        return REMOVAL_COLORS.copy()


def calculate_facet_grid_layout(
    n_items: int,
    max_cols: int = 4
) -> tuple:
    """
    Calculate optimal grid layout for faceted plots.

    Args:
        n_items: Number of items to plot (e.g., clusters)
        max_cols: Maximum columns in grid (default 4)

    Returns:
        Tuple of (nrows, ncols)

    Example:
        >>> nrows, ncols = calculate_facet_grid_layout(10, max_cols=4)
        >>> # nrows=3, ncols=4 (for 10 items)
    """

    ncols = min(n_items, max_cols)
    nrows = int(np.ceil(n_items / ncols))

    return nrows, ncols


def cap_values_at_percentile(
    values: np.ndarray,
    percentile: float = 95.0
) -> np.ndarray:
    """
    Cap extreme values at specified percentile.

    Useful for preventing outliers from dominating color scales.

    Args:
        values: Array of values to cap
        percentile: Percentile threshold (default 95.0)

    Returns:
        Capped array (values above threshold set to threshold)

    Example:
        >>> capped = cap_values_at_percentile(umi_counts, percentile=95)
        >>> # Plot with capped values for better color scale
    """

    threshold = np.percentile(values, percentile)
    capped = np.where(values > threshold, threshold, values)

    n_capped = (values > threshold).sum()
    if n_capped > 0:
        print(f"    Capped {n_capped} values ({n_capped/len(values)*100:.1f}%) at {percentile}th percentile ({threshold:.1f})")

    return capped


def add_mad_threshold_lines(
    ax: plt.Axes,
    thresholds: Dict[str, float],
    axis: str = 'both'
):
    """
    Add MAD threshold reference lines to plot.

    Args:
        ax: Matplotlib axes object
        thresholds: Dict with 'lower' and/or 'upper' threshold values
        axis: Which axis to add lines ('x', 'y', or 'both')

    Example:
        >>> fig, ax = plt.subplots()
        >>> ax.scatter(x_mad, y_mad)
        >>> add_mad_threshold_lines(
        ...     ax,
        ...     {'lower': -3, 'upper': 3},
        ...     axis='both'
        ... )
    """

    line_style = {
        'color': 'red',
        'linestyle': '--',
        'linewidth': 1,
        'alpha': 0.5
    }

    if axis in ['x', 'both']:
        if 'lower' in thresholds:
            ax.axvline(thresholds['lower'], **line_style)
        if 'upper' in thresholds:
            ax.axvline(thresholds['upper'], **line_style)

    if axis in ['y', 'both']:
        if 'lower' in thresholds:
            ax.axhline(thresholds['lower'], **line_style)
        if 'upper' in thresholds:
            ax.axhline(thresholds['upper'], **line_style)


def create_removal_legend(
    ax: plt.Axes,
    counts: Dict[str, int] = None,
    location: str = 'upper right'
):
    """
    Create legend for removal status colors.

    Args:
        ax: Matplotlib axes object
        counts: Optional dict with counts per category
        location: Legend location (matplotlib location string)

    Example:
        >>> fig, ax = plt.subplots()
        >>> # ... scatter plot with removal status colors ...
        >>> create_removal_legend(ax, counts={'Pass': 1000, 'Removed by Step 02': 50})
    """

    import matplotlib.patches as mpatches

    # Create legend handles
    handles = []
    for status in ['Pass', 'Removed by Step 02', 'Removed by Step 03']:
        if counts is not None and status in counts:
            label = f"{status} (n={counts[status]})"
        else:
            label = status

        patch = mpatches.Patch(
            color=REMOVAL_COLORS[status],
            label=label
        )
        handles.append(patch)

    ax.legend(
        handles=handles,
        loc=location,
        frameon=False,
        fontsize=9
    )


def get_cluster_colors(
    n_clusters: int,
    cmap_name: str = 'tab20'
) -> List[str]:
    """
    Get distinct colors for clusters.

    Args:
        n_clusters: Number of clusters
        cmap_name: Matplotlib colormap name (default 'tab20')

    Returns:
        List of hex color strings

    Example:
        >>> colors = get_cluster_colors(15)
        >>> # Use colors[cluster_id] for plotting
    """

    if n_clusters <= 20:
        cmap = mpl.colormaps.get_cmap('tab20')
    else:
        # For many clusters, use continuous colormap
        cmap = mpl.colormaps.get_cmap('hsv')

    colors = [mpl.colors.rgb2hex(cmap(i / n_clusters)) for i in range(n_clusters)]

    return colors


def format_axis_labels(
    ax: plt.Axes,
    x_label: str = None,
    y_label: str = None,
    title: str = None
):
    """
    Apply consistent axis label formatting.

    Args:
        ax: Matplotlib axes object
        x_label: X-axis label (optional)
        y_label: Y-axis label (optional)
        title: Plot title (optional)

    Example:
        >>> format_axis_labels(ax, x_label='UMI count', y_label='Gene count')
    """

    if x_label:
        ax.set_xlabel(x_label, fontsize=10)
    if y_label:
        ax.set_ylabel(y_label, fontsize=10)
    if title:
        ax.set_title(title, fontsize=12, fontweight='bold')


def adjust_plot_limits(
    ax: plt.Axes,
    data: np.ndarray,
    axis: str = 'both',
    padding: float = 0.05
):
    """
    Adjust plot limits based on data range with padding.

    Args:
        ax: Matplotlib axes object
        data: Data array to determine limits from
        axis: Which axis to adjust ('x', 'y', or 'both')
        padding: Fraction of range to add as padding (default 0.05 = 5%)

    Example:
        >>> adjust_plot_limits(ax, mad_values, axis='y', padding=0.1)
    """

    data_min = np.nanmin(data)
    data_max = np.nanmax(data)
    data_range = data_max - data_min

    lower = data_min - (padding * data_range)
    upper = data_max + (padding * data_range)

    if axis in ['x', 'both']:
        ax.set_xlim(lower, upper)
    if axis in ['y', 'both']:
        ax.set_ylim(lower, upper)


def save_figure_to_dict(
    fig: plt.Figure,
    name: str,
    figures_dict: Dict[str, plt.Figure]
) -> Dict[str, plt.Figure]:
    """
    Save figure to dictionary (for returning multiple plots).

    Helper function for returning multiple Figure objects from plotting functions.
    Does NOT save to file.

    Args:
        fig: Matplotlib Figure object
        name: Name/key for this figure
        figures_dict: Dictionary to add figure to

    Returns:
        Updated figures_dict

    Example:
        >>> figures = {}
        >>> fig, ax = plt.subplots()
        >>> # ... create plot ...
        >>> figures = save_figure_to_dict(fig, 'plot1', figures)
        >>> return figures  # Return dict of Figure objects (NO file I/O)
    """

    figures_dict[name] = fig
    return figures_dict


def validate_plot_params(params: dict) -> bool:
    """
    Validate plot parameter dictionary.

    Args:
        params: Dictionary of plot parameters

    Returns:
        True if valid

    Raises:
        ValueError: If invalid parameters

    Example:
        >>> params = {'plot_dpi': 300, 'removal_colors': REMOVAL_COLORS}
        >>> validate_plot_params(params)
    """

    # Check DPI is reasonable
    if 'plot_dpi' in params:
        dpi = params['plot_dpi']
        if not isinstance(dpi, int) or dpi < 50 or dpi > 600:
            raise ValueError(f"Invalid plot_dpi: {dpi} (must be 50-600)")

    # Check color palette exists
    if 'removal_colors' in params:
        colors = params['removal_colors']
        required_keys = ['Pass', 'Removed by Step 02', 'Removed by Step 03']
        missing = [k for k in required_keys if k not in colors]
        if missing:
            raise ValueError(f"Missing removal_colors keys: {missing}")

    return True
