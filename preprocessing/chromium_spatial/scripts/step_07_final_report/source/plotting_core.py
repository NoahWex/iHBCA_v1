"""
Plotting Core Module - Pure Visualization Functions (Shared utilities)

This module contains PURE plotting functions for UMAP and feature visualizations.
Functions return matplotlib Figure objects - they do NOT save files.

Following wrapper/algorithm separation:
- Plotting functions accept in-memory data
- Return Figure objects (wrappers handle saving)
- Accept save_path as optional parameter for direct saving

Input Requirements:
- AnnData object with .obsm['X_umap']
- Color/feature columns in .obs

Output Format:
- matplotlib.Figure objects
- Or None if save_path provided (saves directly)

Reference:
- Agent 1 spec: Section 2 (Plot Specifications)
- Step 02 visualizations.py for patterns
"""

import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from pathlib import Path


def plot_umap_by_category(
    adata,
    color_by,
    title,
    save_path=None,
    legend_loc='right margin',
    palette=None,
    size=5,
    alpha=1.0,
    figsize=(8, 8)
):
    """
    Create UMAP plot colored by a categorical variable.

    Parameters:
    -----------
    adata : AnnData
        Object with .obsm['X_umap'] and color_by column in .obs

    color_by : str
        Column name in .obs to color by (categorical)

    title : str
        Plot title

    save_path : str or Path, optional
        Path to save PNG (if None, returns fig)

    legend_loc : str, optional
        Legend location (scanpy convention: 'right margin', 'on data', None)
        Default: 'right margin'

    palette : str or list, optional
        Colormap/palette name (e.g., 'tab10', 'tab20')
        If None, uses scanpy default

    size : float, optional
        Point size (default: 5)

    alpha : float, optional
        Point transparency (default: 1.0)

    figsize : tuple, optional
        Figure size in inches (default: (8, 8))

    Returns:
    --------
    matplotlib.Figure or None
        Figure object if save_path is None, else saves and returns None

    Notes:
    ------
    - Uses scanpy.pl.umap for consistency with Step 02
    - Handles missing UMAP gracefully (raises error)

    Example:
    --------
    >>> fig = plot_umap_by_category(adata, 'cluster_coarse', 'Clusters', palette='tab10')
    >>> fig.savefig('output.png', dpi=150, bbox_inches='tight')

    >>> # Or save directly
    >>> plot_umap_by_category(adata, 'cluster_coarse', 'Clusters', save_path='output.png')
    """
    # Validate UMAP presence
    if 'X_umap' not in adata.obsm:
        raise ValueError(
            f"UMAP coordinates missing: adata.obsm['X_umap'] not found. "
            f"Cannot generate plot '{title}'"
        )

    # Validate color column
    if color_by not in adata.obs.columns:
        raise ValueError(f"Color column '{color_by}' not found in adata.obs")

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Plot UMAP using scanpy
    sc.pl.umap(
        adata,
        color=color_by,
        palette=palette,
        legend_loc=legend_loc,
        title=title,
        size=size,
        alpha=alpha,
        ax=ax,
        show=False
    )

    # Handle saving or return
    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        return None
    else:
        return fig


def plot_umap_faceted_by_sample(
    adata,
    color_by,
    title_prefix,
    save_dir=None,
    palette=None,
    ncols=None,
    size=5,
    sample_col='sample_id'
):
    """
    Create faceted UMAP plots split by sample.

    Parameters:
    -----------
    adata : AnnData
        Object with .obsm['X_umap'] and sample_col column

    color_by : str
        Column to color by (e.g., 'cluster_coarse')

    title_prefix : str
        Prefix for overall title

    save_dir : str or Path, optional
        Directory to save individual PNGs (one per sample)
        If None, returns dict of figures

    palette : str or list, optional
        Colormap/palette

    ncols : int, optional
        Number of columns (default: sqrt(n_samples))

    size : float, optional
        Point size (default: 5)

    sample_col : str, optional
        Column for sample ID (default: 'sample_id')

    Returns:
    --------
    dict or None
        {sample_id: figure} if save_dir is None, else saves and returns None

    Notes:
    ------
    - Creates one subplot per sample
    - Calculates ncols from sqrt(n_samples) if not provided
    - Shares colormap across facets

    Example:
    --------
    >>> figs = plot_umap_faceted_by_sample(adata, 'cluster_coarse', 'Pat1 Samples')
    >>> for sample_id, fig in figs.items():
    >>>     fig.savefig(f'{sample_id}_clusters.png')

    >>> # Or save directly
    >>> plot_umap_faceted_by_sample(adata, 'cluster_coarse', 'Pat1', save_dir='plots/')
    """
    # Validate UMAP presence
    if 'X_umap' not in adata.obsm:
        raise ValueError("UMAP coordinates missing")

    if color_by not in adata.obs.columns:
        raise ValueError(f"Color column '{color_by}' not found")

    if sample_col not in adata.obs.columns:
        raise ValueError(f"Sample column '{sample_col}' not found")

    # Get unique samples (sorted for consistent ordering)
    samples = sorted(adata.obs[sample_col].unique())
    n_samples = len(samples)

    # Calculate grid dimensions
    if ncols is None:
        ncols = int(np.ceil(np.sqrt(n_samples)))
    nrows = int(np.ceil(n_samples / ncols))

    # Create figure with subplots
    fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 4*nrows))
    axes = np.atleast_1d(axes).flatten()  # Ensure 1D array

    # Plot each sample
    for idx, sample_id in enumerate(samples):
        ax = axes[idx]

        # Subset to this sample
        adata_sample = adata[adata.obs[sample_col] == sample_id]

        # Plot UMAP
        sc.pl.umap(
            adata_sample,
            color=color_by,
            palette=palette,
            title=sample_id,
            size=size,
            ax=ax,
            show=False,
            legend_loc=None  # No legend per facet (saves space)
        )

    # Hide unused subplots
    for idx in range(n_samples, len(axes)):
        axes[idx].axis('off')

    # Overall title
    fig.suptitle(f'{title_prefix} - Faceted by Sample', fontsize=16, y=0.995)

    plt.tight_layout()

    # Handle saving or return
    if save_dir is not None:
        save_dir = Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)
        save_path = save_dir / f"{title_prefix}_umap_facet_sample_by_{color_by}.png"
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        return None
    else:
        return {'combined': fig}


def plot_umap_faceted_by_cluster(
    adata,
    color_by,
    title_prefix,
    save_dir=None,
    palette=None,
    ncols=None,
    size=5,
    cluster_col='cluster_coarse'
):
    """
    Create faceted UMAP plots split by cluster.

    Parameters:
    -----------
    adata : AnnData
        Object with .obsm['X_umap'] and cluster column

    color_by : str
        Column to color by (e.g., 'sample_id')

    title_prefix : str
        Prefix for plot titles

    save_dir : str or Path, optional
        Directory to save PNG
        If None, returns figure

    palette : str or list, optional
        Colormap/palette

    ncols : int, optional
        Number of columns (default: sqrt(n_clusters))

    size : float, optional
        Point size (default: 5)

    cluster_col : str, optional
        Column for clusters (default: 'cluster_coarse')

    Returns:
    --------
    dict or None
        {'combined': figure} if save_dir is None, else saves and returns None

    Notes:
    ------
    - Creates one subplot per cluster
    - DOES NOT assume contiguous cluster IDs (uses .unique())

    Example:
    --------
    >>> fig_dict = plot_umap_faceted_by_cluster(adata, 'sample_id', 'Pat1')
    >>> fig_dict['combined'].savefig('clusters_faceted.png')
    """
    # Validate UMAP presence
    if 'X_umap' not in adata.obsm:
        raise ValueError("UMAP coordinates missing")

    if color_by not in adata.obs.columns:
        raise ValueError(f"Color column '{color_by}' not found")

    if cluster_col not in adata.obs.columns:
        raise ValueError(f"Cluster column '{cluster_col}' not found")

    # Get unique clusters (sorted, don't assume contiguous)
    clusters = sorted(adata.obs[cluster_col].unique())
    n_clusters = len(clusters)

    # Calculate grid dimensions
    if ncols is None:
        ncols = int(np.ceil(np.sqrt(n_clusters)))
    nrows = int(np.ceil(n_clusters / ncols))

    # Create figure with subplots
    fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 4*nrows))
    axes = np.atleast_1d(axes).flatten()

    # Plot each cluster
    for idx, cluster_id in enumerate(clusters):
        ax = axes[idx]

        # Subset to this cluster
        adata_cluster = adata[adata.obs[cluster_col] == cluster_id]

        # Plot UMAP
        sc.pl.umap(
            adata_cluster,
            color=color_by,
            palette=palette,
            title=f"Cluster {cluster_id}",
            size=size,
            ax=ax,
            show=False,
            legend_loc=None
        )

    # Hide unused subplots
    for idx in range(n_clusters, len(axes)):
        axes[idx].axis('off')

    # Overall title
    fig.suptitle(f'{title_prefix} - Faceted by Cluster', fontsize=16, y=0.995)

    plt.tight_layout()

    # Handle saving or return
    if save_dir is not None:
        save_dir = Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)
        save_path = save_dir / f"{title_prefix}_umap_facet_cluster_by_{color_by}.png"
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        return None
    else:
        return {'combined': fig}


def plot_feature_umap(
    adata,
    feature_col,
    title,
    save_path=None,
    vmin=None,
    vmax=None,
    vcenter=None,
    cmap='viridis',
    size=5,
    figsize=(8, 8)
):
    """
    Create UMAP feature plot (continuous values).

    Parameters:
    -----------
    adata : AnnData
        Object with .obsm['X_umap']

    feature_col : str
        Column in .obs with continuous values

    title : str
        Plot title

    save_path : str or Path, optional
        Path to save PNG

    vmin, vmax : float, optional
        Color scale limits

    vcenter : float, optional
        Center value for diverging colormaps (e.g., 0 for MAD plots)

    cmap : str, optional
        Colormap name (default: 'viridis')

    size : float, optional
        Point size (default: 5)

    figsize : tuple, optional
        Figure size (default: (8, 8))

    Returns:
    --------
    matplotlib.Figure or None

    Notes:
    ------
    - For MAD plots: use cmap='RdYlBu_r', vcenter=0, vmin=-3, vmax=3
    - For QC metrics: use cmap='viridis' or 'Reds', cap vmax at 95th percentile

    Example:
    --------
    >>> # Raw metric plot
    >>> vmax_95 = np.percentile(adata.obs['total_counts'], 95)
    >>> fig = plot_feature_umap(adata, 'total_counts', 'UMI Count', vmax=vmax_95)

    >>> # MAD standardized plot
    >>> fig = plot_feature_umap(adata, 'total_counts_mad', 'UMI (MAD units)',
    >>>                         cmap='RdYlBu_r', vcenter=0, vmin=-3, vmax=3)
    """
    # Validate UMAP presence
    if 'X_umap' not in adata.obsm:
        raise ValueError("UMAP coordinates missing")

    if feature_col not in adata.obs.columns:
        raise ValueError(f"Feature column '{feature_col}' not found")

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Plot UMAP with feature coloring
    sc.pl.umap(
        adata,
        color=feature_col,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        vcenter=vcenter,
        title=title,
        size=size,
        ax=ax,
        show=False,
        colorbar_loc='right'
    )

    # Handle saving or return
    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        return None
    else:
        return fig


def plot_removal_reasons_summary(
    removal_series,
    patient_id,
    save_path=None,
    figsize=(12, 6)
):
    """
    Create summary bar plot of removal reasons.

    Parameters:
    -----------
    removal_series : pd.Series
        Categorical series with removal reasons:
        - 'Passing'
        - 'Step 02 Removal (QC)'
        - 'Step 03 Removal (Doublet)'

    patient_id : str
        Patient ID for title

    save_path : str or Path, optional
        Path to save PNG

    figsize : tuple, optional
        Figure size (default: (12, 6))

    Returns:
    --------
    matplotlib.Figure or None

    Notes:
    ------
    - Shows counts and percentages
    - Uses Agent 1's specified color scheme (green/orange/red)
    - Per Agent 2 Section 7 (verified palette)

    Example:
    --------
    >>> fig = plot_removal_reasons_summary(removal_reasons, 'Pat1')
    >>> fig.savefig('removal_summary.png')
    """
    # Define color palette (per Agent 1 spec, verified by Agent 2)
    colors = {
        'Passing': '#2ecc71',
        'Step 02 Removal (QC)': '#f39c12',
        'Step 03 Removal (Doublet)': '#e74c3c'
    }

    # Count each category
    counts = removal_series.value_counts()
    total = len(removal_series)

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Create bar plot
    categories = ['Passing', 'Step 02 Removal (QC)', 'Step 03 Removal (Doublet)']
    bar_colors = [colors[cat] for cat in categories]
    bar_counts = [counts.get(cat, 0) for cat in categories]

    bars = ax.bar(categories, bar_counts, color=bar_colors, alpha=0.8, edgecolor='black')

    # Add count and percentage labels on bars
    for bar, count in zip(bars, bar_counts):
        height = bar.get_height()
        pct = count / total * 100
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            height,
            f'{count:,}\n({pct:.1f}%)',
            ha='center',
            va='bottom',
            fontsize=11,
            fontweight='bold'
        )

    # Formatting
    ax.set_ylabel('Cell Count', fontsize=12)
    ax.set_title(f'{patient_id} - Filtering Summary', fontsize=14, fontweight='bold')
    ax.set_ylim(0, max(bar_counts) * 1.15)  # Space for labels
    ax.grid(axis='y', alpha=0.3, linestyle='--')

    # Rotate x-axis labels if needed
    plt.xticks(rotation=15, ha='right')

    plt.tight_layout()

    # Handle saving or return
    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        return None
    else:
        return fig


def plot_removal_reasons_by_sample(
    adata,
    patient_id,
    save_path=None,
    figsize=None,
    sample_col='sample_id',
    removal_col='removal_reason'
):
    """
    Create stacked bar chart of removal reasons per sample.

    Parameters:
    -----------
    adata : AnnData
        Object with removal_reason column in .obs

    patient_id : str
        Patient ID for title

    save_path : str or Path, optional
        Path to save PNG

    figsize : tuple, optional
        Figure size (auto-calculated if None)

    sample_col : str, optional
        Sample column (default: 'sample_id')

    removal_col : str, optional
        Removal reason column (default: 'removal_reason')

    Returns:
    --------
    matplotlib.Figure or None

    Notes:
    ------
    - Stacked bars showing removal breakdown per sample
    - Uses same color scheme as summary plot

    Example:
    --------
    >>> fig = plot_removal_reasons_by_sample(adata, 'Pat1')
    """
    # Define colors
    colors = {
        'Passing': '#2ecc71',
        'Step 02 Removal (QC)': '#f39c12',
        'Step 03 Removal (Doublet)': '#e74c3c'
    }

    # Create crosstab (samples × removal reasons)
    crosstab = pd.crosstab(adata.obs[sample_col], adata.obs[removal_col])

    # Ensure all categories present (even if zero)
    for cat in ['Passing', 'Step 02 Removal (QC)', 'Step 03 Removal (Doublet)']:
        if cat not in crosstab.columns:
            crosstab[cat] = 0

    # Reorder columns
    crosstab = crosstab[['Passing', 'Step 02 Removal (QC)', 'Step 03 Removal (Doublet)']]

    # Auto-calculate figsize if not provided
    if figsize is None:
        n_samples = len(crosstab)
        width = max(12, n_samples * 0.8)
        figsize = (width, 6)

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Create stacked bar plot
    crosstab.plot(
        kind='bar',
        stacked=True,
        ax=ax,
        color=[colors[cat] for cat in crosstab.columns],
        edgecolor='black',
        alpha=0.8
    )

    # Formatting
    ax.set_xlabel('Sample ID', fontsize=12)
    ax.set_ylabel('Cell Count', fontsize=12)
    ax.set_title(f'{patient_id} - Filtering Breakdown by Sample', fontsize=14, fontweight='bold')
    ax.legend(title='Status', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(axis='y', alpha=0.3, linestyle='--')

    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    # Handle saving or return
    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        return None
    else:
        return fig


# ============================================================================
# HELPER FUNCTIONS (for plot configuration)
# ============================================================================

def get_color_palette(palette_name, n_colors=None):
    """
    Get color palette for categorical plotting.

    Parameters:
    -----------
    palette_name : str
        Palette name ('tab10', 'tab20', etc.)

    n_colors : int, optional
        Number of colors needed

    Returns:
    --------
    list of str
        Hex color codes

    Example:
    --------
    >>> colors = get_color_palette('tab10', n_colors=5)
    """
    import matplotlib.cm as cm

    if palette_name.startswith('tab'):
        cmap = cm.get_cmap(palette_name)
        max_colors = int(palette_name[3:])  # e.g., 'tab10' → 10
        if n_colors is None:
            n_colors = max_colors
        colors = [cmap(i) for i in range(min(n_colors, max_colors))]
        return colors
    else:
        # For other colormaps
        cmap = cm.get_cmap(palette_name)
        if n_colors is None:
            n_colors = 10
        colors = [cmap(i / n_colors) for i in range(n_colors)]
        return colors


def calculate_percentile_cap(adata, metric_col, percentile=95):
    """
    Calculate percentile cap for feature plots.

    Parameters:
    -----------
    adata : AnnData
        Object with metric in .obs

    metric_col : str
        Metric column name

    percentile : float, optional
        Percentile to cap at (default: 95)

    Returns:
    --------
    float
        Value at specified percentile

    Example:
    --------
    >>> vmax = calculate_percentile_cap(adata, 'total_counts', 95)
    >>> fig = plot_feature_umap(adata, 'total_counts', 'UMI', vmax=vmax)
    """
    return np.percentile(adata.obs[metric_col], percentile)
