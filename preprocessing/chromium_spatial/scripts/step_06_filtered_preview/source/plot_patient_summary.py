"""
Step 05 Filtered Preview - Patient Summary Plots

Generate 5 patient-level summary plots showing before/after QC filtering.
NO file I/O - all functions return matplotlib Figure objects.

Author: Agent 3 (Algorithm Developer)
Date: 2025-11-03
"""

from typing import Dict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from plotting_utils import (
    REMOVAL_COLORS,
    set_plotting_style,
    validate_adata_for_plotting
)


def generate_patient_summary_plots(
    adata_full: sc.AnnData,
    adata_filtered: sc.AnnData,
    patient_id: str,
    plot_params: dict
) -> Dict[str, plt.Figure]:
    """
    Generate all 13 patient summary plots.

    Args:
        adata_full: Full dataset (all cells with removal_status)
        adata_filtered: Filtered dataset (passing cells only, reprocessed UMAP)
        patient_id: Patient identifier
        plot_params: Dict with color palette, sizes, etc.

    Returns:
        Dict mapping plot names to matplotlib Figure objects
        Summary plots (5): 'qc_distributions', 'removal_by_cluster', 'removal_by_sample',
                          'qc_scatter_matrix', 'doublet_scatter'
        UMAPs before filtering (4): 'umap_before_removal_status', 'umap_before_clusters',
                                    'umap_before_doublet_score', 'umap_before_sample_id'
        UMAPs after filtering (4): 'umap_after_clusters', 'umap_after_umi',
                                   'umap_after_features', 'umap_after_sample_id'

    Example:
        >>> params = {
        ...     'removal_colors': REMOVAL_COLORS,
        ...     'plot_dpi': 300,
        ...     'scope1_plot_size': [12, 8]
        ... }
        >>> figures = generate_patient_summary_plots(
        ...     adata_full, adata_filtered, 'Pat1', params
        ... )
        >>> # figures['qc_distributions'] is a matplotlib Figure object
    """

    print(f"  Generating patient summary plots for {patient_id}...")

    # Validate inputs
    validate_adata_for_plotting(
        adata_full,
        ['removal_status', 'cluster_coarse', 'sample_id']
    )
    validate_adata_for_plotting(
        adata_filtered,
        ['cluster_coarse', 'sample_id']
    )

    # Initialize dict to collect figures
    figures = {}

    # Generate each plot type
    print(f"    Plots 1-3: Summary plots...")
    figures['qc_distributions'] = plot_qc_distributions(
        adata_full, patient_id, plot_params
    )
    figures['removal_by_cluster'] = plot_removal_by_cluster(
        adata_full, patient_id, plot_params
    )
    figures['removal_by_sample'] = plot_removal_by_sample(
        adata_full, patient_id, plot_params
    )

    print(f"    Plots 4-5: Scatter plots...")
    figures['qc_scatter_matrix'] = plot_qc_scatter_matrix(
        adata_full, patient_id, plot_params
    )
    figures['doublet_scatter'] = plot_doublet_scatter(
        adata_full, patient_id, plot_params
    )

    print(f"    Plots 6-9: UMAPs (before filtering)...")
    # Before filtering (Step 02 UMAP - all cells)
    figures['umap_before_removal_status'] = plot_umap_single(
        adata_full, 'removal_status', 'Before Filtering: Removal Status', patient_id, plot_params, categorical=True
    )
    figures['umap_before_clusters'] = plot_umap_single(
        adata_full, 'cluster_coarse', 'Before Filtering: Clusters', patient_id, plot_params, categorical=True
    )
    figures['umap_before_doublet_score'] = plot_umap_single(
        adata_full, 'scDblFinder_score', 'Before Filtering: Doublet Score', patient_id, plot_params, categorical=False
    )
    figures['umap_before_sample_id'] = plot_umap_single(
        adata_full, 'sample_id', 'Before Filtering: Sample ID', patient_id, plot_params, categorical=True
    )

    print(f"    Plots 10-13: UMAPs (after filtering)...")
    # After filtering (Step 06 UMAP - passing cells only)
    figures['umap_after_clusters'] = plot_umap_single(
        adata_filtered, 'cluster_coarse', 'After Filtering: Clusters', patient_id, plot_params, categorical=True, use_step05_umap=True
    )
    figures['umap_after_umi'] = plot_umap_single(
        adata_filtered, 'total_counts', 'After Filtering: Total UMI', patient_id, plot_params, categorical=False, use_step05_umap=True
    )
    figures['umap_after_features'] = plot_umap_single(
        adata_filtered, 'n_genes_by_counts', 'After Filtering: Detected Features', patient_id, plot_params, categorical=False, use_step05_umap=True
    )
    figures['umap_after_sample_id'] = plot_umap_single(
        adata_filtered, 'sample_id', 'After Filtering: Sample ID', patient_id, plot_params, categorical=True, use_step05_umap=True
    )

    # Cell type UMAPs (if Step 04 data available)
    if 'HBCATransferredLabels.Kumar_2023' in adata_full.obs.columns:
        print(f"    Plots 14-15: Cell type UMAPs (Step 04 labels)...")
        figures['umap_before_celltype'] = plot_umap_single(
            adata_full, 'HBCATransferredLabels.Kumar_2023', 'Before Filtering: Cell Type', patient_id, plot_params, categorical=True
        )
        if 'HBCATransferredLabels.Kumar_2023' in adata_filtered.obs.columns:
            figures['umap_after_celltype'] = plot_umap_single(
                adata_filtered, 'HBCATransferredLabels.Kumar_2023', 'After Filtering: Cell Type', patient_id, plot_params, categorical=True, use_step05_umap=True
            )

    print(f"  All {len(figures)} patient summary plots generated")

    return figures


def plot_qc_distributions(
    adata: sc.AnnData,
    patient_id: str,
    plot_params: dict
) -> plt.Figure:
    """
    Generate QC metric distributions by removal status.

    Faceted violin plot showing MAD-standardized QC metrics split by removal status.

    Args:
        adata: AnnData with removal_status and MAD metrics
        patient_id: Patient identifier
        plot_params: Plot parameters dict

    Returns:
        matplotlib Figure object

    Example:
        >>> fig = plot_qc_distributions(adata_full, 'Pat1', params)
    """

    # Get plot parameters
    figsize = plot_params.get('scope1_plot_size', [12, 8])
    colors = plot_params.get('removal_colors', REMOVAL_COLORS)

    # Create figure
    fig, axes = plt.subplots(1, 3, figsize=figsize)

    # Metrics to plot (MAD-standardized)
    metrics = [
        ('total_counts_mad', 'UMI (MAD)'),
        ('n_genes_by_counts_mad', 'Genes (MAD)'),
        ('pct_counts_mt_mad', 'Mito % (MAD)')
    ]

    # Create violin plots
    for ax, (metric, label) in zip(axes, metrics):
        # Check if MAD column exists
        if metric not in adata.obs.columns:
            ax.text(0.5, 0.5, f'{metric}\nnot available',
                    ha='center', va='center', transform=ax.transAxes)
            ax.set_title(label)
            continue

        # Create violin plot
        sns.violinplot(
            data=adata.obs,
            x='removal_status',
            y=metric,
            order=['Pass', 'Removed by Step 02', 'Removed by Step 03'],
            palette=colors,
            ax=ax,
            cut=0
        )

        ax.set_title(label, fontsize=12, fontweight='bold')
        ax.set_xlabel('Removal Status', fontsize=10)
        ax.set_ylabel('MAD Score', fontsize=10)
        ax.tick_params(axis='x', rotation=45)

        # Add reference line at 0 (median)
        ax.axhline(0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)

    # Overall title
    fig.suptitle(f'{patient_id}: QC Metric Distributions by Removal Status',
                 fontsize=14, fontweight='bold')

    plt.tight_layout()

    return fig


def plot_removal_by_cluster(
    adata: sc.AnnData,
    patient_id: str,
    plot_params: dict
) -> plt.Figure:
    """
    Generate stacked bar chart of removal fractions by cluster.

    Args:
        adata: AnnData with removal_status and cluster_coarse
        patient_id: Patient identifier
        plot_params: Plot parameters dict

    Returns:
        matplotlib Figure object

    Example:
        >>> fig = plot_removal_by_cluster(adata_full, 'Pat1', params)
    """

    # Get plot parameters
    figsize = plot_params.get('scope1_plot_size', [12, 6])
    colors = plot_params.get('removal_colors', REMOVAL_COLORS)

    # Calculate fractions
    removal_by_cluster = adata.obs.groupby(
        ['cluster_coarse', 'removal_status']
    ).size().unstack(fill_value=0)

    # Normalize to fractions
    removal_by_cluster = removal_by_cluster.div(
        removal_by_cluster.sum(axis=1), axis=0
    )

    # Ensure all categories present
    for cat in ['Pass', 'Removed by Step 02', 'Removed by Step 03']:
        if cat not in removal_by_cluster.columns:
            removal_by_cluster[cat] = 0

    # Reorder columns
    removal_by_cluster = removal_by_cluster[
        ['Pass', 'Removed by Step 02', 'Removed by Step 03']
    ]

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Stacked bar plot
    removal_by_cluster.plot(
        kind='bar',
        stacked=True,
        color=[colors[cat] for cat in removal_by_cluster.columns],
        ax=ax,
        width=0.8
    )

    ax.set_xlabel('Cluster ID', fontsize=10)
    ax.set_ylabel('Fraction of Cells', fontsize=10)
    ax.set_title(f'{patient_id}: Removal Fractions by Cluster',
                 fontsize=12, fontweight='bold')
    ax.set_ylim(0, 1)
    ax.legend(title='Removal Status', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)

    plt.tight_layout()

    return fig


def plot_removal_by_sample(
    adata: sc.AnnData,
    patient_id: str,
    plot_params: dict
) -> plt.Figure:
    """
    Generate stacked bar chart of removal fractions by sample.

    Args:
        adata: AnnData with removal_status and sample_id
        patient_id: Patient identifier
        plot_params: Plot parameters dict

    Returns:
        matplotlib Figure object

    Example:
        >>> fig = plot_removal_by_sample(adata_full, 'Pat1', params)
    """

    # Get plot parameters
    figsize = plot_params.get('scope1_plot_size', [12, 6])
    colors = plot_params.get('removal_colors', REMOVAL_COLORS)

    # Calculate fractions
    removal_by_sample = adata.obs.groupby(
        ['sample_id', 'removal_status']
    ).size().unstack(fill_value=0)

    # Normalize to fractions
    removal_by_sample = removal_by_sample.div(
        removal_by_sample.sum(axis=1), axis=0
    )

    # Ensure all categories present
    for cat in ['Pass', 'Removed by Step 02', 'Removed by Step 03']:
        if cat not in removal_by_sample.columns:
            removal_by_sample[cat] = 0

    # Reorder columns
    removal_by_sample = removal_by_sample[
        ['Pass', 'Removed by Step 02', 'Removed by Step 03']
    ]

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Stacked bar plot
    removal_by_sample.plot(
        kind='bar',
        stacked=True,
        color=[colors[cat] for cat in removal_by_sample.columns],
        ax=ax,
        width=0.8
    )

    ax.set_xlabel('Sample ID', fontsize=10)
    ax.set_ylabel('Fraction of Cells', fontsize=10)
    ax.set_title(f'{patient_id}: Removal Fractions by Sample',
                 fontsize=12, fontweight='bold')
    ax.set_ylim(0, 1)
    ax.legend(title='Removal Status', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

    plt.tight_layout()

    return fig


def plot_umap_comparison(
    adata_full: sc.AnnData,
    adata_filtered: sc.AnnData,
    patient_id: str,
    plot_params: dict
) -> plt.Figure:
    """
    Generate side-by-side UMAP comparison (before/after filtering).

    Left: Step 02 UMAP (all cells, colored by removal status)
    Right: Step 05 UMAP (passing cells only, colored by cluster)

    Args:
        adata_full: Full dataset with Step 02 UMAP and removal_status
        adata_filtered: Filtered dataset with Step 05 UMAP
        patient_id: Patient identifier
        plot_params: Plot parameters dict

    Returns:
        matplotlib Figure object

    Example:
        >>> fig = plot_umap_comparison(adata_full, adata_filtered, 'Pat1', params)
    """

    # Get plot parameters
    figsize = plot_params.get('scope1_plot_size', [16, 8])
    colors = plot_params.get('removal_colors', REMOVAL_COLORS)

    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=figsize)

    # Left plot: Before filtering (Step 02 UMAP, colored by removal status)
    # Check if Step 02 UMAP exists
    if 'X_umap' in adata_full.obsm:
        umap_coords = adata_full.obsm['X_umap']

        # Color by removal status
        status_colors = adata_full.obs['removal_status'].map(colors)

        axes[0].scatter(
            umap_coords[:, 0],
            umap_coords[:, 1],
            c=status_colors,
            s=1,
            alpha=0.5,
            rasterized=True
        )

        axes[0].set_xlabel('UMAP 1', fontsize=10)
        axes[0].set_ylabel('UMAP 2', fontsize=10)
        axes[0].set_title('Before Filtering\n(Step 02 UMAP, all cells)',
                          fontsize=12, fontweight='bold')
        axes[0].axis('off')

        # Add cluster labels as text annotations at centroids
        clusters_unique = adata_full.obs['cluster_coarse'].unique()
        for cluster_id in clusters_unique:
            cluster_mask = adata_full.obs['cluster_coarse'] == cluster_id
            cluster_coords = umap_coords[cluster_mask]
            centroid_x = cluster_coords[:, 0].mean()
            centroid_y = cluster_coords[:, 1].mean()
            axes[0].text(centroid_x, centroid_y, str(int(cluster_id)),
                        fontsize=8, fontweight='bold', ha='center', va='center',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='black', alpha=0.7))

        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor=colors['Pass'], label='Pass'),
            Patch(facecolor=colors['Removed by Step 02'], label='Removed by Step 02'),
            Patch(facecolor=colors['Removed by Step 03'], label='Removed by Step 03')
        ]
        axes[0].legend(handles=legend_elements, loc='upper right', frameon=False)

    else:
        axes[0].text(0.5, 0.5, 'Step 02 UMAP\nnot available',
                     ha='center', va='center', transform=axes[0].transAxes)

    # Right plot: After filtering (Step 05 UMAP, colored by cluster)
    if 'X_umap' in adata_filtered.obsm:
        umap_coords = adata_filtered.obsm['X_umap']

        # Color by cluster
        clusters = adata_filtered.obs['cluster_coarse'].astype(int)
        n_clusters = clusters.nunique()

        # Use tab20 colormap for up to 20 clusters
        from matplotlib import cm
        if n_clusters <= 20:
            cmap = cm.get_cmap('tab20')
        else:
            cmap = cm.get_cmap('hsv')

        scatter = axes[1].scatter(
            umap_coords[:, 0],
            umap_coords[:, 1],
            c=clusters,
            s=1,
            alpha=0.5,
            cmap=cmap,
            rasterized=True
        )

        axes[1].set_xlabel('UMAP 1', fontsize=10)
        axes[1].set_ylabel('UMAP 2', fontsize=10)
        axes[1].set_title('After Filtering\n(Step 05 UMAP, passing cells)',
                          fontsize=12, fontweight='bold')
        axes[1].axis('off')

        # Add cluster labels as text annotations at centroids
        clusters_unique_filtered = adata_filtered.obs['cluster_coarse'].unique()
        for cluster_id in clusters_unique_filtered:
            cluster_mask = adata_filtered.obs['cluster_coarse'] == cluster_id
            cluster_coords = umap_coords[cluster_mask]
            centroid_x = cluster_coords[:, 0].mean()
            centroid_y = cluster_coords[:, 1].mean()
            axes[1].text(centroid_x, centroid_y, str(int(cluster_id)),
                        fontsize=8, fontweight='bold', ha='center', va='center',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='black', alpha=0.7))

        # Add colorbar
        cbar = plt.colorbar(scatter, ax=axes[1], fraction=0.046, pad=0.04)
        cbar.set_label('Cluster', fontsize=10)

    else:
        axes[1].text(0.5, 0.5, 'Step 05 UMAP\nnot available',
                     ha='center', va='center', transform=axes[1].transAxes)

    # Overall title
    fig.suptitle(f'{patient_id}: Before/After UMAP Comparison',
                 fontsize=14, fontweight='bold')

    plt.tight_layout()

    return fig


def plot_qc_scatter_matrix(
    adata: sc.AnnData,
    patient_id: str,
    plot_params: dict
) -> plt.Figure:
    """
    Generate pairwise scatter matrix for QC metrics.

    3x3 matrix showing UMI vs Genes, UMI vs Mito, Genes vs Mito
    (MAD-standardized values), colored by removal status.

    Args:
        adata: AnnData with removal_status and MAD metrics
        patient_id: Patient identifier
        plot_params: Plot parameters dict

    Returns:
        matplotlib Figure object

    Example:
        >>> fig = plot_qc_scatter_matrix(adata_full, 'Pat1', params)
    """

    # Get plot parameters
    figsize = plot_params.get('scope1_plot_size', [12, 12])
    colors = plot_params.get('removal_colors', REMOVAL_COLORS)

    # MAD metrics
    metrics = {
        'total_counts_mad': 'UMI (MAD)',
        'n_genes_by_counts_mad': 'Genes (MAD)',
        'pct_counts_mt_mad': 'Mito % (MAD)'
    }

    # Check all metrics exist
    missing = [m for m in metrics.keys() if m not in adata.obs.columns]
    if missing:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, f'Missing MAD metrics:\n{missing}',
                ha='center', va='center', transform=ax.transAxes)
        return fig

    # Prepare data
    plot_data = adata.obs[list(metrics.keys()) + ['removal_status']].copy()
    plot_data = plot_data.dropna()

    # Map colors
    plot_data['color'] = plot_data['removal_status'].map(colors)

    # Create figure
    fig, axes = plt.subplots(3, 3, figsize=figsize)

    metric_names = list(metrics.keys())

    for i, metric_y in enumerate(metric_names):
        for j, metric_x in enumerate(metric_names):

            ax = axes[i, j]

            if i == j:
                # Diagonal: histogram
                for status in ['Pass', 'Removed by Step 02', 'Removed by Step 03']:
                    status_data = plot_data[plot_data['removal_status'] == status]
                    ax.hist(
                        status_data[metric_x].values,
                        bins=30,
                        alpha=0.5,
                        color=colors[status],
                        label=status
                    )
                ax.set_ylabel('Count', fontsize=8)
                if i == 0:
                    ax.legend(fontsize=6, frameon=False)

            elif i > j:
                # Lower triangle: scatter plot
                ax.scatter(
                    plot_data[metric_x].values,
                    plot_data[metric_y].values,
                    c=plot_data['color'].values,
                    s=1,
                    alpha=0.3,
                    rasterized=True
                )

            else:
                # Upper triangle: leave empty or mirror
                ax.axis('off')

            # Set labels
            if i == len(metric_names) - 1:
                ax.set_xlabel(metrics[metric_x], fontsize=8)
            if j == 0 and i != j:
                ax.set_ylabel(metrics[metric_y], fontsize=8)

            # Add grid
            if i != j:
                ax.grid(True, alpha=0.3)

    # Overall title
    fig.suptitle(f'{patient_id}: QC Metric Relationships (MAD-Scaled)',
                 fontsize=14, fontweight='bold')

    plt.tight_layout()

    return fig


def plot_doublet_scatter(
    adata: sc.AnnData,
    patient_id: str,
    plot_params: dict
) -> plt.Figure:
    """
    Generate scatter plot of Doublet Score vs. Detected Features.

    Shows relationship between doublet prediction and feature count,
    colored by removal status.

    Args:
        adata: AnnData with scDblFinder_score, n_genes_by_counts, removal_status
        patient_id: Patient identifier
        plot_params: Plot parameters dict

    Returns:
        matplotlib Figure object
    """

    # Get plot parameters
    figsize = plot_params.get('plot_size', [10, 6])
    colors = plot_params.get('removal_colors', REMOVAL_COLORS)

    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    # Check required columns exist
    required = ['scDblFinder_score', 'n_genes_by_counts', 'removal_status']
    if not all(col in adata.obs.columns for col in required):
        ax.text(0.5, 0.5, f'Required columns missing:\n{required}',
                ha='center', va='center', transform=ax.transAxes)
        return fig

    # Plot by removal status
    for status in ['Pass', 'Removed by Step 02', 'Removed by Step 03']:
        mask = adata.obs['removal_status'] == status
        if mask.sum() > 0:
            ax.scatter(
                adata.obs.loc[mask, 'n_genes_by_counts'],
                adata.obs.loc[mask, 'scDblFinder_score'],
                c=colors[status],
                s=5,
                alpha=0.4,
                label=status,
                rasterized=True
            )

    # Labels and title
    ax.set_xlabel('Detected Features (nGenes)', fontsize=10)
    ax.set_ylabel('Doublet Score (scDblFinder)', fontsize=10)
    ax.set_title(f'{patient_id}: Doublet Score vs. Feature Count',
                 fontsize=12, fontweight='bold')

    # Legend
    ax.legend(loc='upper right', frameon=False, fontsize=9)

    # Grid
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    return fig


def plot_umap_single(
    adata: sc.AnnData,
    color_by: str,
    title: str,
    patient_id: str,
    plot_params: dict,
    categorical: bool = False,
    use_step05_umap: bool = False
) -> plt.Figure:
    """
    Generate single UMAP colored by specified variable.

    Args:
        adata: AnnData with UMAP coordinates
        color_by: Column name in .obs to color points by
        title: Plot title
        patient_id: Patient identifier
        plot_params: Plot parameters dict
        categorical: If True, use categorical coloring; if False, use continuous
        use_step05_umap: If True, use filtered UMAP; if False, use Step 02 UMAP

    Returns:
        matplotlib Figure object
    """

    # Get plot parameters
    figsize = plot_params.get('plot_size', [8, 6])

    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    # Check UMAP exists
    if 'X_umap' not in adata.obsm:
        ax.text(0.5, 0.5, 'UMAP coordinates not available',
                ha='center', va='center', transform=ax.transAxes)
        ax.set_title(f'{patient_id}: {title}', fontsize=12, fontweight='bold')
        return fig

    umap_coords = adata.obsm['X_umap']

    # Check color column exists
    if color_by not in adata.obs.columns:
        ax.text(0.5, 0.5, f'Column not available:\n{color_by}',
                ha='center', va='center', transform=ax.transAxes)
        ax.set_title(f'{patient_id}: {title}', fontsize=12, fontweight='bold')
        return fig

    if categorical:
        # Categorical coloring
        if color_by == 'removal_status':
            # Use predefined removal colors
            colors = plot_params.get('removal_colors', REMOVAL_COLORS)
            unique_vals = adata.obs[color_by].unique()
            for val in unique_vals:
                mask = adata.obs[color_by] == val
                ax.scatter(
                    umap_coords[mask, 0],
                    umap_coords[mask, 1],
                    c=colors.get(val, '#888888'),
                    s=1,
                    alpha=0.5,
                    label=val,
                    rasterized=True
                )
            ax.legend(loc='upper right', frameon=False, fontsize=8, markerscale=3)

        elif color_by == 'cluster_coarse':
            # Cluster coloring with labels
            clusters = adata.obs[color_by].astype(int)
            n_clusters = clusters.nunique()

            # Choose colormap
            from matplotlib import cm
            if n_clusters <= 20:
                cmap = cm.get_cmap('tab20')
            else:
                cmap = cm.get_cmap('hsv')

            scatter = ax.scatter(
                umap_coords[:, 0],
                umap_coords[:, 1],
                c=clusters,
                s=1,
                alpha=0.5,
                cmap=cmap,
                rasterized=True
            )

            # Add cluster centroids as text labels
            for cluster_id in clusters.unique():
                cluster_mask = clusters == cluster_id
                cluster_coords = umap_coords[cluster_mask]
                centroid_x = cluster_coords[:, 0].mean()
                centroid_y = cluster_coords[:, 1].mean()
                ax.text(centroid_x, centroid_y, str(int(cluster_id)),
                        fontsize=8, fontweight='bold', ha='center', va='center',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='black', alpha=0.7))

            # Colorbar
            cbar = plt.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_label('Cluster', fontsize=10)

        else:
            # Generic categorical (e.g., sample_id)
            unique_vals = adata.obs[color_by].unique()
            n_categories = len(unique_vals)

            # Choose colormap
            from matplotlib import cm
            if n_categories <= 20:
                cmap = cm.get_cmap('tab20')
            else:
                cmap = cm.get_cmap('hsv')

            # Map categories to colors
            color_mapping = {val: cmap(i / n_categories) for i, val in enumerate(unique_vals)}

            for val in unique_vals:
                mask = adata.obs[color_by] == val
                ax.scatter(
                    umap_coords[mask, 0],
                    umap_coords[mask, 1],
                    c=[color_mapping[val]],
                    s=1,
                    alpha=0.5,
                    label=val,
                    rasterized=True
                )

            # Always show legend for sample_id, otherwise show if ≤12 categories
            if color_by == 'sample_id':
                # For sample_id, always show legend (critical for batch assessment)
                # Use multiple columns if many samples
                ncols = 2 if n_categories > 8 else 1
                ax.legend(loc='upper right', frameon=False, fontsize=6, markerscale=2, ncol=ncols)
            elif n_categories <= 12:
                ax.legend(loc='upper right', frameon=False, fontsize=7, markerscale=3, ncol=2)

    else:
        # Continuous coloring
        scatter = ax.scatter(
            umap_coords[:, 0],
            umap_coords[:, 1],
            c=adata.obs[color_by],
            s=1,
            alpha=0.5,
            cmap='viridis',
            rasterized=True
        )

        # Colorbar
        cbar = plt.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(color_by.replace('_', ' ').title(), fontsize=10)

    # Labels and title
    ax.set_xlabel('UMAP 1', fontsize=10)
    ax.set_ylabel('UMAP 2', fontsize=10)
    ax.set_title(f'{patient_id}: {title}', fontsize=12, fontweight='bold')
    ax.axis('off')

    plt.tight_layout()

    return fig
