#!/usr/bin/env python3
"""
Step 02 Visualization Module

Generates diagnostic plots for cell filtering QC.
Extracted from archived workflows/qc/scripts/module_2_cell_filter.py.

Visualizations include:
- UMAP plots (cluster, sample, QC status, UMI)
- Low-UMI feature plots (MAD-standardized, diverging colormap)
- MAD ECDF plots (per-cluster, faceted)
- Summary statistics (outlier breakdown, distributions)

Author: Noah Wechter
Date: 2025-10-29
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from pathlib import Path
from scipy.stats import median_abs_deviation


def plot_umap_overview(adata, output_dir, patient_id):
    """
    Generate 2x2 UMAP overview panel.

    Panels:
    1. Coarse cluster assignment
    2. Sample ID
    3. Cell filter pass/fail status
    4. n_UMI gradient (vmax at 95th percentile)

    Args:
        adata: Patient-level AnnData object
        output_dir: Path to visualizations directory
        patient_id: Patient identifier
    """
    print(f"{patient_id}: Generating UMAP overview...")

    # Check if UMAP coordinates exist
    if 'X_umap' not in adata.obsm:
        print(f"{patient_id}: WARNING - UMAP coordinates not found, skipping UMAP overview")
        return

    try:
        fig, axes = plt.subplots(2, 2, figsize=(12, 12))

        # Panel 1: Cluster
        sc.pl.umap(adata, color='cluster_coarse', ax=axes[0,0], show=False, title='Coarse Cluster')

        # Panel 2: Sample
        sc.pl.umap(adata, color='sample_id', ax=axes[0,1], show=False, title='Sample')

        # Panel 3: Pass/Fail
        sc.pl.umap(adata, color='cell_filter_pass_final', ax=axes[1,0], show=False,
                   title='Cell Filter Pass')

        # Panel 4: n_UMI gradient (vmax at 95th percentile)
        umi_q95 = np.percentile(adata.obs['n_umi'], 95)
        sc.pl.umap(adata, color='n_umi', ax=axes[1,1], show=False, title='n_UMI',
                   cmap='viridis', vmax=umi_q95)

        plt.tight_layout()
        output_file = Path(output_dir) / f"{patient_id}_umap.pdf"
        fig.savefig(output_file, dpi=100, bbox_inches='tight')
        plt.close(fig)

        print(f"{patient_id}:   ✓ UMAP overview: {output_file}")

    except Exception as e:
        print(f"{patient_id}: ERROR generating UMAP overview: {e}")


def plot_low_umi_features(adata, output_dir, patient_id):
    """
    Generate 2x2 feature panel highlighting low-UMI debris cells.

    Panels:
    1. UMI (MAD-standardized, diverging colormap) - highlights low UMI debris
    2. Mitochondrial % (vmax at 95th percentile)
    3. n_genes (vmax at 95th percentile)
    4. Any MAD outlier (binary pass/fail)

    Args:
        adata: Patient-level AnnData object
        output_dir: Path to visualizations directory
        patient_id: Patient identifier
    """
    print(f"{patient_id}: Generating low-UMI feature plots...")

    # Check if UMAP coordinates exist
    if 'X_umap' not in adata.obsm:
        print(f"{patient_id}: WARNING - UMAP coordinates not found, skipping low-UMI features")
        return

    try:
        fig, axes = plt.subplots(2, 2, figsize=(12, 12))

        # Panel 1: MAD-standardized UMI per cluster (highlights low UMI debris)
        # Calculate MAD-standardized UMI within each cluster
        umi_mad_standardized = np.zeros(len(adata.obs))
        for cluster in adata.obs['cluster_coarse'].unique():
            cluster_mask = adata.obs['cluster_coarse'] == cluster
            cluster_umi = adata.obs.loc[cluster_mask, 'n_umi'].values

            if len(cluster_umi) > 0:
                cluster_median = np.median(cluster_umi)
                cluster_mad = median_abs_deviation(cluster_umi)

                if cluster_mad > 0:
                    # Standardize: (value - median) / MAD
                    umi_mad_standardized[cluster_mask] = (cluster_umi - cluster_median) / cluster_mad
                else:
                    umi_mad_standardized[cluster_mask] = 0

        adata.obs['umi_mad_std'] = umi_mad_standardized

        # Color scale: -6 to +6 MADs (q95 for symmetry)
        mad_q95 = np.percentile(np.abs(umi_mad_standardized), 95)
        mad_limit = min(mad_q95, 6.0)  # Cap at 6 MADs

        sc.pl.umap(adata, color='umi_mad_std', ax=axes[0,0], show=False,
                   title='UMI (MADs from cluster median)\nYellow=LOW, Purple=HIGH',
                   cmap='RdYlBu_r', vcenter=0, vmin=-mad_limit, vmax=mad_limit)

        # Panel 2: Mito % (q95 limit)
        mito_q95 = np.percentile(adata.obs['mito_percent'], 95)
        sc.pl.umap(adata, color='mito_percent', ax=axes[0,1], show=False,
                   title='Mitochondrial %', cmap='Reds', vmax=mito_q95)

        # Panel 3: Genes (q95 limit)
        genes_q95 = np.percentile(adata.obs['n_genes'], 95)
        sc.pl.umap(adata, color='n_genes', ax=axes[1,0], show=False,
                   title='n_genes', cmap='viridis', vmax=genes_q95)

        # Panel 4: Combined outlier flags
        adata.obs['any_outlier'] = (
            adata.obs['mad_outlier_genes_low'] |
            adata.obs['mad_outlier_genes_high'] |
            adata.obs['mad_outlier_umi_low'] |
            adata.obs['mad_outlier_mito_high']
        )
        sc.pl.umap(adata, color='any_outlier', ax=axes[1,1], show=False,
                   title='Any MAD Outlier', palette=['lightgray', 'red'])

        plt.tight_layout()
        output_file = Path(output_dir) / f"{patient_id}_features_lowUMI.pdf"
        fig.savefig(output_file, dpi=100, bbox_inches='tight')
        plt.close(fig)

        # Clean up temporary columns
        adata.obs.drop(columns=['umi_mad_std', 'any_outlier'], inplace=True)

        print(f"{patient_id}:   ✓ Low-UMI features: {output_file}")

    except Exception as e:
        print(f"{patient_id}: ERROR generating low-UMI features: {e}")


def plot_mad_ecdf(adata, output_dir, patient_id):
    """
    Generate per-patient ECDF plots faceted by cluster, MAD-standardized.

    Creates 3 separate plots:
    1. UMI ECDF (lower threshold only, -3 MAD)
    2. Genes ECDF (both thresholds, -3 and +3 MAD)
    3. Mito ECDF (upper threshold only, +3 MAD)

    Each plot is faceted by cluster, showing:
    - Y-axis: MADs from cluster median
    - X-axis: Empirical cumulative probability
    - Threshold lines and shaded rejection regions

    Args:
        adata: Patient-level AnnData object
        output_dir: Path to visualizations directory
        patient_id: Patient identifier
    """
    print(f"{patient_id}: Generating MAD ECDF plots...")

    try:
        clusters = sorted(adata.obs['cluster_coarse'].unique())
        n_clusters = len(clusters)

        # Calculate grid dimensions
        ncols = min(4, n_clusters)
        nrows = int(np.ceil(n_clusters / ncols))

        # === UMI ECDF (lower threshold only) ===
        try:
            fig, axes = plt.subplots(nrows, ncols, figsize=(ncols*4, nrows*3), squeeze=False)

            for idx, cluster in enumerate(clusters):
                row = idx // ncols
                col = idx % ncols
                ax = axes[row, col]

                # Get cells in this cluster
                cluster_mask = adata.obs['cluster_coarse'] == cluster
                umi_vals = adata.obs.loc[cluster_mask, 'n_umi'].values

                if len(umi_vals) == 0:
                    ax.text(0.5, 0.5, 'No cells', ha='center', va='center')
                    ax.set_title(f'Cluster {cluster}')
                    continue

                # MAD standardization
                umi_median = np.median(umi_vals)
                umi_mad = median_abs_deviation(umi_vals)

                if umi_mad > 0:
                    z_scores = (umi_vals - umi_median) / umi_mad
                else:
                    z_scores = np.zeros_like(umi_vals)

                # ECDF
                sorted_z = np.sort(z_scores)
                ecdf_y = np.arange(1, len(sorted_z)+1) / len(sorted_z)

                ax.plot(ecdf_y, sorted_z, linewidth=2, color='steelblue')

                # Threshold line (lower only)
                ax.axhline(y=-3, color='red', linestyle='--', linewidth=1.5, label='Lower threshold (-3 MAD)')
                ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.5, alpha=0.5, label='Median')

                # Shading
                ax.fill_between([0, 1], -10, -3, color='red', alpha=0.1)

                ax.set_xlabel('Cumulative Probability')
                ax.set_ylabel('MADs from Median')
                ax.set_title(f'Cluster {cluster} (n={len(umi_vals)})')
                ax.set_xlim(0, 1)
                ax.set_ylim(-6, 6)
                ax.grid(True, alpha=0.3)
                if idx == 0:
                    ax.legend(fontsize=8)

            # Hide empty subplots
            for idx in range(n_clusters, nrows*ncols):
                row = idx // ncols
                col = idx % ncols
                axes[row, col].axis('off')

            plt.suptitle(f'{patient_id}: UMI Distribution (MAD-standardized)', fontsize=14, y=0.995)
            plt.tight_layout()

            umi_file = Path(output_dir) / f"{patient_id}_ecdf_umi.pdf"
            fig.savefig(umi_file, dpi=100, bbox_inches='tight')
            plt.close(fig)
            print(f"{patient_id}:   ✓ UMI ECDF: {umi_file}")

        except Exception as e:
            print(f"{patient_id}: ERROR generating UMI ECDF: {e}")

        # === GENES ECDF (both thresholds) ===
        try:
            fig, axes = plt.subplots(nrows, ncols, figsize=(ncols*4, nrows*3), squeeze=False)

            for idx, cluster in enumerate(clusters):
                row = idx // ncols
                col = idx % ncols
                ax = axes[row, col]

                cluster_mask = adata.obs['cluster_coarse'] == cluster
                gene_vals = adata.obs.loc[cluster_mask, 'n_genes'].values

                if len(gene_vals) == 0:
                    ax.text(0.5, 0.5, 'No cells', ha='center', va='center')
                    ax.set_title(f'Cluster {cluster}')
                    continue

                gene_median = np.median(gene_vals)
                gene_mad = median_abs_deviation(gene_vals)

                if gene_mad > 0:
                    z_scores = (gene_vals - gene_median) / gene_mad
                else:
                    z_scores = np.zeros_like(gene_vals)

                sorted_z = np.sort(z_scores)
                ecdf_y = np.arange(1, len(sorted_z)+1) / len(sorted_z)

                ax.plot(ecdf_y, sorted_z, linewidth=2, color='darkgreen')

                # Both thresholds
                ax.axhline(y=-3, color='red', linestyle='--', linewidth=1.5, label='Lower (-3 MAD)')
                ax.axhline(y=3, color='orange', linestyle='--', linewidth=1.5, label='Upper (+3 MAD)')
                ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.5, alpha=0.5)

                # Shading
                ax.fill_between([0, 1], -10, -3, color='red', alpha=0.1)
                ax.fill_between([0, 1], 3, 10, color='orange', alpha=0.1)

                ax.set_xlabel('Cumulative Probability')
                ax.set_ylabel('MADs from Median')
                ax.set_title(f'Cluster {cluster} (n={len(gene_vals)})')
                ax.set_xlim(0, 1)
                ax.set_ylim(-6, 6)
                ax.grid(True, alpha=0.3)
                if idx == 0:
                    ax.legend(fontsize=8)

            for idx in range(n_clusters, nrows*ncols):
                row = idx // ncols
                col = idx % ncols
                axes[row, col].axis('off')

            plt.suptitle(f'{patient_id}: Genes Distribution (MAD-standardized)', fontsize=14, y=0.995)
            plt.tight_layout()

            genes_file = Path(output_dir) / f"{patient_id}_ecdf_genes.pdf"
            fig.savefig(genes_file, dpi=100, bbox_inches='tight')
            plt.close(fig)
            print(f"{patient_id}:   ✓ Genes ECDF: {genes_file}")

        except Exception as e:
            print(f"{patient_id}: ERROR generating genes ECDF: {e}")

        # === MITO ECDF (upper threshold only) ===
        try:
            fig, axes = plt.subplots(nrows, ncols, figsize=(ncols*4, nrows*3), squeeze=False)

            for idx, cluster in enumerate(clusters):
                row = idx // ncols
                col = idx % ncols
                ax = axes[row, col]

                cluster_mask = adata.obs['cluster_coarse'] == cluster
                mito_vals = adata.obs.loc[cluster_mask, 'mito_percent'].values

                if len(mito_vals) == 0:
                    ax.text(0.5, 0.5, 'No cells', ha='center', va='center')
                    ax.set_title(f'Cluster {cluster}')
                    continue

                mito_median = np.median(mito_vals)
                mito_mad = median_abs_deviation(mito_vals)

                if mito_mad > 0:
                    z_scores = (mito_vals - mito_median) / mito_mad
                else:
                    z_scores = np.zeros_like(mito_vals)

                sorted_z = np.sort(z_scores)
                ecdf_y = np.arange(1, len(sorted_z)+1) / len(sorted_z)

                ax.plot(ecdf_y, sorted_z, linewidth=2, color='purple')

                # Upper threshold only
                ax.axhline(y=3, color='orange', linestyle='--', linewidth=1.5, label='Upper threshold (+3 MAD)')
                ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.5, alpha=0.5)

                # Shading
                ax.fill_between([0, 1], 3, 10, color='orange', alpha=0.1)

                ax.set_xlabel('Cumulative Probability')
                ax.set_ylabel('MADs from Median')
                ax.set_title(f'Cluster {cluster} (n={len(mito_vals)})')
                ax.set_xlim(0, 1)
                ax.set_ylim(-6, 6)
                ax.grid(True, alpha=0.3)
                if idx == 0:
                    ax.legend(fontsize=8)

            for idx in range(n_clusters, nrows*ncols):
                row = idx // ncols
                col = idx % ncols
                axes[row, col].axis('off')

            plt.suptitle(f'{patient_id}: Mito % Distribution (MAD-standardized)', fontsize=14, y=0.995)
            plt.tight_layout()

            mito_file = Path(output_dir) / f"{patient_id}_ecdf_mito.pdf"
            fig.savefig(mito_file, dpi=100, bbox_inches='tight')
            plt.close(fig)
            print(f"{patient_id}:   ✓ Mito ECDF: {mito_file}")

        except Exception as e:
            print(f"{patient_id}: ERROR generating mito ECDF: {e}")

    except Exception as e:
        print(f"{patient_id}: ERROR in MAD ECDF plotting: {e}")


def plot_summary_statistics(adata, output_dir, patient_id):
    """
    Generate 2x2 summary statistics panel.

    Panels:
    1. Stacked bar chart - outlier breakdown by cluster (% cells removed per cause)
    2. Violin plot - UMI distribution by cluster (ylim at 95th percentile)
    3. Violin plot - genes distribution by cluster (ylim at 95th percentile)
    4. Bar chart - cluster sizes (cell counts)

    Args:
        adata: Patient-level AnnData object
        output_dir: Path to visualizations directory
        patient_id: Patient identifier
    """
    print(f"{patient_id}: Generating summary statistics...")

    try:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        clusters = sorted(adata.obs['cluster_coarse'].unique())
        n_clusters = len(clusters)

        # Panel 1: Stacked bar chart of removal causes by cluster
        removal_data = {
            'genes_low': [],
            'genes_high': [],
            'umi_low': [],
            'mito_high': []
        }

        for cluster in clusters:
            cluster_mask = adata.obs['cluster_coarse'] == cluster
            n_cluster = cluster_mask.sum()

            if n_cluster > 0:
                removal_data['genes_low'].append(
                    100 * adata.obs.loc[cluster_mask, 'mad_outlier_genes_low'].sum() / n_cluster
                )
                removal_data['genes_high'].append(
                    100 * adata.obs.loc[cluster_mask, 'mad_outlier_genes_high'].sum() / n_cluster
                )
                removal_data['umi_low'].append(
                    100 * adata.obs.loc[cluster_mask, 'mad_outlier_umi_low'].sum() / n_cluster
                )
                removal_data['mito_high'].append(
                    100 * adata.obs.loc[cluster_mask, 'mad_outlier_mito_high'].sum() / n_cluster
                )
            else:
                removal_data['genes_low'].append(0)
                removal_data['genes_high'].append(0)
                removal_data['umi_low'].append(0)
                removal_data['mito_high'].append(0)

        # Stacked bar chart
        x_pos = np.arange(n_clusters)
        width = 0.6

        axes[0,0].bar(x_pos, removal_data['genes_low'], width, label='Genes (low)', color='#d62728')
        axes[0,0].bar(x_pos, removal_data['genes_high'], width, bottom=removal_data['genes_low'],
                     label='Genes (high)', color='#ff7f0e')
        axes[0,0].bar(x_pos, removal_data['umi_low'], width,
                     bottom=np.array(removal_data['genes_low']) + np.array(removal_data['genes_high']),
                     label='UMI (low)', color='#1f77b4')
        axes[0,0].bar(x_pos, removal_data['mito_high'], width,
                     bottom=np.array(removal_data['genes_low']) + np.array(removal_data['genes_high']) + np.array(removal_data['umi_low']),
                     label='Mito (high)', color='#9467bd')

        axes[0,0].set_xlabel('Cluster')
        axes[0,0].set_ylabel('% Cells Removed')
        axes[0,0].set_title('Cell Removal Breakdown by Cluster')
        axes[0,0].set_xticks(x_pos)
        axes[0,0].set_xticklabels([str(c) for c in clusters])
        axes[0,0].legend(fontsize=8, loc='upper right')

        # Panel 2: Violin plot of UMI by cluster (ylim at 95th percentile)
        adata.obs['cluster_coarse_cat'] = adata.obs['cluster_coarse'].astype(str)
        umi_q95 = np.percentile(adata.obs['n_umi'], 95)
        sns.violinplot(data=adata.obs, x='cluster_coarse_cat', y='n_umi', ax=axes[0,1])
        axes[0,1].set_xlabel('Cluster')
        axes[0,1].set_ylabel('n_UMI')
        axes[0,1].set_title('UMI Distribution by Cluster (ylim=q95)')
        axes[0,1].set_ylim(0, umi_q95)
        axes[0,1].tick_params(axis='x', rotation=45)

        # Panel 3: Violin plot of genes by cluster (ylim at 95th percentile)
        genes_q95 = np.percentile(adata.obs['n_genes'], 95)
        sns.violinplot(data=adata.obs, x='cluster_coarse_cat', y='n_genes', ax=axes[1,0])
        axes[1,0].set_xlabel('Cluster')
        axes[1,0].set_ylabel('n_genes')
        axes[1,0].set_title('Gene Distribution by Cluster (ylim=q95)')
        axes[1,0].set_ylim(0, genes_q95)
        axes[1,0].tick_params(axis='x', rotation=45)

        # Panel 4: Cluster sizes
        cluster_sizes = adata.obs['cluster_coarse'].value_counts().sort_index()
        axes[1,1].bar(range(len(cluster_sizes)), cluster_sizes.values, color='steelblue')
        axes[1,1].set_xlabel('Cluster')
        axes[1,1].set_ylabel('Cell Count')
        axes[1,1].set_title('Cluster Sizes')
        axes[1,1].set_xticks(range(len(cluster_sizes)))
        axes[1,1].set_xticklabels(cluster_sizes.index)

        plt.tight_layout()
        output_file = Path(output_dir) / f"{patient_id}_summary.pdf"
        fig.savefig(output_file, dpi=100, bbox_inches='tight')
        plt.close(fig)

        # Clean up temporary column
        adata.obs.drop(columns=['cluster_coarse_cat'], inplace=True)

        print(f"{patient_id}:   ✓ Summary statistics: {output_file}")

    except Exception as e:
        print(f"{patient_id}: ERROR generating summary statistics: {e}")


def generate_all_visualizations(adata, output_dir, patient_id):
    """
    Generate all Step 02 diagnostic plots.

    This is the main entry point called by the wrapper. Generates:
    1. UMAP overview (2x2 panel: cluster, sample, pass/fail, UMI)
    2. Low-UMI features (2x2 panel: MAD-standardized UMI, mito, genes, outliers)
    3. MAD ECDF plots (3 separate faceted plots: UMI, genes, mito)
    4. Summary statistics (2x2 panel: outlier breakdown, distributions, sizes)

    Defensive coding:
    - Skips UMAP plots if coordinates missing
    - Continues on individual plot failures (doesn't crash entire step)
    - Prints clear status messages for each plot

    Args:
        adata: Patient-level AnnData object with filtering results
        output_dir: Path to visualizations directory
        patient_id: Patient identifier
    """
    print(f"{patient_id}: Generating visualizations...")

    # Ensure output directory exists
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Generate each visualization type
    plot_umap_overview(adata, output_dir, patient_id)
    plot_low_umi_features(adata, output_dir, patient_id)
    plot_mad_ecdf(adata, output_dir, patient_id)
    plot_summary_statistics(adata, output_dir, patient_id)

    print(f"{patient_id}: Visualizations complete")
