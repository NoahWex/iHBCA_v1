"""
Step 07 Final Report - Plot Generation Module

Create matplotlib visualizations for Section 3.
Pure computation module - NO file I/O (returns Figure objects, does not save).

Author: Agent 3 (Algorithm Developer)
Date: 2025-11-03
"""

import matplotlib.pyplot as plt
import matplotlib.figure as mpl_figure
import matplotlib.patches as mpatches
import pandas as pd
import numpy as np
import seaborn as sns
from typing import Tuple, List, Dict


def set_plotting_style():
    """
    Set consistent matplotlib style for all plots.

    Configures default plot parameters for publication-quality figures.

    Example:
        >>> set_plotting_style()
        >>> # Now all subsequent plots use this style
    """
    sns.set_style("whitegrid")
    plt.rcParams.update({
        'figure.dpi': 100,
        'savefig.dpi': 300,
        'font.size': 10,
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'legend.fontsize': 9,
        'figure.titlesize': 14,
        'font.family': 'sans-serif'
    })


def plot_distribution(
    df: pd.DataFrame,
    metric_col: str,
    flagged_col: str,
    xlabel: str,
    title: str,
    color: str = 'blue'
) -> mpl_figure.Figure:
    """
    Create histogram with flagged samples highlighted.

    Args:
        df: DataFrame with metrics
        metric_col: Column to plot distribution
        flagged_col: Column indicating flagged status
        xlabel: X-axis label
        title: Plot title
        color: Histogram color (default 'blue')

    Returns:
        matplotlib Figure object (NOT saved)

    Plot elements:
        - Histogram of all samples (colored)
        - Vertical line at median (black dashed)
        - Red markers for flagged samples (top of histogram)
        - Legend showing flagged count

    Notes:
        - Return Figure object, do NOT call fig.savefig()
        - Use consistent color scheme (gray histogram, red flags)

    Example:
        >>> fig = plot_distribution(
        ...     df,
        ...     'cell_removal_rate',
        ...     'cell_removal_rate_flagged',
        ...     'Cell Removal Rate',
        ...     'Distribution of Cell Removal Rates'
        ... )
    """

    fig, ax = plt.subplots(figsize=(10, 6))

    # Get values (drop NaN)
    values = df[metric_col].dropna()

    if len(values) == 0:
        ax.text(0.5, 0.5, 'No data available', ha='center', va='center',
                transform=ax.transAxes, fontsize=14)
        return fig

    # Histogram
    ax.hist(values, bins=20, alpha=0.6, color=color, edgecolor='black')

    # Median line
    median_val = values.median()
    ax.axvline(median_val, color='black', linestyle='--', linewidth=2,
               label=f'Median: {median_val:.3f}')

    # Flagged samples (red vertical lines)
    if flagged_col in df.columns:
        flagged = df[df[flagged_col] == True]
        n_flagged = len(flagged)

        if n_flagged > 0:
            for _, row in flagged.iterrows():
                ax.axvline(row[metric_col], color='red', alpha=0.7, linewidth=1.5)

            # Add flagged count to legend
            ax.plot([], [], 'r-', linewidth=1.5, label=f'Flagged: {n_flagged}')

    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel('Frequency', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(loc='upper right')
    plt.tight_layout()

    return fig


def plot_patient_bar(
    df: pd.DataFrame,
    patient_col: str,
    value_cols: List[str],
    ylabel: str,
    title: str,
    colors: List[str]
) -> mpl_figure.Figure:
    """
    Create stacked bar chart per patient.

    Args:
        df: DataFrame with metrics
        patient_col: Column with patient IDs
        value_cols: Columns to stack (e.g., ['n_cells_removed_step02', 'n_doublets'])
        ylabel: Y-axis label
        title: Plot title
        colors: List of colors for each value_col

    Returns:
        matplotlib Figure object

    Plot elements:
        - X-axis: Patient IDs
        - Y-axis: Summed values
        - Stacked bars colored by value_col
        - Legend identifying each stack

    Notes:
        - Aggregate df by patient (groupby + sum)
        - Use plt.bar() with bottom parameter for stacking

    Example:
        >>> fig = plot_patient_bar(
        ...     df,
        ...     'patient_id',
        ...     ['n_cells_removed_step02', 'n_doublets'],
        ...     'Total Cells Removed',
        ...     'Cell Removal by Patient',
        ...     ['#f39c12', '#e74c3c']
        ... )
    """

    fig, ax = plt.subplots(figsize=(10, 6))

    # Aggregate by patient
    agg_df = df.groupby(patient_col)[value_cols].sum().reset_index()

    if len(agg_df) == 0:
        ax.text(0.5, 0.5, 'No data available', ha='center', va='center',
                transform=ax.transAxes, fontsize=14)
        return fig

    # X positions
    x = np.arange(len(agg_df))
    width = 0.6

    # Create stacked bars
    bottom = np.zeros(len(agg_df))

    for i, (col, color) in enumerate(zip(value_cols, colors)):
        values = agg_df[col].values
        ax.bar(x, values, width, bottom=bottom, label=col, color=color)
        bottom += values

    # Set x-axis labels
    ax.set_xticks(x)
    ax.set_xticklabels(agg_df[patient_col], rotation=45, ha='right')

    ax.set_xlabel('Patient', fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(loc='upper right')
    plt.tight_layout()

    return fig


def plot_vf_retention_scatter(
    df: pd.DataFrame,
    baseline_col: str,
    postqc_col: str,
    flagged_col: str,
    patient_col: str
) -> mpl_figure.Figure:
    """
    Create scatter plot: baseline VFs vs post-QC VFs.

    Args:
        df: DataFrame with VF counts
        baseline_col: Column with baseline VF counts (Step 01)
        postqc_col: Column with post-QC VF counts (Step 05)
        flagged_col: Column indicating flagged samples
        patient_col: Column with patient IDs (for coloring)

    Returns:
        matplotlib Figure object

    Plot elements:
        - Scatter plot: x=baseline, y=post-QC
        - Diagonal reference line (y=x, perfect retention)
        - Flagged samples annotated with sample_id
        - Color samples by patient_id

    Notes:
        - Add diagonal line: plt.plot([xmin, xmax], [xmin, xmax], 'k--', alpha=0.5)
        - Annotate flagged samples: plt.text(x, y, sample_id)

    Example:
        >>> fig = plot_vf_retention_scatter(
        ...     df,
        ...     'n_baseline_vfs',
        ...     'n_post_qc_vfs',
        ...     'vf_retention_rate_flagged',
        ...     'patient_id'
        ... )
    """

    fig, ax = plt.subplots(figsize=(10, 8))

    if len(df) == 0:
        ax.text(0.5, 0.5, 'No data available', ha='center', va='center',
                transform=ax.transAxes, fontsize=14)
        return fig

    # Color palette for patients
    patients = df[patient_col].unique()
    colors = plt.cm.Set2(np.linspace(0, 1, len(patients)))
    patient_colors = {patient: colors[i] for i, patient in enumerate(patients)}

    # Scatter plot for each patient
    for patient_id in patients:
        patient_df = df[df[patient_col] == patient_id]
        ax.scatter(
            patient_df[baseline_col],
            patient_df[postqc_col],
            label=patient_id,
            alpha=0.6,
            s=80,
            color=patient_colors[patient_id]
        )

    # Diagonal reference line (y=x, perfect retention)
    max_vf = max(df[baseline_col].max(), df[postqc_col].max())
    min_vf = min(df[baseline_col].min(), df[postqc_col].min())
    ax.plot([min_vf, max_vf], [min_vf, max_vf], 'k--', alpha=0.3, linewidth=2,
            label='Perfect retention (y=x)')

    # Annotate flagged samples
    if flagged_col in df.columns:
        flagged = df[df[flagged_col] == True]
        for _, row in flagged.iterrows():
            ax.annotate(
                row['sample_id'],
                (row[baseline_col], row[postqc_col]),
                fontsize=8,
                color='red',
                ha='right',
                va='bottom',
                xytext=(-5, 5),
                textcoords='offset points'
            )

    ax.set_xlabel('Baseline Variable Features (Step 01)', fontsize=12)
    ax.set_ylabel('Post-QC Variable Features (Step 05)', fontsize=12)
    ax.set_title('Variable Feature Retention: Baseline vs Post-QC', fontsize=14, fontweight='bold')
    ax.legend(loc='upper left')
    plt.tight_layout()

    return fig


def generate_section3_plots(
    df: pd.DataFrame,
    metrics_config: dict
) -> Dict[str, mpl_figure.Figure]:
    """
    Generate all Section 3 plots.

    Args:
        df: DataFrame with metrics and flags
        metrics_config: Dict with plot parameters (titles, labels, colors, etc.)

    Returns:
        Dict of {plot_name: Figure} for all Section 3 plots:
            - 'cell_removal_distribution': distribution plot
            - 'cell_removal_by_patient': stacked bar
            - 'vf_retention_distribution': distribution plot
            - 'vf_retention_scatter': scatter plot

    Notes:
        - Calls individual plot functions
        - Returns Figure objects (NOT saved)
        - Wrapper will save these to plots/ directory

    Example:
        >>> config = {
        ...     'color_removal': '#3498db',
        ...     'color_vf': '#2ecc71',
        ...     'color_step02': '#f39c12',
        ...     'color_step03': '#e74c3c'
        ... }
        >>> plots = generate_section3_plots(df, config)
        >>> # Returns dict of Figure objects
    """

    # Set plotting style
    set_plotting_style()

    plots = {}

    # Get colors from config
    color_removal = metrics_config.get('color_removal', '#3498db')
    color_vf = metrics_config.get('color_vf', '#2ecc71')
    color_step02 = metrics_config.get('color_step02', '#f39c12')
    color_step03 = metrics_config.get('color_step03', '#e74c3c')

    # 1. Cell Removal Distribution
    plots['cell_removal_distribution'] = plot_distribution(
        df=df,
        metric_col='cell_removal_rate',
        flagged_col='cell_removal_rate_flagged',
        xlabel='Cell Removal Rate',
        title='Distribution of Cell Removal Rates Across All Samples',
        color=color_removal
    )

    # 2. Cell Removal by Patient (stacked bar)
    plots['cell_removal_by_patient'] = plot_patient_bar(
        df=df,
        patient_col='patient_id',
        value_cols=['n_cells_removed_step02', 'n_doublets'],
        ylabel='Total Cells Removed',
        title='Cell Removal Breakdown by Patient',
        colors=[color_step02, color_step03]
    )

    # 3. VF Retention Distribution
    plots['vf_retention_distribution'] = plot_distribution(
        df=df,
        metric_col='vf_retention_rate',
        flagged_col='vf_retention_rate_flagged',
        xlabel='VF Retention Rate',
        title='Distribution of VF Retention Rates Across All Samples',
        color=color_vf
    )

    # 4. VF Retention Scatter
    plots['vf_retention_scatter'] = plot_vf_retention_scatter(
        df=df,
        baseline_col='n_baseline_vfs',
        postqc_col='n_post_qc_vfs',
        flagged_col='vf_retention_rate_flagged',
        patient_col='patient_id'
    )

    return plots
