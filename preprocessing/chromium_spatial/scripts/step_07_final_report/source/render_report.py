"""
Step 07 Final Report - HTML Report Rendering Module

Populate HTML report template and convert to self-contained HTML.
Pure computation module - NO file I/O for reading images (takes Figure objects).

Author: Agent 3 (Algorithm Developer)
Date: 2025-11-03
"""

from jinja2 import Template
import pandas as pd
import numpy as np
from typing import Dict
import base64
import io
import matplotlib.figure as mpl_figure
import matplotlib.pyplot as plt


def create_html_template() -> str:
    """
    Create Jinja2 HTML template for the report.

    Returns:
        HTML template string with placeholders for:
            - {{ summary_tables }} (Section 1)
            - {{ flagged_samples_table }} (Section 2)
            - {{ section3_plots }} (Section 3, base64 embedded)
            - {{ appendix_plots }} (Appendix, base64 embedded)

    Template structure:
        <!DOCTYPE html>
        <html>
        <head>
            <title>Preprocessing QC Report</title>
            <style>
                /* CSS for tables, layout, navigation */
            </style>
        </head>
        <body>
            <h1>Preprocessing QC Report</h1>

            <nav>
                <!-- Table of contents with links -->
            </nav>

            <section id="summary">
                <h2>Section 1: Dataset Summary</h2>
                {{ summary_tables }}
            </section>

            <section id="flagged">
                <h2>Section 2: Flagged Items</h2>
                {{ flagged_samples_table }}
            </section>

            <section id="evidence">
                <h2>Section 3: Evidence by Metric</h2>
                {{ section3_plots }}
            </section>

            <section id="appendix">
                <h2>Appendix: Per-Patient Detailed Plots</h2>
                {{ appendix_plots }}
            </section>
        </body>
        </html>

    Notes:
        - Include CSS for table styling (borders, alternating rows)
        - Add navigation links (jump to sections)
        - Responsive layout (readable on desktop/tablet)

    Example:
        >>> template_str = create_html_template()
        >>> template = Template(template_str)
    """

    html_template = """
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>{{ title }}</title>
    <style>
        body {
            font-family: Arial, Helvetica, sans-serif;
            margin: 40px;
            max-width: 1400px;
            margin-left: auto;
            margin-right: auto;
            background-color: #f8f9fa;
        }

        h1 {
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 15px;
            margin-bottom: 20px;
        }

        h2 {
            color: #34495e;
            margin-top: 40px;
            border-bottom: 2px solid #95a5a6;
            padding-bottom: 10px;
        }

        h3 {
            color: #7f8c8d;
            margin-top: 30px;
        }

        .metadata {
            background-color: #ecf0f1;
            padding: 15px;
            margin: 20px 0;
            border-radius: 5px;
        }

        .toc {
            background-color: #e8f4f8;
            padding: 20px;
            margin: 20px 0;
            border-radius: 5px;
            border-left: 5px solid #3498db;
        }

        .toc h3 {
            margin-top: 0;
            color: #2c3e50;
        }

        .toc ul {
            list-style-type: none;
            padding-left: 0;
        }

        .toc li {
            margin: 10px 0;
        }

        .toc a {
            text-decoration: none;
            color: #2980b9;
            font-size: 16px;
        }

        .toc a:hover {
            color: #3498db;
            text-decoration: underline;
        }

        table {
            border-collapse: collapse;
            width: 100%;
            margin: 20px 0;
            background-color: white;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }

        th, td {
            border: 1px solid #ddd;
            padding: 12px;
            text-align: left;
        }

        th {
            background-color: #3498db;
            color: white;
            font-weight: bold;
        }

        tr:nth-child(even) {
            background-color: #f2f2f2;
        }

        tr:hover {
            background-color: #e8f4f8;
        }

        .table-caption {
            font-weight: bold;
            margin: 20px 0 10px 0;
            font-size: 16px;
            color: #2c3e50;
        }

        img {
            max-width: 100%;
            height: auto;
            margin: 20px 0;
            border: 1px solid #ddd;
            border-radius: 5px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }

        .flagged {
            color: #e74c3c;
            font-weight: bold;
        }

        .plot-container {
            margin: 30px 0;
            text-align: center;
        }

        .plot-title {
            font-weight: bold;
            margin: 15px 0;
            font-size: 14px;
            color: #34495e;
        }

        section {
            background-color: white;
            padding: 30px;
            margin: 30px 0;
            border-radius: 5px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }

        .patient-section {
            margin: 40px 0;
            padding: 20px;
            border-left: 4px solid #3498db;
            background-color: #f8f9fa;
        }

        .warning {
            background-color: #fff3cd;
            border-left: 5px solid #ffc107;
            padding: 15px;
            margin: 20px 0;
        }

        .info {
            background-color: #d1ecf1;
            border-left: 5px solid #17a2b8;
            padding: 15px;
            margin: 20px 0;
        }
    </style>
</head>
<body>
    <h1>{{ title }}</h1>

    <div class="metadata">
        <p><strong>Generated:</strong> {{ timestamp }}</p>
        <p><strong>MAD Threshold:</strong> {{ mad_threshold }}</p>
        <p><strong>Total Samples:</strong> {{ n_samples }}</p>
        <p><strong>Flagged Samples:</strong> {{ n_flagged }}</p>
    </div>

    <div class="toc">
        <h3>Table of Contents</h3>
        <ul>
            <li><a href="#section1">Section 1: Dataset Summary</a></li>
            <li><a href="#section2">Section 2: Flagged Items</a></li>
            <li><a href="#section3">Section 3: Evidence by Metric</a></li>
            <li><a href="#appendix">Appendix: Per-Patient Detailed Plots</a></li>
        </ul>
    </div>

    <section id="section1">
        <h2>Section 1: Dataset Summary</h2>
        {{ summary_tables | safe }}
    </section>

    <section id="section2">
        <h2>Section 2: Flagged Items</h2>
        <div class="info">
            <p><strong>Methodology:</strong> Samples flagged using MAD threshold = {{ mad_threshold }} on dataset-wide distributions</p>
            <p><strong>Priority Metrics:</strong> {{ priority_metrics }}</p>
        </div>
        {{ flagged_samples_table | safe }}
    </section>

    <section id="section3">
        <h2>Section 3: Evidence by Metric</h2>
        {{ section3_plots | safe }}
    </section>

    <section id="appendix">
        <h2>Appendix: Per-Patient Detailed Plots</h2>
        {{ appendix_plots | safe }}
    </section>

</body>
</html>
"""

    return html_template


def dataframe_to_html_table(df: pd.DataFrame, title: str) -> str:
    """
    Convert DataFrame to HTML table.

    Args:
        df: DataFrame to render
        title: Table title

    Returns:
        HTML string with styled table

    Notes:
        - Use df.to_html() with index=False, classes=['table']
        - Add <caption> with title
        - Format numbers (e.g., percentages to 1 decimal)

    Example:
        >>> html = dataframe_to_html_table(df, 'Table 1.1: Input Counts')
    """

    if len(df) == 0:
        return f'<div class="warning"><p><strong>{title}</strong></p><p>No data available</p></div>'

    # Create copy for formatting
    df_display = df.copy()

    # Format numeric columns
    for col in df_display.columns:
        if df_display[col].dtype in [np.float64, np.float32]:
            # Check if values look like percentages (0-1 range)
            if df_display[col].max() <= 1.0 and df_display[col].min() >= 0:
                df_display[col] = df_display[col].apply(lambda x: f'{x:.3f}' if pd.notna(x) else '')
            else:
                df_display[col] = df_display[col].apply(lambda x: f'{x:.2f}' if pd.notna(x) else '')
        elif df_display[col].dtype in [np.int64, np.int32]:
            df_display[col] = df_display[col].apply(lambda x: f'{x:,}' if pd.notna(x) else '')

    # Convert to HTML
    html = df_display.to_html(index=False, escape=False, border=0)

    # Add caption
    table_html = f'<div class="table-caption">{title}</div>\n{html}'

    return table_html


def embed_plot_as_base64(
    fig: mpl_figure.Figure,
    format: str = 'png',
    dpi: int = 300
) -> str:
    """
    Convert matplotlib Figure to base64-encoded image string.

    Args:
        fig: matplotlib Figure object
        format: Image format ('png' or 'svg')
        dpi: Resolution for PNG

    Returns:
        Base64 string: "data:image/png;base64,iVBORw0KG..."

    Implementation:
        1. Save figure to BytesIO buffer
        2. Encode buffer as base64
        3. Return data URI string

    Notes:
        - Use fig.savefig(buf, format=format, dpi=dpi, bbox_inches='tight')
        - Close figure after encoding: plt.close(fig)

    Example:
        >>> fig, ax = plt.subplots()
        >>> ax.plot([1, 2, 3], [1, 2, 3])
        >>> b64_str = embed_plot_as_base64(fig)
        >>> # Use in HTML: <img src="{b64_str}">
    """

    # Create BytesIO buffer
    buf = io.BytesIO()

    # Save figure to buffer
    fig.savefig(buf, format=format, dpi=dpi, bbox_inches='tight')

    # Get buffer contents
    buf.seek(0)
    image_bytes = buf.read()

    # Encode as base64
    b64_encoded = base64.b64encode(image_bytes).decode('utf-8')

    # Close buffer and figure
    buf.close()
    plt.close(fig)

    # Return data URI
    mime_type = f'image/{format}'
    data_uri = f'data:{mime_type};base64,{b64_encoded}'

    return data_uri


def embed_image_file_as_base64(image_path: str) -> str:
    """
    Load image file and convert to base64-encoded string.

    Args:
        image_path: Path to image file

    Returns:
        Base64 data URI string

    Notes:
        - This function DOES perform file I/O (reading image files)
        - Only used for embedding Step 06 plots (appendix)
        - Wrapper will provide file paths

    Example:
        >>> b64_str = embed_image_file_as_base64('/path/to/plot.png')
    """

    with open(image_path, 'rb') as f:
        image_bytes = f.read()

    b64_encoded = base64.b64encode(image_bytes).decode('utf-8')

    # Detect format from extension
    if image_path.lower().endswith('.png'):
        mime_type = 'image/png'
    elif image_path.lower().endswith('.jpg') or image_path.lower().endswith('.jpeg'):
        mime_type = 'image/jpeg'
    elif image_path.lower().endswith('.svg'):
        mime_type = 'image/svg+xml'
    else:
        mime_type = 'image/png'  # Default

    data_uri = f'data:{mime_type};base64,{b64_encoded}'

    return data_uri


def render_html_report(
    df_metrics: pd.DataFrame,
    df_flagged: pd.DataFrame,
    section3_figs: Dict[str, mpl_figure.Figure],
    appendix_plot_paths: Dict[str, Dict[str, str]],
    config: dict
) -> str:
    """
    Render complete HTML report.

    Args:
        df_metrics: Full metrics DataFrame
        df_flagged: Flagged samples DataFrame (long format with metric column)
        section3_figs: Dict of {plot_name: Figure} for Section 3
        appendix_plot_paths: Dict of {patient: {plot_type: file_path}} for Step 06 plots
        config: Report configuration (MAD threshold, title, etc.)

    Returns:
        Complete HTML string (self-contained)

    Process:
        1. Create summary tables (Section 1) from df_metrics
        2. Create flagged items table (Section 2) from df_flagged
        3. Embed Section 3 plots as base64 images
        4. Load and embed appendix plots as base64
        5. Populate Jinja2 template
        6. Return rendered HTML

    Notes:
        - All images embedded as base64 (no external files)
        - Include metadata: generation date, MAD threshold
        - Add clickable links from flagged table to appendix

    Example:
        >>> html = render_html_report(
        ...     df_metrics,
        ...     df_flagged,
        ...     section3_figs,
        ...     appendix_plot_paths,
        ...     config
        ... )
    """

    from datetime import datetime

    # Extract config
    title = config.get('report_title', 'Preprocessing QC Report')
    mad_threshold = config.get('mad_threshold', 2.0)
    priority_metrics = config.get('priority_metrics', ['cell_removal_rate', 'vf_retention_rate'])

    # Create template
    template = Template(create_html_template())

    # Section 1: Summary Tables
    summary_tables_html = ''

    # Table 1.1: Input Counts
    if len(df_metrics) > 0:
        table_1_1 = df_metrics[['patient_id', 'sample_id', 'n_cells_baseline', 'n_baseline_vfs']].copy()
        table_1_1.columns = ['Patient', 'Sample ID', 'Baseline Cells', 'Baseline VFs']
        summary_tables_html += dataframe_to_html_table(table_1_1, 'Table 1.1: Input Counts')

        # Table 1.2: Output Counts
        table_1_2 = df_metrics[['patient_id', 'sample_id', 'n_cells_final', 'n_post_qc_vfs']].copy()
        cell_retention = (df_metrics['n_cells_final'] / df_metrics['n_cells_baseline'] * 100)
        vf_retention = (df_metrics['vf_retention_rate'] * 100)
        table_1_2['Cell Retention %'] = cell_retention
        table_1_2['VF Retention %'] = vf_retention
        table_1_2.columns = ['Patient', 'Sample ID', 'Final Cells', 'Post-QC VFs', 'Cell Retention %', 'VF Retention %']
        summary_tables_html += dataframe_to_html_table(table_1_2, 'Table 1.2: Output Counts')

        # Table 1.3: Removal Accounting
        table_1_3 = df_metrics[['patient_id', 'sample_id', 'n_cells_removed_step02', 'n_doublets']].copy()
        table_1_3['Total Removed'] = table_1_3['n_cells_removed_step02'] + table_1_3['n_doublets']
        table_1_3['Removal Rate %'] = (df_metrics['cell_removal_rate'] * 100)
        table_1_3.columns = ['Patient', 'Sample ID', 'Removed by Step 02', 'Removed by Step 03',
                             'Total Removed', 'Removal Rate %']
        summary_tables_html += dataframe_to_html_table(table_1_3, 'Table 1.3: Removal Accounting')

    # Section 2: Flagged Items Table
    if len(df_flagged) > 0:
        table_2_1 = df_flagged.copy()
        # Add link to appendix
        table_2_1['Link'] = table_2_1['patient_id'].apply(lambda x: f'<a href="#appendix_{x}">View Details</a>')
        table_2_1 = table_2_1[['sample_id', 'patient_id', 'metric', 'value', 'mad_score', 'Link']]
        table_2_1.columns = ['Sample ID', 'Patient', 'Metric', 'Value', 'MAD Score', 'Link to Appendix']
        flagged_table_html = dataframe_to_html_table(table_2_1, 'Table 2.1: Flagged Samples')
    else:
        flagged_table_html = '<div class="info"><p><strong>No samples flagged</strong></p></div>'

    # Section 3: Evidence Plots
    section3_html = ''

    section3_html += '<h3>3A. Cell Removal (Steps 02 + 03)</h3>\n'
    if 'cell_removal_distribution' in section3_figs:
        b64_img = embed_plot_as_base64(section3_figs['cell_removal_distribution'])
        section3_html += f'<div class="plot-container">\n'
        section3_html += f'  <div class="plot-title">Cell Removal Rate Distribution</div>\n'
        section3_html += f'  <img src="{b64_img}" alt="Cell Removal Distribution">\n'
        section3_html += f'</div>\n'

    if 'cell_removal_by_patient' in section3_figs:
        b64_img = embed_plot_as_base64(section3_figs['cell_removal_by_patient'])
        section3_html += f'<div class="plot-container">\n'
        section3_html += f'  <div class="plot-title">Cell Removal by Patient</div>\n'
        section3_html += f'  <img src="{b64_img}" alt="Cell Removal by Patient">\n'
        section3_html += f'</div>\n'

    section3_html += '<h3>3B. Variable Feature Retention (Steps 01 → 04)</h3>\n'
    if 'vf_retention_distribution' in section3_figs:
        b64_img = embed_plot_as_base64(section3_figs['vf_retention_distribution'])
        section3_html += f'<div class="plot-container">\n'
        section3_html += f'  <div class="plot-title">VF Retention Rate Distribution</div>\n'
        section3_html += f'  <img src="{b64_img}" alt="VF Retention Distribution">\n'
        section3_html += f'</div>\n'

    if 'vf_retention_scatter' in section3_figs:
        b64_img = embed_plot_as_base64(section3_figs['vf_retention_scatter'])
        section3_html += f'<div class="plot-container">\n'
        section3_html += f'  <div class="plot-title">VF Retention Scatter</div>\n'
        section3_html += f'  <img src="{b64_img}" alt="VF Retention Scatter">\n'
        section3_html += f'</div>\n'

    # Appendix: Per-Patient Plots
    appendix_html = ''

    for patient_id in sorted(appendix_plot_paths.keys()):
        patient_plots = appendix_plot_paths[patient_id]

        appendix_html += f'<div class="patient-section" id="appendix_{patient_id}">\n'
        appendix_html += f'<h3>Patient: {patient_id}</h3>\n'

        # Patient summary plots
        if 'patient_summary' in patient_plots:
            appendix_html += '<h4>Overview</h4>\n'
            for plot_name, plot_path in patient_plots['patient_summary'].items():
                try:
                    b64_img = embed_image_file_as_base64(plot_path)
                    appendix_html += f'<div class="plot-container">\n'
                    appendix_html += f'  <div class="plot-title">{plot_name.replace("_", " ").title()}</div>\n'
                    appendix_html += f'  <img src="{b64_img}" alt="{plot_name}">\n'
                    appendix_html += f'</div>\n'
                except Exception as e:
                    appendix_html += f'<div class="warning"><p>Error loading {plot_name}: {e}</p></div>\n'

        # Sample detail plots
        if 'sample_detail' in patient_plots:
            appendix_html += '<h4>Sample Detail</h4>\n'
            # sample_detail is a flat dict: {'plot_name': '/path/to/plot.png'}
            for plot_name, plot_path in patient_plots['sample_detail'].items():
                try:
                    b64_img = embed_image_file_as_base64(plot_path)
                    appendix_html += f'<div class="plot-container">\n'
                    appendix_html += f'  <div class="plot-title">{plot_name.replace("_", " ").title()}</div>\n'
                    appendix_html += f'  <img src="{b64_img}" alt="{plot_name}">\n'
                    appendix_html += f'</div>\n'
                except Exception as e:
                    appendix_html += f'<div class="warning"><p>Error loading {plot_name}: {e}</p></div>\n'

        appendix_html += '</div>\n'

    # Render template
    html_content = template.render(
        title=title,
        timestamp=datetime.now().isoformat(),
        mad_threshold=mad_threshold,
        priority_metrics=', '.join(priority_metrics),
        n_samples=len(df_metrics),
        n_flagged=len(df_flagged),
        summary_tables=summary_tables_html,
        flagged_samples_table=flagged_table_html,
        section3_plots=section3_html,
        appendix_plots=appendix_html
    )

    return html_content
