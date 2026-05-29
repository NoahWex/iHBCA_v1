#!/usr/bin/env python3
"""
load_counts.py - Load Xenium count matrices with three compartment layers

Builds sparse count matrices from transcripts.parquet, separating transcripts
by overlaps_nucleus to create whole/nuclear/cytoplasmic layers.

Usage:
    from shared.loaders import load_counts

    # Create new AnnData with counts
    adata = load_counts(sample_path)

    # Or add to existing AnnData
    adata = load_counts(sample_path, adata=existing_adata)

Output:
    adata.X             - Whole cell counts (cells x genes)
    adata.layers['whole']       - Same as X
    adata.layers['nuclear']     - Transcripts with overlaps_nucleus == 1
    adata.layers['cytoplasmic'] - Transcripts with overlaps_nucleus == 0
    adata.obsm['spatial']       - Cell centroids (x, y)
    adata.obs                   - Cell metadata including QC metrics
"""

from pathlib import Path
from typing import Dict, Any, Optional, Union
import warnings

import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad


def build_count_matrix(
    transcripts: pd.DataFrame,
    cell_ids: np.ndarray,
    genes: np.ndarray
) -> sparse.csr_matrix:
    """
    Build sparse count matrix from transcripts dataframe.

    Args:
        transcripts: DataFrame with feature_name and cell_id columns
        cell_ids: Array of cell IDs (defines column order)
        genes: Array of gene names (defines row order)

    Returns:
        Sparse CSR matrix (genes x cells), transposed to (cells x genes) by caller
    """
    gene_to_idx = {g: i for i, g in enumerate(genes)}
    cell_to_idx = {c: i for i, c in enumerate(cell_ids)}

    # Filter to known genes and cells
    valid_mask = (
        transcripts['feature_name'].isin(gene_to_idx) &
        transcripts['cell_id'].isin(cell_to_idx)
    )
    valid_transcripts = transcripts[valid_mask]

    rows = valid_transcripts['feature_name'].map(gene_to_idx).values
    cols = valid_transcripts['cell_id'].map(cell_to_idx).values
    data = np.ones(len(rows), dtype=np.float32)

    matrix = sparse.coo_matrix(
        (data, (rows, cols)),
        shape=(len(genes), len(cell_ids))
    ).tocsr()

    return matrix


def load_counts(
    sample_path: Union[str, Path],
    adata: Optional[ad.AnnData] = None,
    qv_cutoff: int = 20,
    metadata: Optional[Dict[str, Any]] = None,
) -> ad.AnnData:
    """
    Load Xenium count data with three compartment layers.

    Args:
        sample_path: Path to Xenium sample directory
        adata: Existing AnnData to add counts to (creates new if None)
        qv_cutoff: Quality value cutoff for transcripts (default: 20)
        metadata: Optional dict with sample metadata (xenium_id, patient_id, etc.)

    Returns:
        AnnData with:
            - X: whole cell counts
            - layers['whole', 'nuclear', 'cytoplasmic']
            - obsm['spatial']: cell centroids
            - obs: cell metadata and QC metrics
    """
    sample_path = Path(sample_path)

    # Load cells metadata
    cells_path = sample_path / "cells.parquet"
    if not cells_path.exists():
        cells_path = sample_path / "cells.csv.gz"
    if not cells_path.exists():
        raise FileNotFoundError(f"No cells file found in {sample_path}")

    if str(cells_path).endswith('.parquet'):
        try:
            cells_df = pd.read_parquet(cells_path)
        except Exception as e:
            csv_fallback = sample_path / "cells.csv.gz"
            if csv_fallback.exists():
                warnings.warn(f"Corrupt parquet {cells_path}, falling back to CSV: {e}")
                cells_df = pd.read_csv(csv_fallback)
            else:
                raise
    else:
        cells_df = pd.read_csv(cells_path)

    # Decode bytes if needed
    if cells_df['cell_id'].dtype == object and isinstance(cells_df['cell_id'].iloc[0], bytes):
        cells_df['cell_id'] = cells_df['cell_id'].str.decode('utf-8')

    # Load transcripts
    transcripts_path = sample_path / "transcripts.parquet"
    if not transcripts_path.exists():
        transcripts_path = sample_path / "transcripts.csv.gz"
    if not transcripts_path.exists():
        raise FileNotFoundError(f"No transcripts file found in {sample_path}")

    if str(transcripts_path).endswith('.parquet'):
        try:
            transcripts = pd.read_parquet(transcripts_path)
        except Exception as e:
            # Parquet files truncated at block boundaries are surprisingly common
            csv_fallback = sample_path / "transcripts.csv.gz"
            if csv_fallback.exists():
                warnings.warn(f"Corrupt parquet {transcripts_path}, falling back to CSV: {e}")
                transcripts = pd.read_csv(csv_fallback)
            else:
                raise
    else:
        transcripts = pd.read_csv(transcripts_path)

    # Decode bytes if needed
    for col in ['feature_name', 'cell_id']:
        if transcripts[col].dtype == object and isinstance(transcripts[col].iloc[0], bytes):
            transcripts[col] = transcripts[col].str.decode('utf-8')

    # Filter transcripts: QV cutoff and remove controls
    transcripts = transcripts[
        (transcripts['qv'] >= qv_cutoff) &
        (~transcripts['feature_name'].str.startswith('BLANK_')) &
        (~transcripts['feature_name'].str.startswith('NegControl'))
    ]

    # Get gene list
    genes = np.sort(transcripts['feature_name'].unique())
    cell_ids = cells_df['cell_id'].values
    n_genes = len(genes)
    n_cells = len(cells_df)

    # Split by overlaps_nucleus
    nuclear_transcripts = transcripts[transcripts['overlaps_nucleus'] == 1]
    cytoplasmic_transcripts = transcripts[transcripts['overlaps_nucleus'] == 0]

    # Build three matrices
    matrix_whole = build_count_matrix(transcripts, cell_ids, genes)
    matrix_nuclear = build_count_matrix(nuclear_transcripts, cell_ids, genes)
    matrix_cytoplasmic = build_count_matrix(cytoplasmic_transcripts, cell_ids, genes)

    # Create or update AnnData
    if adata is None:
        adata = ad.AnnData(
            X=matrix_whole.T,  # cells x genes
            obs=pd.DataFrame(index=[f"cell_{i}" for i in range(n_cells)]),
            var=pd.DataFrame(index=genes)
        )
    else:
        # Validate dimensions match
        if adata.n_obs != n_cells:
            raise ValueError(f"Cell count mismatch: adata has {adata.n_obs}, counts has {n_cells}")
        adata.X = matrix_whole.T
        adata.var = pd.DataFrame(index=genes)

    # Add layers
    adata.layers['whole'] = adata.X.copy()
    adata.layers['nuclear'] = matrix_nuclear.T
    adata.layers['cytoplasmic'] = matrix_cytoplasmic.T

    # Add spatial coordinates
    adata.obs['x_centroid'] = cells_df['x_centroid'].values
    adata.obs['y_centroid'] = cells_df['y_centroid'].values
    adata.obsm['spatial'] = np.column_stack([
        adata.obs['x_centroid'].values,
        adata.obs['y_centroid'].values
    ])

    # Add cell metadata from cells.parquet
    for col in ['cell_area', 'nucleus_area', 'transcript_counts',
                'control_probe_counts', 'total_counts']:
        if col in cells_df.columns:
            adata.obs[col] = cells_df[col].values

    # Add sample metadata
    if metadata:
        for key in ['xenium_id', 'patient_id', 'position_id']:
            if key in metadata:
                adata.obs[key] = metadata[key]

        # Add demographics
        demo = metadata.get('demographics', {})
        for key in ['age', 'menopausal_status', 'brca_mutation_status']:
            if key in demo:
                adata.obs[key] = demo[key]

    # Compute per-layer QC metrics
    adata.obs['nCount_Whole'] = np.array(adata.layers['whole'].sum(axis=1)).flatten()
    adata.obs['nCount_Nuclear'] = np.array(adata.layers['nuclear'].sum(axis=1)).flatten()
    adata.obs['nCount_Cytoplasmic'] = np.array(adata.layers['cytoplasmic'].sum(axis=1)).flatten()

    adata.obs['nFeature_Whole'] = np.array((adata.layers['whole'] > 0).sum(axis=1)).flatten()
    adata.obs['nFeature_Nuclear'] = np.array((adata.layers['nuclear'] > 0).sum(axis=1)).flatten()
    adata.obs['nFeature_Cytoplasmic'] = np.array((adata.layers['cytoplasmic'] > 0).sum(axis=1)).flatten()

    # Derived metrics
    adata.obs['nuclear_fraction'] = adata.obs['nCount_Nuclear'] / adata.obs['nCount_Whole'].clip(lower=1)
    adata.obs['cytoplasmic_fraction'] = adata.obs['nCount_Cytoplasmic'] / adata.obs['nCount_Whole'].clip(lower=1)
    adata.obs['transcripts_per_gene'] = adata.obs['nCount_Whole'] / adata.obs['nFeature_Whole'].clip(lower=1)

    # Nucleus to cell area ratio
    if 'nucleus_area' in adata.obs.columns and 'cell_area' in adata.obs.columns:
        adata.obs['nucleus_ratio'] = adata.obs['nucleus_area'] / adata.obs['cell_area'].clip(lower=1)

    # QC filter: (transcripts/gene > 1) OR (total_counts > 10)
    adata.obs['qc_pass'] = (adata.obs['transcripts_per_gene'] > 1) | (adata.obs['nCount_Whole'] > 10)

    # Store transcript counts in uns for reference
    adata.uns['transcript_counts'] = {
        'total': int(len(transcripts)),
        'nuclear': int(len(nuclear_transcripts)),
        'cytoplasmic': int(len(cytoplasmic_transcripts)),
        'nuclear_fraction': float(len(nuclear_transcripts) / len(transcripts)) if len(transcripts) > 0 else 0.0,
    }
    adata.uns['qv_cutoff'] = qv_cutoff

    return adata


# CLI for standalone use
if __name__ == '__main__':
    import argparse
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent))
    from manifest_utils import get_sample_by_id, get_sample_by_index

    parser = argparse.ArgumentParser(description="Load Xenium counts")
    parser.add_argument('--id', type=str, help='Xenium ID')
    parser.add_argument('--index', '-i', type=int, help='Sample index')
    parser.add_argument('--path', '-p', type=str, help='Direct path to sample')
    parser.add_argument('--output', '-o', required=True, help='Output h5ad path')
    parser.add_argument('--qv-cutoff', type=int, default=20)

    args = parser.parse_args()

    if args.path:
        sample_path = Path(args.path)
        metadata = None
    elif args.id:
        sample = get_sample_by_id(args.id)
        if not sample:
            print(f"Sample not found: {args.id}")
            sys.exit(1)
        sample_path = Path(sample['path'])
        metadata = sample
    elif args.index is not None:
        sample = get_sample_by_index(args.index)
        sample_path = Path(sample['path'])
        metadata = sample
    else:
        parser.print_help()
        sys.exit(1)

    print(f"Loading counts from {sample_path}")
    adata = load_counts(sample_path, qv_cutoff=args.qv_cutoff, metadata=metadata)

    print(f"Saving to {args.output}")
    adata.write(args.output)

    print(f"Done: {adata.n_obs} cells, {adata.n_vars} genes")
    print(f"Layers: {list(adata.layers.keys())}")
