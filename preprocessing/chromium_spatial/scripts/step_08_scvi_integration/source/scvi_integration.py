"""
Step 08 scVI Integration - Pure Algorithm Functions

This module implements batch correction and integration using scVI
(Single-Cell Variational Inference). All functions are pure computation
with ZERO file I/O operations (except manifest update function).

Functions (Proven - Copied from archived functional implementation):
    create_gene_union: Standardize AnnData objects to union gene space
    create_vf_subset: Extract variable feature subset for scVI training
    train_scvi_model: Train scVI model with patient-level batch correction
    generate_embeddings_and_clusters: Generate latent, UMAP, clusters
    prepare_integrated_adata: Assemble complete AnnData with integration results
    prepare_seurat_export_data: Extract data for R Seurat object conversion

Functions (Novel - Implemented for GAP-003):
    update_preprocessing_manifest_with_integration: Atomic manifest update

Author: Pattern-Oriented Development Cycle - Phase 4A
Date: 2025-11-07
Version: 2.0 (Regenerated from archived functional v1.0)
Provenance:
    - Proven functions: .scratch/08/archive/scripts/scvi_integration.py
    - Novel function: Implemented per INTERFACE_CONTRACT.yaml INT-NOVEL-001
"""

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
from typing import List, Dict, Tuple
from scipy import sparse
from pathlib import Path
from datetime import datetime
import logging
import yaml
import shutil

logger = logging.getLogger(__name__)


# =============================================================================
# PROVEN FUNCTIONS - COPIED VERBATIM FROM ARCHIVE
# =============================================================================

def create_gene_union(adatas: List[ad.AnnData]) -> List[ad.AnnData]:
    """
    Standardize all patient AnnData objects to union gene space.

    Takes per-patient AnnData objects (each with potentially different genes)
    and returns standardized objects where all have the same gene space (union).
    Missing genes for a patient are filled with zeros.

    Args:
        adatas: List of AnnData objects, one per patient

    Returns:
        List of standardized AnnData objects with union gene space

    Raises:
        ValueError: If input list is empty or contains invalid objects

    Notes:
        - Union gene space typically ~20K genes
        - Missing genes filled with zeros (not NaN)
        - Gene order standardized (sorted alphabetically)
        - Based on hbca_analysis lines 198-237
        - Preserves all obs (cell metadata) and obsm (embeddings if present)
        - Var columns merged (takes first non-null value)

    Example:
        >>> pat1 = ad.AnnData(np.random.rand(100, 15000))
        >>> pat2 = ad.AnnData(np.random.rand(200, 18000))
        >>> pat1.var.index = [f'GENE{i}' for i in range(15000)]
        >>> pat2.var.index = [f'GENE{i}' for i in range(18000)]
        >>> standardized = create_gene_union([pat1, pat2])
        >>> print(standardized[0].shape, standardized[1].shape)
        (100, 18000) (200, 18000)  # Both now have same gene count
    """
    if len(adatas) == 0:
        raise ValueError("Input list is empty")

    if not all(isinstance(adata, ad.AnnData) for adata in adatas):
        raise ValueError("All inputs must be AnnData objects")

    logger.info(f"Creating gene union from {len(adatas)} AnnData objects")

    # Collect all unique genes across all objects
    all_genes = set()
    for adata in adatas:
        all_genes.update(adata.var.index.tolist())

    # Sort genes alphabetically for consistency
    union_genes = sorted(list(all_genes))
    logger.info(f"Union gene space: {len(union_genes)} genes")

    # Standardize each AnnData to union gene space
    standardized_adatas = []

    for i, adata in enumerate(adatas):
        logger.info(f"Processing object {i+1}/{len(adatas)}: {adata.n_obs} cells × {adata.n_vars} genes")

        # Identify missing genes
        current_genes = set(adata.var.index)
        missing_genes = [g for g in union_genes if g not in current_genes]

        if len(missing_genes) > 0:
            logger.info(f"  Adding {len(missing_genes)} missing genes (filled with zeros)")

            # Create zero matrix for missing genes
            if sparse.issparse(adata.X):
                # Preserve sparsity
                missing_data = sparse.csr_matrix((adata.n_obs, len(missing_genes)), dtype=adata.X.dtype)
            else:
                missing_data = np.zeros((adata.n_obs, len(missing_genes)), dtype=adata.X.dtype)

            # Create AnnData for missing genes
            missing_var = pd.DataFrame(index=missing_genes)
            missing_adata = ad.AnnData(X=missing_data, obs=adata.obs, var=missing_var)

            # Concatenate existing + missing genes
            adata_expanded = ad.concat([adata, missing_adata], axis=1, merge='first')
        else:
            adata_expanded = adata.copy()

        # Reorder to match union gene order
        adata_standardized = adata_expanded[:, union_genes].copy()

        # Validate
        assert adata_standardized.n_vars == len(union_genes), \
            f"Gene count mismatch: {adata_standardized.n_vars} != {len(union_genes)}"
        assert all(adata_standardized.var.index == union_genes), \
            "Gene order doesn't match union"

        standardized_adatas.append(adata_standardized)
        logger.info(f"  Standardized to {adata_standardized.n_obs} cells × {adata_standardized.n_vars} genes")

    logger.info(f"Successfully standardized {len(adatas)} objects to union gene space")
    return standardized_adatas


def create_vf_subset(adata_union: ad.AnnData, vf_genes: List[str]) -> ad.AnnData:
    """
    Extract variable feature subset for scVI training.

    Args:
        adata_union: AnnData with full gene space
        vf_genes: List of variable feature gene names (~13K)

    Returns:
        AnnData subset containing only VF genes

    Raises:
        ValueError: If VF gene list is empty or no VF genes found in adata

    Notes:
        - VF genes must exist in adata_union.var.index
        - Preserves all obs (cell metadata)
        - Only subsets var (genes)
        - Missing VF genes logged as warning (not error)
        - Based on hbca_analysis lines 239-250
        - Returns view (not copy) for memory efficiency

    Example:
        >>> adata_full = ad.AnnData(np.random.rand(1000, 20000))
        >>> adata_full.var.index = [f'GENE{i}' for i in range(20000)]
        >>> vf_genes = [f'GENE{i}' for i in range(0, 13000)]
        >>> adata_vf = create_vf_subset(adata_full, vf_genes)
        >>> print(adata_vf.shape)
        (1000, 13000)  # Same cells, subset genes
    """
    if len(vf_genes) == 0:
        raise ValueError("VF gene list is empty")

    if not isinstance(adata_union, ad.AnnData):
        raise ValueError("Input must be AnnData object")

    logger.info(f"Creating VF subset: {len(vf_genes)} VF genes from {adata_union.n_vars} total genes")

    # Find which VF genes are present in adata
    available_genes = set(adata_union.var.index)
    present_vf_genes = [g for g in vf_genes if g in available_genes]
    missing_vf_genes = [g for g in vf_genes if g not in available_genes]

    if len(missing_vf_genes) > 0:
        logger.warning(f"  {len(missing_vf_genes)}/{len(vf_genes)} VF genes not found in adata")
        logger.warning(f"  First 10 missing: {missing_vf_genes[:10]}")

    if len(present_vf_genes) == 0:
        raise ValueError("No VF genes found in adata")

    logger.info(f"  Found {len(present_vf_genes)}/{len(vf_genes)} VF genes ({100*len(present_vf_genes)/len(vf_genes):.1f}%)")

    # Subset to VF genes
    adata_vf = adata_union[:, present_vf_genes].copy()

    # Mark genes as variable features
    adata_vf.var['highly_variable'] = True

    logger.info(f"VF subset created: {adata_vf.n_obs} cells × {adata_vf.n_vars} VF genes")

    return adata_vf


def train_scvi_model(adata_vf: ad.AnnData, config: dict) -> scvi.model.SCVI:
    """
    Train scVI model with patient-level batch correction.

    Args:
        adata_vf: AnnData with VF genes only (~13K genes)
        config: Dict with parameters:
            - n_latent: 50 (latent dimension)
            - n_layers: 3 (VAE layers)
            - dropout_rate: 0.1
            - max_epochs: 300
            - train_size: 0.9 (validation split)
            - early_stopping: true
            - early_stopping_patience: 15
            - batch_key: "patient_id"

    Returns:
        Trained scVI model (scvi.model.SCVI object)

    Raises:
        ValueError: If required config parameters missing or invalid
        KeyError: If batch_key not found in adata.obs

    Notes:
        - Requires adata_vf.obs['patient_id'] column
        - GPU-accelerated if available
        - Saves training history in model.history
        - Based on hbca_analysis lines 453-457 (with improvements)
        - Uses validation split for early stopping
        - Trains on raw counts (not normalized)

    Example:
        >>> config = {
        ...     'n_latent': 50,
        ...     'batch_key': 'patient_id',
        ...     'max_epochs': 300,
        ...     'train_size': 0.9,
        ...     'early_stopping': True,
        ...     'early_stopping_patience': 15
        ... }
        >>> vae = train_scvi_model(adata_vf, config)
        >>> print(f"Trained for {len(vae.history['elbo_train'])} epochs")
    """
    # Validate config
    required_params = ['n_latent', 'batch_key', 'max_epochs']
    for param in required_params:
        if param not in config:
            raise ValueError(f"Missing required config parameter: {param}")

    batch_key = config['batch_key']
    if batch_key not in adata_vf.obs.columns:
        raise KeyError(f"Batch key '{batch_key}' not found in adata.obs columns")

    logger.info("Training scVI model")
    logger.info(f"  Data: {adata_vf.n_obs} cells × {adata_vf.n_vars} genes")
    logger.info(f"  Batch key: {batch_key} ({adata_vf.obs[batch_key].nunique()} batches)")
    logger.info(f"  Latent dimension: {config.get('n_latent', 50)}")
    logger.info(f"  Max epochs: {config.get('max_epochs', 300)}")

    # Setup scVI model
    scvi.model.SCVI.setup_anndata(
        adata_vf,
        batch_key=batch_key,
        layer=None  # Use .X (raw counts)
    )

    # Create model
    model = scvi.model.SCVI(
        adata_vf,
        n_latent=config.get('n_latent', 50),
        n_layers=config.get('n_layers', 3),
        dropout_rate=config.get('dropout_rate', 0.1),
        gene_likelihood='nb'  # Negative binomial for count data
    )

    logger.info("Model created, beginning training...")

    # Train model
    train_kwargs = {
        'max_epochs': config.get('max_epochs', 300),
        'train_size': config.get('train_size', 0.9),
        'batch_size': config.get('batch_size', 128)
    }

    # Add early stopping if enabled
    if config.get('early_stopping', True):
        train_kwargs['early_stopping'] = True
        train_kwargs['early_stopping_patience'] = config.get('early_stopping_patience', 15)
        logger.info(f"  Early stopping enabled (patience={train_kwargs['early_stopping_patience']})")

    # Check for GPU
    use_gpu = config.get('use_gpu', True)
    if use_gpu:
        import torch
        if torch.cuda.is_available():
            logger.info("  Using GPU acceleration")
            train_kwargs['accelerator'] = 'gpu'
        else:
            logger.warning("  GPU requested but not available, using CPU")
            train_kwargs['accelerator'] = 'cpu'
    else:
        train_kwargs['accelerator'] = 'cpu'

    # Train
    model.train(**train_kwargs)

    # Log training summary (avoid accessing history to prevent Series formatting errors)
    logger.info("Training complete")

    return model


def generate_embeddings_and_clusters(
    vae: scvi.model.SCVI,
    adata_full: ad.AnnData,
    config: dict
) -> dict:
    """
    Generate latent embeddings, UMAP, and multi-resolution clusters.

    Args:
        vae: Trained scVI model
        adata_full: AnnData with full gene space (for embedding generation)
        config: Dict with parameters:
            - n_neighbors: 30
            - umap_min_dist: 0.3
            - clustering_resolutions: [0.3, 0.5, 0.8, 1.0, 5.0]

    Returns:
        Dict with:
            'latent': np.ndarray (cells × 50)
            'umap': np.ndarray (cells × 2)
            'clusters': dict of {resolution: np.ndarray}
                e.g., {'0.3': array([0, 1, 0, ...]), '0.5': array([...]), ...}

    Raises:
        ValueError: If config parameters invalid or missing

    Notes:
        - Latent embeddings: 50D scVI latent space
        - UMAP: Computed on latent space (not raw counts)
        - Clusters: Leiden algorithm on scVI neighbors
        - All resolutions computed in single pass
        - Cluster labels are strings (for consistency with R)
        - Based on scanpy workflow

    Example:
        >>> results = generate_embeddings_and_clusters(vae, adata, config)
        >>> print(results['latent'].shape)  # (276426, 50)
        >>> print(results['umap'].shape)    # (276426, 2)
        >>> print(results['clusters'].keys())  # dict_keys(['0.3', '0.5', '0.8', '1.0', '5.0'])
    """
    logger.info("Generating embeddings and clusters")

    # Extract latent embeddings
    logger.info("  Extracting scVI latent embeddings...")
    latent = vae.get_latent_representation()
    logger.info(f"  Latent embeddings: {latent.shape}")

    # Create temporary AnnData for clustering (using latent space)
    adata_latent = ad.AnnData(X=latent, obs=adata_full.obs.copy())

    # Compute neighbors on latent space
    n_neighbors = config.get('n_neighbors', 30)
    logger.info(f"  Computing neighbors (n_neighbors={n_neighbors})...")
    sc.pp.neighbors(adata_latent, n_neighbors=n_neighbors, use_rep='X')

    # Compute UMAP on latent space
    logger.info("  Computing UMAP...")
    umap_min_dist = config.get('umap_min_dist', 0.3)
    sc.tl.umap(adata_latent, min_dist=umap_min_dist)
    umap_coords = adata_latent.obsm['X_umap']
    logger.info(f"  UMAP coordinates: {umap_coords.shape}")

    # Multi-resolution clustering
    resolutions = config.get('clustering_resolutions', [0.3, 0.5, 0.8, 1.0, 5.0])
    logger.info(f"  Computing Leiden clustering at {len(resolutions)} resolutions...")

    clusters = {}
    for res in resolutions:
        logger.info(f"    Resolution {res}...")
        sc.tl.leiden(adata_latent, resolution=res, key_added=f'leiden_{res}')

        # Get cluster labels as array
        cluster_labels = adata_latent.obs[f'leiden_{res}'].values

        # Convert to numeric strings (consistent with R factor levels)
        clusters[str(res)] = cluster_labels.astype(str)

        n_clusters = len(np.unique(cluster_labels))
        logger.info(f"      {n_clusters} clusters identified")

    logger.info(f"Generated embeddings and {len(clusters)} cluster resolutions")

    return {
        'latent': latent,
        'umap': umap_coords,
        'clusters': clusters
    }


def prepare_integrated_adata(
    adata_full: ad.AnnData,
    latent_embeddings: np.ndarray,
    umap_coords: np.ndarray,
    clusters: dict,
    config: dict
) -> ad.AnnData:
    """
    Assemble complete AnnData object with all integration results.

    Args:
        adata_full: Original AnnData with raw counts and metadata
        latent_embeddings: scVI 50D embeddings (cells × 50)
        umap_coords: UMAP 2D coordinates (cells × 2)
        clusters: Dict of cluster assignments per resolution
        config: Dict (for any additional metadata)

    Returns:
        Complete AnnData with:
            .X: Raw counts (sparse)
            .obs: Original metadata + cluster columns
            .obsm['X_scvi']: Latent embeddings
            .obsm['X_umap']: UMAP coordinates
            .var: Gene metadata
            .uns: Integration metadata

    Raises:
        ValueError: If dimensions don't match

    Notes:
        - Preserves original counts in .X
        - Adds 'leiden_scvi_0.3', 'leiden_scvi_0.5', etc. to .obs
        - Embeddings stored in .obsm for scanpy compatibility
        - Ready to save as .h5ad
        - Gene names in .var.index

    Example:
        >>> adata_integrated = prepare_integrated_adata(
        ...     adata_full, latent, umap, clusters, config
        ... )
        >>> print(adata_integrated.obsm.keys())
        dict_keys(['X_scvi', 'X_umap'])
        >>> print('leiden_scvi_0.3' in adata_integrated.obs.columns)
        True
    """
    logger.info("Assembling integrated AnnData object")

    # Validate dimensions
    n_cells = adata_full.n_obs
    if latent_embeddings.shape[0] != n_cells:
        raise ValueError(f"Latent embeddings ({latent_embeddings.shape[0]}) doesn't match cells ({n_cells})")
    if umap_coords.shape[0] != n_cells:
        raise ValueError(f"UMAP coords ({umap_coords.shape[0]}) doesn't match cells ({n_cells})")

    logger.info(f"  Cells: {n_cells}")
    logger.info(f"  Genes: {adata_full.n_vars}")
    logger.info(f"  Latent dims: {latent_embeddings.shape[1]}")
    logger.info(f"  Cluster resolutions: {len(clusters)}")

    # Create copy of adata to preserve original
    adata_integrated = adata_full.copy()

    # Add embeddings to obsm
    adata_integrated.obsm['X_scvi'] = latent_embeddings
    adata_integrated.obsm['X_umap'] = umap_coords
    logger.info("  Added embeddings to .obsm")

    # Add cluster assignments to obs
    for res, cluster_labels in clusters.items():
        if len(cluster_labels) != n_cells:
            raise ValueError(f"Cluster labels for res={res} ({len(cluster_labels)}) doesn't match cells ({n_cells})")

        col_name = f'leiden_scvi_{res}'
        adata_integrated.obs[col_name] = cluster_labels

    logger.info(f"  Added {len(clusters)} cluster resolution columns to .obs")

    # Add integration metadata to uns
    adata_integrated.uns['scvi_integration'] = {
        'n_latent': latent_embeddings.shape[1],
        'clustering_resolutions': list(clusters.keys()),
        'n_cells': n_cells,
        'n_genes': adata_full.n_vars,
        'config': config
    }
    logger.info("  Added integration metadata to .uns")

    logger.info("Integrated AnnData object assembled successfully")

    return adata_integrated


def prepare_seurat_export_data(adata_integrated: ad.AnnData) -> dict:
    """
    Extract data for R Seurat object conversion.

    Args:
        adata_integrated: Complete AnnData from prepare_integrated_adata()

    Returns:
        Dict with:
            'counts': scipy.sparse matrix or np.ndarray (genes × cells, transposed!)
            'gene_names': list of gene names
            'cell_ids': list of cell barcodes
            'scvi_embeddings': np.ndarray (cells × 50)
            'umap_coords': np.ndarray (cells × 2)
            'metadata': pd.DataFrame with all obs columns

    Raises:
        ValueError: If required embeddings not found in adata

    Notes:
        - Counts matrix TRANSPOSED for R (genes × cells)
        - Gene names must match row order of transposed matrix
        - Metadata includes patient_id, sample_id, clusters
        - R script will create Seurat @reductions from embeddings
        - Preserves sparsity in counts matrix

    Example:
        >>> export_data = prepare_seurat_export_data(adata_integrated)
        >>> print(export_data['counts'].shape)  # (20000, 276426) - genes × cells
        >>> print(export_data['scvi_embeddings'].shape)  # (276426, 50)
        >>> print(export_data['metadata'].columns)
        Index(['patient_id', 'sample_id', ..., 'leiden_scvi_0.3', ...])
    """
    logger.info("Preparing data for Seurat export")

    # Validate required embeddings present
    if 'X_scvi' not in adata_integrated.obsm:
        raise ValueError("scVI embeddings (X_scvi) not found in adata.obsm")
    if 'X_umap' not in adata_integrated.obsm:
        raise ValueError("UMAP coordinates (X_umap) not found in adata.obsm")

    # Transpose counts matrix for R (genes × cells)
    logger.info("  Transposing counts matrix for R...")
    counts_transposed = adata_integrated.X.T
    if sparse.issparse(counts_transposed):
        # Convert to CSC for efficient column operations in R
        counts_transposed = counts_transposed.tocsc()
        logger.info(f"  Counts: {counts_transposed.shape} (sparse CSC)")
    else:
        logger.info(f"  Counts: {counts_transposed.shape} (dense)")

    # Extract gene names (must match transposed matrix row order)
    gene_names = adata_integrated.var.index.tolist()
    logger.info(f"  Gene names: {len(gene_names)}")

    # Extract cell IDs
    cell_ids = adata_integrated.obs.index.tolist()
    logger.info(f"  Cell IDs: {len(cell_ids)}")

    # Extract embeddings
    scvi_embeddings = adata_integrated.obsm['X_scvi']
    umap_coords = adata_integrated.obsm['X_umap']
    logger.info(f"  scVI embeddings: {scvi_embeddings.shape}")
    logger.info(f"  UMAP coords: {umap_coords.shape}")

    # Extract metadata (all obs columns)
    metadata = adata_integrated.obs.copy()
    logger.info(f"  Metadata: {metadata.shape} ({len(metadata.columns)} columns)")

    # Validate dimensions
    assert counts_transposed.shape[0] == len(gene_names), \
        f"Gene count mismatch: {counts_transposed.shape[0]} != {len(gene_names)}"
    assert counts_transposed.shape[1] == len(cell_ids), \
        f"Cell count mismatch: {counts_transposed.shape[1]} != {len(cell_ids)}"
    assert scvi_embeddings.shape[0] == len(cell_ids), \
        f"scVI embedding cell count mismatch: {scvi_embeddings.shape[0]} != {len(cell_ids)}"
    assert umap_coords.shape[0] == len(cell_ids), \
        f"UMAP cell count mismatch: {umap_coords.shape[0]} != {len(cell_ids)}"

    export_data = {
        'counts': counts_transposed,
        'gene_names': gene_names,
        'cell_ids': cell_ids,
        'scvi_embeddings': scvi_embeddings,
        'umap_coords': umap_coords,
        'metadata': metadata
    }

    logger.info("Seurat export data prepared successfully")

    return export_data


# =============================================================================
# NOVEL FUNCTION - IMPLEMENTED FOR GAP-003
# =============================================================================

def update_preprocessing_manifest_with_integration(
    manifest_path: Path,
    integration_metadata: Dict,
    output_paths: Dict,
    timestamp: str = None
) -> Dict:
    """
    Atomically update preprocessing_manifest.yaml with integrated_objects section.

    Parameters
    ----------
    manifest_path : Path
        Path to preprocessing_manifest.yaml
    integration_metadata : Dict
        Metadata with keys: patients (list), total_cells (int), total_genes (int),
        n_latent (int), batch_key (str), vf_source (str), algorithm (str)
    output_paths : Dict
        Nested dict with keys: embeddings, clusters, objects (each with file paths)
    timestamp : str, optional
        ISO format timestamp, defaults to current datetime

    Returns
    -------
    Dict
        Status dict with keys: success (bool), backup_path (str or None),
        error (str or None), summary (dict or None)

    Notes
    -----
    Implements atomic update pattern:
    1. Load existing manifest
    2. Create timestamped backup
    3. Build integrated_objects structure
    4. Write updated manifest
    5. Validate YAML parseability
    6. Restore backup on any failure

    The function creates a three-tier structure in the manifest:
    - filtering_status: unfiltered | filtered | final
    - scope: all_cells | compartment_specific
    - Integration details: metadata and output paths

    Example
    -------
    >>> metadata = {
    ...     'patients': ['Pat1', 'Pat2', 'UCI604', 'UCI220228'],
    ...     'total_cells': 276543,
    ...     'total_genes': 12847,
    ...     'n_latent': 50,
    ...     'batch_key': 'patient_id',
    ...     'vf_source': 'step_04_post_qc_vfs',
    ...     'algorithm': 'scvi-tools'
    ... }
    >>> output_paths = {
    ...     'embeddings': {
    ...         'latent': 'outputs/08_ScviIntegration/embeddings/latent_50d.csv',
    ...         'umap': 'outputs/08_ScviIntegration/embeddings/umap_2d.csv'
    ...     },
    ...     'clusters': {
    ...         'leiden_res_0_3': 'outputs/08_ScviIntegration/clusters/leiden_res_0.3.csv'
    ...     },
    ...     'objects': {
    ...         'anndata': 'outputs/08_ScviIntegration/integrated_objects/integrated.h5ad',
    ...         'scvi_model': 'outputs/08_ScviIntegration/model/scvi_model'
    ...     }
    ... }
    >>> result = update_preprocessing_manifest_with_integration(
    ...     Path('preprocessing_manifest.yaml'),
    ...     metadata,
    ...     output_paths
    ... )
    >>> print(result['success'])
    True
    """
    if timestamp is None:
        timestamp = datetime.now().isoformat()

    result = {
        'success': False,
        'backup_path': None,
        'error': None,
        'summary': None
    }

    try:
        # Step 1: Load existing manifest
        logger.info(f"Loading existing manifest: {manifest_path}")
        with open(manifest_path, 'r') as f:
            manifest = yaml.safe_load(f)

        # Step 2: Create timestamped backup
        backup_timestamp = timestamp.replace(':', '-').replace('.', '-')
        backup_path = str(manifest_path) + f".backup_{backup_timestamp}"
        logger.info(f"Creating backup: {backup_path}")
        shutil.copy2(manifest_path, backup_path)
        result['backup_path'] = backup_path

        # Step 3: Build integrated_objects structure
        logger.info("Building integrated_objects structure...")

        # Initialize integrated_objects if not present
        if 'integrated_objects' not in manifest:
            manifest['integrated_objects'] = {}

        # Initialize unfiltered tier if not present
        if 'unfiltered' not in manifest['integrated_objects']:
            manifest['integrated_objects']['unfiltered'] = {}

        # Build the complete structure for unfiltered.all_cells
        manifest['integrated_objects']['unfiltered']['all_cells'] = {
            'filtering_status': 'unfiltered',
            'scope': 'all_cells',
            'type': 'scvi-integration',
            'status': 'complete',
            'timestamp': timestamp,
            'created_by': 'step_08_scvi_integration',
            'metadata': integration_metadata,
            'outputs': output_paths
        }

        # Create placeholder sections for downstream steps
        if 'filtered' not in manifest['integrated_objects']:
            manifest['integrated_objects']['filtered'] = {
                'all_cells': {'status': None}
            }

        if 'final' not in manifest['integrated_objects']:
            manifest['integrated_objects']['final'] = {
                'epithelial': {'status': None},
                'stromal': {'status': None},
                'immune': {'status': None}
            }

        logger.info("  Added integrated_objects.unfiltered.all_cells")
        logger.info(f"  Patients: {integration_metadata.get('patients', [])}")
        logger.info(f"  Total cells: {integration_metadata.get('total_cells', 0)}")
        logger.info(f"  Total genes: {integration_metadata.get('total_genes', 0)}")

        # Step 4: Write updated manifest
        logger.info(f"Writing updated manifest: {manifest_path}")
        with open(manifest_path, 'w') as f:
            yaml.dump(manifest, f, default_flow_style=False, sort_keys=False, allow_unicode=True)

        # Step 5: Validate YAML parseability
        logger.info("Validating updated manifest...")
        with open(manifest_path, 'r') as f:
            yaml.safe_load(f)

        # Success
        logger.info("Manifest update successful")
        result['success'] = True
        result['summary'] = {
            'section_added': 'integrated_objects.unfiltered.all_cells',
            'patients': integration_metadata.get('patients', []),
            'total_cells': integration_metadata.get('total_cells', 0),
            'output_count': len(output_paths)
        }

    except Exception as e:
        # Step 6: Restore backup on failure
        error_msg = str(e)
        logger.error(f"Manifest update failed: {error_msg}")
        result['error'] = error_msg

        if result['backup_path'] and Path(result['backup_path']).exists():
            try:
                logger.warning(f"Restoring backup from: {result['backup_path']}")
                shutil.copy2(result['backup_path'], manifest_path)
                logger.info("Backup restored successfully")
            except Exception as restore_error:
                restore_msg = f"Backup restore failed: {str(restore_error)}"
                logger.error(restore_msg)
                result['error'] += f" | {restore_msg}"

    return result
