"""
Step 12 Filtered scVI Integration - Pure Algorithm Functions

This module reuses Step 08 scVI integration algorithms.
The filtering difference (FOUR-WAY instead of THREE-WAY) is handled
in the wrapper, not in the pure algorithms.

Imports all functions from Step 08 source module.

Author: Pattern-Oriented Development Cycle
Date: 2025-12-01
Version: 1.0
Provenance: Step 08 scVI Integration algorithms
"""

# Import all functions from Step 08 source
import sys
from pathlib import Path

# Add Step 08 source to path for imports
STEP_08_SOURCE = Path(__file__).parent.parent.parent / "step_08_scvi_integration" / "source"
sys.path.insert(0, str(STEP_08_SOURCE))

from scvi_integration import (
    create_gene_union,
    create_vf_subset,
    train_scvi_model,
    generate_embeddings_and_clusters,
    prepare_integrated_adata,
    prepare_seurat_export_data,
    update_preprocessing_manifest_with_integration
)

__all__ = [
    'create_gene_union',
    'create_vf_subset',
    'train_scvi_model',
    'generate_embeddings_and_clusters',
    'prepare_integrated_adata',
    'prepare_seurat_export_data',
    'update_preprocessing_manifest_with_integration'
]
