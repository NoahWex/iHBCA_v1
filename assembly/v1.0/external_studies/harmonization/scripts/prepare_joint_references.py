#!/usr/bin/env python3
"""
prepare_joint_references.py
============================
Convert Austin's author_share files to parquet format for efficient access.

Usage:
    python prepare_joint_references.py

Inputs (from author_share/):
    - X_scVI100.csv (2.7GB) - Joint 100-dimensional scVI embedding
    - ihbca_level1.5_annotations.csv (850MB) - Cell annotations

Outputs (to harmonization/outputs/joint_objects/):
    - X_scVI100_full.parquet - Joint embedding in parquet format
    - cell_annotations.parquet - Cell annotations in parquet format
    - joint_references.yaml - Manifest with paths and counts npz reference

Note: The counts matrix (8.5GB npz) is NOT copied - it's referenced in place.
"""

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import yaml
from pathlib import Path
from datetime import datetime

# =============================================================================
# Configuration
# =============================================================================

BASE_PATH = Path("${SOURCE_IHBCAV1_EXTERNAL_STUDIES%/}")
AUTHOR_SHARE = Path("${SOURCE_AUTHOR_SHARE}")
OUTPUT_DIR = BASE_PATH / "harmonization" / "outputs" / "joint_objects"

# Source files
JOINT_EMBEDDING_CSV = AUTHOR_SHARE / "X_scVI100.csv"
CELL_ANNOTATIONS_CSV = AUTHOR_SHARE / "ihbca_level1.5_annotations.csv"
COUNTS_NPZ = AUTHOR_SHARE / "preintegration_HBCA_inner_ENSEMBL_simple_counts_matrix.npz"

# =============================================================================
# Functions
# =============================================================================

def convert_embedding_to_parquet():
    """Convert X_scVI100.csv to parquet format."""
    print(f"\n[1/3] Converting joint embedding to parquet...")
    print(f"  Source: {JOINT_EMBEDDING_CSV}")

    if not JOINT_EMBEDDING_CSV.exists():
        raise FileNotFoundError(f"Joint embedding not found: {JOINT_EMBEDDING_CSV}")

    # Read CSV (first column is cell ID as index)
    print("  Reading CSV (this may take a while for 2.7GB)...")
    df = pd.read_csv(JOINT_EMBEDDING_CSV, index_col=0)

    print(f"  Shape: {df.shape[0]:,} cells × {df.shape[1]} dimensions")

    # Reset index to make cell_id a column
    df = df.reset_index()
    df.columns = ['cell_id'] + [f'scVI_{i}' for i in range(100)]

    # Write to parquet
    output_file = OUTPUT_DIR / "X_scVI100_full.parquet"
    print(f"  Writing parquet: {output_file.name}")
    df.to_parquet(output_file, index=False, compression='snappy')

    # Get file size
    size_mb = output_file.stat().st_size / (1024 * 1024)
    print(f"  Output size: {size_mb:.1f} MB (vs ~2700 MB CSV)")

    return {
        'file': str(output_file),
        'cells': df.shape[0],
        'dimensions': 100,
        'size_mb': round(size_mb, 1)
    }


def convert_annotations_to_parquet():
    """Convert ihbca_level1.5_annotations.csv to parquet format."""
    print(f"\n[2/3] Converting cell annotations to parquet...")
    print(f"  Source: {CELL_ANNOTATIONS_CSV}")

    if not CELL_ANNOTATIONS_CSV.exists():
        raise FileNotFoundError(f"Cell annotations not found: {CELL_ANNOTATIONS_CSV}")

    # Read CSV
    print("  Reading CSV (this may take a while for 850MB)...")
    df = pd.read_csv(CELL_ANNOTATIONS_CSV, low_memory=False)

    print(f"  Shape: {df.shape[0]:,} cells × {df.shape[1]} columns")

    # Write to parquet
    output_file = OUTPUT_DIR / "cell_annotations.parquet"
    print(f"  Writing parquet: {output_file.name}")
    df.to_parquet(output_file, index=False, compression='snappy')

    # Get file size
    size_mb = output_file.stat().st_size / (1024 * 1024)
    print(f"  Output size: {size_mb:.1f} MB (vs ~850 MB CSV)")

    return {
        'file': str(output_file),
        'cells': df.shape[0],
        'columns': df.shape[1],
        'size_mb': round(size_mb, 1)
    }


def create_manifest(embedding_info, annotations_info):
    """Create joint_references.yaml manifest."""
    print(f"\n[3/3] Creating manifest...")

    manifest = {
        'version': '1.0',
        'created': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'description': 'Joint embedding and annotation references for Phase B',

        'joint_embedding': {
            'file': 'X_scVI100_full.parquet',
            'source': str(JOINT_EMBEDDING_CSV),
            'cells': embedding_info['cells'],
            'dimensions': embedding_info['dimensions'],
            'size_mb': embedding_info['size_mb'],
            'format': 'parquet (snappy compression)',
            'notes': 'scVI latent space (100 dimensions) from Reed et al. iHBCA integration'
        },

        'cell_annotations': {
            'file': 'cell_annotations.parquet',
            'source': str(CELL_ANNOTATIONS_CSV),
            'cells': annotations_info['cells'],
            'columns': annotations_info['columns'],
            'size_mb': annotations_info['size_mb'],
            'format': 'parquet (snappy compression)',
            'notes': 'Level 1.5 cell type annotations from iHBCA'
        },

        'counts_matrix': {
            'file': str(COUNTS_NPZ),
            'format': 'scipy sparse npz',
            'size_gb': round(COUNTS_NPZ.stat().st_size / (1024**3), 1) if COUNTS_NPZ.exists() else None,
            'notes': 'Pre-integration counts (ENSEMBL IDs). Reference in place, not copied.'
        }
    }

    output_file = OUTPUT_DIR / "joint_references.yaml"
    with open(output_file, 'w') as f:
        yaml.dump(manifest, f, default_flow_style=False, sort_keys=False)

    print(f"  Saved: {output_file.name}")

    return manifest


def main():
    print("=" * 60)
    print("Prepare Joint References")
    print("=" * 60)
    print(f"Output directory: {OUTPUT_DIR}")

    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Convert files
    embedding_info = convert_embedding_to_parquet()
    annotations_info = convert_annotations_to_parquet()

    # Create manifest
    manifest = create_manifest(embedding_info, annotations_info)

    print("\n" + "=" * 60)
    print("SUCCESS: Joint references prepared")
    print(f"  Embedding: {embedding_info['cells']:,} cells × {embedding_info['dimensions']} dims")
    print(f"  Annotations: {annotations_info['cells']:,} cells × {annotations_info['columns']} cols")
    print(f"  Output: {OUTPUT_DIR}")
    print("=" * 60)


if __name__ == "__main__":
    main()
