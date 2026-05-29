#!/usr/bin/env python3
"""
manifest_utils.py - Xenium data access utilities

Loads Xenium sample metadata from raw_data_manifest.yaml.
Follows patterns from project/03_Integration/.../utils/config.py

Usage:
    from shared.manifest_utils import get_xenium_samples, get_sample_path

    samples = get_xenium_samples()
    for s in samples:
        print(f"{s['xenium_id']}: {s['path']}")
"""

import os
from pathlib import Path
from typing import Dict, List, Any, Optional

import yaml


def get_project_root() -> Path:
    """Get project root from environment or derive from script location."""
    if 'PROJECT_ROOT' in os.environ:
        return Path(os.environ['PROJECT_ROOT'])
    # Default: assume we're in project/06_Xenium/shared/
    return Path(__file__).parents[3]


def load_yaml(path: Path) -> Dict[str, Any]:
    """Load YAML file."""
    with open(path) as f:
        return yaml.safe_load(f)


def load_raw_data_manifest() -> Dict[str, Any]:
    """Load raw_data_manifest.yaml."""
    manifest_path = get_project_root() / "project/config/raw_data_manifest.yaml"
    return load_yaml(manifest_path)


def get_xenium_samples(manifest: Optional[Dict] = None) -> List[Dict[str, Any]]:
    """
    Extract all Xenium samples from raw_data_manifest.

    Args:
        manifest: Pre-loaded manifest, or loads fresh if None

    Returns:
        List of dicts with keys:
            - xenium_id: Unique identifier (e.g., Pat1_P1_xenium_1)
            - patient_id: Patient (e.g., Pat1)
            - position_id: Position (e.g., P1)
            - path: Full path to Xenium data directory
            - original_folder: Original folder name from instrument
            - demographics: Patient demographics (age, menopausal_status, etc.)
    """
    if manifest is None:
        manifest = load_raw_data_manifest()

    samples = []

    for patient_id, patient in manifest.get('patients', {}).items():
        demographics = patient.get('demographics', {})

        for position_id, position in patient.get('positions', {}).items():
            anatomical = position.get('anatomical', {})

            for xenium in position.get('xenium', []):
                samples.append({
                    'xenium_id': xenium['xenium_id'],
                    'patient_id': patient_id,
                    'position_id': position_id,
                    'path': xenium['path'],
                    'original_folder': xenium.get('original_folder', ''),
                    'demographics': demographics,
                    'anatomical': anatomical,
                })

    return samples


def get_sample_by_id(xenium_id: str, manifest: Optional[Dict] = None) -> Optional[Dict[str, Any]]:
    """
    Get a specific sample by xenium_id.

    Args:
        xenium_id: Sample ID (e.g., Pat1_P1_xenium_1)
        manifest: Pre-loaded manifest

    Returns:
        Sample dict or None if not found
    """
    samples = get_xenium_samples(manifest)
    for s in samples:
        if s['xenium_id'] == xenium_id:
            return s
    return None


def get_sample_by_index(index: int, manifest: Optional[Dict] = None) -> Dict[str, Any]:
    """
    Get sample by index (for SLURM array jobs).

    Args:
        index: 0-based index
        manifest: Pre-loaded manifest

    Returns:
        Sample dict

    Raises:
        IndexError if index out of range
    """
    samples = get_xenium_samples(manifest)
    if index < 0 or index >= len(samples):
        raise IndexError(f"Sample index {index} out of range (0-{len(samples)-1})")
    return samples[index]


def validate_sample_path(sample: Dict[str, Any]) -> bool:
    """
    Check if required files exist in sample directory.

    Required files:
        - cells.parquet (or cells.csv.gz)
        - transcripts.parquet (or transcripts.csv.gz)
        - cell_feature_matrix.h5

    Returns:
        True if all required files exist
    """
    sample_path = Path(sample['path'])

    # Check cells file
    cells_ok = (
        (sample_path / "cells.parquet").exists() or
        (sample_path / "cells.csv.gz").exists()
    )

    # Check transcripts file
    transcripts_ok = (
        (sample_path / "transcripts.parquet").exists() or
        (sample_path / "transcripts.csv.gz").exists()
    )

    # Check expression matrix
    matrix_ok = (sample_path / "cell_feature_matrix.h5").exists()

    return cells_ok and transcripts_ok and matrix_ok


def get_sample_files(sample: Dict[str, Any]) -> Dict[str, Path]:
    """
    Get paths to standard Xenium files for a sample.

    Returns dict with keys:
        - cells: Path to cells file
        - transcripts: Path to transcripts file
        - matrix: Path to cell_feature_matrix.h5
        - cell_boundaries: Path to cell boundaries (if exists)
        - nucleus_boundaries: Path to nucleus boundaries (if exists)
        - morphology: Path to morphology image (if exists)
    """
    sample_path = Path(sample['path'])

    files = {}

    # Cells - prefer parquet
    if (sample_path / "cells.parquet").exists():
        files['cells'] = sample_path / "cells.parquet"
    elif (sample_path / "cells.csv.gz").exists():
        files['cells'] = sample_path / "cells.csv.gz"

    # Transcripts - prefer parquet
    if (sample_path / "transcripts.parquet").exists():
        files['transcripts'] = sample_path / "transcripts.parquet"
    elif (sample_path / "transcripts.csv.gz").exists():
        files['transcripts'] = sample_path / "transcripts.csv.gz"

    # Matrix
    if (sample_path / "cell_feature_matrix.h5").exists():
        files['matrix'] = sample_path / "cell_feature_matrix.h5"

    # Optional files
    for name, filename in [
        ('cell_boundaries', 'cell_boundaries.parquet'),
        ('nucleus_boundaries', 'nucleus_boundaries.parquet'),
        ('morphology', 'morphology.ome.tif'),
    ]:
        path = sample_path / filename
        if path.exists():
            files[name] = path

    return files


def print_manifest_summary(manifest: Optional[Dict] = None):
    """Print summary of Xenium samples in manifest."""
    samples = get_xenium_samples(manifest)

    print(f"Total Xenium samples: {len(samples)}")
    print()

    # Group by patient
    by_patient = {}
    for s in samples:
        pid = s['patient_id']
        if pid not in by_patient:
            by_patient[pid] = []
        by_patient[pid].append(s)

    for patient_id in sorted(by_patient.keys()):
        patient_samples = by_patient[patient_id]
        demo = patient_samples[0]['demographics']
        print(f"{patient_id}: {len(patient_samples)} samples")
        print(f"  Age: {demo.get('age', 'N/A')}, Menopausal: {demo.get('menopausal_status', 'N/A')}")

        # List positions
        positions = sorted(set(s['position_id'] for s in patient_samples))
        print(f"  Positions: {', '.join(positions)}")
        print()


# CLI for testing
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Xenium manifest utilities")
    parser.add_argument('--list', '-l', action='store_true', help='List all samples')
    parser.add_argument('--validate', '-v', action='store_true', help='Validate sample paths')
    parser.add_argument('--index', '-i', type=int, help='Get sample by index')
    parser.add_argument('--id', type=str, help='Get sample by xenium_id')

    args = parser.parse_args()

    if args.index is not None:
        sample = get_sample_by_index(args.index)
        print(f"Sample [{args.index}]:")
        for k, v in sample.items():
            if k != 'demographics':
                print(f"  {k}: {v}")

    elif args.id:
        sample = get_sample_by_id(args.id)
        if sample:
            print(f"Sample {args.id}:")
            for k, v in sample.items():
                if k != 'demographics':
                    print(f"  {k}: {v}")
        else:
            print(f"Sample not found: {args.id}")

    elif args.validate:
        samples = get_xenium_samples()
        valid = 0
        invalid = []
        for s in samples:
            if validate_sample_path(s):
                valid += 1
            else:
                invalid.append(s['xenium_id'])

        print(f"Valid: {valid}/{len(samples)}")
        if invalid:
            print(f"Invalid samples: {', '.join(invalid)}")

    else:
        print_manifest_summary()
