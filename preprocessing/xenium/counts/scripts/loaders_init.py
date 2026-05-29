"""
Modular Xenium data loaders.

Each loader is independent and can be called separately.
Loaders add data to existing AnnData objects (except load_images which returns separate ImageContainer).

Usage:
    from shared.loaders import load_counts, load_molecules, load_boundaries, load_images

    # Start with counts (creates adata)
    adata = load_counts(sample_path)

    # Add molecules to existing adata
    adata = load_molecules(adata, sample_path)

    # Add boundaries to existing adata
    adata = load_boundaries(adata, sample_path)

    # Images return separate container (squidpy pattern)
    img = load_images(sample_path)
"""

from .load_counts import load_counts
from .load_molecules import load_molecules
from .load_boundaries import load_boundaries
from .load_images import load_images, load_images_lazy

__all__ = ['load_counts', 'load_molecules', 'load_boundaries', 'load_images', 'load_images_lazy']
