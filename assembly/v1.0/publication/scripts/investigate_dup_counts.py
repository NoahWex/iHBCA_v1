#!/usr/bin/env python3
"""Investigate HCA validator 'duplicated raw counts' in integrated objects.

The HCA validator reports N cells with duplicated raw counts in adata.X.
This script determines:
1. What the validator actually checks (duplicate rows in X)
2. Which cells are affected
3. Which studies they come from
4. Whether they're exact duplicates or near-duplicates
"""

import argparse
import time
import sys
import numpy as np
import scipy.sparse as sp
import anndata as ad
from collections import Counter


def find_duplicate_rows_sparse(X, batch_size=10000):
    """Find rows in sparse matrix X that are exact duplicates of other rows.

    Uses row hashing for O(n) detection instead of O(n^2) pairwise comparison.
    Returns dict mapping canonical row index -> list of duplicate row indices.
    """
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    elif not sp.isspmatrix_csr(X):
        X = X.tocsr()

    n_rows = X.shape[0]
    print(f"  Scanning {n_rows:,} rows for duplicate count vectors...")

    # Hash each row by its data content
    # For sparse CSR: each row is defined by indices[indptr[i]:indptr[i+1]]
    # and data[indptr[i]:indptr[i+1]]
    row_hashes = {}  # hash -> list of row indices

    t0 = time.time()
    for i in range(n_rows):
        start, end = X.indptr[i], X.indptr[i + 1]
        indices = X.indices[start:end]
        data = X.data[start:end]
        # Create hashable representation
        row_key = (tuple(indices), tuple(data))
        h = hash(row_key)

        if h not in row_hashes:
            row_hashes[h] = [i]
        else:
            row_hashes[h].append(i)

        if (i + 1) % 500000 == 0:
            elapsed = time.time() - t0
            print(f"    Processed {i + 1:,}/{n_rows:,} rows ({elapsed:.1f}s)")

    elapsed = time.time() - t0
    print(f"  Hashing complete in {elapsed:.1f}s")

    # Filter to groups with duplicates
    dup_groups = {k: v for k, v in row_hashes.items() if len(v) > 1}

    # Verify hash collisions by exact comparison
    verified_groups = []
    false_collisions = 0
    for h, indices_list in dup_groups.items():
        if len(indices_list) == 2:
            i, j = indices_list
            row_i = X[i]
            row_j = X[j]
            if (row_i - row_j).nnz == 0:
                verified_groups.append(indices_list)
            else:
                false_collisions += 1
        else:
            # For groups > 2, do pairwise verification
            # Group by actual content
            sub_groups = []
            used = set()
            for idx_a in range(len(indices_list)):
                if idx_a in used:
                    continue
                group = [indices_list[idx_a]]
                for idx_b in range(idx_a + 1, len(indices_list)):
                    if idx_b in used:
                        continue
                    row_a = X[indices_list[idx_a]]
                    row_b = X[indices_list[idx_b]]
                    if (row_a - row_b).nnz == 0:
                        group.append(indices_list[idx_b])
                        used.add(idx_b)
                if len(group) > 1:
                    sub_groups.append(group)
                    used.add(idx_a)
            verified_groups.extend(sub_groups)

    if false_collisions > 0:
        print(f"  Hash collisions (not true duplicates): {false_collisions}")

    return verified_groups


def main():
    parser = argparse.ArgumentParser(description="Investigate duplicate count vectors")
    parser.add_argument("h5ad", help="Path to h5ad file")
    parser.add_argument("--max-groups", type=int, default=50,
                        help="Max duplicate groups to print details for")
    args = parser.parse_args()

    print(f"Loading {args.h5ad}...")
    t0 = time.time()
    adata = ad.read_h5ad(args.h5ad, backed="r")
    print(f"  Shape: {adata.shape} in {time.time() - t0:.1f}s")

    # Check for duplicate obs_names (cell IDs)
    dup_ids = adata.obs_names.duplicated()
    n_dup_ids = dup_ids.sum()
    print(f"\n--- Duplicate cell IDs (obs_names) ---")
    print(f"  {n_dup_ids:,} duplicate cell IDs out of {adata.n_obs:,}")
    if n_dup_ids > 0:
        dup_examples = adata.obs_names[dup_ids][:10]
        print(f"  Examples: {list(dup_examples)}")

    # Load X into memory for duplicate row detection
    print(f"\n--- Loading X matrix into memory ---")
    t0 = time.time()
    X = adata.X[:]  # Load from backed mode
    if sp.issparse(X):
        X = X.tocsr()
    print(f"  Loaded: {X.shape}, nnz={X.nnz:,} in {time.time() - t0:.1f}s")

    # Check for zero rows
    row_nnz = np.diff(X.indptr)
    n_zero = (row_nnz == 0).sum()
    if n_zero > 0:
        print(f"\n  WARNING: {n_zero:,} all-zero rows (empty cells)")

    # Find duplicate rows
    print(f"\n--- Finding duplicate count vectors ---")
    dup_groups = find_duplicate_rows_sparse(X)

    n_groups = len(dup_groups)
    n_dup_cells = sum(len(g) for g in dup_groups)
    # The validator likely counts all cells in dup groups minus one per group
    # (i.e., the "extra" duplicates)
    n_extra = sum(len(g) - 1 for g in dup_groups)

    print(f"\n--- Results ---")
    print(f"  Duplicate groups: {n_groups:,}")
    print(f"  Total cells in dup groups: {n_dup_cells:,}")
    print(f"  Extra duplicates (n - 1 per group): {n_extra:,}")

    if n_groups == 0:
        print("\n  No duplicate count vectors found.")
        # Check if validator means something else
        print("\n  The HCA validator may be checking something other than")
        print("  exact row equality. Consider checking the validator source.")
        return

    # Analyze by study
    print(f"\n--- Study distribution of duplicate cells ---")
    has_dataset = "dataset" in adata.obs.columns
    has_donor = "donor_id" in adata.obs.columns

    obs = adata.obs

    if has_dataset:
        study_col = "dataset"
    elif has_donor:
        study_col = "donor_id"  # Will extract study prefix
    else:
        study_col = None

    # Collect all cells in dup groups
    all_dup_indices = []
    for g in dup_groups:
        all_dup_indices.extend(g)

    if study_col:
        dup_studies = obs.iloc[all_dup_indices][study_col].values
        if study_col == "donor_id":
            # Extract study prefix (e.g., "Gray" from "Gray_HBCA_Donor_1")
            dup_studies = [str(s).split("_")[0] for s in dup_studies]

        study_counts = Counter(dup_studies)
        for study, count in sorted(study_counts.items(), key=lambda x: -x[1]):
            print(f"  {study}: {count:,} cells")

    # Analyze group sizes
    group_sizes = [len(g) for g in dup_groups]
    size_counts = Counter(group_sizes)
    print(f"\n--- Group size distribution ---")
    for size, count in sorted(size_counts.items()):
        print(f"  Size {size}: {count:,} groups")

    # Cross-study analysis: are duplicates within same study or across studies?
    if study_col:
        print(f"\n--- Cross-study duplicate analysis ---")
        n_cross = 0
        n_within = 0
        cross_pairs = Counter()

        for g in dup_groups:
            studies_in_group = set()
            for idx in g:
                s = str(obs.iloc[idx][study_col])
                if study_col == "donor_id":
                    s = s.split("_")[0]
                studies_in_group.add(s)

            if len(studies_in_group) > 1:
                n_cross += 1
                pair = tuple(sorted(studies_in_group))
                cross_pairs[pair] += 1
            else:
                n_within += 1

        print(f"  Within-study duplicate groups: {n_within:,}")
        print(f"  Cross-study duplicate groups: {n_cross:,}")
        if cross_pairs:
            print(f"  Cross-study pairs:")
            for pair, count in sorted(cross_pairs.items(), key=lambda x: -x[1]):
                print(f"    {' × '.join(pair)}: {count:,} groups")

    # Print example groups
    print(f"\n--- Example duplicate groups (first {min(args.max_groups, n_groups)}) ---")
    for i, g in enumerate(dup_groups[:args.max_groups]):
        cells = [adata.obs_names[idx] for idx in g]
        row = X[g[0]]
        nnz = row.nnz

        if study_col:
            studies = [str(obs.iloc[idx][study_col]) for idx in g]
            if study_col == "donor_id":
                studies = [s.split("_")[0] for s in studies]
            study_info = f" studies={set(studies)}"
        else:
            study_info = ""

        print(f"  Group {i + 1}: {len(g)} cells, nnz={nnz}{study_info}")
        for cell in cells[:5]:
            print(f"    {cell}")
        if len(cells) > 5:
            print(f"    ... and {len(cells) - 5} more")

    # Check specifically for zero-count duplicates
    n_zero_dup_groups = sum(1 for g in dup_groups if X[g[0]].nnz == 0)
    if n_zero_dup_groups > 0:
        n_zero_dup_cells = sum(len(g) for g in dup_groups if X[g[0]].nnz == 0)
        print(f"\n--- Zero-count duplicates ---")
        print(f"  {n_zero_dup_groups:,} groups are all-zero rows")
        print(f"  {n_zero_dup_cells:,} cells total in zero-count groups")
        print(f"  These are empty cells with identical (zero) count vectors")


if __name__ == "__main__":
    main()
