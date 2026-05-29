#!/usr/bin/env python3
"""Quick diagnostic: inspect obs group structure in h5ad."""
import h5py
import sys

h5ad = sys.argv[1]
with h5py.File(h5ad, "r") as f:
    obs = f["obs"]
    print("obs attrs:", dict(obs.attrs))
    print()
    for key in list(obs.keys())[:30]:
        item = obs[key]
        if isinstance(item, h5py.Dataset):
            print(f"  {key}: Dataset, shape={item.shape}, dtype={item.dtype}")
            if item.shape[0] > 0 and item.shape[0] < 10:
                print(f"    values: {item[()]}")
            elif item.shape[0] > 0:
                print(f"    first 3: {item[:3]}")
        elif isinstance(item, h5py.Group):
            print(f"  {key}: Group, keys={list(item.keys())}")
            for sk in item.keys():
                si = item[sk]
                if isinstance(si, h5py.Dataset):
                    print(f"    {sk}: shape={si.shape}, dtype={si.dtype}, attrs={dict(si.attrs)}")
                    if si.shape[0] <= 5:
                        print(f"      values: {si[()]}")
                    else:
                        print(f"      first 3: {si[:3]}")
            print(f"    group attrs: {dict(item.attrs)}")
    print(f"\n  Total keys: {len(list(obs.keys()))}")
