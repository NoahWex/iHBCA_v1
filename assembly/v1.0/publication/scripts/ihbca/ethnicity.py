"""Shared ethnicity normalization for iHBCA assembly pipeline.

Functions for normalizing multi-value ethnicity fields and applying
HANCESTRO term downgrades between HCA-native and CxG 5.3.2 schemas.
"""

import re

import pandas as pd

from ihbca.constants import HCA_TO_CXG_HANCESTRO


def normalize_multi_ethnicity(series, sep=" || "):
    """Normalize multi-value ethnicity: sorted, ascending order.

    sep: ' || ' for HCA (default), ',' for CxG.
    """
    def _normalize(val):
        if pd.isna(val) or val in ("unknown", "nan", ""):
            return val
        # Split on either ' || ' or ',' (handle both inputs)
        parts = [p.strip() for p in re.split(r'\s*\|\|\s*|,', str(val))]
        parts = [p for p in parts if p]
        parts = sorted(set(parts))
        return sep.join(parts)
    return series.map(_normalize)


def downgrade_hancestro_terms(series):
    """Downgrade HCA-native HANCESTRO terms to CxG 5.3.2 :0004 branch.

    Replaces :0601/:0602 terms with their :0004 equivalents per
    HCA_TO_CXG_HANCESTRO mapping. Only call this for --target cxg builds.
    """
    result = series
    for hca_term, cxg_term in HCA_TO_CXG_HANCESTRO.items():
        result = result.str.replace(hca_term, cxg_term, regex=False)
    return result
