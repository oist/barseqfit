"""Multiple-testing adjustment."""

from __future__ import annotations

import numpy as np


def benjamini_hochberg(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values over the non-NaN entries of ``p``.

    ``n`` is the number of non-NaN entries; NaN entries pass through as NaN.
    Matches ``stats::p.adjust(p, "BH")``: sort descending, cumulative minimum of
    ``p * n / rank``, clipped at 1. Empty input returns an empty array.
    """
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    valid = np.flatnonzero(~np.isnan(p))
    n = valid.size
    if n == 0:
        return out
    order = valid[np.argsort(-p[valid], kind="stable")]
    ranks = np.arange(n, 0, -1)
    adjusted = np.minimum.accumulate((n / ranks) * p[order])
    out[order] = np.minimum(adjusted, 1.0)
    return out
