"""Shared comparison fit behind ``fit_wt`` and ``fit_pair``.

A comparison is the grouped-binomial fit of one pair ``(from, to)``. Both
public fits reduce to a list of pairs and call ``fit_comparisons`` once; this
module owns the count-matrix pivot, design matrix, chunk loop, batch-to-result
conversion and the scalar-oracle fallback. The IRLS numerics live in
``_fit`` and are looked up at call time so that tests can patch them there.
"""

from __future__ import annotations

import logging
import math

import numpy as np
import polars as pl
from scipy.stats import norm

# Import the module, not names: `from ._fit import <name>` at module level would
# fail on the circular import; look attributes up at call time.
from . import _fit

_log = logging.getLogger(__name__)


def _failed(message: str, converged=None, iterations=None) -> dict:
    return {
        "m": None,
        "se": None,
        "p": None,
        "status": "failed",
        "dispersion": None,
        "df_residual": None,
        "diagnostics": _fit._diagnostics(converged, iterations, error=message),
    }


def _comparison_result(
    m, se, converged, iterations, warning, pearson, df_residual
) -> dict:
    """Result dict for the interaction estimate ``m`` with standard error ``se``.

    ``se`` is the raw Wald SE (dispersion fixed at 1) and ``p`` a z test.
    ``dispersion`` is the raw Pearson ``phi = pearson / df_residual``, reported
    as a diagnostic. With no residual df (one replicate per time point) ``phi``
    is NaN and reported as null; the fit is unaffected.
    """
    m, se = float(m), float(se)
    phi = float(pearson) / df_residual if df_residual > 0 else float("nan")
    p = float(2 * norm.sf(abs(m / se))) if se > 0 else float("nan")
    if not all(math.isfinite(value) for value in (m, se, p)):
        return _failed("non-finite estimate", bool(converged), int(iterations))
    return {
        "m": m,
        "se": se,
        "p": p,
        "status": "converged" if converged else "not_converged",
        "dispersion": phi if math.isfinite(phi) else None,
        "df_residual": int(df_residual),
        "diagnostics": _fit._diagnostics(
            bool(converged), int(iterations), warning=warning
        ),
    }


def _fit_comparison(
    joined: pl.DataFrame,
    reference: str,
    test: str,
    subset: pl.DataFrame | None = None,
) -> dict:
    """Fit ``cbind(Nreads, Diff) ~ 1 + type*time`` for one comparison.

    ``type`` is 1 for ``test``; ``time`` is the numeric sample time. Returns
    ``m, se, p, status, dispersion, df_residual, diagnostics``. ``subset``
    optionally supplies the reference and test rows already sliced from
    ``joined`` in its row order.
    """
    if subset is None:
        subset = joined.filter(pl.col("sequence_id").is_in([reference, test]))
    is_test = np.asarray(
        [value == test for value in subset.get_column("sequence_id").to_list()],
        dtype=float,
    )
    time = subset.get_column("time").to_numpy().astype(float)
    later = time == time.max()
    for variant in (0.0, 1.0):
        for is_later in (False, True):
            if not np.any((is_test == variant) & (later == is_later)):
                raise ValueError(f"comparison {test!r} is missing a type/time cell")
    design = np.column_stack([np.ones(subset.height), is_test, time, is_test * time])
    success = subset.get_column("Nreads").to_numpy().astype(float)
    failure = subset.get_column("Diff").to_numpy().astype(float)

    try:
        with np.errstate(all="ignore"):
            fit = _fit._logit_fit(success, failure, design)
    except Exception as exc:  # noqa: BLE001 - any fit error yields a failed row
        return _failed(str(exc))
    with np.errstate(all="ignore"):
        return _comparison_result(
            fit.beta[3],
            fit.se[3],
            fit.converged,
            fit.iterations,
            " | ".join(dict.fromkeys(fit.warnings)) or None,
            fit.pearson,
            fit.df_residual,
        )


def _needs_oracle(result: dict) -> bool:
    """True unless the batch result is a clean convergence.

    Near separation the IRLS iterates are chaotic, so only clean fits are
    trusted from the batch solver; everything else is re-run through the scalar
    ``_logit_fit``.
    """
    diagnostics = result["diagnostics"]
    return (
        result["status"] != "converged"
        or diagnostics["warning"] is not None
        or diagnostics["error"] is not None
    )


def _count_matrix(
    joined: pl.DataFrame, columns: list[str], ids: list, values: str
) -> np.ndarray:
    """Pivot ``values`` to an ``(len(ids), len(columns))`` matrix in ``ids`` order."""
    # Pivot row order is unspecified, so reindex explicitly by sequence_id.
    wide = joined.pivot(on="sample_id", index="sequence_id", values=values)
    position = {key: row for row, key in enumerate(wide["sequence_id"].to_list())}
    rows = [position[key] for key in ids]
    return wide.select(columns).to_numpy().astype(float)[rows]


def _fit_comparisons_vectorised(
    joined: pl.DataFrame,
    samples: pl.DataFrame,
    sources: list,
    targets: list,
    chunk_size: int,
) -> list[dict]:
    """Fit every pair in batches; one result dict per pair, in order.

    Counts are pivoted once to a sequence-by-sample matrix and each pair's two
    rows are gathered into its design block. Pairs without a clean convergence
    are re-fitted by the scalar oracle.
    """
    columns = [str(value) for value in samples.get_column("sample_id").to_list()]
    ids = list(dict.fromkeys([*sources, *targets]))
    index = {key: row for row, key in enumerate(ids)}
    nreads = _count_matrix(joined, columns, ids, "Nreads")
    diff = _count_matrix(joined, columns, ids, "Diff")
    time = samples.get_column("time").to_numpy().astype(float)
    n = len(columns)
    is_test = np.concatenate([np.zeros(n), np.ones(n)])
    both_time = np.concatenate([time, time])
    design = np.column_stack([np.ones(2 * n), is_test, both_time, is_test * both_time])

    from_idx = np.array([index[value] for value in sources])
    to_idx = np.array([index[value] for value in targets])
    results: list[dict] = []
    for start in range(0, len(sources), chunk_size):
        block = slice(start, start + chunk_size)
        size = len(from_idx[block])
        success = np.hstack([nreads[from_idx[block]], nreads[to_idx[block]]])
        failure = np.hstack([diff[from_idx[block]], diff[to_idx[block]]])
        stacked = np.broadcast_to(design, (size, 2 * n, design.shape[1]))
        fit = _fit._logit_fit_batch(success, failure, stacked)
        for k in range(size):
            if fit.errors[k] is not None:
                results.append(_failed(fit.errors[k]))
                continue
            results.append(
                _comparison_result(
                    fit.beta[k, 3],
                    fit.se[k, 3],
                    fit.converged[k],
                    fit.iterations[k],
                    " | ".join(dict.fromkeys(fit.warnings[k])) or None,
                    fit.pearson[k],
                    fit.df_residual,
                )
            )
    unclean = [k for k, r in enumerate(results) if _needs_oracle(r)]
    _log.debug(
        "vectorised: %d of %d comparisons re-fitted by the scalar oracle",
        len(unclean),
        len(results),
    )
    if unclean:
        # Slice each sequence's rows once; restoring the joined row order keeps
        # the scalar fit identical to filtering the full table.
        wanted = {value for k in unclean for value in (sources[k], targets[k])}
        slices = (
            joined.with_row_index("_row")
            .filter(pl.col("sequence_id").is_in(list(wanted)))
            .partition_by("sequence_id", as_dict=True, maintain_order=True)
        )
        for k in unclean:
            subset = pl.concat([slices[(sources[k],)], slices[(targets[k],)]]).sort(
                "_row"
            )
            results[k] = _fit_comparison(
                joined, sources[k], targets[k], subset.drop("_row")
            )
    return results


def fit_comparisons(
    joined: pl.DataFrame,
    samples: pl.DataFrame,
    sources: list,
    targets: list,
    engine: str,
    chunk_size: int,
) -> list[dict]:
    """Fit the comparison ``(sources[k], targets[k])`` for every ``k``.

    ``joined`` is the validated counts joined to ``samples`` (``_validate``).
    Returns one result dict (``m, se, p, status, dispersion, df_residual,
    diagnostics, flags``) per pair, in
    order. ``engine`` is ``"vectorised"`` or ``"scalar"``; both agree within
    1e-9 and ``chunk_size`` does not change it.
    """
    if engine == "vectorised":
        results = _fit_comparisons_vectorised(
            joined, samples, sources, targets, chunk_size
        )
    else:
        results = [
            _fit_comparison(joined, source, target)
            for source, target in zip(sources, targets, strict=True)
        ]
    # Flags are per genotype from the raw counts, unioned over both genotypes
    # so a pair and its reversal agree; they sit outside the engines.
    by_genotype = _fit._flags(joined)
    return [
        {
            **result,
            "flags": _fit._join_flags(by_genotype[source], by_genotype[target]),
        }
        for result, source, target in zip(results, sources, targets, strict=True)
    ]
