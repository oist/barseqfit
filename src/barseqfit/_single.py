"""Per-sequence growth rates: ``fit_sequences``.

The growth rate ``r`` of one sequence is the slope per ``time_unit`` of
``logit(count / sample_total)`` against ``time`` from that sequence's own rows
(grouped binomial GLM, ``cbind(Nreads, Diff) ~ 1 + time``, design ``[ones,
time]``). Because the ``type*time`` design of ``fit_wt`` is block-separable,
``fit_wt``'s ``m`` equals ``r_i - r_wt`` at convergence.

This is a sibling of ``_comparison`` rather than part of it: a comparison is a
two-genotype block with four coefficients, a growth rate a one-genotype block
with two, and there is no reference, no ``type`` column and no dispersion
option. The IRLS numerics (``_logit_fit``, ``_logit_fit_batch``), the result
dict, the scalar-oracle test and the count pivot are the ones ``_comparison``
uses, looked up at call time so that tests can patch them in ``_fit``.
"""

from __future__ import annotations

import logging

import numpy as np
import polars as pl

from . import _comparison, _fit

_log = logging.getLogger(__name__)

ZERO_COUNT_FLAG = "zero_count_replicate"

_REPLICATE_SCHEMA = {
    "sequence_id": pl.String,
    "replicate_id": pl.String,
    "r": pl.Float64,
    "se": pl.Float64,
    "status": pl.String,
    "flags": pl.String,
    "diagnostics": _fit._DIAGNOSTICS,
    "elapsed_time": pl.Float64,
    "time_unit": pl.String,
}

_SEQUENCE_SCHEMA = {
    "sequence_id": pl.String,
    "r": pl.Float64,
    "se_pooled": pl.Float64,
    "p_pooled": pl.Float64,
    "n_replicates": pl.Int64,
    "n_replicates_converged": pl.Int64,
    "r_replicate_mean": pl.Float64,
    "r_replicate_sd": pl.Float64,
    "se_replicate": pl.Float64,
    "status": pl.String,
    "dispersion": pl.Float64,
    "df_residual": pl.Int64,
    "flags": pl.String,
    "diagnostics": _fit._DIAGNOSTICS,
    "elapsed_time": pl.Float64,
    "time_unit": pl.String,
    "pseudocount": pl.String,
}


def _result(fit, k) -> dict:
    """Result dict (``r, se, p, status, ...``) of fit ``k`` of a batch."""
    if fit.errors[k] is not None:
        return _comparison._failed(fit.errors[k])
    with np.errstate(all="ignore"):
        return _comparison._comparison_result(
            fit.beta[k, 1],
            fit.se[k, 1],
            fit.converged[k],
            fit.iterations[k],
            " | ".join(dict.fromkeys(fit.warnings[k])) or None,
            fit.pearson[k],
            fit.df_residual,
        )


def _fit_sequence(subset: pl.DataFrame) -> dict:
    """Fit ``cbind(Nreads, Diff) ~ 1 + time`` to one sequence's rows (scalar)."""
    time = subset.get_column("time").to_numpy().astype(float)
    design = np.column_stack([np.ones(subset.height), time])
    success = subset.get_column("Nreads").to_numpy().astype(float)
    failure = subset.get_column("Diff").to_numpy().astype(float)
    try:
        with np.errstate(all="ignore"):
            fit = _fit._logit_fit(success, failure, design)
    except Exception as exc:  # noqa: BLE001 - any fit error yields a failed row
        return _comparison._failed(str(exc))
    with np.errstate(all="ignore"):
        return _comparison._comparison_result(
            fit.beta[1],
            fit.se[1],
            fit.converged,
            fit.iterations,
            " | ".join(dict.fromkeys(fit.warnings)) or None,
            fit.pearson,
            fit.df_residual,
        )


def _fit_vectorised(
    joined: pl.DataFrame,
    samples: pl.DataFrame,
    ids: list,
    chunk_size: int,
) -> list[dict]:
    """Fit every sequence in batches; one result dict per id, in order.

    Counts are pivoted once to a sequence-by-sample matrix and stacked in
    chunks of ``chunk_size``. Sequences without a clean convergence are
    re-fitted by the scalar oracle.
    """
    columns = [str(value) for value in samples.get_column("sample_id").to_list()]
    nreads = _comparison._count_matrix(joined, columns, ids, "Nreads")
    diff = _comparison._count_matrix(joined, columns, ids, "Diff")
    time = samples.get_column("time").to_numpy().astype(float)
    design = np.column_stack([np.ones(len(columns)), time])
    results: list[dict] = []
    for start in range(0, len(ids), chunk_size):
        block = slice(start, start + chunk_size)
        success, failure = nreads[block], diff[block]
        stacked = np.broadcast_to(design, (len(success), *design.shape))
        fit = _fit._logit_fit_batch(success, failure, stacked)
        results.extend(_result(fit, k) for k in range(len(success)))
    unclean = [
        k for k, result in enumerate(results) if _comparison._needs_oracle(result)
    ]
    _log.debug(
        "vectorised: %d of %d sequences re-fitted by the scalar oracle",
        len(unclean),
        len(results),
    )
    if unclean:
        slices = joined.filter(
            pl.col("sequence_id").is_in([ids[k] for k in unclean])
        ).partition_by("sequence_id", as_dict=True, maintain_order=True)
        for k in unclean:
            results[k] = _fit_sequence(slices[(ids[k],)])
    return results


def _fit_all(
    joined: pl.DataFrame,
    samples: pl.DataFrame,
    ids: list,
    engine: str,
    chunk_size: int,
) -> list[dict]:
    """One result dict per id from the chosen engine (``vectorised`` or ``scalar``)."""
    if engine == "vectorised":
        return _fit_vectorised(joined, samples, ids, chunk_size)
    slices = joined.partition_by("sequence_id", as_dict=True, maintain_order=True)
    return [_fit_sequence(slices[(value,)]) for value in ids]


def _replicate_samples(samples: pl.DataFrame) -> list[tuple[object, pl.DataFrame]]:
    """Split ``samples`` into one before/after pair per ``replicate_id``.

    Validates exactly one sample per (replicate_id, time) and that every
    replicate has both time points; returns replicates in order of first
    appearance in ``samples``.
    """
    samples = _fit._as_frame(samples, "samples", {"sample_id", "replicate_id", "time"})
    if samples.get_column("replicate_id").null_count():
        raise ValueError("samples.replicate_id must not contain nulls")
    if samples.get_column("time").null_count():
        raise ValueError("samples.time must be finite")
    times = samples.get_column("time").unique().sort().to_list()
    if len(times) != 2:
        raise ValueError("fits require exactly two sample time points")
    per_cell = samples.group_by("replicate_id", "time").len()
    if per_cell.filter(pl.col("len") > 1).height:
        raise ValueError(
            "samples must have exactly one sample per (replicate_id, time)"
        )
    replicates = samples.get_column("replicate_id").unique(maintain_order=True)
    if per_cell.height != 2 * len(replicates):
        raise ValueError("every replicate_id must have a sample at both time points")
    return [(rep, samples.filter(pl.col("replicate_id") == rep)) for rep in replicates]


def _zero_count_ids(joined: pl.DataFrame) -> set:
    """Sequence ids with a raw zero count in ``joined`` (before any pseudocount)."""
    zero = joined.filter(pl.col("count") == 0)
    return set(zero.get_column("sequence_id").unique().to_list())


def _replicate_table(
    counts: pl.DataFrame,
    samples: pl.DataFrame,
    sequences: pl.DataFrame,
    *,
    elapsed_time: float,
    time_unit: str,
    engine: str,
    chunk_size: int,
) -> pl.DataFrame:
    """Per-replicate fits for every sequence (inputs already validated in full)."""
    ids = sequences.get_column("sequence_id").to_list()
    counts = _fit._as_frame(counts, "counts", {"sequence_id", "sample_id", "count"})
    by_replicate = []
    replicate_ids = []
    for replicate, subset in _replicate_samples(samples):
        # Each replicate is a stand-alone 2-row fit: the pseudocount and the
        # flags come from that replicate's own raw counts.
        joined = _fit._validate(
            counts.filter(pl.col("sample_id").is_in(subset["sample_id"].to_list())),
            subset,
            sequences,
            elapsed_time,
            time_unit,
        )
        results = _fit_all(joined, subset, ids, engine, chunk_size)
        flags = _fit._flags(joined)
        zeros = _zero_count_ids(joined)
        by_replicate.append(
            [
                {
                    **result,
                    "flags": _fit._join_flags(
                        flags[value], [ZERO_COUNT_FLAG] if value in zeros else []
                    ),
                }
                for value, result in zip(ids, results, strict=True)
            ]
        )
        replicate_ids.append(str(replicate))
    rows = [
        {
            "sequence_id": value,
            "replicate_id": replicate_ids[j],
            "r": by_replicate[j][i]["m"],
            "se": by_replicate[j][i]["se"],
            "status": by_replicate[j][i]["status"],
            "flags": by_replicate[j][i]["flags"],
            "diagnostics": by_replicate[j][i]["diagnostics"],
            "elapsed_time": float(elapsed_time),
            "time_unit": time_unit,
        }
        for i, value in enumerate(ids)
        for j in range(len(by_replicate))
    ]
    return pl.DataFrame(rows, schema=_REPLICATE_SCHEMA)


def _aggregate(replicates: pl.DataFrame, ids: list) -> pl.DataFrame:
    """Replicate columns per sequence, from the ``_replicate_table`` output."""
    converged = pl.col("status") == "converged"
    grouped = replicates.group_by("sequence_id").agg(
        pl.len().cast(pl.Int64).alias("n_replicates"),
        converged.sum().cast(pl.Int64).alias("n_replicates_converged"),
        pl.col("r").filter(converged).mean().alias("r_replicate_mean"),
        pl.col("r").filter(converged).std(ddof=1).alias("r_replicate_sd"),
        pl.col("flags")
        .str.contains(ZERO_COUNT_FLAG, literal=True)
        .any()
        .alias("_zero_count"),
    )
    return (
        pl.DataFrame({"sequence_id": ids}, schema={"sequence_id": pl.String})
        .join(grouped, on="sequence_id", how="left", maintain_order="left")
        .with_columns(
            pl.col("n_replicates", "n_replicates_converged").fill_null(0),
            pl.col("_zero_count").fill_null(False),
        )
        .with_columns(
            pl.when(pl.col("n_replicates_converged") >= 2)
            .then(pl.col("r_replicate_sd"))
            .otherwise(None)
            .alias("r_replicate_sd"),
        )
        .with_columns(
            (pl.col("r_replicate_sd") / pl.col("n_replicates_converged").sqrt()).alias(
                "se_replicate"
            )
        )
    )


def fit_sequences(
    counts: pl.DataFrame,
    samples: pl.DataFrame,
    sequences: pl.DataFrame,
    *,
    elapsed_time: float,
    time_unit: str,
    pseudocount: str,
    engine: str = "vectorised",
    chunk_size: int = 2000,
) -> pl.DataFrame:
    """Estimate each sequence's own growth rate ``r``.

    ``r`` is the slope per ``time_unit`` of ``logit(count / sample_total)``
    against ``samples.time``, a grouped binomial GLM with logit link
    ``cbind(Nreads, Diff) ~ 1 + time`` fitted on the sequence's own rows (all
    samples, pooled over replicates). This is the per-genotype half of
    ``fit_wt``'s ``type*time`` model: the design is block-separable, so at
    convergence ``fit_wt``'s ``m`` equals ``r_i - r_wt`` and, under fixed
    dispersion, its ``se`` equals ``sqrt(se_pooled_i**2 + se_pooled_wt**2)``.

    Inputs, validation, ``engine`` and ``chunk_size`` are as in ``fit_wt``
    (there is no ``reference_id`` and no ``dispersion`` argument; the pooled
    standard error is the fixed-dispersion one), with one stricter requirement
    on the replicate structure: ``samples`` must have exactly one sample per
    ``(replicate_id, time)`` and every replicate must have both time points,
    otherwise ``ValueError``. ``pseudocount`` is a required keyword argument
    with no default; the only accepted value is ``"all_zeros"`` (every zero
    count is replaced by 1), and any other value raises ``ValueError``.

    One row per sequence of ``sequences``, in order, with columns:

    - ``r``, ``se_pooled``, ``p_pooled``: the pooled estimate, its
      fixed-dispersion Wald SE (comparable with ``fit_wt``'s ``se`` under the
      default) and the two-sided z test of ``r = 0``. Null when ``status`` is
      ``failed``.
    - ``n_replicates``: replicates (``replicate_id``), each with both time
      points; ``n_replicates_converged``: those whose replicate fit (the same
      model on that replicate's two samples) has status ``converged``.
    - ``r_replicate_mean``: mean of the converged replicate slopes;
      ``r_replicate_sd``: their sample SD (ddof 1), null below two converged
      replicates; ``se_replicate = r_replicate_sd / sqrt(n_replicates_converged)``,
      null when the SD is. ``se_replicate`` is the SE that captures
      between-replicate overdispersion; ``se_pooled`` does not.
    - ``status``: of the pooled IRLS only: ``converged``, ``not_converged``
      (estimate kept) or ``failed`` (error or non-finite estimate). A failed
      replicate fit never changes it; it only lowers
      ``n_replicates_converged``.
    - ``dispersion``: the pooled fit's Pearson ``phi`` (diagnostic only; not
      used by ``se_pooled``); ``df_residual``: rows minus 2. Both null when the
      pooled fit fails and, for ``dispersion``, when there is no residual df
      (one replicate).
    - ``flags``: ``censored`` as in ``fit_wt`` (from the
      pooled raw counts) plus ``zero_count_replicate``, sorted and
      ``|``-joined, null when none. ``zero_count_replicate`` is set when any
      replicate has a raw zero count at either time. It propagates to the
      pooled row because ``r_replicate_mean`` and ``se_replicate`` are
      computed from pseudocount-patched slopes of those replicates, so the
      caller can see which sequences' replicate spread rests on a patched
      zero. Flags never change any estimate or status.
    - ``diagnostics``: of the pooled fit.
    - ``elapsed_time``, ``time_unit``, ``pseudocount``: the arguments, echoed
      on every row (``pseudocount`` is a new String column).

    The pseudocount applies to the pooled fit over all of a sequence's samples
    and, separately, to each replicate over its own two samples.
    """
    _fit._check_pseudocount(pseudocount)
    joined = _fit._check_common(
        counts, samples, sequences, elapsed_time, time_unit, engine, chunk_size
    )
    ids = sequences.get_column("sequence_id").to_list()
    results = _fit_all(joined, samples, ids, engine, chunk_size)
    flags = _fit._flags(joined)
    replicates = _replicate_table(
        counts,
        samples,
        sequences,
        elapsed_time=elapsed_time,
        time_unit=time_unit,
        engine=engine,
        chunk_size=chunk_size,
    )
    aggregated = _aggregate(replicates, ids)
    zeros = set(
        aggregated.filter(pl.col("_zero_count")).get_column("sequence_id").to_list()
    )
    pooled = pl.DataFrame(
        [
            {
                "sequence_id": value,
                "r": result["m"],
                "se_pooled": result["se"],
                "p_pooled": result["p"],
                "status": result["status"],
                "dispersion": result["dispersion"],
                "df_residual": result["df_residual"],
                "flags": _fit._join_flags(
                    flags[value], [ZERO_COUNT_FLAG] if value in zeros else []
                ),
                "diagnostics": result["diagnostics"],
                "elapsed_time": float(elapsed_time),
                "time_unit": time_unit,
                "pseudocount": pseudocount,
            }
            for value, result in zip(ids, results, strict=True)
        ],
        schema={
            key: dtype
            for key, dtype in _SEQUENCE_SCHEMA.items()
            if key
            not in {
                "n_replicates",
                "n_replicates_converged",
                "r_replicate_mean",
                "r_replicate_sd",
                "se_replicate",
            }
        },
    )
    return pooled.join(
        aggregated.drop("_zero_count"),
        on="sequence_id",
        how="left",
        maintain_order="left",
    ).select(list(_SEQUENCE_SCHEMA))
