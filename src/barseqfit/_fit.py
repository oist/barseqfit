"""Grouped-binomial wild-type fitness estimation."""

from __future__ import annotations

import math
from typing import Any, NamedTuple

import numpy as np
import polars as pl
import scipy.linalg

from . import _comparison
from ._adjust import benjamini_hochberg


def _as_frame(value: Any, name: str, columns: set[str]) -> pl.DataFrame:
    if not isinstance(value, pl.DataFrame):
        raise TypeError(f"{name} must be a Polars DataFrame")
    missing = columns.difference(value.columns)
    if missing:
        raise ValueError(f"{name} is missing required columns: {sorted(missing)}")
    return value


def _require_finite_integers(frame: pl.DataFrame, column: str, label: str) -> None:
    numeric = pl.col(column).cast(pl.Float64, strict=False)
    invalid = frame.filter(
        numeric.is_null() | ~numeric.is_finite() | (numeric != numeric.floor())
    )
    if invalid.height:
        raise ValueError(f"{label} must contain finite integer values")


def _y_log_y(y: np.ndarray, mu: np.ndarray) -> np.ndarray:
    """R's ``y_log_y``: ``y * log(y / mu)``, zero where ``y == 0``."""
    return np.where(y != 0.0, y * np.log(y / mu), 0.0)


def _deviance(success: np.ndarray, failure: np.ndarray, mu: np.ndarray) -> float:
    """Binomial deviance as R's ``binomial()$dev.resids``.

    Computed on the proportion scale, ``2 * n * (y*log(y/mu) +
    (1-y)*log((1-y)/(1-mu)))``, with ``1 - y`` and ``1 - mu`` formed first.
    The algebraically equal count-scale form loses about 1e-8 absolute
    accuracy when ``n`` is large, which flips the ``1e-8`` convergence test
    against R near the threshold.
    """
    total = success + failure
    y = success / total
    with np.errstate(divide="ignore", invalid="ignore"):
        terms = 2.0 * total * (_y_log_y(y, mu) + _y_log_y(1.0 - y, 1.0 - mu))
    return float(np.sum(terms))


_EPS = float(np.finfo(float).eps)


class IrlsResult(NamedTuple):
    """Result of ``_logit_fit`` (one fit) or ``_logit_fit_batch`` (``G`` fits).

    ``beta`` and ``se`` are ``(p,)`` or ``(G, p)``; ``warnings`` is the list of
    glm.fit messages (one list per fit in the batch); ``errors`` is ``None``
    for the scalar fit (which raises instead) and ``None`` or the failure
    message per fit in the batch. ``pearson`` is the Pearson chi-square of the
    final iteration (``_pearson_statistic``, NaN for failed batch fits) and
    ``df_residual`` the rows minus coefficients.
    """

    beta: np.ndarray
    se: np.ndarray
    converged: Any
    iterations: Any
    warnings: list
    errors: Any
    pearson: Any
    df_residual: int


def _linkinv(eta: np.ndarray) -> np.ndarray:
    """R's logit linkinv: exp(eta) clamped to [eps, 1/eps] beyond |eta| > 30."""
    tmp = np.where(eta < -30, _EPS, np.where(eta > 30, 1.0 / _EPS, np.exp(eta)))
    return tmp / (1.0 + tmp)


def _mu_eta(eta: np.ndarray) -> np.ndarray:
    opexp = 1.0 + np.exp(eta)
    return np.where(np.abs(eta) > 30, _EPS, np.exp(eta) / (opexp * opexp))


def _numeric_rank(matrix: np.ndarray, tol: float = 1e-11) -> int:
    """Rank from pivoted QR: diagonal of R relative to its largest entry."""
    diagonal = np.abs(np.diag(scipy.linalg.qr(matrix, mode="r", pivoting=True)[0]))
    if diagonal.size == 0 or diagonal.max() == 0:
        return 0
    return int(np.sum(diagonal > tol * diagonal.max()))


def _logit_fit(
    success: np.ndarray,
    failure: np.ndarray,
    design: np.ndarray,
) -> IrlsResult:
    """Fit grouped binomial logistic regression with Fisher scoring.

    Mirrors ``stats::glm.fit`` with ``epsilon=1e-8, maxit=25``: starting
    values ``mu = (y*n + 0.5) / (n + 1)``, initial deviance taken at the
    starting values, convergence when ``|dev - devold| / (|dev| + 0.1) <
    epsilon``, R's clamped logit inverse link and step-halving when the
    deviance is not finite. Returns an ``IrlsResult`` (``errors`` is ``None``;
    ``warnings`` uses glm.fit's message text). When the
    iteration limit is hit the last estimate is still returned with
    ``converged=False``. Raises ``ValueError`` for a rank-deficient design
    (QR tolerance 1e-11) or when no valid step can be found.
    """
    maxit = 25
    total = success + failure
    observed = success / total
    mu = (success + 0.5) / (total + 1.0)
    eta = np.log(mu / (1 - mu))
    mu = _linkinv(eta)
    previous_deviance = _deviance(success, failure, mu)
    converged = False
    boundary = False
    beta = np.zeros(design.shape[1], dtype=float)
    coefold = None
    weighted_design = design
    iteration = 0
    messages: list[str] = []

    for iteration in range(1, maxit + 1):
        mu_eta = _mu_eta(eta)
        variance = mu * (1.0 - mu)
        working = eta + (observed - mu) / mu_eta
        root_weights = np.sqrt(total * mu_eta**2 / variance)
        weighted_design = design * root_weights[:, None]
        weighted_response = working * root_weights
        if _numeric_rank(weighted_design) < design.shape[1]:
            raise ValueError("rank-deficient design")
        beta, solved = _qr_solve(weighted_design[None], weighted_response[None])
        if not solved:
            raise ValueError("rank-deficient design")
        eta = design @ beta
        mu = _linkinv(eta)
        deviance = _deviance(success, failure, mu)
        if not math.isfinite(deviance):
            if coefold is None:
                raise ValueError(
                    "no valid set of coefficients has been found: "
                    "please supply starting values"
                )
            messages.append("step size truncated due to divergence")
            step = 1
            while not math.isfinite(deviance):
                if step > maxit:
                    raise ValueError("inner loop 1; cannot correct step size")
                step += 1
                beta = (beta + coefold) / 2.0
                eta = design @ beta
                mu = _linkinv(eta)
                deviance = _deviance(success, failure, mu)
            boundary = True
        if abs(deviance - previous_deviance) / (abs(deviance) + 0.1) < 1e-8:
            converged = True
            break
        coefold = beta
        previous_deviance = deviance

    if not converged:
        messages.append("glm.fit: algorithm did not converge")
    if boundary:
        messages.append("glm.fit: algorithm stopped at boundary value")
    eps = 10 * _EPS
    if np.any(mu > 1 - eps) or np.any(mu < eps):
        messages.append("glm.fit: fitted probabilities numerically 0 or 1 occurred")

    # R's summary uses the weighted model matrix of the final iteration,
    # built immediately before the last coefficient update.
    covariance = np.linalg.inv(weighted_design.T @ weighted_design)
    pearson = float(_pearson_statistic(success, failure, root_weights**2, eta, mu))
    return IrlsResult(
        beta,
        np.sqrt(np.diag(covariance)),
        converged,
        iteration,
        messages,
        None,
        pearson,
        len(success) - design.shape[1],
    )


def _deviance_batch(
    success: np.ndarray, failure: np.ndarray, mu: np.ndarray
) -> np.ndarray:
    """Row-wise binomial deviance for ``(G, n)`` arrays, as ``_deviance``."""
    total = success + failure
    y = success / total
    with np.errstate(divide="ignore", invalid="ignore"):
        terms = 2.0 * total * (_y_log_y(y, mu) + _y_log_y(1.0 - y, 1.0 - mu))
    return np.sum(terms, axis=-1)


def _lstsq_batch(a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Least squares ``a[g] x ~ b[g]`` for every ``g`` by Householder QR.

    R's ``glm.fit`` solves each IRLS step with a Householder QR (LINPACK
    ``dqrdc2``). A LAPACK QR followed by back-substitution tracks that to
    about 1e-10 in the deviance, where an SVD solve is about 30 times noisier;
    that noise decides the ``1e-8`` convergence test at the threshold. Returns
    ``(x, ok)``; ``ok`` is False where the inputs are not finite, the
    factorisation fails, or ``R`` is numerically rank deficient (diagonal below
    ``1e-11`` of its largest entry, as ``_numeric_rank``).
    """
    count, _, p = a.shape
    out = np.full((count, p), np.nan)
    ok = np.isfinite(a).all(axis=(1, 2)) & np.isfinite(b).all(axis=1)
    if ok.any():
        sel = np.flatnonzero(ok)
        try:
            q, r = np.linalg.qr(a[sel])
        except np.linalg.LinAlgError:
            q = r = None
        if q is None:
            # Retry one at a time so only the offending genotype is marked failed.
            for i in sel:
                try:
                    out[i], ok[i] = _qr_solve(a[i : i + 1], b[i : i + 1])
                except np.linalg.LinAlgError:
                    ok[i] = False
        else:
            out[sel], ok[sel] = _back_substitute(r, np.einsum("gnp,gn->gp", q, b[sel]))
    return out, ok


def _qr_solve(a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, bool]:
    q, r = np.linalg.qr(a)
    x, ok = _back_substitute(r, np.einsum("gnp,gn->gp", q, b))
    return x[0], bool(ok[0])


def _back_substitute(r: np.ndarray, rhs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Solve upper-triangular ``r[g] x = rhs[g]``; ok False if rank deficient."""
    p = r.shape[1]
    diagonal = np.abs(np.diagonal(r, axis1=1, axis2=2))
    peak = diagonal.max(axis=1) if p else np.zeros(len(r))
    ok = (peak > 0) & np.all(diagonal > 1e-11 * peak[:, None], axis=1)
    x = np.full(rhs.shape, np.nan)
    with np.errstate(all="ignore"):
        for j in range(p - 1, -1, -1):
            partial = np.einsum("gk,gk->g", r[:, j, j + 1 :], x[:, j + 1 :])
            x[:, j] = (rhs[:, j] - partial) / r[:, j, j]
    return x, ok


def _invert_batch(a: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Invert every ``a[g]``; singular matrices give NaN and ``ok=False``."""
    try:
        return np.linalg.inv(a), np.ones(len(a), dtype=bool)
    except np.linalg.LinAlgError:
        out = np.full(a.shape, np.nan)
        ok = np.zeros(len(a), dtype=bool)
        for i in range(len(a)):
            try:
                out[i] = np.linalg.inv(a[i])
                ok[i] = True
            except np.linalg.LinAlgError:
                pass
        return out, ok


def _logit_fit_batch(
    success: np.ndarray,
    failure: np.ndarray,
    design: np.ndarray,
) -> IrlsResult:
    """Fit many grouped binomial logistic regressions at once.

    Vectorised counterpart of ``_logit_fit`` with the same glm.fit semantics
    (starting values, deviance convergence test, ``maxit=25``, R's clamped
    inverse link, step-halving, warning text). ``success`` and ``failure`` are
    ``(G, n)``, ``design`` is ``(G, n, p)``. Each iteration updates only the
    unconverged genotypes, so converged genotypes freeze and ``iterations``
    matches the scalar fit per genotype. Coefficients come from a batched
    Householder-QR least-squares solve (``_lstsq_batch``), the same solver as
    the scalar fit; the normal matrix is formed only for the standard errors.
    Returns an ``IrlsResult`` where ``warnings`` is a list of message lists and
    ``errors`` holds ``None`` or the failure message per genotype (beta, se and
    pearson are NaN for failed genotypes); ``pearson`` is the ``(G,)`` array of
    Pearson chi-squares of each genotype's final iteration.
    """
    maxit = 25
    g, _, p = design.shape
    total = success + failure
    observed = success / total
    mu = (success + 0.5) / (total + 1.0)
    eta = np.log(mu / (1 - mu))
    mu = _linkinv(eta)
    previous = _deviance_batch(success, failure, mu)
    beta = np.zeros((g, p))
    coefold = np.zeros((g, p))
    normal = np.full((g, p, p), np.nan)
    weight_used = np.full(mu.shape, np.nan)
    converged = np.zeros(g, dtype=bool)
    iterations = np.zeros(g, dtype=np.int64)
    truncated = np.zeros(g, dtype=bool)
    boundary = np.zeros(g, dtype=bool)
    errors: list[str | None] = [None] * g

    # Rank check once on the unweighted design. Rank loss induced by the IRLS
    # weights is not detected here: the QR solve flags a near-zero
    # diagonal of R as not ok and such genotypes end non-clean, so the scalar oracle re-fits them.
    if g and design.strides[0] == 0:
        # Broadcast of one shared design: factorise it once.
        r = np.broadcast_to(np.linalg.qr(design[0], mode="r"), (g, p, p))
    else:
        r = np.linalg.qr(design, mode="r")
    diagonal = np.abs(np.diagonal(r, axis1=1, axis2=2))
    peak = diagonal.max(axis=1) if diagonal.shape[1] else np.zeros(g)
    rank = np.sum(diagonal > 1e-11 * peak[:, None], axis=1)
    deficient = (peak == 0) | (rank < p)
    for i in np.flatnonzero(deficient):
        errors[i] = "rank-deficient design"
    active = ~deficient

    with np.errstate(all="ignore"):
        for iteration in range(1, maxit + 1):
            idx = np.flatnonzero(active)
            if idx.size == 0:
                break
            x = design[idx]
            e, m = eta[idx], mu[idx]
            mu_eta = _mu_eta(e)
            root = np.sqrt(total[idx] * mu_eta**2 / (m * (1.0 - m)))
            working = e + (observed[idx] - m) / mu_eta
            xw = x * root[:, :, None]
            normal[idx] = xw.transpose(0, 2, 1) @ xw
            weight_used[idx] = root**2
            b_new, ok = _lstsq_batch(xw, working * root)
            if not ok.all():
                for i in idx[~ok]:
                    errors[i] = "rank-deficient design"
                active[idx[~ok]] = False
                idx, x, b_new = idx[ok], x[ok], b_new[ok]
            s, f = success[idx], failure[idx]
            e_new = np.einsum("gnp,gp->gn", x, b_new)
            m_new = _linkinv(e_new)
            dev = _deviance_batch(s, f, m_new)
            alive = np.ones(idx.size, dtype=bool)
            bad = ~np.isfinite(dev)
            if bad.any():
                if iteration == 1:
                    for i in idx[bad]:
                        errors[i] = (
                            "no valid set of coefficients has been found: "
                            "please supply starting values"
                        )
                    alive &= ~bad
                else:
                    truncated[idx[bad]] = True
                    pos = np.flatnonzero(bad)
                    step = 1
                    while pos.size:
                        if step > maxit:
                            for i in idx[pos]:
                                errors[i] = "inner loop 1; cannot correct step size"
                            alive[pos] = False
                            break
                        step += 1
                        b_new[pos] = (b_new[pos] + coefold[idx[pos]]) / 2.0
                        e_new[pos] = np.einsum("gnp,gp->gn", x[pos], b_new[pos])
                        m_new[pos] = _linkinv(e_new[pos])
                        dev[pos] = _deviance_batch(s[pos], f[pos], m_new[pos])
                        pos = pos[~np.isfinite(dev[pos])]
                    boundary[idx[bad & alive]] = True
            active[idx[~alive]] = False
            sel = idx[alive]
            beta[sel], eta[sel], mu[sel] = b_new[alive], e_new[alive], m_new[alive]
            dev = dev[alive]
            iterations[sel] = iteration
            done = np.abs(dev - previous[sel]) / (np.abs(dev) + 0.1) < 1e-8
            converged[sel[done]] = True
            active[sel[done]] = False
            go = sel[~done]
            coefold[go] = beta[go]
            previous[go] = dev[~done]

        eps = 10 * _EPS
        edge = np.any(mu > 1 - eps, axis=1) | np.any(mu < eps, axis=1)
        good = np.array([e is None for e in errors], dtype=bool)
        se = np.full((g, p), np.nan)
        if good.any():
            inverse, ok = _invert_batch(normal[good])
            se[good] = np.sqrt(np.diagonal(inverse, axis1=1, axis2=2))
            for i in np.flatnonzero(good)[~ok]:
                errors[i] = "Singular matrix"
    beta[~np.array([e is None for e in errors])] = np.nan

    warnings: list[list[str]] = []
    for i in range(g):
        messages: list[str] = []
        if truncated[i]:
            messages.append("step size truncated due to divergence")
        if not converged[i]:
            messages.append("glm.fit: algorithm did not converge")
        if boundary[i]:
            messages.append("glm.fit: algorithm stopped at boundary value")
        if edge[i]:
            messages.append("glm.fit: fitted probabilities numerically 0 or 1 occurred")
        warnings.append(messages)
    with np.errstate(all="ignore"):
        pearson = _pearson_statistic(success, failure, weight_used, eta, mu)
    pearson[~np.array([e is None for e in errors], dtype=bool)] = np.nan
    return IrlsResult(
        beta,
        se,
        converged,
        iterations,
        warnings,
        errors,
        pearson,
        success.shape[-1] - p,
    )


_DIAGNOSTICS = pl.Struct(
    {
        "converged": pl.Boolean,
        "iterations": pl.Int64,
        "warning": pl.String,
        "error": pl.String,
    }
)


def _diagnostics(converged=None, iterations=None, warning=None, error=None) -> dict:
    return {
        "converged": converged,
        "iterations": iterations,
        "warning": warning,
        "error": error,
    }


def _pearson_statistic(success, failure, weight, eta, mu) -> np.ndarray:
    """Pearson statistic ``sum(w * r^2)`` as in R's ``summary.glm`` dispersion.

    ``r = (y - mu) / mu.eta(eta)`` are the working residuals at the final
    ``eta, mu`` and ``w = n * mu.eta^2 / var(mu)`` the working weights the
    final IRLS iteration started from (glm.fit keeps those stale weights, as it
    does for the standard errors). At exact convergence this is the Pearson
    chi-square ``sum(n * (y - mu)^2 / (mu * (1 - mu)))``; at the 1e-8 deviance
    stopping point it differs from it by ~1e-5 relative. Divided by the
    residual df it equals ``summary(glm(..., quasibinomial))$dispersion``.
    """
    observed = success / (success + failure)
    resid = (observed - mu) / _mu_eta(eta)
    return np.sum(weight * resid**2, axis=-1)


_PSEUDOCOUNT = "all_zeros"
_DISPERSION = "fixed"


def _check_pseudocount(pseudocount: str) -> None:
    if not isinstance(pseudocount, str) or pseudocount != _PSEUDOCOUNT:
        raise ValueError(
            f"pseudocount must be {_PSEUDOCOUNT!r} (every zero count is replaced "
            f"by 1); got {pseudocount!r}"
        )


def _check_dispersion(dispersion: str) -> None:
    if not isinstance(dispersion, str) or dispersion != _DISPERSION:
        raise ValueError(
            f"dispersion must be {_DISPERSION!r} (binomial dispersion fixed at 1); "
            f"got {dispersion!r}"
        )


def _flags(joined: pl.DataFrame) -> dict[str, list[str]]:
    """Data-quality flags per genotype, from the raw counts.

    Returns ``sequence_id -> sorted flag names`` (empty list when clean).
    ``censored``: every raw count at the later time is 0. Flags never alter any
    fit; they use ``count``, not the pseudocount-patched ``Nreads``.
    """
    # Validation enforces two time points and complete coverage, so the
    # per-group max time equals the global later time; take it once anyway.
    later_time = joined["time"].max()
    later = pl.col("time") == later_time
    totals = joined.group_by("sequence_id", maintain_order=True).agg(
        pl.col("count").filter(later).sum().alias("after"),
    )
    return {
        sequence_id: ["censored"] if after == 0 else []
        for sequence_id, after in totals.iter_rows()
    }


def _add_p_adj(rows: list[dict]) -> list[dict]:
    """Return new rows with ``p_adj``: BH over the non-null ``p`` of the call."""
    p = np.array(
        [np.nan if row["p"] is None else row["p"] for row in rows], dtype=float
    )
    adjusted = benjamini_hochberg(p)
    return [
        {**row, "p_adj": None if np.isnan(value) else float(value)}
        for row, value in zip(rows, adjusted, strict=True)
    ]


def _join_flags(*names: list[str]) -> str | None:
    """Union of flag lists as a sorted ``|``-joined string, ``None`` when empty."""
    union = sorted(set().union(*names))
    return "|".join(union) or None


def _validate_optional_sequence_columns(sequences: pl.DataFrame) -> None:
    """Check the optional ``sequence`` and ``aa_id`` columns; they never affect a fit."""
    if "sequence" in sequences.columns:
        column = sequences.get_column("sequence")
        if column.dtype != pl.String or column.null_count():
            raise ValueError("sequences.sequence must be non-null strings")
    if "aa_id" in sequences.columns:
        column = sequences.get_column("aa_id")
        if column.dtype != pl.String or column.null_count():
            raise ValueError("sequences.aa_id must be non-null strings")
        if (column == "").any():
            raise ValueError("sequences.aa_id must be nonempty")


def _validate(
    counts: pl.DataFrame,
    samples: pl.DataFrame,
    sequences: pl.DataFrame,
    elapsed_time: float,
    time_unit: str,
):
    """Validate inputs shared by every fit; return joined.

    Every zero count is replaced by 1 in ``Nreads`` (the published rule);
    ``Diff`` keeps the raw count.
    """
    counts = _as_frame(counts, "counts", {"sequence_id", "sample_id", "count"})
    samples = _as_frame(
        samples,
        "samples",
        {"sample_id", "replicate_id", "condition", "time", "sample_total"},
    )
    sequences = _as_frame(sequences, "sequences", {"sequence_id"})
    _validate_optional_sequence_columns(sequences)
    if (
        isinstance(elapsed_time, bool)
        or not isinstance(elapsed_time, int | float | np.integer | np.floating)
        or not math.isfinite(float(elapsed_time))
        or elapsed_time <= 0
    ):
        raise ValueError("elapsed_time must be a finite positive number")
    if not isinstance(time_unit, str) or not time_unit.strip():
        raise ValueError("time_unit must be a nonempty string")
    if counts.is_empty() or samples.is_empty() or sequences.is_empty():
        raise ValueError("counts, samples, and sequences must be nonempty")
    for frame, name, key in (
        (samples, "samples", "sample_id"),
        (sequences, "sequences", "sequence_id"),
    ):
        if frame.get_column(key).n_unique() != frame.height:
            raise ValueError(f"{name}.{key} must be unique")
    if (
        counts.select(pl.struct("sequence_id", "sample_id").n_unique()).item()
        != counts.height
    ):
        raise ValueError("counts must have one row per sequence_id and sample_id")

    # Reject broken foreign keys rather than silently losing records in joins.
    if counts.join(samples.select("sample_id"), on="sample_id", how="anti").height:
        raise ValueError("counts contains sample_id values absent from samples")
    if counts.join(
        sequences.select("sequence_id"), on="sequence_id", how="anti"
    ).height:
        raise ValueError("counts contains sequence_id values absent from sequences")
    expected_pairs = sequences.select("sequence_id").join(
        samples.select("sample_id"), how="cross"
    )
    missing_pairs = expected_pairs.join(
        counts.select("sequence_id", "sample_id"),
        on=["sequence_id", "sample_id"],
        how="anti",
    )
    if missing_pairs.height:
        example = missing_pairs.row(0, named=True)
        raise ValueError(
            "counts must provide complete sequence/sample coverage; "
            f"missing {example['sequence_id']!r} in sample {example['sample_id']!r}"
        )
    if samples.get_column("condition").n_unique() != 1:
        raise ValueError("fits currently require one condition")
    if not samples.schema["time"].is_numeric():
        raise ValueError("samples.time must be numeric")
    times = samples.get_column("time")
    if times.null_count() or not times.cast(pl.Float64).is_finite().all():
        raise ValueError("samples.time must be finite")
    time_values = times.unique().sort().to_list()
    if len(time_values) != 2:
        raise ValueError("fits require exactly two sample time points")
    span = float(time_values[-1]) - float(time_values[0])
    if abs(span - float(elapsed_time)) > 1e-9:
        raise ValueError(
            f"elapsed_time ({elapsed_time!r}) must equal max(time) - min(time) "
            f"({span!r}) within 1e-9"
        )
    _require_finite_integers(counts, "count", "count")
    _require_finite_integers(samples, "sample_total", "sample_total")
    if samples.filter(pl.col("sample_total") <= 0).height:
        raise ValueError("sample_total must be positive")
    if counts.filter((pl.col("count") < 0) | (pl.col("count") > 2**53)).height:
        raise ValueError("count must be a nonnegative integer")
    if samples.filter(pl.col("sample_total") > 2**53).height:
        raise ValueError("sample_total exceeds the supported exact integer range")

    joined = counts.join(samples, on="sample_id", how="left", validate="m:1")
    if joined.filter(pl.col("count") > pl.col("sample_total")).height:
        raise ValueError("count cannot exceed sample_total")
    # Keep the published Diff based on the raw count before zero replacement.
    joined = joined.with_columns(
        (pl.col("sample_total") - pl.col("count")).alias("Diff"),
        pl.when(pl.col("count") == 0)
        .then(1)
        .otherwise(pl.col("count"))
        .alias("Nreads"),
    )
    return joined


def _check_common(
    counts, samples, sequences, elapsed_time, time_unit, engine, chunk_size
):
    """Input checks shared by ``fit_wt`` and ``fit_pair``; returns the joined table."""
    if engine not in ("vectorised", "scalar"):
        raise ValueError("engine must be 'vectorised' or 'scalar'")
    if (
        isinstance(chunk_size, bool)
        or not isinstance(chunk_size, int | np.integer)
        or chunk_size < 1
    ):
        raise ValueError("chunk_size must be a positive integer")
    return _validate(counts, samples, sequences, elapsed_time, time_unit)


def fit_wt(
    counts: pl.DataFrame,
    samples: pl.DataFrame,
    sequences: pl.DataFrame,
    *,
    reference_id: str,
    elapsed_time: float,
    time_unit: str,
    pseudocount: str,
    dispersion: str,
    engine: str = "vectorised",
    chunk_size: int = 2000,
) -> pl.DataFrame:
    """Estimate each sequence's fitness difference from a wild-type reference.

    ``counts`` has one row per sequence and sample. ``samples`` supplies each
    sample's full denominator in ``sample_total``; ``sequences`` maps stable
    IDs to DNA. The interaction in the grouped-binomial type-by-time model is
    reported as ``m`` (per ``time_unit``, using ``samples.time`` directly) with
    its raw Wald standard error and two-sided p-value. ``elapsed_time`` must
    equal ``max(time) - min(time)`` within 1e-9. The reference row has
    ``m = 0``, null ``se``/``p`` and status ``"reference"``. Status is
    ``converged``, ``not_converged`` (estimate kept) or ``failed`` (error or
    non-finite estimate, null results). ``engine`` selects the batched
    ``"vectorised"`` IRLS (default) or the per-sequence ``"scalar"`` loop kept
    as the oracle; both give the same output. ``chunk_size`` is the number of
    sequences per vectorised batch (it does not change the output).

    ``pseudocount`` and ``dispersion`` are required keyword arguments with no
    default, so that the modelling choices are explicit in the calling code.
    The only accepted values are ``pseudocount="all_zeros"`` (every zero count
    is replaced by 1 at both time points before fitting, the published rule)
    and ``dispersion="fixed"`` (binomial dispersion fixed at 1, Wald z test).
    Any other value raises ``ValueError``.

    The ``flags`` column holds data-quality flags, sorted and ``|``-joined,
    null when none: ``censored`` (every raw later-time count of the genotype
    is 0). A comparison carries the union over both genotypes; the reference
    row carries the reference's own flags. Flags never change ``m``, ``se``,
    ``p`` or ``status``.

    ``p_adj`` (right after ``p``) is the Benjamini-Hochberg adjusted ``p`` over
    the non-null ``p`` of the call; null for the reference and failed rows.
    ``dispersion`` is the raw Pearson ``phi`` (Pearson chi-square over the
    residual degrees of freedom; a diagnostic, not used by ``se``) and
    ``df_residual`` the residual degrees of freedom; both are null for the
    reference row and failed rows.

    Returns one row per sequence of ``sequences``, in order, with columns
    ``sequence_id, reference_id, m, se, p, p_adj, status, dispersion,
    df_residual, flags, diagnostics, elapsed_time, time_unit, pseudocount,
    dispersion_model``. The last four echo the arguments: ``elapsed_time``,
    ``time_unit`` and the string columns ``pseudocount`` and
    ``dispersion_model`` (the value of ``dispersion``), repeated on every row.
    """
    _check_pseudocount(pseudocount)
    _check_dispersion(dispersion)
    joined = _check_common(
        counts, samples, sequences, elapsed_time, time_unit, engine, chunk_size
    )
    ids = sequences.get_column("sequence_id").to_list()
    if reference_id not in ids:
        raise ValueError(f"reference_id {reference_id!r} is absent from sequences")

    reference_flags = _flags(joined.filter(pl.col("sequence_id") == reference_id))[
        reference_id
    ]
    reference_result = {
        "m": 0.0,
        "se": None,
        "p": None,
        "status": "reference",
        "dispersion": None,
        "df_residual": None,
        "flags": _join_flags(reference_flags),
        "diagnostics": _diagnostics(),
    }
    tests = [value for value in ids if value != reference_id]
    results = _comparison.fit_comparisons(
        joined, samples, [reference_id] * len(tests), tests, engine, chunk_size
    )
    by_position = list(zip(tests, results, strict=True))
    by_position.insert(ids.index(reference_id), (reference_id, reference_result))
    adjusted = _add_p_adj([result for _, result in by_position])
    by_position = [
        (sequence_id, result)
        for (sequence_id, _), result in zip(by_position, adjusted, strict=True)
    ]
    rows = [
        {
            "sequence_id": sequence_id,
            "reference_id": reference_id,
            **result,
            "elapsed_time": float(elapsed_time),
            "time_unit": time_unit,
            "pseudocount": pseudocount,
            "dispersion_model": dispersion,
        }
        for sequence_id, result in by_position
    ]
    return pl.DataFrame(
        rows,
        schema={
            "sequence_id": pl.String,
            "reference_id": pl.String,
            "m": pl.Float64,
            "se": pl.Float64,
            "p": pl.Float64,
            "p_adj": pl.Float64,
            "status": pl.String,
            "dispersion": pl.Float64,
            "df_residual": pl.Int64,
            "flags": pl.String,
            "diagnostics": _DIAGNOSTICS,
            "elapsed_time": pl.Float64,
            "time_unit": pl.String,
            "pseudocount": pl.String,
            "dispersion_model": pl.String,
        },
    )


def fit_pair(
    counts: pl.DataFrame,
    samples: pl.DataFrame,
    sequences: pl.DataFrame,
    *,
    pairs: pl.DataFrame,
    elapsed_time: float,
    time_unit: str,
    pseudocount: str,
    dispersion: str,
    engine: str = "vectorised",
    chunk_size: int = 2000,
) -> pl.DataFrame:
    """Estimate pairwise fitness differences, one binomial GLM per pair.

    ``pairs`` has columns ``from`` and ``to`` (sequence IDs). ``from`` is the
    reference level and ``to`` the test level, so ``m`` is the fitness of
    ``to`` relative to ``from`` per ``time_unit``. Reversing a pair negates
    ``m`` and leaves ``se`` and ``p`` unchanged. One row per pair, in input
    order. ``engine`` selects the batched ``"vectorised"`` IRLS (default) or
    the per-pair ``"scalar"`` loop kept as the oracle; both give the same
    output. ``chunk_size`` is the number of pairs per vectorised batch (it does
    not change the output).

    ``pseudocount`` and ``dispersion`` are required keyword arguments with no
    default. The only accepted values are ``pseudocount="all_zeros"`` (every
    zero count is replaced by 1 at both time points before fitting, the
    published rule) and ``dispersion="fixed"`` (binomial dispersion fixed at 1,
    Wald z test). Any other value raises ``ValueError``.

    ``se`` and ``p`` are the raw Wald standard error and z-test p-value;
    ``p_adj`` is the Benjamini-Hochberg adjusted ``p`` over the non-null ``p``
    of the call, placed right after ``p``. ``dispersion`` is the raw Pearson
    ``phi`` (a diagnostic, not used by ``se``) and ``df_residual`` the residual
    degrees of freedom. ``flags`` is as in ``fit_wt`` (``censored``, the union
    over both genotypes).

    Returns one row per pair with columns ``from, to, m, se, p, p_adj, status,
    dispersion, df_residual, flags, diagnostics, elapsed_time, time_unit,
    pseudocount, dispersion_model``. The last four echo the arguments:
    ``elapsed_time``, ``time_unit`` and the string columns ``pseudocount`` and
    ``dispersion_model`` (the value of ``dispersion``), repeated on every row.
    """
    _check_pseudocount(pseudocount)
    _check_dispersion(dispersion)
    joined = _check_common(
        counts, samples, sequences, elapsed_time, time_unit, engine, chunk_size
    )
    pairs = _as_frame(pairs, "pairs", {"from", "to"})
    if pairs.is_empty():
        raise ValueError("pairs must be nonempty")
    sources = pairs.get_column("from").to_list()
    targets = pairs.get_column("to").to_list()
    if any(a == b for a, b in zip(sources, targets, strict=True)):
        raise ValueError("pairs must not compare a sequence with itself (from == to)")
    if len(set(zip(sources, targets, strict=True))) != len(sources):
        raise ValueError("pairs must not contain duplicate (from, to) rows")
    known = set(sequences.get_column("sequence_id").to_list())
    if any(value not in known for value in sources + targets):
        raise ValueError("pairs contains sequence_id values absent from sequences")

    fitted = _add_p_adj(
        _comparison.fit_comparisons(
            joined, samples, sources, targets, engine, chunk_size
        )
    )
    rows = [
        {
            "from": source,
            "to": target,
            **result,
            "elapsed_time": float(elapsed_time),
            "time_unit": time_unit,
            "pseudocount": pseudocount,
            "dispersion_model": dispersion,
        }
        for source, target, result in zip(sources, targets, fitted, strict=True)
    ]
    return pl.DataFrame(
        rows,
        schema={
            "from": pl.String,
            "to": pl.String,
            "m": pl.Float64,
            "se": pl.Float64,
            "p": pl.Float64,
            "p_adj": pl.Float64,
            "status": pl.String,
            "dispersion": pl.Float64,
            "df_residual": pl.Int64,
            "flags": pl.String,
            "diagnostics": _DIAGNOSTICS,
            "elapsed_time": pl.Float64,
            "time_unit": pl.String,
            "pseudocount": pl.String,
            "dispersion_model": pl.String,
        },
    )
