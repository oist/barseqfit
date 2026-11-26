"""The ``fit_sequences`` / ``fit_wt`` identity as a pure comparison.

``fit_wt`` uses a type-by-time design that is block-separable, so at convergence
``m == r_i - r_ref`` and, under fixed dispersion,
``se == sqrt(se_pooled_i**2 + se_pooled_ref**2)``, where ``r`` and ``se_pooled``
come from ``fit_sequences``. This module takes the two fitted tables and an
``IdentityRule`` built by the caller; it performs no I/O and reads no config.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass, field

import polars as pl


@dataclass(frozen=True)
class IdentityRule:
    """Absolute tolerances on ``m`` and ``se`` for the identity."""

    m_tol: float
    se_tol: float

    @classmethod
    def from_run_config(cls, config: Mapping) -> IdentityRule:
        identity = config["identity"]
        return cls(m_tol=identity["m_tol"], se_tol=identity["se_tol"])


@dataclass(frozen=True)
class IdentityVerdict:
    checked: int
    skipped: int
    max_m: float
    max_se: float
    outside_m: list[str] = field(default_factory=list)
    outside_se: list[str] = field(default_factory=list)
    failures: list[str] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        return not self.failures

    def lines(self, rule: IdentityRule) -> list[str]:
        """Markdown bullet lines for a verdict report."""
        return [
            f"- Rows checked: {self.checked}; skipped (either fit not converged): {self.skipped}",
            f"- max |(r_i - r_ref) - m|: {self.max_m:.3g} (tolerance {rule.m_tol:g}); rows outside: {len(self.outside_m)}",
            f"- max |sqrt(se_pooled_i^2 + se_pooled_ref^2) - se|: {self.max_se:.3g} (tolerance {rule.se_tol:g}); rows outside: {len(self.outside_se)}",
        ]


def sequences_identity(
    seq: pl.DataFrame, wt: pl.DataFrame, reference_id: str, rule: IdentityRule
) -> IdentityVerdict:
    """Compare ``fit_wt`` output ``wt`` with the ``fit_sequences`` output ``seq``.

    Rows where either fit has a status other than ``converged`` are counted as
    skipped. An id set that differs between the two tables, or a reference
    missing from ``seq``, is a failure.
    """
    failures: list[str] = []
    ref_rows = seq.filter(pl.col("sequence_id") == reference_id)
    if ref_rows.height != 1:
        failures.append(
            f"reference {reference_id} not found in the fit_sequences table"
        )
        return IdentityVerdict(0, 0, math.nan, math.nan, failures=failures)
    ref = ref_rows.row(0, named=True)
    seq_ids = set(seq["sequence_id"].to_list()) - {reference_id}
    wt_ids = set(wt["sequence_id"].to_list()) - {reference_id}
    if seq_ids != wt_ids:
        failures.append(
            f"id sets differ ({len(seq_ids - wt_ids)} only in fit_sequences, "
            f"{len(wt_ids - seq_ids)} only in fit_wt)"
        )
    both = seq.join(wt, on="sequence_id", suffix="_wt").filter(
        pl.col("sequence_id") != reference_id
    )
    ok = (pl.col("status") == "converged") & (pl.col("status_wt") == "converged")
    use = both.filter(ok).with_columns(
        d_m=((pl.col("r") - ref["r"]) - pl.col("m")).abs(),
        d_se=(
            ((pl.col("se_pooled") ** 2 + ref["se_pooled"] ** 2).sqrt()) - pl.col("se")
        ).abs(),
    )
    outside_m = use.filter(pl.col("d_m") >= rule.m_tol)["sequence_id"].to_list()
    outside_se = use.filter(pl.col("d_se") >= rule.se_tol)["sequence_id"].to_list()
    if outside_m:
        failures.append(f"{len(outside_m)} rows with m outside {rule.m_tol:g}")
    if outside_se:
        failures.append(f"{len(outside_se)} rows with se outside {rule.se_tol:g}")
    return IdentityVerdict(
        checked=use.height,
        skipped=both.height - use.height,
        max_m=float(use["d_m"].max()) if use.height else math.nan,
        max_se=float(use["d_se"].max()) if use.height else math.nan,
        outside_m=outside_m,
        outside_se=outside_se,
        failures=failures,
    )
