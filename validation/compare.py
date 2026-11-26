"""Private tolerance rule shared by the Papkou2023 tests and the full-data gate.

``Tolerance`` is the single owner of the comparison rule (absolute,
relative, ``m``-only noise floor, ``p`` underflow). Fixture tests use
``Tolerance.fixture()``; the gates build the same object from the Run config
with ``Tolerance.from_run_config()``. Nothing here reads ``validation.toml``;
values arrive as arguments.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import NamedTuple

import numpy as np


class Residuals(NamedTuple):
    """Per-entry outcome of ``Tolerance.compare``."""

    abs_residual: np.ndarray
    rel_residual: np.ndarray
    outside: np.ndarray
    underflow: np.ndarray


@dataclass(frozen=True)
class Tolerance:
    """The comparison rule as one value."""

    abs_tol: float
    rel_tol: float
    noise_floor: float
    # p below this on both sides is IEEE underflow territory (scipy returns 0 or
    # a denormal where R's pnorm returns another denormal); treated as equal.
    p_underflow: float = 1e-300

    @classmethod
    def fixture(cls) -> Tolerance:
        """The defaults used by the fixture tests."""
        return cls(abs_tol=1e-6, rel_tol=1e-6, noise_floor=1e-13)

    @classmethod
    def from_run_config(cls, config: Mapping) -> Tolerance:
        """Build the rule from a parsed Run config (``[gate]`` section)."""
        gate = config["gate"]
        return cls(
            abs_tol=gate["abs_tol"],
            rel_tol=gate["rel_tol"],
            noise_floor=gate["noise_floor"],
        )

    def floor_for(self, column: str | None) -> float:
        """The noise floor applies to ``m`` only; se and p get no floor."""
        return self.noise_floor if column == "m" else 0.0

    def compare(
        self,
        actual,
        expected,
        column: str | None = None,
        *,
        allow_both_nan: bool = True,
    ) -> Residuals:
        """Residuals and out-of-tolerance mask; nulls must already be NaN."""
        a = np.asarray(actual, dtype=float)
        e = np.asarray(expected, dtype=float)
        with np.errstate(invalid="ignore", divide="ignore"):
            abs_residual = np.abs(a - e)
            rel_residual = abs_residual / np.abs(e)
        outside = outside_tolerance(
            a,
            e,
            allow_both_nan=allow_both_nan,
            abs_tol=self.abs_tol,
            rel_tol=self.rel_tol,
            noise_floor=self.floor_for(column),
        )
        underflow = np.zeros(a.shape, dtype=bool)
        if column == "p":
            with np.errstate(invalid="ignore"):
                both_underflow = (a < self.p_underflow) & (e < self.p_underflow)
            outside = outside & ~both_underflow
            underflow = both_underflow & (a != e)
        return Residuals(abs_residual, rel_residual, outside, underflow)


_FIXTURE = Tolerance.fixture()


def outside_tolerance(
    actual,
    expected,
    *,
    allow_both_nan: bool = True,
    abs_tol: float = _FIXTURE.abs_tol,
    rel_tol: float = _FIXTURE.rel_tol,
    noise_floor: float = 0.0,
) -> np.ndarray:
    """Boolean mask of entries that fail the Papkou2023 comparison rule.

    Absolute 1e-6 always; additionally relative 1e-6 when ``|expected| < 1``
    (small values such as p = 4e-6 pass the absolute check too easily).
    The relative bound is floored at ``noise_floor`` (default 0: no floor; only
    ``m`` takes the floor, see ``Tolerance.floor_for``): a difference below it is
    indistinguishable from double rounding for O(1) quantities (m = -6.2257e-12
    vs -6.2255e-12 differ by 1.3e-16 yet the bare relative bound is 6e-18), so
    the relative rule only bites above the floor.
    Non-finite entries fail unless both are NaN (when ``allow_both_nan``) or equal.
    """
    a = np.asarray(actual, dtype=float)
    e = np.asarray(expected, dtype=float)
    with np.errstate(invalid="ignore"):
        diff = np.abs(a - e)
        bad = ~(diff <= abs_tol)
        bad |= (np.abs(e) < 1) & ~(diff <= np.maximum(rel_tol * np.abs(e), noise_floor))
    ok = a == e
    if allow_both_nan:
        ok = ok | (np.isnan(a) & np.isnan(e))
    return bad & ~ok


def assert_close(
    actual,
    expected,
    label: str,
    *,
    column: str | None = None,
    allow_both_nan: bool = False,
    tolerance: Tolerance = _FIXTURE,
) -> None:
    """Assert agreement under the Papkou2023 rule; both-NaN fails unless allowed.

    The noise floor applies only when ``column == "m"``; se and p are
    relative-rounded quantities and get no floor. ``column == "p"`` also treats
    both-sides underflow as equal.
    """
    residuals = tolerance.compare(
        actual, expected, column, allow_both_nan=allow_both_nan
    )
    assert not residuals.outside.any(), label
