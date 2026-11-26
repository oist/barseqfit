"""Unit tests for the Papkou2023 comparison rule."""

from __future__ import annotations

import numpy as np
import pytest

from validation.compare import Tolerance, assert_close, outside_tolerance

FIXTURE = Tolerance.fixture()


def _bad(actual: float, expected: float) -> bool:
    return bool(
        outside_tolerance([actual], [expected], noise_floor=FIXTURE.noise_floor)[0]
    )


def test_noise_floor_passes_rounding_noise():
    assert not _bad(6e-12 + 1e-16, 6e-12)


def test_diff_above_floor_fails_at_tiny_expected():
    assert _bad(6e-12 + 5e-13, 6e-12)


def test_relative_rule_intact_above_floor():
    assert _bad(1e-3 + 1e-7, 1e-3)


def test_p_like_values_fail_without_noise_floor():
    # p is relative-rounded (norm.sf), so the floor must not apply to it.
    assert bool(outside_tolerance([5e-16], [1e-20])[0])  # no floor by default
    assert not outside_tolerance([5e-16], [1e-20], noise_floor=FIXTURE.noise_floor)[0]


def test_floor_applies_to_m_only():
    assert_close(6e-12 + 1e-16, 6e-12, "m floor", column="m")
    with pytest.raises(AssertionError):
        assert_close(6e-12 + 1e-16, 6e-12, "se no floor", column="se")


def test_floor_for():
    assert FIXTURE.floor_for("m") == FIXTURE.noise_floor
    assert FIXTURE.floor_for("se") == 0.0
    assert FIXTURE.floor_for("p") == 0.0
    assert FIXTURE.floor_for(None) == 0.0


def test_p_underflow_is_equal_but_close_values_still_fail():
    assert_close(1e-310, 0.0, "p underflow", column="p")
    with pytest.raises(AssertionError):
        assert_close(5e-16, 1e-20, "p tiny", column="p")


def test_underflow_only_flagged_for_p():
    res = FIXTURE.compare([1e-310], [0.0], "p")
    assert res.underflow.tolist() == [True]
    assert not res.outside.any()
    assert not FIXTURE.compare([1e-310], [0.0], "se").underflow.any()
    assert not FIXTURE.compare([0.0], [0.0], "p").underflow.any()


def test_both_nan_behaviour():
    nan = float("nan")
    with pytest.raises(AssertionError):
        assert_close(nan, nan, "nan strict")
    assert_close(nan, nan, "nan allowed", allow_both_nan=True)
    assert FIXTURE.compare([nan], [nan]).outside.tolist() == [False]
    assert FIXTURE.compare([nan], [nan], allow_both_nan=False).outside.tolist() == [
        True
    ]
    with pytest.raises(AssertionError):
        assert_close(nan, 1.0, "one nan", allow_both_nan=True)


def test_compare_residuals():
    res = FIXTURE.compare([1.0, 0.5], [1.0, 0.0], "se")
    assert res.abs_residual.tolist() == [0.0, 0.5]
    assert res.rel_residual[0] == 0.0 and np.isinf(res.rel_residual[1])


def test_from_run_config_round_trip():
    config = {"gate": {"abs_tol": 1e-5, "rel_tol": 2e-6, "noise_floor": 3e-12}}
    tol = Tolerance.from_run_config(config)
    assert tol == Tolerance(1e-5, 2e-6, 3e-12)
    assert tol.p_underflow == 1e-300
    assert (
        Tolerance.from_run_config(
            {"gate": {"abs_tol": 1e-6, "rel_tol": 1e-6, "noise_floor": 1e-13}}
        )
        == Tolerance.fixture()
    )
