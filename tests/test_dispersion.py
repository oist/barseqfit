"""Dispersion diagnostics of fit_wt and fit_pair (fixed dispersion)."""

import math

import numpy as np
import polars as pl
import pytest
from scipy.stats import norm

import barseqfit as bsf

from ._papkou import REFERENCE_ID, TIME_UNIT, load
from .test_fit_wt import ENGINES, _random_inputs

DF = 20  # 24 rows - 4 coefficients


def _wt(d, engine="vectorised", **kwargs):
    return bsf.fit_wt(
        d["counts"],
        d["samples"],
        d["sequences"],
        reference_id=REFERENCE_ID,
        elapsed_time=d["elapsed"],
        time_unit=TIME_UNIT,
        engine=engine,
        **kwargs,
        pseudocount="all_zeros",
        dispersion="fixed",
    )


def _pair(d, engine="vectorised", **kwargs):
    return bsf.fit_pair(
        d["counts"],
        d["samples"],
        d["sequences"],
        pairs=d["pairs"],
        elapsed_time=d["elapsed"],
        time_unit=TIME_UNIT,
        engine=engine,
        **kwargs,
        pseudocount="all_zeros",
        dispersion="fixed",
    )


@pytest.mark.parametrize("engine", ENGINES)
def test_dispersion_column_reports_raw_pearson_phi(data_dir, engine):
    d = load(data_dir)
    default = _wt(d, engine)
    tests = default.filter(pl.col("sequence_id") != REFERENCE_ID)
    # The raw phi is reported whatever the mode; fixed leaves se and p alone.
    assert tests["dispersion"].null_count() == 0
    assert set(tests["df_residual"]) == {DF}
    ref = default.filter(pl.col("sequence_id") == REFERENCE_ID).row(0, named=True)
    assert ref["dispersion"] is None and ref["df_residual"] is None
    assert default.columns[-9:-5] == ["status", "dispersion", "df_residual", "flags"]
    assert default.schema["dispersion"] == pl.Float64
    assert default.schema["df_residual"] == pl.Int64


@pytest.mark.parametrize("fit", [_wt, _pair])
def test_engines_agree(data_dir, fit):
    d = load(data_dir)
    fast, slow = fit(d, "vectorised"), fit(d, "scalar")
    for key in ("m", "se", "p", "dispersion"):
        a, b = fast[key].drop_nulls().to_numpy(), slow[key].drop_nulls().to_numpy()
        np.testing.assert_allclose(a, b, rtol=1e-9, atol=0)
    assert fast["df_residual"].to_list() == slow["df_residual"].to_list()


def test_engines_agree_on_random_inputs():
    inputs = _random_inputs(0, genotypes=40, total=1000, high=500, edges=False)
    kwargs = {
        "reference_id": "wt",
        "elapsed_time": 1.0,
        "time_unit": "g",
        "pseudocount": "all_zeros",
        "dispersion": "fixed",
    }
    out = {engine: bsf.fit_wt(*inputs, engine=engine, **kwargs) for engine in ENGINES}
    fast, slow = out["vectorised"], out["scalar"]
    assert fast["status"].to_list() == slow["status"].to_list()
    assert fast["dispersion"].drop_nulls().len() == 40
    for key in ("m", "se", "p"):
        a = fast[key].drop_nulls().to_numpy()
        b = slow[key].drop_nulls().to_numpy()
        np.testing.assert_allclose(a, b, rtol=1e-9, atol=1e-15)


def test_se_and_p_are_raw_wald_whatever_the_dispersion_estimate(data_dir):
    d = load(data_dir)
    for row in _pair(d).iter_rows(named=True):
        assert row["dispersion"] > 1  # over-dispersed fixture, yet se is not scaled
        assert row["df_residual"] == DF
        assert row["p"] == pytest.approx(
            2 * norm.sf(abs(row["m"] / row["se"])), abs=1e-12, rel=1e-12
        )


def _saturated_inputs():
    """One replicate per time point: 4 rows per comparison, zero residual df."""
    samples = pl.DataFrame(
        {
            "sample_id": ["b", "a"],
            "replicate_id": ["r1", "r1"],
            "condition": ["c", "c"],
            "time": [0, 2],
            "sample_total": [1000, 1000],
        }
    )
    table = {"wt": (100, 150), "x": (80, 260)}
    counts = pl.DataFrame(
        [
            {"sequence_id": g, "sample_id": s, "count": c}
            for g, row in table.items()
            for s, c in zip(["b", "a"], row, strict=True)
        ]
    )
    sequences = pl.DataFrame({"sequence_id": list(table), "sequence": ["AAAA"] * 2})
    return counts, samples, sequences, table


def _closed_form(table):
    def log_odds(n):
        return math.log(n / (1000 - n))

    m = (
        (log_odds(table["x"][1]) - log_odds(table["x"][0]))
        - (log_odds(table["wt"][1]) - log_odds(table["wt"][0]))
    ) / 2.0
    se = (
        math.sqrt(sum(1 / n + 1 / (1000 - n) for row in table.values() for n in row))
        / 2.0
    )
    return m, se


@pytest.mark.parametrize("engine", ENGINES)
def test_zero_residual_df_keeps_fixed_fit(engine):
    counts, samples, sequences, table = _saturated_inputs()
    out = bsf.fit_wt(
        counts,
        samples,
        sequences,
        reference_id="wt",
        elapsed_time=2.0,
        time_unit="g",
        engine=engine,
        pseudocount="all_zeros",
        dispersion="fixed",
    ).filter(pl.col("sequence_id") == "x")
    m, se = _closed_form(table)
    row = out.row(0, named=True)
    assert row["status"] == "converged"
    assert row["m"] == pytest.approx(m, abs=1e-8)
    assert row["se"] == pytest.approx(se, abs=1e-8)
    assert row["p"] == pytest.approx(2 * norm.sf(abs(m / se)), abs=1e-8)
    assert row["dispersion"] is None
    assert row["df_residual"] == 0
