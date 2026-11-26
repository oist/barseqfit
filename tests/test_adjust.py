"""Benjamini-Hochberg adjusted p."""

import math

import numpy as np
import polars as pl
import pytest

import barseqfit as bsf
from barseqfit import _fit as fit_module
from barseqfit._adjust import benjamini_hochberg

from ._papkou import REFERENCE_ID, TIME_UNIT, load
from .test_fit_wt import ENGINES, _random_inputs

# Generated with R 4.1.3:
#   set.seed(42); p <- round(runif(20)^2, 6); dput(p); dput(p.adjust(p, "BH"))
R_P = [
    0.83687, 0.87811, 0.081876, 0.689643, 0.411837, 0.269461, 0.542562,
    0.018135, 0.431639, 0.497116, 0.209528, 0.517122, 0.873612, 0.065244,
    0.213715, 0.883627, 0.956927, 0.013803, 0.225622, 0.313973,
]  # fmt: skip
R_P_ADJ = [
    0.930133684210526, 0.930133684210526, 0.40938, 0.919524, 0.775088571428571,
    0.6736525, 0.775088571428571, 0.18135, 0.775088571428571, 0.775088571428571,
    0.644634285714286, 0.775088571428571, 0.930133684210526, 0.40938,
    0.644634285714286, 0.930133684210526, 0.956927, 0.18135, 0.644634285714286,
    0.697717777777778,
]  # fmt: skip


def test_hand_values():
    # sorted 0.01, 0.03, 0.04, 0.2 with n = 4: 0.04, 0.06, 0.0533.., 0.2,
    # then the cumulative minimum from the largest p gives 0.04, 0.0533.., 0.0533.., 0.2.
    got = benjamini_hochberg(np.array([0.01, 0.04, 0.03, 0.2]))
    np.testing.assert_allclose(got, [0.04, 0.16 / 3, 0.16 / 3, 0.2], rtol=1e-15)


def test_nan_passthrough_uses_non_nan_count():
    got = benjamini_hochberg(np.array([0.01, np.nan, 0.04, 0.03, 0.2]))
    assert math.isnan(got[1])
    # n = 4 (the NaN is not counted), same as without it.
    np.testing.assert_allclose(
        got[[0, 2, 3, 4]], [0.04, 0.16 / 3, 0.16 / 3, 0.2], rtol=1e-15
    )


def test_empty_and_all_nan():
    assert benjamini_hochberg(np.array([])).shape == (0,)
    assert np.isnan(benjamini_hochberg(np.array([np.nan, np.nan]))).all()


def test_ties():
    got = benjamini_hochberg(np.array([0.02, 0.02, 0.02, 0.5]))
    np.testing.assert_allclose(got, [0.02 * 4 / 3] * 3 + [0.5], rtol=1e-15)


def test_monotone_and_clipped_at_one():
    rng = np.random.default_rng(0)
    p = rng.uniform(size=500) ** 0.5  # skewed towards 1
    got = benjamini_hochberg(p)
    assert got.max() <= 1.0
    assert got.max() == p.max()  # the largest p is its own adjusted value
    assert (got >= p).all()
    order = np.argsort(p, kind="stable")
    assert (np.diff(got[order]) >= 0).all()
    np.testing.assert_array_equal(
        benjamini_hochberg(np.array([0.9, 0.95])), [0.95, 0.95]
    )


def test_matches_r_p_adjust():
    got = benjamini_hochberg(np.array(R_P))
    np.testing.assert_allclose(got, R_P_ADJ, rtol=1e-15, atol=0)


@pytest.mark.parametrize("engine", ENGINES)
def test_p_adj_wt_equals_pair_on_reference_pairs(engine):
    inputs = _random_inputs(0, genotypes=30, total=1000, high=500, edges=False)
    ids = inputs[2]["sequence_id"].to_list()
    tests = [value for value in ids if value != "wt"]
    options = {"elapsed_time": 1.0, "time_unit": "g", "engine": engine}
    wt = bsf.fit_wt(
        *inputs,
        reference_id="wt",
        **options,
        pseudocount="all_zeros",
        dispersion="fixed",
    )
    pair = bsf.fit_pair(
        *inputs,
        pairs=pl.DataFrame({"from": ["wt"] * len(tests), "to": tests}),
        **options,
        pseudocount="all_zeros",
        dispersion="fixed",
    )
    wt = wt.filter(pl.col("sequence_id") != "wt")
    assert wt["p_adj"].to_list() == pair["p_adj"].to_list()
    assert wt["p_adj"].null_count() == 0


@pytest.mark.parametrize("engine", ENGINES)
def test_p_adj_column_position_and_values(engine):
    inputs = _random_inputs(0, genotypes=30, total=1000, high=500, edges=False)
    wt = bsf.fit_wt(
        *inputs,
        reference_id="wt",
        elapsed_time=1.0,
        time_unit="g",
        engine=engine,
        pseudocount="all_zeros",
        dispersion="fixed",
    )
    assert wt.columns[wt.columns.index("p") + 1] == "p_adj"
    assert wt.schema["p_adj"] == pl.Float64
    ref = wt.filter(pl.col("status") == "reference")
    assert ref["p_adj"].null_count() == 1
    rest = wt.filter(pl.col("status") != "reference")
    expected = benjamini_hochberg(rest["p"].to_numpy())
    np.testing.assert_allclose(rest["p_adj"].to_numpy(), expected, rtol=1e-15)


@pytest.mark.parametrize("engine", ENGINES)
def test_p_adj_null_on_failed_rows(monkeypatch, engine):
    counts, samples, sequences = _random_inputs(
        3, genotypes=6, total=400, high=300, edges=False
    )
    options = {
        "reference_id": "wt",
        "elapsed_time": 1.0,
        "time_unit": "g",
        "engine": engine,
        "pseudocount": "all_zeros",
        "dispersion": "fixed",
    }
    real_scalar = fit_module._logit_fit
    real_batch = fit_module._logit_fit_batch

    # Record each genotype's response in fit order, then fail two of them.
    seen: list[tuple] = []

    def spy(success, failure, design):
        seen.append(tuple(success))
        return real_scalar(success, failure, design)

    monkeypatch.setattr(fit_module, "_logit_fit", spy)
    bsf.fit_wt(counts, samples, sequences, **{**options, "engine": "scalar"})
    assert len(seen) == 6
    poisoned = {seen[1], seen[4]}

    def broken_scalar(success, failure, design):
        if tuple(success) in poisoned:
            raise ValueError("rank-deficient design")
        return real_scalar(success, failure, design)

    def broken_batch(success, failure, design):
        fit = real_batch(success, failure, design)
        errors = [
            "rank-deficient design" if tuple(row) in poisoned else None
            for row in success
        ]
        return fit._replace(errors=errors)

    monkeypatch.setattr(fit_module, "_logit_fit", broken_scalar)
    monkeypatch.setattr(fit_module, "_logit_fit_batch", broken_batch)
    out = bsf.fit_wt(counts, samples, sequences, **options)
    failed = out.filter(pl.col("status") == "failed")
    assert failed.height == 2
    assert failed["p_adj"].null_count() == failed.height
    assert out.filter(pl.col("status") == "reference")["p_adj"][0] is None
    rest = out.filter(~pl.col("status").is_in(["failed", "reference"]))
    assert rest.height == 4
    np.testing.assert_allclose(
        rest["p_adj"].to_numpy(),
        benjamini_hochberg(rest["p"].to_numpy()),
        rtol=1e-15,
    )


@pytest.mark.parametrize("engine", ENGINES)
def test_fixture_p_adj_matches_r(data_dir, engine):
    d = load(data_dir)
    wt = bsf.fit_wt(
        d["counts"],
        d["samples"],
        d["sequences"],
        reference_id=REFERENCE_ID,
        elapsed_time=d["elapsed"],
        time_unit=TIME_UNIT,
        engine=engine,
        pseudocount="all_zeros",
        dispersion="fixed",
    )
    exp = pl.read_csv(
        data_dir / "expected_wt_padj_r.csv", schema_overrides={"p_adj": pl.Float64}
    )
    assert wt["sequence_id"].to_list() == exp["sequence_id"].to_list()
    _assert_close(wt["p_adj"].to_list(), exp["p_adj"].to_list())
    pair = bsf.fit_pair(
        d["counts"],
        d["samples"],
        d["sequences"],
        pairs=d["pairs"],
        elapsed_time=d["elapsed"],
        time_unit=TIME_UNIT,
        engine=engine,
        pseudocount="all_zeros",
        dispersion="fixed",
    )
    exp = pl.read_csv(
        data_dir / "expected_pair_padj_r.csv", schema_overrides={"p_adj": pl.Float64}
    )
    assert pair["to"].to_list() == exp["to"].to_list()
    _assert_close(pair["p_adj"].to_list(), exp["p_adj"].to_list())


def _assert_close(got, expected):
    for a, b in zip(got, expected, strict=True):
        if b is None:
            assert a is None
        else:
            # p are tiny (1e-134): relative, not absolute. Python and R agree on
            # the fixture p to about 1e-12 relative (the 1e-9 engine rule).
            assert a == pytest.approx(b, rel=1e-9, abs=0)


def test_defaults_leave_other_columns_unchanged(data_dir):
    d = load(data_dir)
    out = bsf.fit_pair(
        d["counts"],
        d["samples"],
        d["sequences"],
        pairs=d["pairs"],
        elapsed_time=d["elapsed"],
        time_unit=TIME_UNIT,
        pseudocount="all_zeros",
        dispersion="fixed",
    )
    exp = pl.read_csv(data_dir / "expected_pair.csv")
    assert out["from"].to_list() == exp["from"].to_list()
    for key in ("m", "se", "p"):
        _assert_close(out[key].to_list(), exp[key].to_list())
