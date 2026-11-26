import math

import polars as pl
import pytest
from polars.testing import assert_frame_equal

import barseqfit as bsf

from .test_fit_wt import (
    ENGINES,
    _close,
    _count_logit_fit,
    _random_inputs,
    _tiny_inputs,
)


def _pairs(inputs, count=60):
    ids = inputs[2]["sequence_id"].to_list()
    rows = [(ids[i], ids[(i * 7 + 3) % len(ids)]) for i in range(len(ids))]
    rows += [(ids[i], ids[(i + 1) % len(ids)]) for i in range(len(ids))]
    rows = [r for r in dict.fromkeys(rows) if r[0] != r[1]][:count]
    return pl.DataFrame({"from": [a for a, _ in rows], "to": [b for _, b in rows]})


def _fit(inputs, pairs, engine="vectorised", **kwargs):
    return bsf.fit_pair(
        *inputs,
        pairs=pairs,
        elapsed_time=1.0,
        time_unit="generations",
        engine=engine,
        **{"pseudocount": "all_zeros", "dispersion": "fixed", **kwargs},
    )


def _assert_engines_agree(inputs, pairs, select=lambda row: True):
    fast_rows = _fit(inputs, pairs, "vectorised").iter_rows(named=True)
    slow_rows = _fit(inputs, pairs, "scalar").iter_rows(named=True)
    checked = 0
    for fast, slow in zip(fast_rows, slow_rows, strict=True):
        if not select(slow):
            continue
        checked += 1
        label = (slow["from"], slow["to"])
        assert (fast["from"], fast["to"]) == label
        assert fast["status"] == slow["status"], label
        assert fast["diagnostics"] == slow["diagnostics"], label
        for key in ("m", "se", "p"):
            assert _close(fast[key], slow[key]), f"{label} {key}"
    assert checked


@pytest.mark.parametrize("engine", ENGINES)
def test_fit_pair_matches_analytic_contrast(engine):
    pairs = pl.DataFrame({"from": ["wt"], "to": ["variant"]})
    row = _fit(_tiny_inputs(), pairs, engine).row(0, named=True)
    assert row["status"] == "converged"
    assert abs(row["m"] - math.log(3)) < 1e-8


def test_fit_pair_rejects_unknown_engine():
    pairs = pl.DataFrame({"from": ["wt"], "to": ["variant"]})
    with pytest.raises(ValueError, match="engine"):
        _fit(_tiny_inputs(), pairs, "gpu")


@pytest.mark.parametrize(("total", "high"), [(50, 20), (1000, 500)])
def test_engines_agree_on_well_conditioned_data(total, high, monkeypatch):
    inputs = _random_inputs(0, genotypes=40, total=total, high=high, edges=False)
    pairs = _pairs(inputs)
    scalar = _fit(inputs, pairs, "scalar")
    calls = _count_logit_fit(monkeypatch)
    _fit(inputs, pairs, "vectorised")
    assert calls == []
    assert (scalar["status"] == "converged").all()
    assert scalar["diagnostics"].struct.field("warning").is_null().all()
    _assert_engines_agree(inputs, pairs)


def test_vectorised_is_independent_of_input_row_order(monkeypatch):
    inputs = _random_inputs(1, genotypes=300, total=1000, high=500, edges=False)
    pairs = _pairs(inputs, count=400)
    counts, samples, sequences = (
        frame.sample(fraction=1.0, shuffle=True, seed=seed)
        for frame, seed in zip(inputs, (11, 12, 13), strict=True)
    )
    calls = _count_logit_fit(monkeypatch)
    shuffled = _fit((counts, samples, sequences), pairs, "vectorised")
    assert calls == []
    monkeypatch.undo()
    expected = _fit(inputs, pairs, "scalar")
    assert shuffled["from"].to_list() == expected["from"].to_list()
    assert shuffled["to"].to_list() == expected["to"].to_list()
    for fast, slow in zip(
        shuffled.iter_rows(named=True), expected.iter_rows(named=True), strict=True
    ):
        label = (slow["from"], slow["to"])
        assert fast["status"] == slow["status"], label
        for key in ("m", "se", "p"):
            assert _close(fast[key], slow[key]), f"{label} {key}"


def test_engines_agree_on_edge_data():
    inputs = _random_inputs(0, genotypes=40, total=1000, high=10)
    _assert_engines_agree(inputs, _pairs(inputs))


def test_engines_agree_when_not_converging():
    inputs = _random_inputs(0, genotypes=40, total=1000, high=10)
    _assert_engines_agree(
        inputs, _pairs(inputs), lambda row: row["status"] == "not_converged"
    )


def test_engines_both_flag_non_convergence():
    inputs = _random_inputs(0, genotypes=40, total=1000, high=10)
    pairs = _pairs(inputs)
    for engine in ENGINES:
        assert (_fit(inputs, pairs, engine)["status"] == "not_converged").any(), engine


def test_vectorised_falls_back_to_scalar_fit_on_separation_data(monkeypatch):
    inputs = _random_inputs(0, genotypes=40, total=1000, high=10)
    pairs = _pairs(inputs)
    expected = _fit(inputs, pairs, "scalar")
    unclean = expected.filter(
        (pl.col("status") != "converged")
        | pl.col("diagnostics").struct.field("warning").is_not_null()
    ).height
    assert unclean

    calls = _count_logit_fit(monkeypatch)
    _fit(inputs, pairs, "vectorised")
    assert len(calls) == unclean


@pytest.mark.parametrize("engine", ENGINES)
def test_reversing_pairs_negates_m(engine):
    inputs = _random_inputs(2, genotypes=20, total=1000, high=500, edges=False)
    pairs = _pairs(inputs, count=25)
    reverse = pairs.select(pl.col("to").alias("from"), pl.col("from").alias("to"))
    forward, backward = _fit(inputs, pairs, engine), _fit(inputs, reverse, engine)
    for f, r in zip(
        forward.iter_rows(named=True), backward.iter_rows(named=True), strict=True
    ):
        assert (r["from"], r["to"]) == (f["to"], f["from"])
        assert _close(r["m"], -f["m"])
        assert _close(r["se"], f["se"])
        assert _close(r["p"], f["p"])


@pytest.mark.parametrize("chunk_size", [1, 7, 2000])
def test_chunk_size_does_not_change_output(chunk_size):
    inputs = _random_inputs(3, genotypes=25, total=50, high=20, edges=True)
    pairs = _pairs(inputs, count=30)
    assert_frame_equal(_fit(inputs, pairs, chunk_size=chunk_size), _fit(inputs, pairs))


@pytest.mark.parametrize("chunk_size", [0, -1, 1.5, True])
def test_chunk_size_must_be_a_positive_integer(chunk_size):
    inputs = _random_inputs(0, genotypes=5, total=50, high=20)
    with pytest.raises(ValueError, match="chunk_size"):
        _fit(inputs, _pairs(inputs, count=3), chunk_size=chunk_size)
