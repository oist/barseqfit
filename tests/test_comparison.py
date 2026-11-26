"""fit_wt is the pair fit over (reference, x) plus a synthetic reference row."""

import math

import polars as pl
import pytest

import barseqfit as bsf

from .test_fit_wt import ENGINES, _random_inputs

TOLERANCE = 1e-6
NUMERIC = ("m", "se", "p", "p_adj", "dispersion", "elapsed_time")


def _same_number(a, b):
    if a is None or b is None:
        return a is None and b is None
    if math.isnan(a) or math.isnan(b):
        return math.isnan(a) and math.isnan(b)
    return math.isclose(a, b, rel_tol=0, abs_tol=TOLERANCE)


@pytest.mark.parametrize("engine", ENGINES)
@pytest.mark.parametrize(
    ("seed", "edges", "total", "high"),
    [(0, False, 1000, 500), (3, True, 50, 20), (0, True, 1000, 10)],
)
def test_fit_wt_equals_fit_pair_on_reference_pairs(engine, seed, edges, total, high):
    inputs = _random_inputs(seed, genotypes=30, total=total, high=high, edges=edges)
    reference = "wt"
    ids = inputs[2]["sequence_id"].to_list()
    tests = [value for value in ids if value != reference]
    assert reference in ids
    options = {"elapsed_time": 1.0, "time_unit": "generations", "engine": engine}

    wt = bsf.fit_wt(
        *inputs,
        reference_id=reference,
        **options,
        pseudocount="all_zeros",
        dispersion="fixed",
    )
    pair = bsf.fit_pair(
        *inputs,
        pairs=pl.DataFrame({"from": [reference] * len(tests), "to": tests}),
        **options,
        pseudocount="all_zeros",
        dispersion="fixed",
    )

    wt = wt.filter(pl.col("sequence_id") != reference).rename(
        {"sequence_id": "to", "reference_id": "from"}
    )
    assert set(wt.columns) == set(pair.columns)
    wt = wt.select(pair.columns)
    assert wt.height == pair.height == len(tests)
    assert wt["to"].to_list() == pair["to"].to_list() == tests
    if (seed, edges, total, high) == (0, True, 1000, 10):
        # Exercises the scalar-oracle hand-off; must not revert to all converged.
        assert "not_converged" in pair["status"].to_list()
    assert wt["from"].to_list() == pair["from"].to_list()
    for column in pair.columns:
        if column in NUMERIC:
            for a, b in zip(wt[column].to_list(), pair[column].to_list(), strict=True):
                assert _same_number(a, b), (column, a, b)
        else:
            assert wt[column].to_list() == pair[column].to_list(), column
