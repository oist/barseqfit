"""The gate rule run on the fixture: one test per branch of the chain."""

from __future__ import annotations

import dataclasses
import functools
from pathlib import Path

import polars as pl
import pytest

import barseqfit as bsf
from validation import gate as _gate
from validation.compare import Tolerance
from validation.gate import (
    COLUMNS,
    GateRule,
    compare_chain,
    keyset,
    residuals,
    scalar_vs_vectorised,
    subset_residual_lines,
)

from ._papkou import REFERENCE_ID, TIME_UNIT, load

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "tests" / "data"
KEYS = {"wt": ["sequence_id"], "pair": ["from", "to"]}
ROW = {"wt": 5, "pair": 3}

RULE = GateRule(
    tolerance=Tolerance.fixture(),
    tie_fraction=0.0005,
    tie_se_abs=1e-4,
    tie_p_abs=1e-4,
    scalar_tol=1e-9,
)
# Fixture tables are tiny, so the default fraction gives a tie cap of 0.
LOOSE = dataclasses.replace(RULE, tie_fraction=0.1)


@functools.cache
def _frames(kind: str) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """(fitted, published, synthetic R table) from the committed fixture."""
    d = load(DATA)
    floats = {"m": pl.Float64, "se": pl.Float64, "p": pl.Float64}
    if kind == "wt":
        fitted = bsf.fit_wt(
            d["counts"],
            d["samples"],
            d["sequences"],
            reference_id=REFERENCE_ID,
            elapsed_time=d["elapsed"],
            time_unit=TIME_UNIT,
            engine="vectorised",
            pseudocount="all_zeros",
            dispersion="fixed",
        )
        published = pl.read_csv(DATA / "expected_wt.csv", schema_overrides=floats)
        not_ref = pl.col("sequence_id") != REFERENCE_ID
        fitted, published = fitted.filter(not_ref), published.filter(not_ref)
    else:
        fitted = bsf.fit_pair(
            d["counts"],
            d["samples"],
            d["sequences"],
            pairs=d["pairs"],
            elapsed_time=d["elapsed"],
            time_unit=TIME_UNIT,
            engine="vectorised",
            pseudocount="all_zeros",
            dispersion="fixed",
        )
        published = pl.read_csv(DATA / "expected_pair.csv", schema_overrides=floats)
    keys = KEYS[kind]
    r_table = fitted.select(
        *keys,
        "m",
        "se",
        "p",
        "status",
        pl.format(
            "iterations={}", pl.col("diagnostics").struct.field("iterations")
        ).alias("diagnostics"),
    )
    return fitted, published, r_table


def _shift(frame: pl.DataFrame, row: int, column: str, delta: float) -> pl.DataFrame:
    idx = pl.int_range(pl.len())
    return frame.with_columns(
        pl.when(idx == row)
        .then(pl.col(column) + delta)
        .otherwise(pl.col(column))
        .alias(column)
    )


def _set(frame: pl.DataFrame, row: int, column: str, value) -> pl.DataFrame:
    idx = pl.int_range(pl.len())
    return frame.with_columns(
        pl.when(idx == row)
        .then(pl.lit(value, dtype=frame.schema[column]))
        .otherwise(pl.col(column))
        .alias(column)
    )


def _make_tie(r_table: pl.DataFrame, row: int, se_delta: float) -> pl.DataFrame:
    """Give one R row a different iteration count and a shifted se."""
    r_table = _set(r_table, row, "diagnostics", "iterations=99")
    return _shift(r_table, row, "se", se_delta)


@pytest.fixture(params=["wt", "pair"])
def kind(request) -> str:
    return request.param


def _chain(kind, rule=RULE, *, fitted=None, published=None, r_table="default"):
    base_fitted, base_pub, base_r = _frames(kind)
    return compare_chain(
        base_fitted if fitted is None else fitted,
        base_pub if published is None else published,
        base_r if isinstance(r_table, str) else r_table,
        KEYS[kind],
        rule,
    )


def test_rule_from_run_config():
    config = {
        "gate": {
            "abs_tol": 1e-6,
            "rel_tol": 1e-6,
            "noise_floor": 1e-13,
            "tie_fraction": 0.0005,
            "tie_se_abs": 1e-4,
            "tie_p_abs": 1e-4,
            "scalar_tol": 1e-9,
        }
    }
    assert GateRule.from_run_config(config) == RULE


def test_clean_frames_pass(kind):
    chain = _chain(kind)
    assert chain.passed and chain.failures == []
    assert chain.res_r is not None
    assert not any("e.g." in line for line in chain.r_lines)


def test_residual_columns_in_order(kind):
    fitted, published, _ = _frames(kind)
    res = residuals(fitted, published, "pub", KEYS[kind], RULE.tolerance)
    expected = list(KEYS[kind])
    for c in COLUMNS:
        if c == "p":
            expected.append("p_uf_pub")
        expected += [f"{c}_abs_pub", f"{c}_rel_pub", f"{c}_out_pub"]
    assert res.columns == expected


def test_subset_residual_lines_filters_to_subset(kind):
    fitted, published, _ = _frames(kind)
    subset = fitted.select(KEYS[kind]).head(3)
    lines = subset_residual_lines(fitted, published, KEYS[kind], subset, RULE.tolerance)
    assert lines[0] == "- Rows compared with the published table: 3"
    assert len(lines) == 1 + 2 + len(COLUMNS)
    empty = subset_residual_lines(
        fitted, published, KEYS[kind], subset.clear(), RULE.tolerance
    )
    assert empty[0] == "- Rows compared with the published table: 0"


def test_threshold_tie_is_counted_not_failed(kind):
    _, _, r_table = _frames(kind)
    chain = _chain(kind, LOOSE, r_table=_make_tie(r_table, ROW[kind], 5e-6))
    assert chain.failures == []
    assert any("threshold ties (iterations differ): 1 " in x for x in chain.rule_lines)
    assert sum("| 99 |" in x for x in chain.tie_lines) == 1


def test_tie_cap_exceeded(kind):
    _, _, r_table = _frames(kind)
    chain = _chain(kind, RULE, r_table=_make_tie(r_table, ROW[kind], 5e-6))
    assert "threshold ties 1 exceed cap 0" in chain.failures
    assert not chain.passed


def test_tie_se_residual_over_cap(kind):
    _, _, r_table = _frames(kind)
    chain = _chain(kind, LOOSE, r_table=_make_tie(r_table, ROW[kind], 2e-4))
    assert any(
        f.startswith("tie se residual") and "exceeds" in f for f in chain.failures
    )


def test_r_vs_published_outlier_is_excused(kind):
    # Published differs from both fits (which agree): R's outlier excuses Python's.
    _, published, _ = _frames(kind)
    published = _shift(published, ROW[kind], "se", 1e-4)
    chain = _chain(kind, LOOSE, published=published)
    assert chain.failures == []
    assert any(
        "excused by a tie or an R outlier (reported): 1" in x for x in chain.rule_lines
    )


def test_python_outliers_exceed_r_outliers_plus_ties(kind):
    # R and published agree with each other but not with Python.
    _, published, r_table = _frames(kind)
    published = _shift(published, ROW[kind], "se", 1e-4)
    r_table = _shift(r_table, ROW[kind], "se", 1e-4)
    chain = _chain(kind, LOOSE, published=published, r_table=r_table)
    assert "se/p outside tolerance vs published on gated rows: 1" in chain.failures
    assert any("exceed R-vs-published" in f for f in chain.failures)


def test_p_underflow_is_treated_as_equal(kind):
    fitted, published, r_table = _frames(kind)
    row = ROW[kind]
    chain = _chain(
        kind,
        fitted=_set(fitted, row, "p", 1e-310),
        published=_set(published, row, "p", 0.0),
        r_table=_set(r_table, row, "p", 0.0),
    )
    assert chain.failures == []
    assert any("1 rows vs R, 1 vs published" in x for x in chain.rule_lines)


def test_id_set_mismatch(kind):
    _, _, r_table = _frames(kind)
    chain = _chain(kind, r_table=r_table.slice(1))
    assert any(f.startswith("R ids != published ids") for f in chain.failures)


def test_r_table_absent(kind):
    _, published, _ = _frames(kind)
    chain = _chain(
        kind, published=_shift(published, ROW[kind], "se", 1e-4), r_table=None
    )
    assert chain.res_r is None
    assert "se/p outside tolerance vs published: 1" in chain.failures
    assert chain.rule_lines == [
        "- R table absent: se, p asserted strictly vs published."
    ]


def test_m_outside_tolerance_vs_published(kind):
    _, published, _ = _frames(kind)
    chain = _chain(kind, published=_shift(published, ROW[kind], "m", 1e-3))
    assert "m outside tolerance vs published: 1" in chain.failures


def test_m_outside_tolerance_vs_r(kind):
    _, _, r_table = _frames(kind)
    chain = _chain(kind, r_table=_shift(r_table, ROW[kind], "m", 1e-3))
    assert any(f.startswith("m outside tolerance vs R:") for f in chain.failures)


def test_se_p_outside_tolerance_vs_r_on_matched_rows(kind):
    _, _, r_table = _frames(kind)
    chain = _chain(kind, r_table=_shift(r_table, ROW[kind], "se", 1e-3))
    assert any(
        f.startswith("se/p outside tolerance vs R on matched rows:")
        for f in chain.failures
    )


def test_tie_p_residual_over_cap(kind):
    _, _, r_table = _frames(kind)
    r_table = _set(r_table, ROW[kind], "diagnostics", "iterations=99")
    r_table = _shift(r_table, ROW[kind], "p", 2e-4)
    chain = _chain(kind, LOOSE, r_table=r_table)
    assert any(
        f.startswith("tie p residual") and "exceeds" in f for f in chain.failures
    )


def test_r_vs_published_outliers_exceed_cap(kind):
    _, published, _ = _frames(kind)
    published = _shift(published, ROW[kind], "se", 1e-4)
    chain = _chain(kind, RULE, published=published)
    assert any(f.startswith("R-vs-published outliers") for f in chain.failures)


def _sv(se: float, se_v: float) -> pl.DataFrame:
    diag = {"converged": True, "iterations": 5, "warning": None, "error": None}
    return pl.DataFrame(
        {
            "m": [1.0],
            "m_v": [1.0],
            "se": [se],
            "se_v": [se_v],
            "p": [0.5],
            "p_v": [0.5],
            "status": ["converged"],
            "status_v": ["converged"],
            "diagnostics": [diag],
            "diagnostics_v": [diag],
        }
    )


def test_scalar_vs_vectorised_both_nan_agree():
    nan = float("nan")
    stats = scalar_vs_vectorised(_sv(nan, nan), RULE.scalar_tol)
    assert stats["bad"] == 0 and stats["ok"]
    stats = scalar_vs_vectorised(_sv(nan, 1.0), RULE.scalar_tol)
    assert stats["bad"] == 1 and not stats["ok"]


def test_keyset():
    frame = pl.DataFrame({"a": ["x", "y"], "b": ["u", "v"]})
    assert keyset(frame, ["a"]) == {"x", "y"}
    assert keyset(frame, ["a", "b"]) == {("x", "u"), ("y", "v")}


def test_gate_module_is_pure():
    names = vars(_gate)
    assert "tomllib" not in names and "load_run_config" not in names
    assert not hasattr(_gate, "logging") and not hasattr(_gate, "os")
    text = Path(_gate.__file__).read_text()
    for forbidden in (
        "import logging",
        "from logging",
        "print(",
        "read_parquet",
        "read_csv",
        "write_csv",
        "open(",
        "tomllib",
        "validation.toml",
    ):
        assert forbidden not in text, forbidden


def test_script_adapter_holds_no_rule_literals():
    text = (REPO / "scripts" / "_gate_common.py").read_text()
    for forbidden in ("1e-13", "1e-300", "FLOOR_FOR"):
        assert forbidden not in text
