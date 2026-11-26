import polars as pl
import pytest

import barseqfit as bsf
from barseqfit import _fit

from ._papkou import TIME_UNIT, assert_close, load
from .test_fit_wt import _tiny_inputs

ENGINES = ["vectorised", "scalar"]


def _run(d, pairs, engine="vectorised"):
    return bsf.fit_pair(
        d["counts"],
        d["samples"],
        d["sequences"],
        pairs=pairs,
        elapsed_time=d["elapsed"],
        time_unit=TIME_UNIT,
        engine=engine,
        pseudocount="all_zeros",
        dispersion="fixed",
    )


@pytest.mark.parametrize("engine", ENGINES)
def test_fit_pair_reproduces_published_table(data_dir, engine):
    d = load(data_dir)
    expected = pl.read_csv(data_dir / "expected_pair.csv")
    result = _run(d, d["pairs"], engine)
    assert result.height == expected.height == 10
    for got, want in zip(
        result.iter_rows(named=True), expected.iter_rows(named=True), strict=True
    ):
        assert (got["from"], got["to"]) == (want["from"], want["to"])
        assert got["status"] == "converged"
        for key in ("m", "se", "p"):
            assert_close(
                got[key], want[key], f"{want['from']}->{want['to']} {key}", column=key
            )


@pytest.mark.parametrize("engine", ENGINES)
def test_fit_pair_reversal_negates_m(data_dir, engine):
    d = load(data_dir)
    forward = _run(d, d["pairs"], engine)
    reverse = _run(
        d,
        d["pairs"].select(pl.col("to").alias("from"), pl.col("from").alias("to")),
        engine,
    )
    for f, r in zip(
        forward.iter_rows(named=True), reverse.iter_rows(named=True), strict=True
    ):
        assert (r["from"], r["to"]) == (f["to"], f["from"])
        assert r["m"] == pytest.approx(-f["m"], abs=1e-9, rel=0)
        assert r["se"] == pytest.approx(f["se"], abs=1e-9, rel=0)
        assert r["p"] == pytest.approx(f["p"], abs=1e-9, rel=0)


def test_fit_pair_rejects_self_pair(data_dir):
    d = load(data_dir)
    a = d["pairs"]["from"][0]
    with pytest.raises(ValueError, match="from == to"):
        _run(d, pl.DataFrame({"from": [a], "to": [a]}))


def test_fit_pair_rejects_duplicate_pairs(data_dir):
    d = load(data_dir)
    with pytest.raises(ValueError, match="duplicate"):
        _run(d, pl.concat([d["pairs"][:1], d["pairs"][:1]]))


def test_fit_pair_rejects_empty_pairs(data_dir):
    d = load(data_dir)
    with pytest.raises(ValueError, match="pairs must be nonempty"):
        _run(d, d["pairs"].clear())


def test_fit_pair_rejects_unknown_sequence():
    pairs = pl.DataFrame({"from": ["wt"], "to": ["nope"]})
    with pytest.raises(ValueError, match="absent from sequences"):
        bsf.fit_pair(
            *_tiny_inputs(),
            pairs=pairs,
            elapsed_time=1.0,
            time_unit="generations",
            pseudocount="all_zeros",
            dispersion="fixed",
        )


def test_fit_pair_validates_self_pair_before_unknown_ids():
    pairs = pl.DataFrame({"from": ["nope"], "to": ["nope"]})
    with pytest.raises(ValueError, match="from == to"):
        bsf.fit_pair(
            *_tiny_inputs(),
            pairs=pairs,
            elapsed_time=1.0,
            time_unit="generations",
            pseudocount="all_zeros",
            dispersion="fixed",
        )


def test_fit_pair_rejects_non_numeric_time():
    counts, samples, sequences = _tiny_inputs()
    samples = samples.with_columns(pl.col("time").cast(pl.String))
    pairs = pl.DataFrame({"from": ["wt"], "to": ["variant"]})
    with pytest.raises(ValueError, match="samples.time must be numeric"):
        bsf.fit_pair(
            counts,
            samples,
            sequences,
            pairs=pairs,
            elapsed_time=1.0,
            time_unit="generations",
            pseudocount="all_zeros",
            dispersion="fixed",
        )


def test_fit_pair_failed_status_for_rank_deficient_design(monkeypatch):
    original = _fit._logit_fit

    def broken(success, failure, design, **kwargs):
        return original(success, failure, design * [1.0, 1.0, 1.0, 0.0], **kwargs)

    monkeypatch.setattr(_fit, "_logit_fit", broken)
    pairs = pl.DataFrame({"from": ["wt"], "to": ["variant"]})
    row = bsf.fit_pair(
        *_tiny_inputs(),
        pairs=pairs,
        elapsed_time=1.0,
        time_unit="generations",
        engine="scalar",
        pseudocount="all_zeros",
        dispersion="fixed",
    ).row(0, named=True)
    assert row["status"] == "failed"
    assert (row["m"], row["se"], row["p"]) == (None, None, None)
    assert "rank-deficient" in row["diagnostics"]["error"]
    assert (row["dispersion"], row["df_residual"]) == (None, None)


def test_fit_pair_records_fitted_probabilities_warning():
    # Complete separation with huge totals drives mu to within 10 * eps of 1.
    big = 10**9
    samples = pl.DataFrame(
        {
            "sample_id": ["a0", "a1", "b0", "b1"],
            "replicate_id": ["a", "a", "b", "b"],
            "condition": ["c"] * 4,
            "time": [0.0, 1.0, 0.0, 1.0],
            "sample_total": [big] * 4,
        }
    )
    counts = pl.DataFrame(
        {
            "sequence_id": ["wt"] * 4 + ["variant"] * 4,
            "sample_id": ["a0", "a1", "b0", "b1"] * 2,
            "count": [big // 2, big // 2, big // 2, big // 2, 0, big, 0, big],
        }
    )
    sequences = pl.DataFrame({"sequence_id": ["wt", "variant"], "sequence": ["A", "C"]})
    pairs = pl.DataFrame({"from": ["wt"], "to": ["variant"]})
    row = bsf.fit_pair(
        counts,
        samples,
        sequences,
        pairs=pairs,
        elapsed_time=1.0,
        time_unit="g",
        pseudocount="all_zeros",
        dispersion="fixed",
    ).row(0, named=True)
    assert "fitted probabilities numerically 0 or 1 occurred" in (
        row["diagnostics"]["warning"] or ""
    )
