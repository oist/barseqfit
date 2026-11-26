"""Data-quality flag: censored."""

import polars as pl
import pytest

import barseqfit as bsf

from ._papkou import REFERENCE_ID, TIME_UNIT, load
from .test_fit_wt import ENGINES

SAMPLES = ["b1", "b2", "a1", "a2"]
# Raw counts at (b1, b2, a1, a2); total 1000 per sample.
COUNTS = {
    "wt": [100, 110, 120, 130],
    "ok": [90, 95, 200, 210],
    "cens": [50, 60, 0, 0],
    "low": [2, 3, 4, 5],
    "lowafter": [90, 95, 1, 2],
    "sporadic": [0, 80, 150, 170],
}


def _inputs():
    samples = pl.DataFrame(
        {
            "sample_id": SAMPLES,
            "replicate_id": ["r1", "r2", "r1", "r2"],
            "condition": ["c"] * 4,
            "time": [0, 0, 2, 2],
            "sample_total": [1000] * 4,
        }
    )
    counts = pl.DataFrame(
        [
            {"sequence_id": g, "sample_id": s, "count": c}
            for g, row in COUNTS.items()
            for s, c in zip(SAMPLES, row, strict=True)
        ]
    )
    sequences = pl.DataFrame(
        {"sequence_id": list(COUNTS), "sequence": ["AAAA"] * len(COUNTS)}
    )
    return counts, samples, sequences


def _wt(engine="vectorised", **kwargs):
    return bsf.fit_wt(
        *_inputs(),
        reference_id="wt",
        elapsed_time=2.0,
        time_unit="g",
        engine=engine,
        **kwargs,
        pseudocount="all_zeros",
        dispersion="fixed",
    )


def _flags(frame, key="sequence_id"):
    return dict(zip(frame[key], frame["flags"], strict=True))


def _pair(pairs, engine="vectorised", **kwargs):
    return bsf.fit_pair(
        *_inputs(),
        pairs=pl.DataFrame(
            pairs, schema={"from": pl.String, "to": pl.String}, orient="row"
        ),
        elapsed_time=2.0,
        time_unit="g",
        engine=engine,
        **kwargs,
        pseudocount="all_zeros",
        dispersion="fixed",
    )


@pytest.mark.parametrize("engine", ENGINES)
def test_defaults_flag_only_censored(engine):
    out = _wt(engine)
    assert _flags(out) == {
        "wt": None,
        "ok": None,
        "cens": "censored",
        "low": None,
        "lowafter": None,
        "sporadic": None,
    }
    assert out.schema["flags"] == pl.String
    assert out.columns[-7:] == [
        "df_residual",
        "flags",
        "diagnostics",
        "elapsed_time",
        "time_unit",
        "pseudocount",
        "dispersion_model",
    ]


@pytest.mark.parametrize("engine", ENGINES)
def test_reference_row_carries_its_own_flags(engine):
    counts, samples, sequences = _inputs()
    counts = counts.with_columns(
        pl.when(
            (pl.col("sequence_id") == "wt") & pl.col("sample_id").is_in(["a1", "a2"])
        )
        .then(0)
        .otherwise(pl.col("count"))
        .alias("count")
    )
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
    )
    row = out.filter(pl.col("sequence_id") == "wt").row(0, named=True)
    assert row["status"] == "reference"
    assert row["flags"] == "censored"
    # A comparison against a censored reference carries the reference's flag too.
    assert _flags(out)["ok"] == "censored"
    assert _wt(engine).filter(pl.col("sequence_id") == "wt")["flags"].item() is None


@pytest.mark.parametrize("engine", ENGINES)
def test_pair_flags_symmetric_under_reversal(engine):
    pairs = [("ok", "cens"), ("cens", "ok"), ("low", "ok"), ("ok", "sporadic")]
    out = _pair(pairs, engine)
    got = {
        (f, t): fl
        for f, t, fl in zip(out["from"], out["to"], out["flags"], strict=True)
    }
    assert got[("ok", "cens")] == got[("cens", "ok")] == "censored"
    assert got[("low", "ok")] is None
    assert got[("ok", "sporadic")] is None
    both = _pair([("cens", "low"), ("low", "cens")], engine)
    assert set(both["flags"]) == {"censored"}


@pytest.mark.parametrize("engine", ENGINES)
def test_censored_flag_follows_raw_counts_not_the_pseudocount(engine):
    """The zero later-time counts are patched to 1 for the fit, yet still flagged."""
    wt = _wt(engine)
    assert _flags(wt)["cens"] == "censored"
    assert _flags(wt)["sporadic"] is None
    pair = _pair([("ok", "cens"), ("cens", "ok")], engine)
    assert pair["flags"].to_list() == ["censored", "censored"]


def test_engines_identical_flags():
    assert _wt("vectorised")["flags"].to_list() == _wt("scalar")["flags"].to_list()


@pytest.mark.parametrize("engine", ENGINES)
def test_fit_wt_and_fit_pair_agree_on_flags(engine):
    wt = _wt(engine).filter(pl.col("sequence_id") != "wt")
    tests = wt["sequence_id"].to_list()
    pair = _pair([("wt", t) for t in tests], engine)
    assert pair["flags"].to_list() == wt["flags"].to_list()


@pytest.mark.parametrize("engine", ENGINES)
def test_fixture_flags_null_under_defaults(data_dir, engine):
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
    assert wt["flags"].null_count() == wt.height
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
    assert pair["flags"].null_count() == pair.height


def test_fixture_has_no_censored_genotype(data_dir):
    d = load(data_dir)
    later = d["samples"].filter(pl.col("time") == pl.col("time").max())["sample_id"]
    after = (
        d["counts"]
        .filter(pl.col("sample_id").is_in(later.implode()))
        .group_by("sequence_id")
        .agg(pl.col("count").sum())
    )
    assert after.filter(pl.col("count") == 0).height == 0


def _many_inputs(n=50):
    """Reference plus `n` genotypes; flags follow the index (known by construction)."""
    samples = _inputs()[1]
    rows, expected = [], {}
    tab = {"wt": [100, 110, 120, 130]}
    for i in range(n):
        gid = f"g{i:02d}"
        if i % 5 == 0:
            tab[gid] = [50 + i, 60, 0, 0]
            expected[gid] = "censored"
        elif i % 7 == 0:
            tab[gid] = [2, 3, 4, 5]
            expected[gid] = None
        else:
            tab[gid] = [90 + i, 95, 200, 210]
            expected[gid] = None
    for g, row in tab.items():
        rows += [
            {"sequence_id": g, "sample_id": s, "count": c}
            for s, c in zip(SAMPLES, row, strict=True)
        ]
    sequences = pl.DataFrame(
        {"sequence_id": list(tab), "sequence": ["AAAA"] * len(tab)}
    )
    return pl.DataFrame(rows), samples, sequences, expected


@pytest.mark.parametrize("engine", ENGINES)
def test_flags_correct_across_chunks(engine):
    counts, samples, sequences, expected = _many_inputs()
    wt = bsf.fit_wt(
        counts,
        samples,
        sequences,
        reference_id="wt",
        elapsed_time=2.0,
        time_unit="g",
        engine=engine,
        chunk_size=7,
        pseudocount="all_zeros",
        dispersion="fixed",
    )
    got = dict(zip(wt["sequence_id"], wt["flags"], strict=True))
    assert got.pop("wt") is None
    assert got == expected
    pairs = pl.DataFrame({"from": ["wt"] * len(expected), "to": list(expected)})
    pr = bsf.fit_pair(
        counts,
        samples,
        sequences,
        pairs=pairs,
        elapsed_time=2.0,
        time_unit="g",
        engine=engine,
        chunk_size=7,
        pseudocount="all_zeros",
        dispersion="fixed",
    )
    assert dict(zip(pr["to"], pr["flags"], strict=True)) == expected


# Same hand-constructed counts as COUNTS (and as flag_data() in the R test-flags.R).
PARITY_EXPECTED = {
    "wt": None,
    "ok": None,
    "cens": "censored",
    "low": None,
    "lowafter": None,
    "sporadic": None,
}


@pytest.mark.parametrize("engine", ENGINES)
def test_flag_parity_literals(engine):
    """Literals mirrored verbatim in r/tests/testthat/test-flags.R."""
    wt = _wt(engine)
    assert _flags(wt) == PARITY_EXPECTED
    others = [g for g in COUNTS if g != "wt"]
    pr = _pair([("wt", g) for g in others], engine)
    assert _flags(pr, "to") == {g: PARITY_EXPECTED[g] for g in others}
