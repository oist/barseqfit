import polars as pl
import pytest

import barseqfit as bsf

from ._papkou import REFERENCE_ID, TIME_UNIT, assert_close, load


@pytest.mark.parametrize("engine", ["vectorised", "scalar"])
def test_fit_wt_reproduces_published_table(data_dir, engine):
    d = load(data_dir)
    expected = pl.read_csv(
        data_dir / "expected_wt.csv",
        schema_overrides={"se": pl.Float64, "p": pl.Float64},
    )
    result = bsf.fit_wt(
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
    assert result.height == expected.height
    got = {r["sequence_id"]: r for r in result.iter_rows(named=True)}
    for row in expected.iter_rows(named=True):
        actual = got[row["sequence_id"]]
        if row["sequence_id"] == REFERENCE_ID:
            assert actual["status"] == "reference"
            assert actual["m"] == 0
            assert actual["se"] is None and actual["p"] is None
            continue
        assert actual["status"] == "converged"
        for key in ("m", "se", "p"):
            assert_close(
                actual[key], row[key], f"{row['sequence_id']} {key}", column=key
            )


def test_fit_wt_rejects_wrong_elapsed_time(data_dir):
    d = load(data_dir)
    with pytest.raises(ValueError, match="elapsed_time"):
        bsf.fit_wt(
            d["counts"],
            d["samples"],
            d["sequences"],
            reference_id=REFERENCE_ID,
            elapsed_time=1.0,
            time_unit=TIME_UNIT,
            pseudocount="all_zeros",
            dispersion="fixed",
        )
