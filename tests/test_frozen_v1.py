"""Regression against frozen outputs of an earlier version.

``tests/data/frozen_v1/*.parquet`` are earlier outputs on the Papkou et al.
2023 subset fixture (see ``tests/data/frozen_v1/README.md`` for their origin). Current code
with ``pseudocount="all_zeros"`` and ``dispersion="fixed"`` must reproduce them,
after dropping the echo columns added since. The comparison is exact for strings and integers and to relative 1e-9 (absolute 1e-12) for floats:
bitwise float equality is not portable across CPUs and BLAS builds (it failed on
some machines), while 1e-9 (the engine agreement) still catches any change of method.
"""

from pathlib import Path

import polars as pl
import pytest
from polars.testing import assert_frame_equal

from ._frozen_cases import CASES, ECHO_COLUMNS

DATA = Path(__file__).parent / "data"
FROZEN = DATA / "frozen_v1"


@pytest.fixture(scope="module")
def inputs():
    samples = pl.read_csv(DATA / "samples.csv")
    return {
        "counts": pl.read_csv(DATA / "counts.csv"),
        "samples": samples,
        "sequences": pl.read_csv(DATA / "sequences.csv"),
        "pairs": pl.read_csv(DATA / "pairs.csv"),
        "common": dict(
            elapsed_time=samples["time"].max() - samples["time"].min(),
            time_unit="generations",
        ),
    }


@pytest.mark.parametrize("name", list(CASES))
def test_matches_frozen_baseline(inputs, name):
    expected = pl.read_parquet(FROZEN / f"{name}.parquet")
    actual = CASES[name](inputs).drop(ECHO_COLUMNS, strict=False)
    if "diagnostics" in actual.columns:
        # Iteration counts and warning text flip at the 1e-8 deviance threshold on
        # near-separated rows across platforms (seen on some machines); compare the rest.
        for frame in (actual, expected):
            assert frame["diagnostics"].struct.field("error").is_null().to_list() == (
                expected["diagnostics"].struct.field("error").is_null().to_list()
            )
        actual, expected = actual.drop("diagnostics"), expected.drop("diagnostics")
    assert_frame_equal(actual, expected, rtol=1e-9, atol=1e-12)
