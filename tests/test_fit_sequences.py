"""Per-sequence growth rates: fit_sequences and its replicate aggregation."""

import math

import numpy as np
import polars as pl
import pytest
from scipy.stats import norm

import barseqfit as bsf
from barseqfit import _fit as fit_module
from barseqfit import _single as single_module

from ._papkou import REFERENCE_ID, TIME_UNIT, load
from .test_fit_wt import ENGINES, _close, _random_inputs

TOTAL = 10_000
TIMES = (0.0, 2.0)
# sequence -> {replicate: (before, after)}
TABLE = {
    "a": {"r1": (100, 150), "r2": (120, 140), "r3": (90, 170)},
    "b": {"r1": (80, 260), "r2": (100, 220), "r3": (70, 300)},
}
COLUMNS = [
    "sequence_id",
    "r",
    "se_pooled",
    "p_pooled",
    "n_replicates",
    "n_replicates_converged",
    "r_replicate_mean",
    "r_replicate_sd",
    "se_replicate",
    "status",
    "dispersion",
    "df_residual",
    "flags",
    "diagnostics",
    "elapsed_time",
    "time_unit",
    "pseudocount",
]
REPLICATE_COLUMNS = [
    "sequence_id",
    "replicate_id",
    "r",
    "se",
    "status",
    "flags",
    "diagnostics",
    "elapsed_time",
    "time_unit",
]


def _inputs(table=TABLE):
    reps = list(next(iter(table.values())))
    samples = pl.DataFrame(
        {
            "sample_id": [f"b_{r}" for r in reps] + [f"a_{r}" for r in reps],
            "replicate_id": reps + reps,
            "condition": ["c"] * (2 * len(reps)),
            "time": [TIMES[0]] * len(reps) + [TIMES[1]] * len(reps),
            "sample_total": [TOTAL] * (2 * len(reps)),
        }
    )
    counts = pl.DataFrame(
        [
            {"sequence_id": g, "sample_id": f"{kind}_{r}", "count": cell[i]}
            for g, by_rep in table.items()
            for r, cell in by_rep.items()
            for i, kind in enumerate(("b", "a"))
        ]
    )
    sequences = pl.DataFrame(
        {"sequence_id": list(table), "sequence": ["AAAA"] * len(table)}
    )
    return counts, samples, sequences


def _sequences(inputs, engine="vectorised", **kwargs):
    kwargs = {"pseudocount": "all_zeros", **kwargs}
    return bsf.fit_sequences(
        *inputs, elapsed_time=2.0, time_unit="g", engine=engine, **kwargs
    )


def _replicates(inputs, engine="vectorised", chunk_size=2000, elapsed_time=2.0):
    """Internal per-replicate table behind the replicate columns of fit_sequences."""
    counts, samples, sequences = inputs
    fit_module._check_common(
        counts, samples, sequences, elapsed_time, "g", engine, chunk_size
    )
    return single_module._replicate_table(
        counts,
        samples,
        sequences,
        elapsed_time=elapsed_time,
        time_unit="g",
        engine=engine,
        chunk_size=chunk_size,
    )


def _logit(count):
    p = count / TOTAL
    return math.log(p / (1 - p))


@pytest.mark.parametrize("engine", ENGINES)
def test_replicate_slope_is_the_closed_form(engine):
    out = _replicates(_inputs(), engine)
    assert out.columns == REPLICATE_COLUMNS
    assert out.schema["replicate_id"] == pl.String
    assert out["sequence_id"].to_list() == ["a"] * 3 + ["b"] * 3
    assert out["replicate_id"].to_list() == ["r1", "r2", "r3"] * 2
    for row in out.iter_rows(named=True):
        before, after = TABLE[row["sequence_id"]][row["replicate_id"]]
        expected = (_logit(after) - _logit(before)) / (TIMES[1] - TIMES[0])
        assert row["status"] == "converged"
        assert abs(row["r"] - expected) < 1e-12
        assert row["flags"] is None
        assert row["elapsed_time"] == 2.0 and row["time_unit"] == "g"


@pytest.mark.parametrize("engine", ENGINES)
def test_aggregation_matches_numpy_on_the_replicate_table(engine):
    inputs = _inputs()
    replicates = _replicates(inputs, engine)
    out = _sequences(inputs, engine)
    assert out.columns == COLUMNS
    assert out.schema == pl.Schema(
        {
            "sequence_id": pl.String,
            "r": pl.Float64,
            "se_pooled": pl.Float64,
            "p_pooled": pl.Float64,
            "n_replicates": pl.Int64,
            "n_replicates_converged": pl.Int64,
            "r_replicate_mean": pl.Float64,
            "r_replicate_sd": pl.Float64,
            "se_replicate": pl.Float64,
            "status": pl.String,
            "dispersion": pl.Float64,
            "df_residual": pl.Int64,
            "flags": pl.String,
            "diagnostics": out.schema["diagnostics"],
            "elapsed_time": pl.Float64,
            "time_unit": pl.String,
            "pseudocount": pl.String,
        }
    )
    assert out["pseudocount"].to_list() == ["all_zeros"] * 2
    assert out["sequence_id"].to_list() == ["a", "b"]
    for row in out.iter_rows(named=True):
        slopes = replicates.filter(pl.col("sequence_id") == row["sequence_id"])[
            "r"
        ].to_numpy()
        assert row["n_replicates"] == row["n_replicates_converged"] == 3
        assert row["r_replicate_mean"] == pytest.approx(slopes.mean(), abs=1e-12)
        assert row["r_replicate_sd"] == pytest.approx(slopes.std(ddof=1), abs=1e-12)
        assert row["se_replicate"] == pytest.approx(
            slopes.std(ddof=1) / math.sqrt(3), abs=1e-12
        )
        assert row["status"] == "converged"
        assert row["df_residual"] == 4
        assert row["se_pooled"] > 0 and 0 <= row["p_pooled"] <= 1
        assert row["dispersion"] is not None
        assert row["flags"] is None
        assert row["elapsed_time"] == 2.0 and row["time_unit"] == "g"


@pytest.mark.parametrize("engine", ENGINES)
def test_null_replicate_statistics_with_one_replicate(engine):
    table = {"a": {"r1": (100, 150)}, "b": {"r1": (80, 260)}}
    out = _sequences(_inputs(table), engine)
    for row in out.iter_rows(named=True):
        before, after = table[row["sequence_id"]]["r1"]
        expected = (_logit(after) - _logit(before)) / 2.0
        assert row["n_replicates"] == row["n_replicates_converged"] == 1
        assert row["r_replicate_mean"] == pytest.approx(expected, abs=1e-12)
        assert row["r_replicate_sd"] is None
        assert row["se_replicate"] is None
        # One replicate is a saturated pooled fit: still converged, no dispersion.
        assert row["status"] == "converged"
        assert row["r"] == pytest.approx(expected, abs=1e-8)
        assert row["df_residual"] == 0 and row["dispersion"] is None


def _replicate_table(rows):
    """Hand-built replicate table: (sequence, r, status) rows."""
    return pl.DataFrame(
        {
            "sequence_id": [row[0] for row in rows],
            "replicate_id": [f"r{i}" for i in range(len(rows))],
            "r": [row[1] for row in rows],
            "status": [row[2] for row in rows],
            "flags": [row[3] if len(row) > 3 else None for row in rows],
        },
        schema_overrides={"r": pl.Float64, "flags": pl.String},
    )


def test_failed_replicates_reduce_the_converged_count_only():
    rows = [
        # three converged, one failed with a null slope
        ("a", 0.5, "converged"),
        ("a", 0.7, "converged"),
        ("a", 0.9, "converged"),
        ("a", None, "failed"),
        # two converged, one not converged (slope kept but excluded)
        ("b", 0.2, "converged"),
        ("b", 0.4, "converged"),
        ("b", 99.0, "not_converged"),
        # one converged
        ("c", 0.3, "converged"),
        ("c", None, "failed"),
        # none converged
        ("d", None, "failed"),
        ("d", 5.0, "not_converged"),
    ]
    out = single_module._aggregate(_replicate_table(rows), ["a", "b", "c", "d", "e"])
    by = {row["sequence_id"]: row for row in out.iter_rows(named=True)}
    assert [by[k]["n_replicates"] for k in "abcde"] == [4, 3, 2, 2, 0]
    assert [by[k]["n_replicates_converged"] for k in "abcde"] == [3, 2, 1, 0, 0]
    a = np.array([0.5, 0.7, 0.9])
    assert by["a"]["r_replicate_mean"] == pytest.approx(a.mean(), abs=1e-12)
    assert by["a"]["r_replicate_sd"] == pytest.approx(a.std(ddof=1), abs=1e-12)
    assert by["a"]["se_replicate"] == pytest.approx(
        a.std(ddof=1) / math.sqrt(3), abs=1e-12
    )
    assert by["b"]["r_replicate_mean"] == pytest.approx(0.3, abs=1e-12)
    assert by["b"]["r_replicate_sd"] == pytest.approx(np.std([0.2, 0.4], ddof=1))
    assert by["b"]["se_replicate"] == pytest.approx(by["b"]["r_replicate_sd"] / 2**0.5)
    for key in "cde":
        assert by[key]["r_replicate_sd"] is None
        assert by[key]["se_replicate"] is None
    assert by["c"]["r_replicate_mean"] == pytest.approx(0.3, abs=1e-12)
    # No converged replicate (or no replicate at all): null mean.
    assert by["d"]["r_replicate_mean"] is None
    assert by["e"]["r_replicate_mean"] is None
    assert (
        out.schema["n_replicates"] == out.schema["n_replicates_converged"] == pl.Int64
    )


def test_aggregate_reports_zero_count_flag():
    rows = [("a", 0.1, "converged", "zero_count_replicate"), ("b", 0.2, "converged")]
    out = single_module._aggregate(_replicate_table(rows), ["a", "b"])
    assert out["_zero_count"].to_list() == [True, False]


@pytest.mark.parametrize("engine", ENGINES)
def test_zero_count_replicate_flag(engine):
    table = {
        "a": {"r1": (100, 0), "r2": (120, 140), "r3": (90, 170)},
        "b": {"r1": (80, 260), "r2": (100, 220), "r3": (70, 300)},
    }
    inputs = _inputs(table)
    replicates = _replicates(inputs, engine)
    # The zero is still fitted (all_zeros pseudocount) and flagged on its replicate.
    assert replicates["status"].to_list() == ["converged"] * 6
    assert replicates["flags"].to_list() == [
        "censored|zero_count_replicate",
        None,
        None,
        None,
        None,
        None,
    ]
    flags = dict(_sequences(inputs, engine).select("sequence_id", "flags").iter_rows())
    assert flags == {"a": "zero_count_replicate", "b": None}
    # The patched success count is 1; Diff stays at the raw sample_total - 0.
    r1 = replicates.row(0, named=True)
    assert r1["r"] == pytest.approx((math.log(1 / TOTAL) - _logit(100)) / 2.0, abs=1e-9)


def test_reference_row_is_a_real_fit(data_dir):
    d = load(data_dir)
    out = bsf.fit_sequences(
        d["counts"],
        d["samples"],
        d["sequences"],
        elapsed_time=d["elapsed"],
        time_unit=TIME_UNIT,
        pseudocount="all_zeros",
    )
    assert out.height == d["sequences"].height
    assert out["sequence_id"].to_list() == d["sequences"]["sequence_id"].to_list()
    wt = out.filter(pl.col("sequence_id") == REFERENCE_ID).row(0, named=True)
    assert wt["status"] == "converged"
    assert wt["se_pooled"] > 0
    assert wt["r"] != 0
    assert wt["n_replicates"] == wt["n_replicates_converged"] == 6
    assert wt["se_replicate"] > 0
    assert wt["df_residual"] == 10
    assert (out["status"] == "converged").all()


class TestValidation:
    @pytest.mark.parametrize("fit", [_sequences, _replicates])
    def test_missing_columns(self, fit):
        counts, samples, sequences = _inputs()
        with pytest.raises(ValueError, match="counts is missing required columns"):
            fit((counts.drop("count"), samples, sequences))
        with pytest.raises(ValueError, match="samples is missing required columns"):
            fit((counts, samples.drop("sample_total"), sequences))
        with pytest.raises(ValueError, match="sequences is missing required columns"):
            fit((counts, samples, sequences.drop("sequence_id")))

    @pytest.mark.parametrize("fit", [_sequences, _replicates])
    def test_requires_dataframes(self, fit):
        counts, samples, sequences = _inputs()
        with pytest.raises(TypeError, match="counts must be a Polars DataFrame"):
            fit((counts.to_dicts(), samples, sequences))

    @pytest.mark.parametrize("fit", [_sequences, _replicates])
    def test_three_time_points(self, fit):
        counts, samples, sequences = _inputs()
        samples = samples.with_columns(
            pl.when(pl.col("sample_id") == "b_r1")
            .then(1.0)
            .otherwise(pl.col("time"))
            .alias("time")
        )
        with pytest.raises(ValueError, match="exactly two sample time points"):
            fit((counts, samples, sequences))

    @pytest.mark.parametrize("function", ["fit_sequences"])
    @pytest.mark.parametrize(
        ("kwargs", "match"),
        [
            ({"elapsed_time": 1.0}, "elapsed_time"),
            ({"elapsed_time": float("nan")}, "elapsed_time"),
            ({"elapsed_time": 2.0, "time_unit": " "}, "time_unit"),
            ({"elapsed_time": 2.0, "engine": "gpu"}, "engine"),
            ({"elapsed_time": 2.0, "chunk_size": 0}, "chunk_size"),
            ({"elapsed_time": 2.0, "pseudocount": "some"}, "pseudocount"),
        ],
    )
    def test_invalid_arguments(self, function, kwargs, match):
        kwargs = {"time_unit": "g", "pseudocount": "all_zeros", **kwargs}
        with pytest.raises(ValueError, match=match):
            getattr(bsf, function)(*_inputs(), **kwargs)

    @pytest.mark.parametrize("function", ["fit_sequences"])
    def test_missing_sequence_sample_count(self, function):
        counts, samples, sequences = _inputs()
        counts = counts.filter(
            ~((pl.col("sequence_id") == "b") & (pl.col("sample_id") == "b_r1"))
        )
        with pytest.raises(ValueError, match="complete sequence/sample coverage"):
            getattr(bsf, function)(
                counts,
                samples,
                sequences,
                elapsed_time=2.0,
                time_unit="g",
                pseudocount="all_zeros",
            )

    @pytest.mark.parametrize("function", ["fit_sequences"])
    @pytest.mark.parametrize("bad", [2.5, float("nan")])
    def test_noninteger_counts(self, function, bad):
        counts, samples, sequences = _inputs()
        counts = counts.with_columns(pl.col("count").cast(pl.Float64))
        counts[0, "count"] = bad
        with pytest.raises(ValueError, match="must contain finite integer values"):
            getattr(bsf, function)(
                counts,
                samples,
                sequences,
                elapsed_time=2.0,
                time_unit="g",
                pseudocount="all_zeros",
            )

    def test_replicate_at_one_time_only(self):
        counts, samples, sequences = _inputs()
        samples = samples.filter(pl.col("sample_id") != "a_r3")
        counts = counts.filter(pl.col("sample_id") != "a_r3")
        with pytest.raises(ValueError, match="both time points"):
            _sequences((counts, samples, sequences))


def test_replicates_reject_duplicate_replicate_time_sample():
    counts, samples, sequences = _inputs()
    samples = samples.with_columns(
        pl.when(pl.col("sample_id") == "b_r2")
        .then(pl.lit("r1"))
        .otherwise(pl.col("replicate_id"))
        .alias("replicate_id")
    )
    for fit in (_replicates, _sequences):
        with pytest.raises(ValueError, match="exactly one sample per"):
            fit((counts, samples, sequences))


def test_replicates_reject_replicate_missing_a_time_point():
    counts, samples, sequences = _inputs()
    samples = samples.filter(pl.col("sample_id") != "a_r3")
    counts = counts.filter(pl.col("sample_id") != "a_r3")
    with pytest.raises(ValueError, match="both time points"):
        _replicates((counts, samples, sequences))


def test_replicate_id_is_string_for_integer_ids():
    counts, samples, sequences = _inputs()
    mapping = {"r1": 1, "r2": 2, "r3": 3}
    samples = samples.with_columns(
        pl.col("replicate_id").replace_strict(mapping, return_dtype=pl.Int64)
    )
    out = _replicates((counts, samples, sequences))
    assert out.schema["replicate_id"] == pl.String
    assert out["replicate_id"].to_list() == ["1", "2", "3"] * 2


def test_output_follows_sequences_order():
    counts, samples, sequences = _inputs()
    out = _sequences((counts, samples, sequences.reverse()))
    assert out["sequence_id"].to_list() == ["b", "a"]
    reps = _replicates((counts, samples, sequences.reverse()))
    assert reps["sequence_id"].to_list() == ["b"] * 3 + ["a"] * 3


def _count_logit_fit(monkeypatch):
    calls = []
    original = fit_module._logit_fit

    def counting(*args, **kwargs):
        calls.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(fit_module, "_logit_fit", counting)
    return calls


def _assert_tables_agree(fast, slow, keys):
    assert fast.height == slow.height
    for a, b in zip(
        fast.iter_rows(named=True), slow.iter_rows(named=True), strict=True
    ):
        assert a["status"] == b["status"], a
        assert a["diagnostics"] == b["diagnostics"], a
        for key in keys:
            assert _close(a[key], b[key]), (a["sequence_id"], key)


def test_engines_agree_on_clean_data_without_fallback(monkeypatch):
    inputs = _random_inputs(0, genotypes=40, total=1000, high=500, edges=False)
    kwargs = {"elapsed_time": 1.0, "time_unit": "g", "pseudocount": "all_zeros"}
    slow = bsf.fit_sequences(*inputs, engine="scalar", **kwargs)
    calls = _count_logit_fit(monkeypatch)
    fast = bsf.fit_sequences(*inputs, engine="vectorised", **kwargs)
    assert calls == []
    _assert_tables_agree(
        fast,
        slow,
        ("r", "se_pooled", "p_pooled", "r_replicate_mean", "se_replicate"),
    )


@pytest.mark.parametrize("chunk_size", [1, 7, 2000])
def test_engines_agree_with_edge_counts(chunk_size):
    # Zero and full-total counts: separation sends fits to the scalar oracle.
    inputs = _random_inputs(1, genotypes=30, total=1000, high=500, edges=True)
    kwargs = {"elapsed_time": 1.0, "time_unit": "g", "pseudocount": "all_zeros"}
    slow = bsf.fit_sequences(*inputs, engine="scalar", **kwargs)
    fast = bsf.fit_sequences(*inputs, chunk_size=chunk_size, **kwargs)
    _assert_tables_agree(
        fast, slow, ("r", "se_pooled", "p_pooled", "r_replicate_mean", "se_replicate")
    )
    assert fast["n_replicates_converged"].to_list() == (
        slow["n_replicates_converged"].to_list()
    )
    fast_rep = _replicates(inputs, chunk_size=chunk_size, elapsed_time=1.0)
    slow_rep = _replicates(inputs, engine="scalar", elapsed_time=1.0)
    _assert_tables_agree(fast_rep, slow_rep, ("r", "se"))
    assert fast["flags"].to_list() == slow["flags"].to_list()


def test_p_pooled_is_two_sided_z_test():
    out = _sequences(_inputs())
    for row in out.iter_rows(named=True):
        assert row["p_pooled"] == pytest.approx(
            2 * norm.sf(abs(row["r"] / row["se_pooled"])), abs=1e-15
        )
    assert np.isfinite(out["se_pooled"].to_numpy()).all()


@pytest.mark.parametrize("engine", ENGINES)
def test_identity_with_fit_wt(data_dir, engine):
    # fit_wt's type*time design is block-separable, so at convergence
    # m == r_i - r_wt and, under fixed dispersion, se == sqrt(se_i^2 + se_wt^2).
    d = load(data_dir)
    kwargs = {"elapsed_time": d["elapsed"], "time_unit": TIME_UNIT, "engine": engine}
    seq = bsf.fit_sequences(
        d["counts"], d["samples"], d["sequences"], **kwargs, pseudocount="all_zeros"
    )
    wt_fit = bsf.fit_wt(
        d["counts"],
        d["samples"],
        d["sequences"],
        reference_id=REFERENCE_ID,
        **kwargs,
        pseudocount="all_zeros",
        dispersion="fixed",
    )
    wt = seq.filter(pl.col("sequence_id") == REFERENCE_ID).row(0, named=True)
    joined = seq.join(wt_fit, on="sequence_id", suffix="_wt").filter(
        pl.col("status_wt") != "reference"
    )
    assert joined.height == seq.height - 1
    checked = 0
    for row in joined.iter_rows(named=True):
        if row["status"] != "converged" or row["status_wt"] != "converged":
            continue
        checked += 1
        label = row["sequence_id"]
        assert abs((row["r"] - wt["r"]) - row["m"]) < 1e-7, label
        se = math.sqrt(row["se_pooled"] ** 2 + wt["se_pooled"] ** 2)
        # 2e-5, not tighter: glm.fit stops on a relative deviance change of 1e-8,
        # and the joint 4-parameter fit_wt and the 2-parameter per-sequence fits
        # can stop on different IRLS iterations; the se comes from the weights at
        # the start of the final iteration. Observed max 1.27e-5 (CGCACTTGC);
        # where the iteration counts match the gap is ~1e-13.
        assert abs(se - row["se"]) < 2e-5, label
    assert checked == seq.height - 1
