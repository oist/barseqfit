import math

import numpy as np
import polars as pl
import pytest
from polars.testing import assert_frame_equal
from scipy.stats import norm

import barseqfit as bsf
from barseqfit import _fit as fit_module

ENGINES = ["vectorised", "scalar"]


def _tiny_inputs():
    sample_rows = []
    count_rows = []
    for time, variant_count in ((0.0, 2), (1.0, 4)):
        for replicate in range(6):
            sample_id = f"t{int(time)}-r{replicate}"
            sample_rows.append(
                {
                    "sample_id": sample_id,
                    "replicate_id": f"r{replicate}",
                    "condition": "constant",
                    "time": time,
                    "sample_total": 8,
                }
            )
            count_rows.extend(
                [
                    {"sequence_id": "wt", "sample_id": sample_id, "count": 2},
                    {
                        "sequence_id": "variant",
                        "sample_id": sample_id,
                        "count": variant_count,
                    },
                ]
            )

    return (
        pl.DataFrame(count_rows),
        pl.DataFrame(sample_rows),
        pl.DataFrame({"sequence_id": ["wt", "variant"], "sequence": ["AAAA", "AAAT"]}),
    )


# Sample times are 0 and 1, so elapsed_time must be 1.0 (contract: max - min).
def _fit(inputs, engine="vectorised", **kwargs):
    kwargs = {
        "elapsed_time": 1.0,
        "pseudocount": "all_zeros",
        "dispersion": "fixed",
        **kwargs,
    }
    return bsf.fit_wt(
        *inputs,
        reference_id="wt",
        time_unit="generations",
        engine=engine,
        **kwargs,
    )


@pytest.mark.parametrize("engine", ENGINES)
def test_fit_wt_matches_analytic_two_timepoint_binomial_contrast(engine):
    result = _fit(_tiny_inputs(), engine)

    assert isinstance(result, pl.DataFrame)
    # The reference row is included (m=0, null se/p, status "reference").
    assert result.height == 2
    ref = result.filter(pl.col("sequence_id") == "wt").row(0, named=True)
    assert (ref["m"], ref["se"], ref["p"], ref["status"]) == (
        0.0,
        None,
        None,
        "reference",
    )
    row = result.filter(pl.col("sequence_id") == "variant").row(0, named=True)
    assert row["diagnostics"]["converged"] is True
    assert row["sequence_id"] == "variant"
    assert row["reference_id"] == "wt"
    assert row["status"] == "converged"
    assert math.isclose(row["m"], math.log(3), rel_tol=0, abs_tol=1e-8)
    # The analytic Fisher-information SE is sqrt(5/12). The R reference
    # uses the final IRLS weights at its deviance stopping point, which can
    # differ slightly from that limit.
    expected_se = math.sqrt(5 / 12)
    assert math.isclose(row["se"], expected_se, rel_tol=0, abs_tol=2e-6)
    expected_p = 2 * norm.sf(math.log(3) / expected_se)
    assert math.isclose(row["p"], expected_p, rel_tol=0, abs_tol=2e-6)
    assert row["elapsed_time"] == 1.0
    assert row["time_unit"] == "generations"


@pytest.mark.parametrize("engine", ENGINES)
def test_fit_wt_rejects_missing_sequence_sample_count(engine):
    counts, samples, sequences = _tiny_inputs()
    counts = counts.filter(
        ~((pl.col("sequence_id") == "variant") & (pl.col("sample_id") == "t0-r0"))
    )

    with pytest.raises(ValueError, match="complete sequence/sample coverage"):
        bsf.fit_wt(
            counts,
            samples,
            sequences,
            reference_id="wt",
            elapsed_time=1.0,
            time_unit="generations",
            engine=engine,
            pseudocount="all_zeros",
            dispersion="fixed",
        )


@pytest.mark.parametrize("engine", ENGINES)
@pytest.mark.parametrize(
    ("input_name", "replacement"),
    [
        ("count", 2.5),
        ("count", float("nan")),
        ("sample_total", 8.5),
        ("sample_total", float("nan")),
    ],
)
def test_fit_wt_rejects_noninteger_or_nonfinite_counts_and_totals(
    input_name, replacement, engine
):
    counts, samples, sequences = _tiny_inputs()
    if input_name == "count":
        counts = counts.with_columns(
            pl.when(
                (pl.col("sequence_id") == "variant") & (pl.col("sample_id") == "t0-r0")
            )
            .then(replacement)
            .otherwise(pl.col("count"))
            .alias("count")
        )
    else:
        samples = samples.with_columns(
            pl.when(pl.col("sample_id") == "t0-r0")
            .then(replacement)
            .otherwise(pl.col("sample_total"))
            .alias("sample_total")
        )

    with pytest.raises(ValueError, match="must contain finite integer values"):
        bsf.fit_wt(
            counts,
            samples,
            sequences,
            reference_id="wt",
            elapsed_time=1.0,
            time_unit="generations",
            engine=engine,
            pseudocount="all_zeros",
            dispersion="fixed",
        )


@pytest.mark.parametrize("engine", ENGINES)
def test_fit_wt_rejects_non_numeric_time(engine):
    counts, samples, sequences = _tiny_inputs()
    samples = samples.with_columns(pl.col("time").cast(pl.String))
    with pytest.raises(ValueError, match="samples.time must be numeric"):
        _fit((counts, samples, sequences), engine)


@pytest.mark.parametrize("engine", ENGINES)
@pytest.mark.parametrize("elapsed", [True, "1", None])
def test_fit_wt_rejects_non_numeric_elapsed_time(elapsed, engine):
    with pytest.raises(ValueError, match="elapsed_time"):
        bsf.fit_wt(
            *_tiny_inputs(),
            reference_id="wt",
            elapsed_time=elapsed,
            time_unit="generations",
            engine=engine,
            pseudocount="all_zeros",
            dispersion="fixed",
        )


def test_fit_wt_rejects_unknown_engine():
    with pytest.raises(ValueError, match="engine"):
        _fit(_tiny_inputs(), "gpu")


def _random_inputs(seed, *, genotypes, total, high, n_per_time=6, edges=True):
    """Wild type plus random genotypes.

    With ``edges`` every 7th count is 0 and every 11th is ``total``; without,
    counts lie in ``[1, high)`` and the data are clean.
    """
    rng = np.random.default_rng(seed)
    samples, counts = [], []
    for time in (0.0, 1.0):
        for replicate in range(n_per_time):
            samples.append(
                {
                    "sample_id": f"t{int(time)}-r{replicate}",
                    "replicate_id": f"r{replicate}",
                    "condition": "constant",
                    "time": time,
                    "sample_total": total,
                }
            )
    ids = ["wt", *(f"g{i}" for i in range(genotypes))]
    position = 0
    for sequence_id in ids:
        for sample in samples:
            position += 1
            count = int(rng.integers(0, high) if edges else rng.integers(1, high))
            if edges and position % 7 == 0:
                count = 0
            if edges and position % 11 == 0:
                count = total
            counts.append(
                {
                    "sequence_id": sequence_id,
                    "sample_id": sample["sample_id"],
                    "count": count,
                }
            )
    sequences = pl.DataFrame({"sequence_id": ids, "sequence": ["A"] * len(ids)})
    return pl.DataFrame(counts), pl.DataFrame(samples), sequences


def _close(a, b):
    if a is None or b is None:
        return a is None and b is None
    return math.isclose(a, b, rel_tol=0, abs_tol=1e-9)


def _assert_engines_agree(inputs, select=lambda row: True):
    vectorised = _fit(inputs, "vectorised").iter_rows(named=True)
    scalar = _fit(inputs, "scalar").iter_rows(named=True)
    checked = 0
    for fast, slow in zip(vectorised, scalar, strict=True):
        if not select(slow):
            continue
        checked += 1
        label = slow["sequence_id"]
        assert fast["status"] == slow["status"], label
        assert fast["diagnostics"] == slow["diagnostics"], label
        for key in ("m", "se", "p"):
            assert _close(fast[key], slow[key]), f"{label} {key}"
    assert checked


def _count_logit_fit(monkeypatch):
    calls = []
    original = fit_module._logit_fit

    def counting(*args, **kwargs):
        calls.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(fit_module, "_logit_fit", counting)
    return calls


@pytest.mark.parametrize(("total", "high"), [(50, 20), (1000, 500)])
def test_engines_agree_on_well_conditioned_data(total, high, monkeypatch):
    # No zero or full-total counts, so every fit is clean and the batch solver
    # must produce all of them without the scalar fallback.
    inputs = _random_inputs(0, genotypes=40, total=total, high=high, edges=False)
    scalar = _fit(inputs, "scalar")
    calls = _count_logit_fit(monkeypatch)
    _fit(inputs, "vectorised")
    assert calls == []
    rows = scalar.filter(pl.col("status") != "reference")
    assert (rows["status"] == "converged").all()
    assert rows["diagnostics"].struct.field("warning").is_null().all()
    _assert_engines_agree(inputs)


def test_vectorised_is_independent_of_input_row_order(monkeypatch):
    inputs = _random_inputs(1, genotypes=3000, total=1000, high=500, edges=False)
    counts, samples, sequences = (
        frame.sample(fraction=1.0, shuffle=True, seed=seed)
        for frame, seed in zip(inputs, (11, 12, 13), strict=True)
    )
    calls = _count_logit_fit(monkeypatch)
    shuffled = _fit((counts, samples, sequences), "vectorised")
    assert calls == []
    monkeypatch.undo()

    expected = _fit(inputs, "scalar").sort("sequence_id")
    # Output follows the (shuffled) sequences order; compare per sequence_id.
    assert shuffled["sequence_id"].to_list() == sequences["sequence_id"].to_list()
    shuffled = shuffled.sort("sequence_id")
    assert shuffled["sequence_id"].to_list() == expected["sequence_id"].to_list()
    for fast, slow in zip(
        shuffled.iter_rows(named=True), expected.iter_rows(named=True), strict=True
    ):
        label = slow["sequence_id"]
        assert fast["status"] == slow["status"], label
        for key in ("m", "se", "p"):
            assert _close(fast[key], slow[key]), f"{label} {key}"


# Near complete separation the IRLS iterates reach coefficients of ~1e15 and
# every step is dominated by rounding error, so batch and scalar solvers follow
# different orbits. The vectorised engine therefore hands every non-clean fit
# back to the scalar oracle, and agreement holds by construction.
def test_engines_agree_when_not_converging():
    inputs = _random_inputs(0, genotypes=40, total=1000, high=10)
    _assert_engines_agree(inputs, lambda row: row["status"] == "not_converged")


def test_engines_agree_at_fitted_probability_boundary():
    inputs = _random_inputs(0, genotypes=40, total=1000, high=10)
    _assert_engines_agree(
        inputs,
        lambda row: "numerically 0 or 1" in (row["diagnostics"]["warning"] or ""),
    )


def test_vectorised_falls_back_to_scalar_fit_on_separation_data(monkeypatch):
    inputs = _random_inputs(0, genotypes=40, total=1000, high=10)
    expected = _fit(inputs, "scalar")
    unclean = expected.filter(
        (pl.col("status") != "reference")
        & (
            (pl.col("status") != "converged")
            | pl.col("diagnostics").struct.field("warning").is_not_null()
        )
    ).height
    assert unclean

    calls = _count_logit_fit(monkeypatch)
    _fit(inputs, "vectorised")
    assert len(calls) == unclean


def test_engines_both_flag_boundary_and_non_convergence():
    """The separation fixture exercises both boundary and non-convergence paths."""
    inputs = _random_inputs(0, genotypes=40, total=1000, high=10)
    for engine in ENGINES:
        rows = _fit(inputs, engine).filter(pl.col("status") != "reference")
        warnings = rows["diagnostics"].struct.field("warning").drop_nulls().to_list()
        assert any("numerically 0 or 1" in w for w in warnings), engine
        assert (rows["status"] == "not_converged").any(), engine


def test_rank_deficient_design_fails_in_both_engines():
    # fit_wt always builds a full-rank design (two times x two types), so the
    # engines are compared directly with a duplicated column.
    success = np.array([[3.0, 5.0, 2.0, 6.0, 4.0, 4.0]])
    failure = 10.0 - success
    base = np.column_stack([np.ones(6), [0, 0, 0, 1, 1, 1], [0, 1, 0, 1, 0, 1]]).astype(
        float
    )
    design = np.column_stack([base, base[:, 2]])
    with pytest.raises(ValueError, match="rank-deficient design") as scalar_error:
        fit_module._logit_fit(success[0], failure[0], design)
    errors = fit_module._logit_fit_batch(success, failure, design[None]).errors
    assert errors == [str(scalar_error.value)]


@pytest.mark.parametrize("chunk_size", [1, 7, 2000])
def test_chunk_size_does_not_change_output(chunk_size):
    inputs = _random_inputs(3, genotypes=25, total=50, high=20, edges=True)
    assert_frame_equal(
        _fit(inputs, "vectorised", chunk_size=chunk_size), _fit(inputs, "vectorised")
    )


@pytest.mark.parametrize("chunk_size", [0, -1, 1.5, True])
def test_chunk_size_must_be_a_positive_integer(chunk_size):
    with pytest.raises(ValueError, match="chunk_size"):
        _fit(
            _random_inputs(0, genotypes=5, total=50, high=20),
            "vectorised",
            chunk_size=chunk_size,
        )


@pytest.mark.parametrize("engine", ENGINES)
def test_fit_wt_keeps_reference_row_at_its_position(engine):
    inputs = _random_inputs(1, genotypes=6, total=1000, high=500, edges=False)
    result = bsf.fit_wt(
        *inputs,
        reference_id="g4",
        elapsed_time=1.0,
        time_unit="generations",
        engine=engine,
        pseudocount="all_zeros",
        dispersion="fixed",
    )
    assert result["sequence_id"].to_list() == inputs[2]["sequence_id"].to_list()
    assert result.filter(pl.col("sequence_id") == "g4")["status"].to_list() == [
        "reference"
    ]


# Zero counts: every zero cell of Nreads becomes 1; Diff keeps the raw count.
_ZERO_COUNTS = {
    "wt": [10, 12, 11, 13],
    "sporadic": [0, 5, 7, 9],  # one zero before
    "afterzero": [8, 6, 0, 0],  # all-zero after
    "beforezero": [0, 0, 4, 6],  # all-zero before
}
_ZERO_SAMPLES = ["b1", "b2", "a1", "a2"]
_ZERO_EXPECTED = {
    "wt": [10, 12, 11, 13],
    "sporadic": [1, 5, 7, 9],
    "afterzero": [8, 6, 1, 1],
    "beforezero": [1, 1, 4, 6],
}


def _zero_inputs():
    samples = pl.DataFrame(
        {
            "sample_id": _ZERO_SAMPLES,
            "replicate_id": ["r1", "r2", "r1", "r2"],
            "condition": ["c"] * 4,
            "time": [0, 0, 2, 2],
            "sample_total": [100] * 4,
        }
    )
    counts = pl.DataFrame(
        [
            {"sequence_id": g, "sample_id": s, "count": c}
            for g, row in _ZERO_COUNTS.items()
            for s, c in zip(_ZERO_SAMPLES, row, strict=True)
        ]
    )
    sequences = pl.DataFrame(
        {"sequence_id": list(_ZERO_COUNTS), "sequence": ["AAAA"] * 4}
    )
    return counts, samples, sequences


def test_all_zeros_patches_every_zero_count():
    counts, samples, sequences = _zero_inputs()
    joined = fit_module._validate(counts, samples, sequences, 2, "g")
    for genotype, expected in _ZERO_EXPECTED.items():
        rows = joined.filter(pl.col("sequence_id") == genotype)
        by_sample = dict(zip(rows["sample_id"], rows["Nreads"], strict=True))
        assert [by_sample[s] for s in _ZERO_SAMPLES] == expected, genotype
        diff = dict(zip(rows["sample_id"], rows["Diff"], strict=True))
        assert [diff[s] for s in _ZERO_SAMPLES] == [
            100 - c for c in _ZERO_COUNTS[genotype]
        ]


@pytest.mark.parametrize("engine", ENGINES)
def test_zero_counts_reach_the_fit(engine):
    """A genotype with zeros is fitted on the patched counts, not dropped."""
    out = _fit(_zero_inputs(), engine, elapsed_time=2.0)
    assert set(out["status"]) <= {"reference", "converged", "not_converged"}
    assert out["m"].null_count() == 0


@pytest.mark.parametrize("engine", ENGINES)
def test_optional_sequence_columns_do_not_change_fit(engine):
    counts, samples, sequences = _tiny_inputs()
    kw = dict(reference_id="wt", elapsed_time=1.0, time_unit="h", engine=engine)
    base = bsf.fit_wt(
        counts, samples, sequences, **kw, pseudocount="all_zeros", dispersion="fixed"
    )
    for frame in (
        sequences.drop("sequence"),
        sequences.with_columns(aa_id=pl.Series(["P1", "P1"])),
    ):
        assert_frame_equal(
            bsf.fit_wt(
                counts,
                samples,
                frame,
                **kw,
                pseudocount="all_zeros",
                dispersion="fixed",
            ),
            base,
        )


@pytest.mark.parametrize(
    "column, values, match",
    [
        ("sequence", [1, 2], "non-null strings"),
        ("sequence", ["AAAA", None], "non-null strings"),
        ("aa_id", ["P1", None], "non-null strings"),
        ("aa_id", ["P1", ""], "nonempty"),
        ("aa_id", [1, 2], "non-null strings"),
    ],
)
def test_optional_sequence_columns_validated(column, values, match):
    counts, samples, sequences = _tiny_inputs()
    bad = sequences.with_columns(pl.Series(column, values))
    kw = dict(reference_id="wt", elapsed_time=1.0, time_unit="h")
    with pytest.raises(ValueError, match=match):
        bsf.fit_wt(
            counts, samples, bad, **kw, pseudocount="all_zeros", dispersion="fixed"
        )
