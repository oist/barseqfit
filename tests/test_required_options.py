"""Required ``pseudocount`` and ``dispersion``: accepted values, echo columns, API surface."""

import inspect

import polars as pl
import pytest

import barseqfit as bsf

from .test_fit_wt import _tiny_inputs

REJECTED_PSEUDOCOUNTS = ["after_zeros", "separated_only", "none"]
REJECTED_DISPERSIONS = ["quasi"]
COMMON = {"elapsed_time": 1.0, "time_unit": "generations"}
PAIRS = pl.DataFrame({"from": ["wt"], "to": ["variant"]})

CALLS = {
    "fit_wt": lambda **kw: bsf.fit_wt(
        *_tiny_inputs(), reference_id="wt", **COMMON, **kw
    ),
    "fit_pair": lambda **kw: bsf.fit_pair(*_tiny_inputs(), pairs=PAIRS, **COMMON, **kw),
    "fit_sequences": lambda **kw: bsf.fit_sequences(*_tiny_inputs(), **COMMON, **kw),
}
# fit_sequences has no dispersion choice.
CHOICES = {
    "fit_wt": {"pseudocount": "all_zeros", "dispersion": "fixed"},
    "fit_pair": {"pseudocount": "all_zeros", "dispersion": "fixed"},
    "fit_sequences": {"pseudocount": "all_zeros"},
}


def test_public_api_is_the_three_fit_functions():
    assert sorted(bsf.__all__) == ["fit_pair", "fit_sequences", "fit_wt"]
    for removed in (
        "fit_replicates",
        "summarise_replicates",
        "fit_sequence_replicates",
        "neighbour_index",
        "error_neighbour_ratio_from_epsilon",
    ):
        assert not hasattr(bsf, removed)


@pytest.mark.parametrize(
    ("name", "missing"),
    [(name, arg) for name, choices in CHOICES.items() for arg in choices],
)
def test_required_argument_has_no_default(name, missing):
    kwargs = {k: v for k, v in CHOICES[name].items() if k != missing}
    with pytest.raises(TypeError, match=missing):
        CALLS[name](**kwargs)
    parameter = inspect.signature(getattr(bsf, name)).parameters[missing]
    assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
    assert parameter.default is inspect.Parameter.empty


def test_fit_sequences_takes_no_dispersion():
    with pytest.raises(TypeError, match="dispersion"):
        CALLS["fit_sequences"](dispersion="fixed", **CHOICES["fit_sequences"])


@pytest.mark.parametrize("name", CALLS)
@pytest.mark.parametrize("value", REJECTED_PSEUDOCOUNTS)
def test_other_pseudocounts_are_rejected(name, value):
    kwargs = {**CHOICES[name], "pseudocount": value}
    with pytest.raises(ValueError, match="pseudocount must be 'all_zeros'"):
        CALLS[name](**kwargs)


@pytest.mark.parametrize("name", ["fit_wt", "fit_pair"])
@pytest.mark.parametrize("value", REJECTED_DISPERSIONS)
def test_other_dispersions_are_rejected(name, value):
    kwargs = {**CHOICES[name], "dispersion": value}
    with pytest.raises(ValueError, match="dispersion must be 'fixed'"):
        CALLS[name](**kwargs)


@pytest.mark.parametrize("name", CALLS)
@pytest.mark.parametrize("value", [None, 1, True])
def test_non_string_pseudocount_is_rejected(name, value):
    kwargs = {**CHOICES[name], "pseudocount": value}
    with pytest.raises(ValueError, match="pseudocount"):
        CALLS[name](**kwargs)


@pytest.mark.parametrize("name", CALLS)
def test_echo_columns_follow_the_arguments(name):
    out = CALLS[name](**CHOICES[name])
    tail = ["elapsed_time", "time_unit", "pseudocount"]
    if name != "fit_sequences":
        tail.append("dispersion_model")
    assert out.columns[-len(tail) :] == tail
    assert out["pseudocount"].to_list() == ["all_zeros"] * out.height
    assert out.schema["pseudocount"] == pl.String
    if name == "fit_sequences":
        assert "dispersion_model" not in out.columns
    else:
        assert out["dispersion_model"].to_list() == ["fixed"] * out.height
        assert out.schema["dispersion_model"] == pl.String
        # The estimate column keeps its name and meaning.
        assert out.schema["dispersion"] == pl.Float64


@pytest.mark.parametrize("name", CALLS)
@pytest.mark.parametrize(
    "removed",
    [
        "min_before_reads",
        "min_after_reads",
        "family",
        "dispersion_method",
        "error_neighbour_ratio",
    ],
)
def test_removed_options_are_rejected(name, removed):
    with pytest.raises(TypeError, match=removed):
        CALLS[name](**CHOICES[name], **{removed: None})
