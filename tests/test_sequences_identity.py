"""The fit_sequences / fit_wt identity as a pure comparison on the fixture."""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

import barseqfit as bsf
from validation.config import load_run_config
from validation.identity import IdentityRule, sequences_identity

from ._papkou import REFERENCE_ID, TIME_UNIT, load

REPO = Path(__file__).resolve().parent.parent
RULE = IdentityRule.from_run_config(load_run_config(REPO / "validation.toml"))


@pytest.fixture(scope="module")
def fits():
    d = load(REPO / "tests" / "data")
    kwargs = {"elapsed_time": d["elapsed"], "time_unit": TIME_UNIT}
    seq = bsf.fit_sequences(
        d["counts"], d["samples"], d["sequences"], **kwargs, pseudocount="all_zeros"
    )
    wt = bsf.fit_wt(
        d["counts"],
        d["samples"],
        d["sequences"],
        reference_id=REFERENCE_ID,
        **kwargs,
        pseudocount="all_zeros",
        dispersion="fixed",
    )
    return seq, wt


def test_rule_from_run_config():
    assert RULE == IdentityRule(m_tol=1e-7, se_tol=2e-5)


def test_identity_holds_on_the_fixture(fits):
    seq, wt = fits
    verdict = sequences_identity(seq, wt, REFERENCE_ID, RULE)
    assert verdict.passed, verdict.failures
    assert verdict.checked == seq.height - 1
    assert verdict.skipped == 0
    assert verdict.max_m < RULE.m_tol and verdict.max_se < RULE.se_tol


def test_a_perturbed_m_is_reported(fits):
    seq, wt = fits
    victim = wt.filter(pl.col("status") == "converged")["sequence_id"][0]
    bad = wt.with_columns(
        pl.when(pl.col("sequence_id") == victim)
        .then(pl.col("m") + 1e-5)
        .otherwise(pl.col("m"))
        .alias("m")
    )
    verdict = sequences_identity(seq, bad, REFERENCE_ID, RULE)
    assert not verdict.passed
    assert verdict.outside_m == [victim]
    assert verdict.outside_se == []


def test_a_perturbed_se_is_reported(fits):
    seq, wt = fits
    victim = wt.filter(pl.col("status") == "converged")["sequence_id"][0]
    bad = wt.with_columns(
        pl.when(pl.col("sequence_id") == victim)
        .then(pl.col("se") * 1.01)
        .otherwise(pl.col("se"))
        .alias("se")
    )
    verdict = sequences_identity(seq, bad, REFERENCE_ID, RULE)
    assert verdict.outside_se == [victim] and verdict.outside_m == []


def test_unconverged_rows_are_skipped_not_checked(fits):
    seq, wt = fits
    victim = wt.filter(pl.col("status") == "converged")["sequence_id"][0]
    skipped = seq.with_columns(
        pl.when(pl.col("sequence_id") == victim)
        .then(pl.lit("failed"))
        .otherwise(pl.col("status"))
        .alias("status")
    )
    verdict = sequences_identity(skipped, wt, REFERENCE_ID, RULE)
    assert verdict.skipped == 1 and verdict.checked == seq.height - 2


def test_id_set_mismatch_fails(fits):
    seq, wt = fits
    dropped = wt.filter(pl.col("sequence_id") != REFERENCE_ID).slice(1)
    verdict = sequences_identity(seq, dropped, REFERENCE_ID, RULE)
    assert not verdict.passed
    assert any("id" in line for line in verdict.failures)


def test_missing_reference_fails(fits):
    seq, wt = fits
    verdict = sequences_identity(
        seq.filter(pl.col("sequence_id") != REFERENCE_ID), wt, REFERENCE_ID, RULE
    )
    assert not verdict.passed
