"""Case table for the frozen-v1 test.

``name -> callable(inputs) -> DataFrame``. The baseline files in
``tests/data/frozen_v1`` predate the required ``pseudocount`` and ``dispersion``
arguments; the cases pass them explicitly and the test drops the echo columns
they add (see ``ECHO_COLUMNS``) before comparing.
"""

import barseqfit as bsf

REFERENCE_ID = "GCCGATCTC"

# Columns added after the baselines were generated; they echo the arguments.
ECHO_COLUMNS = ("pseudocount", "dispersion_model")

_CHOICES = {"pseudocount": "all_zeros"}


def _wt(**kw):
    return lambda i: bsf.fit_wt(
        i["counts"],
        i["samples"],
        i["sequences"],
        reference_id=REFERENCE_ID,
        **i["common"],
        **_CHOICES,
        dispersion="fixed",
        **kw,
    )


def _pair(**kw):
    return lambda i: bsf.fit_pair(
        i["counts"],
        i["samples"],
        i["sequences"],
        pairs=i["pairs"],
        **i["common"],
        **_CHOICES,
        dispersion="fixed",
        **kw,
    )


def _sequences(**kw):
    return lambda i: bsf.fit_sequences(
        i["counts"], i["samples"], i["sequences"], **i["common"], **_CHOICES, **kw
    )


CASES = {
    "wt_default": _wt(),
    "wt_scalar": _wt(engine="scalar"),
    "pair_default": _pair(),
    "pair_scalar": _pair(engine="scalar"),
    "sequences": _sequences(),
}
