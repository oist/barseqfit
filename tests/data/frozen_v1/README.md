# Frozen baselines

The `*.parquet` files here are outputs of an earlier version, run on the
Papkou et al. 2023 subset fixture in `tests/data/`. They are kept as proof that
the default output has not changed since that version. `tests/test_frozen_v1.py`
runs the current code with `pseudocount="all_zeros"` and `dispersion="fixed"`
and drops the `pseudocount` and `dispersion_model` echo columns, which were
added later. It then compares the frames: floats to relative 1e-9 (absolute
1e-12), everything else exactly. The `diagnostics` column is not compared in
full: only whether its `error` field is null is checked for each row, and the
remaining fields (iteration counts and warning text, which can flip at the
convergence threshold on near-separated rows across platforms) are dropped.

The generator script is not part of the public tree, so the files are not
regenerated: a change that makes this test fail changes the numerical results.
