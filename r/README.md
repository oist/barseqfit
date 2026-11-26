# R package and validation

This directory is the R side of BarSeqFit: the R package `barseqfit` (base R
only) plus the validation material used to check the Python package.

```
DESCRIPTION, NAMESPACE   package metadata; exports fit_wt, fit_pair, papkou_store,
                         read_papkou_store, read_papkou_export, read_test_data,
                         test_data_dir
R/env.R                  machine-environment resolver (.env parser; store, scratch and
                         cores from the shell, a repo-root .env filling gaps)
R/io.R                   papkou_store() (resolves PAPKOU2023RUGGED_DATA through env.R)
                         and test_data_dir() (repo-root tests/data)
R/papkou.R               Papkou2023 constants (reference genotype, Generation time) and the
                         Export format: read_papkou_store() (the only reader of the raw RDS
                         files), read_papkou_export(), the Export writer and .papkou_subset()
R/config.R               Run config loader and the timing-file writer
R/fit_wt.R, R/fit_pair.R fit_wt() and fit_pair() (one binomial GLM per sequence/pair)
R/glm_internal.R         shared data assembly and single-fit helper
tests/testthat/          tests against tests/data/ plus the synthetic 24-row check
renv.lock, renv/         pinned R 4.1.3 environment
validation/               authors_scripts.md5 (md5 of the authors' scripts in the store) and run_authors_pipeline.R; see docs/validation.md
```

`samples` requires `sample_id, condition, replicate_id, time, sample_total`; `time` is numeric
time in `time_unit` and `elapsed_time` must equal `max(time) - min(time)`.

`read_papkou_store()` and `read_papkou_export()` return the same shape: the requested
tables plus `reference_id`, `elapsed_time`, `time_unit` and `manifest`, as does
`read_test_data()` for the fixture. The fit functions take `cores` as an argument
(default 1) and read no environment variable; the scripts in `scripts/` resolve the
machine environment (`.resolve_machine_env`) and pass the core count in.

Fixtures live once at the repo root in `tests/data/` and are shared with the
Python tests. The data store is configured with the `PAPKOU2023RUGGED_DATA`
environment variable (or a repo-root `.env`, see `.env.example`).

## Tests

From this directory, under R 4.1.3:

```sh
Rscript -e 'testthat::test_dir("tests/testthat", load_package = "source")'
```

`test-synthetic.R` holds the deterministic 24-row synthetic fixture (six
replicate samples per time, eight reads each; WT 2/8 at both times, variant 2/8
then 4/8). Analytically `m = log(3)` and `se = sqrt(5/12)`; finite-iteration
IRLS differs in the last digits, so the test uses a loose tolerance.
