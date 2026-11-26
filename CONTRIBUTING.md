# Contributing

Thank you for your interest in BarSeqFit. By taking part you agree to follow the [Code of Conduct](CODE_OF_CONDUCT.md).

## Scope

The public API is deliberately limited to what the accompanying paper uses: `fit_sequences`, `fit_pair` and `fit_wt`, with the model choices inherited from Papkou et al. 2023. Bug fixes, documentation improvements and tests are very welcome. For new features or changes to the model, please open an issue first so we can agree whether it fits the scope.

Coding agents should read [AGENTS.md](AGENTS.md) before making changes.

## Set up

You need [uv](https://docs.astral.sh/uv/) and [just](https://just.systems/).

```sh
git clone https://github.com/oist/barseqfit
cd barseqfit
just sync
```

## Run the tests

```sh
just test
```

This runs `uv run pytest tests`. Pytest also collects the doctests in `src/barseqfit` and `README.md` when run without arguments (`uv run pytest`).

## The R side

The `r/` package is the reference implementation that the Python results are checked against. It needs R 4.1.3 and [renv](https://rstudio.github.io/renv/). Set `RSCRIPT` if `Rscript` on your path is a different version, then run:

```sh
just test-r
```

`just check` runs both the Python and R suites. Some recipes need the Papkou et al. 2023 data store described in `docs/validation.md`; ordinary contributions do not.

## Report a bug

Open an issue using the bug report template. Please include the BarSeqFit version, your Python version and operating system, a minimal example (the fixture in `tests/data` is a good starting point), and the full error message.

## Pull requests

- Keep each pull request focused on one change.
- Add or update tests, and make sure `just test` passes.
- Describe what changed and why in the pull request.
