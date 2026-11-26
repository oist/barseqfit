# AGENTS.md

Guidance for coding agents working in this repository.

## What this is

BarSeqFit estimates the growth rate `r` of each DNA sequence, and the fitness difference `m` between sequences, from barcode sequencing counts at two time points. It uses a grouped-binomial logit regression and reproduces the analysis of Papkou et al. 2023. The Python package is `barseqfit`; an R package in `r/` is the reference it is checked against.

## Layout

- `src/barseqfit/`: the package (`_fit.py`, `_single.py`, `_comparison.py`, `_adjust.py`).
- `tests/`: pytest suite. `tests/data/` holds the small Papkou et al. 2023 fixture, the expected tables and `frozen_v1/` baselines.
- `r/`: the R reference package with its own `renv.lock` and testthat tests.
- `validation/`: Python library for the full-data validation (gate, identity, comparison); `validation.toml` holds its parameters.
- `scripts/`: R and Python scripts that regenerate fixtures and run the full-data validation.
- `docs/`: `model.md` (the model) and `validation.md` (validation records).
- `GLOSSARY.md`: the project vocabulary. `CONTRIBUTING.md`, `RELEASING.md`, `CHANGELOG.md`.

## Setup and commands

You need `uv` and `just`. Run `just --list` to see every recipe.

- `just sync`: install the development environment (`uv sync`).
- `just test`: Python tests (`uv run pytest tests`). Plain `uv run pytest` also runs the doctests in `src/barseqfit` and `README.md`.
- `just test-r`: R testthat suite. It needs R 4.1.3 and renv; set `RSCRIPT` if `Rscript` on your path is another version.
- `just check`: Python and R suites together.

The repository defines no lint or format recipe, so do not invent one.

## Rules

- Use the terms in `GLOSSARY.md` (for example growth rate, fitness difference, reference sequence, pair, comparison) in code, docs and messages.
- The public API is limited to what the paper uses: `fit_sequences`, `fit_pair` and `fit_wt`. `pseudocount` and `dispersion` are required and each accepts one value (`"all_zeros"` and `"fixed"`). Do not add public functions or options; for a new feature or a model change, ask for an issue first.
- Do not edit the generated tables in `tests/data` by hand. Regenerate them with `just make-fixture` (counts, samples, sequences, pairs, expected tables) and `just make-padj-expected` (`expected_*_padj_r.csv`). Both need R 4.1.3, renv and the Papkou et al. 2023 data store (`$PAPKOU2023RUGGED_DATA`).
- Do not edit or regenerate `tests/data/frozen_v1/`. Its generator is not in this repository. If `tests/test_frozen_v1.py` fails, your change altered the numerical results; treat that as a bug unless the maintainer says otherwise.
- Use R only to regenerate ground truth and to run the Papkou et al. 2023 export (`just export-papkou`). Ordinary changes need Python only.
- Add a line to `CHANGELOG.md` for every user-visible change.
- Keep the README examples correct; they run as tests.

## Proposing a change

1. Create a branch off `main`.
2. Keep the change to one topic. Add or update tests.
3. Run `just test` (and `just test-r` if you touched `r/`) until it passes.
4. Open a pull request that says what changed and why. Follow `CONTRIBUTING.md`.
