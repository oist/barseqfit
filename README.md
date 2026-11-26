# BarSeqFit

[![Licence: MIT](https://img.shields.io/badge/licence-MIT-blue.svg)](https://github.com/oist/barseqfit/blob/main/LICENSE)
![Python 3.12 | 3.13](https://img.shields.io/badge/python-3.12%20%7C%203.13-blue.svg)
<!-- TODO: add the DOI badge once the Zenodo record exists -->
<!-- TODO: add the coverage badge once coverage is published -->

BarSeqFit estimates how fast each DNA sequence grows in a pooled barcode sequencing experiment. From read counts at two time points it returns the growth rate `r` of each sequence, and the fitness difference `m` between two sequences or between a sequence and a wild-type reference. Both come from a grouped-binomial logit regression of the counts, with a standard error and a p-value. BarSeqFit reproduces the analysis of Papkou et al. 2023 and is the fitting software for the accompanying paper. The Python package is `barseqfit`.

## Install

Install the package from GitHub:

```bash
pip install "barseqfit @ git+https://github.com/oist/barseqfit.git@v0.2.0"
```

To work on the code, clone the repository and install the development environment with the project's tool (`uv`):

```bash
uv sync
```

Publishing to PyPI is a later step for the maintainer. After that publish, `pip install barseqfit` installs the same package.

### Supported platforms and requirements

- Python 3.12 and 3.13 on Linux, macOS and Windows. A continuous integration workflow for all six combinations is included in `.github/workflows/tests.yml`; it is not enabled on this repository yet, so run `uv run pytest` locally.
- Python dependencies: numpy, polars and scipy.
- R 4.1.3 with [renv](https://rstudio.github.io/renv/) is needed only to regenerate the ground-truth tables and to run the Papkou et al. 2023 export. It is tested on Linux. You do not need R to use the package.

## Quickstart

The example uses the small test data set in `tests/data`, a 21-sequence subset of Papkou et al. 2023. Run it from the repository root. The code below is also run as a test, so it stays correct.

```python
>>> import polars as pl
>>> import barseqfit as bsf
>>> counts = pl.read_csv("tests/data/counts.csv")
>>> samples = pl.read_csv("tests/data/samples.csv")
>>> sequences = pl.read_csv("tests/data/sequences.csv")
>>> pairs = pl.read_csv("tests/data/pairs.csv")
>>> elapsed = samples["time"].max() - samples["time"].min()
>>> elapsed
3.5166197103891137

```

Time is measured in generations here, so `r` and `m` are changes in logit abundance per generation. First, the growth rate of every sequence:

```python
>>> growth = bsf.fit_sequences(
...     counts, samples, sequences,
...     elapsed_time=elapsed, time_unit="generations", pseudocount="all_zeros",
... )
>>> growth.shape
(21, 17)
>>> growth["status"].unique().to_list()
['converged']
>>> growth.filter(pl.col("sequence_id") == "AAGACGAAT")["r"].round(2).item()
-0.67

```

Then the fitness difference of each listed pair, where `m` is the fitness of `to` relative to `from`. The result is compared with the published value for the same pairs:

```python
>>> pair_fit = bsf.fit_pair(
...     counts, samples, sequences, pairs=pairs,
...     elapsed_time=elapsed, time_unit="generations",
...     pseudocount="all_zeros", dispersion="fixed",
... )
>>> pair_fit.shape
(10, 15)
>>> published = pl.read_csv("tests/data/expected_pair.csv")
>>> bool((pair_fit["m"] - published["m"]).abs().max() < 1e-6)
True

```

And the fitness difference of every sequence against the wild type:

```python
>>> wt_fit = bsf.fit_wt(
...     counts, samples, sequences, reference_id="GCCGATCTC",
...     elapsed_time=elapsed, time_unit="generations",
...     pseudocount="all_zeros", dispersion="fixed",
... )
>>> sorted(wt_fit["status"].unique().to_list())
['converged', 'reference']

```

## Input tables

All three functions take Polars data frames.

| Table | Required columns | Meaning |
| :-- | :-- | :-- |
| `counts` | `sequence_id`, `sample_id`, `count` | One row per sequence and sample. Non-negative integers. Zeros are kept. Every sequence must appear in every sample. |
| `samples` | `sample_id`, `condition`, `replicate_id`, `time`, `sample_total` | One row per sample. `time` is numeric, in `time_unit`, with exactly two distinct values (before and after). `sample_total` is the total read count of the sample, a positive integer at least as large as every `count` in it. All samples share one `condition`. |
| `sequences` | `sequence_id` | One row per sequence. Optional `sequence` and `aa_id` columns (strings) are carried for your convenience and never affect a fit. |
| `pairs` | `from`, `to` | `fit_pair` only. Sequence ids. `m` is the fitness of `to` relative to `from`. A pair may not compare a sequence with itself or repeat. |

Further arguments:

- `elapsed_time`: the time between the two time points. It must equal `max(time) - min(time)` within 1e-9.
- `time_unit`: a non-empty label, such as `"generations"`. It is echoed in the output and gives the unit of `r` and `m`.
- `reference_id` (`fit_wt` only): the wild-type sequence id.
- `engine` and `chunk_size`: see [docs/model.md](https://github.com/oist/barseqfit/blob/main/docs/model.md). The defaults are right for most uses.

`fit_sequences` also needs exactly one sample per `replicate_id` and time point, and both time points for every replicate.

## Output

Each function returns one row per sequence (`fit_sequences`, `fit_wt`) or per pair (`fit_pair`), in input order. The columns `elapsed_time` and `time_unit` echo the arguments on every row.

### `fit_sequences`

| Column | Meaning |
| :-- | :-- |
| `r` | Growth rate of the sequence: slope of the logit of its read fraction against time, per `time_unit`, pooled over all samples. |
| `se_pooled`, `p_pooled` | Wald standard error of `r` and two-sided p-value of `r = 0`, assuming binomial variance (dispersion fixed at 1). They do not see disagreement between replicates. |
| `n_replicates`, `n_replicates_converged` | Replicates in `samples`, and how many of their separate fits converged. |
| `r_replicate_mean`, `r_replicate_sd` | Mean and sample standard deviation of the converged per-replicate slopes. The standard deviation is null below two converged replicates. |
| `se_replicate` | `r_replicate_sd / sqrt(n_replicates_converged)`. It reflects the spread between replicates, which `se_pooled` does not. Null when the standard deviation is. |
| `status` | `converged`, `not_converged` (estimate kept) or `failed` (estimate null), for the pooled fit. |
| `dispersion`, `df_residual` | Pearson dispersion of the pooled fit (a diagnostic only) and its residual degrees of freedom. |
| `flags` | Data-quality flags, `|`-joined, null when none. `censored`: every raw count at the later time is zero. `zero_count_replicate`: some replicate has a zero count at either time. Flags never change an estimate. |
| `diagnostics` | Struct with `converged`, `iterations`, `warning` and `error` of the pooled fit. |
| `pseudocount` | The `pseudocount` argument, echoed. |

### `fit_wt` and `fit_pair`

| Column | Meaning |
| :-- | :-- |
| `sequence_id` and `reference_id` (`fit_wt`), or `from` and `to` (`fit_pair`) | The compared sequences. |
| `m` | Fitness difference per `time_unit`: the type-by-time interaction of the logit model. Positive means the sequence (or `to`) grows faster than the wild type (or `from`). |
| `se`, `p` | Wald standard error of `m` and two-sided p-value of `m = 0`, assuming binomial variance. |
| `p_adj` | Benjamini-Hochberg adjusted `p` over the non-null `p` of the call. It depends on which comparisons are fitted together. |
| `status` | `converged`, `not_converged` (estimate kept), `failed` (estimate null) or, for the reference row of `fit_wt` only, `reference` (`m = 0`, null `se` and `p`). |
| `dispersion`, `df_residual` | Pearson dispersion (a diagnostic only) and residual degrees of freedom. |
| `flags` | As above; only `censored`, taken over both sequences of the comparison. |
| `diagnostics` | Struct with `converged`, `iterations`, `warning` and `error`. |
| `pseudocount`, `dispersion_model` | The `pseudocount` and `dispersion` arguments, echoed. |

The model, the status values and the identity between `r` and `m` are described in [docs/model.md](https://github.com/oist/barseqfit/blob/main/docs/model.md).

## Model choices inherited from Papkou et al. 2023

BarSeqFit keeps the choices of the published analysis so that it reproduces the published estimates. Four of them matter when you read a result:

- **Zeros are replaced by 1 at all time points.** Every zero count becomes 1 before fitting (`pseudocount="all_zeros"`). The complement count is computed from the raw count first, so the replacement does not change it.
- **Binomial dispersion is fixed at 1.** Standard errors and p-values come from the plain binomial model, with a normal (Wald) test (`dispersion="fixed"`). The reported `dispersion` is a diagnostic and does not scale the standard error.
- **No read-depth filter.** Every sequence in `sequences` is fitted, however few reads it has. Sequences with no reads at the later time are flagged `censored`.
- **The complement count is taken as given.** For each sample, the count of all other reads is `sample_total - count`, with `sample_total` read from your `samples` table. BarSeqFit does not recompute it from the rows of `counts`, so subsetting `counts` does not change it.

`pseudocount` and `dispersion` are required arguments with no default, so that every analysis states them. Version 0.2.0 accepts only `"all_zeros"` and `"fixed"`; any other value raises an error.

## Validation

The test suite checks the fits against the estimates published with Papkou et al. 2023 for a 21-sequence subset, to a tolerance of 1e-6, and checks the fast engine against a simple scalar reference. An R implementation in `r/` gives an independent check, and the full published data (261,333 sequences and 3,518,696 pairs) can be checked with the recipes below. The repository has no simulation test: the claim is reproduction of published estimates, as described in [docs/validation.md](https://github.com/oist/barseqfit/blob/main/docs/validation.md).

Run the tests:

```sh
just test      # Python tests
just test-r    # R tests; needs R 4.1.3 (set RSCRIPT to its path)
```

`uv run pytest` with no arguments also runs the doctests in this README.

The full-data recipes need the Zenodo data of Papkou et al. 2023 (record 10.5281/zenodo.13837982). Set `PAPKOU2023RUGGED_DATA` to the folder that holds it, in your shell or in a `.env` file (copy `.env.example`). Then run `just validate-full-wt` or `just validate-full-pair`. They are slow; the pair fit in R takes about two hours on four cores.

## Related packages

- [topographer](https://github.com/oist/topographer) — builds merge trees for graphs that carry a scalar on each vertex, for example fitness landscapes built from barseqfit's fitness tables.

## Publications using this package

Dubanevics et al. "Topological characterization of combinatorial fitness landscapes using merge trees." Manuscript in preparation (2026). Code: https://github.com/igorsdub/dubanevics2026topological

## How to cite

If you use BarSeqFit, please cite the software through its Zenodo record. See [`CITATION.cff`](https://github.com/oist/barseqfit/blob/main/CITATION.cff) for the machine-readable form, which GitHub shows as "Cite this repository".

<!-- TODO: add the DOI of the Zenodo record here and in CITATION.cff -->

Please also cite the analysis that BarSeqFit reproduces:

- Papkou A, García-Pastor L, Escudero JA, Wagner A. A rugged yet easily navigable fitness landscape. *Science*. 2023. doi:[10.1126/science.adh3860](https://doi.org/10.1126/science.adh3860)
- The authors' data and code record, from which the test data subset is taken (CC-BY-4.0): doi:[10.5281/zenodo.13837982](https://doi.org/10.5281/zenodo.13837982)

## Licence

MIT. Copyright (c) 2026 Okinawa Institute of Science and Technology Graduate University. See [LICENSE](https://github.com/oist/barseqfit/blob/main/LICENSE).

## Contributing

Contributions are welcome; see [CONTRIBUTING.md](https://github.com/oist/barseqfit/blob/main/CONTRIBUTING.md). Please follow the [Code of Conduct](https://github.com/oist/barseqfit/blob/main/CODE_OF_CONDUCT.md).
