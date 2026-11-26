# The model

This page describes what BarSeqFit computes. Every rule below is the one the code applies. The choices that BarSeqFit inherits from Papkou et al. 2023 are marked as such.

## Data layout

An experiment has a set of sequences and a set of samples. Each sample is one sequencing run of the pool, taken at one of two time points. In the Papkou et al. 2023 data there are six biological replicates, each with a "before" sample and an "after" sample, so 12 samples.

The data come in tables.

- `counts` holds `count`, the number of reads of a sequence in a sample.
- `samples` holds, for each sample, its `replicate_id`, its `time` and its `sample_total`, the total number of reads in the sample.
- `sequences` lists the sequences.

Every sequence must have a count in every sample. For a sequence in a sample, the number of reads of all other sequences is `Diff = sample_total - count`. BarSeqFit takes `sample_total` from your table and does not recompute it from `counts`. This is inherited from Papkou et al. 2023, where the total is the row sum of the full count matrix, taken before any sequence was left out.

There must be exactly two distinct values of `time`. The smaller is "before" and the larger is "after". The number of samples at each time is not restricted, except for `fit_sequences`, which needs exactly one sample per replicate and time point.

## Time and units

`samples.time` is used as it is. It is numeric and measured in `time_unit`, a label that you choose (the Papkou et al. 2023 data use generations). `elapsed_time` is the gap between the two time points and must equal `max(time) - min(time)` to within 1e-9; the check guards against passing a time in the wrong unit. All rates `r` and differences `m` are per `time_unit`. `elapsed_time` and `time_unit` are copied to every row of the output.

## Pseudocount

Before fitting, each count that is zero is replaced by 1, at both time points (`pseudocount="all_zeros"`). Call the result `Nreads`. `Diff` is computed from the original count before the replacement, so a replaced zero has `Nreads = 1` and `Diff = sample_total`. The number of binomial trials is `Nreads + Diff`. It equals `sample_total` except where a zero was replaced, when it is `sample_total + 1`. This follows the published analysis, and `docs/validation.md` records that it reproduces the published numbers.

No read-depth filter is applied. A sequence with all counts zero is still fitted, and is flagged (see Flags).

## A single sequence: `fit_sequences`

For one sequence, let `p_j` be the fraction of reads that belong to it in sample `j`, at time `t_j`. The model is a grouped-binomial regression with a logit link:

```
Nreads_j ~ Binomial(Nreads_j + Diff_j, p_j)
logit(p_j) = a + r * t_j
```

The intercept `a` is the starting abundance on the logit scale. The slope `r` is the growth rate of the sequence: the change of the logit of its read fraction per `time_unit`. Every sample of the sequence is used, so the estimate pools all replicates. Each sequence is fitted separately. The wild type is an ordinary sequence here.

`se_pooled` is the Wald standard error of `r` from the binomial model, and `p_pooled` is the two-sided normal p-value of `r = 0` (see Standard errors).

### Pooled and per-replicate estimates

The binomial model assumes that the counts vary only by sampling. If replicates differ more than that, `se_pooled` is too small. `fit_sequences` therefore also fits each replicate on its own. A replicate has two samples (before and after) and the model has two parameters, so the fit passes exactly through both points and the slope is

```
(logit(p_after) - logit(p_before)) / (t_after - t_before)
```

For each sequence, the slopes of the replicates whose fit converged are summarised by

- `r_replicate_mean`, their mean;
- `r_replicate_sd`, their sample standard deviation (divisor `n - 1`), null when fewer than two replicates converged;
- `se_replicate = r_replicate_sd / sqrt(n_replicates_converged)`, null when the standard deviation is.

`se_replicate` measures the spread of independent replicate slopes, so it includes the differences between replicates. `se_pooled` does not. `n_replicates` is the number of replicates in `samples` (the same for every sequence), and `n_replicates_converged` is how many of them converged for this sequence. A replicate that fails or does not converge is left out of the summary and never changes the `status` of the pooled row.

The pseudocount is applied to each replicate's own two counts, in the same way as for the pooled fit.

## A sequence against the wild type: `fit_wt`

To compare a sequence with a reference (the wild type), the reference rows and the sequence rows are fitted together, with an indicator `type` that is 0 for the reference and 1 for the sequence:

```
Nreads_j ~ Binomial(Nreads_j + Diff_j, p_j)
logit(p_j) = b0 + b1 * type + b2 * t_j + b3 * type * t_j
```

The reference has slope `b2`, and the sequence has slope `b2 + b3`. The interaction `b3` is the fitness difference `m`: how much faster the sequence grows than the reference, per `time_unit`. Positive `m` means enrichment relative to the reference. `fit_wt` runs this fit once for every sequence other than the reference. The reference row of the output has `m = 0`, null `se`, `p` and `p_adj`, and status `reference`.

## A pair: `fit_pair`

`fit_pair` uses the same model for each row of `pairs`, with `from` as the reference (`type = 0`) and `to` as the test sequence (`type = 1`). So `m` is the fitness of `to` relative to `from`. Reversing a pair changes the sign of `m` and leaves `se` and `p` unchanged. The orientation of the pairs in the published data is arbitrary and is used as given.

`fit_wt` is the special case in which every pair has the reference as `from`.

## How `r` and `m` relate

The wild-type model is a sum of two blocks, one for the reference rows and one for the sequence rows, that share no parameter once `type` is included. Therefore the slopes agree with the single-sequence fits. With `r_x` the growth rate of sequence `x` from `fit_sequences`:

```
m = r_to - r_from
se(m) = sqrt(se_pooled_to^2 + se_pooled_from^2)
```

The second line holds under fixed dispersion. The two fits stop their iteration separately, so the identity holds only up to the stopping error. The checks in the test suite use the tolerances in `validation.toml`: a maximum absolute difference of `1e-7` for `m` and `2e-5` for `se`. The same check runs on the full data in the full-data gate described in `docs/validation.md`.

## Fitting, standard errors and p-values

Each model is fitted by iteratively reweighted least squares (Fisher scoring), written to follow R's `glm.fit` for a binomial family with logit link. It starts from `mu = (y * n + 0.5) / (n + 1)`, and stops when the relative change in deviance, `|dev - dev_old| / (|dev| + 0.1)`, is below `1e-8`, or after 25 iterations. The standard error is the square root of the diagonal of the inverse of the weighted cross-product matrix, using the weights that the final iteration started from (R's convention).

The dispersion of the binomial model is fixed at 1 (`dispersion="fixed"`, inherited from Papkou et al. 2023). The test is a Wald test with a normal reference distribution:

```
p = 2 * P(Z > |estimate / se|)
```

`p_adj` is the Benjamini-Hochberg adjustment of `p`, over the non-null `p` of one call of `fit_wt` or `fit_pair`. It is the same as R's `p.adjust(p, "BH")`. Because it depends on the set of comparisons fitted together, it is not comparable between calls with different pair sets. `fit_sequences` reports no adjusted p-value.

`dispersion` is the Pearson statistic of the final fit divided by `df_residual`. It is a diagnostic, and it does not scale `se`. A value well above 1 says that the counts vary more than the binomial model allows. `df_residual` is the number of rows fitted minus the number of coefficients: the number of samples minus 2 for `fit_sequences`, and twice the number of samples minus 4 for `fit_wt` and `fit_pair` (20 for the 12 samples of the Papkou et al. 2023 data). With one replicate per time point `df_residual` is 0 and `dispersion` is null. Both are null on failed rows.

## Status, flags and diagnostics

`status` describes the fit of the row:

| Value | Meaning |
| :-- | :-- |
| `converged` | The iteration stopped on the deviance rule. |
| `not_converged` | 25 iterations were reached. The last estimate is kept. |
| `failed` | The fit raised an error, or `m`, `se` or `p` was not finite. `m`, `se`, `p`, `p_adj`, `dispersion` and `df_residual` are null. |
| `reference` | `fit_wt` only: the reference row. |

`flags` holds data-quality flags, sorted and joined by `|`, or null when there is none. A flag never changes an estimate or a status.

| Flag | Meaning |
| :-- | :-- |
| `censored` | Every raw count of the sequence at the later time is zero. For a comparison, the union over both sequences. |
| `zero_count_replicate` | `fit_sequences` only: some replicate has a raw zero count at either time, so that replicate slope rests on a replaced count. |

`diagnostics` is a struct of `converged` (boolean), `iterations`, `warning` and `error` (strings, null when absent). Several warnings are joined with ` | `. The warnings that can occur are:

- `step size truncated due to divergence`
- `glm.fit: algorithm did not converge`
- `glm.fit: algorithm stopped at boundary value`
- `glm.fit: fitted probabilities numerically 0 or 1 occurred`

The errors that can occur are `rank-deficient design`, `no valid set of coefficients has been found: please supply starting values`, `inner loop 1; cannot correct step size`, `Singular matrix` and `non-finite estimate`. A sequence whose counts are all zero at one time point is a typical source of warnings, because its likelihood has no finite maximum without the pseudocount.

## Engines and `chunk_size`

Two engines run the same iteration and give the same output.

- `engine="vectorised"` (the default) solves many fits at once with numpy. Fits are processed in batches of `chunk_size` sequences or pairs (default 2000). `chunk_size` limits memory and does not change the output.
- `engine="scalar"` fits one sequence or pair at a time. It is slow and simple, and serves as the reference for the vectorised engine.

Both compute the deviance on the proportion scale and solve each step by a Householder QR, as R's `glm.fit` does. Near a boundary the iteration is sensitive to rounding, so the vectorised engine hands every fit that is not a clean convergence (any status other than `converged`, any warning, any error) to the scalar engine and keeps its result. The two engines agree on `m`, `se` and `p` to within 1e-9 in the tests. `docs/validation.md` records the agreement on a sample of the full data.
