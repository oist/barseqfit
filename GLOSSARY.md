# BarSeqFit

BarSeqFit describes barcode sequencing measurements used to estimate sequence fitness across samples and time points.

## Measurements and sequences

**Barcode**:
A DNA sequence used to identify a construct whose abundance is measured by sequencing.

**Sequence**:
A DNA molecule associated with a stable identity and a measured fitness effect. A sequence may be the wild type or a variant.

**Wild type (WT)**:
The reference sequence used as the baseline for interpreting variant fitness.
_Avoid_: control sequence, unless referring to a distinct experimental control.

**Variant**:
A sequence that differs from the wild type and whose fitness is compared with it.

**Sample**:
A sequencing measurement associated with a biological replicate, a condition and a time point.

**Read count**:
The number of sequencing reads assigned to a barcode in a sample.

**Elapsed time**:
The duration between the two time points of an experiment, reported with an explicit unit.

**Generation time**:
The elapsed time between the two time points of Papkou et al. 2023 in generations, 3.5166197103891137. It is a constant recovered from the published analysis so that the authors' code reproduces the published tables, and it is stored as the `time` of the after samples (in `tests/data/samples.csv` and as a constant in the R helpers); it is not derived from the counts. The rounded 3.51662 reproduces the published tables only to about 1e-7.

## Fitness and growth

**Growth rate**:
`r`, the change in the logit of a sequence's read fraction per unit of time, estimated by `fit_sequences` from that sequence's own samples.

**Fitness difference**:
`m`, the estimated difference between the growth rates of two sequences per unit of time. Positive and negative values indicate relative enrichment and depletion. It equals the difference of the two growth rates up to numerical tolerance.

**Reference sequence**:
The sequence whose fitness is the baseline of a fit. In the wild-type fit it is the WT sequence (`GCCGATCTC` in the Papkou et al. 2023 data), reported with `m = 0` and missing `se` and `p`.

**Pair**:
An ordered tuple `(from, to)` of two sequences, fitted by its own binomial regression. `m` is the fitness of `to` relative to `from`. Any two sequences may form a pair. The wild-type fit is the set of pairs (reference sequence, x) for every other sequence x, plus the reference row.

**Comparison**:
The fit of one pair, whichever public function (`fit_wt` or `fit_pair`) asked for it.

**Edge**:
A pair whose sequences lie at a chosen Hamming distance in a sequence network. The published Papkou et al. 2023 edge table uses distance one. Reversing an edge negates `m` and leaves `se` and `p` unchanged.

**Orientation**:
Which sequence of a pair is `from`. The published orientation is arbitrary and is used as given.

**Replicate fit**:
The fit of one sequence within a single replicate, from the replicate's before and after samples, made by `fit_sequences`. With two samples and two parameters it passes through both points, so it has no residual degrees of freedom.

**Replicate SE**:
`se_replicate`, the standard error of a sequence's growth rate from the spread of its converged replicate fits, `r_replicate_sd / sqrt(n_replicates_converged)`. Null below two converged replicates. Unlike the pooled `se_pooled`, it reflects differences between replicates.

## Model choices inherited from Papkou et al. 2023

**Pseudocount**:
The replacement of every zero read count by one, at both time points. The complement count is computed from the raw count before the replacement.

**Dispersion**:
The Pearson statistic of a fit divided by its residual degrees of freedom, reported in the `dispersion` column as a diagnostic. Binomial dispersion is fixed at 1 for the standard error and the p-value (`dispersion="fixed"`). Null when there are no residual degrees of freedom.

**Residual degrees of freedom**:
The number of fitted rows minus the number of coefficients (`df_residual`).

**Censored sequence**:
A sequence whose raw counts are all zero at the later time. It is flagged `censored`.

**Flags**:
The `flags` column: data-quality flag names (`censored`, `zero_count_replicate`), sorted and joined by `|`, null when none. Flags are not `status` values and never change an estimate.

**zero_count_replicate**:
The flag set by `fit_sequences` on a sequence when any of its replicates has a raw zero count at either time, so that its replicate slope rests on a replaced count.

**Adjusted p**:
`p_adj`, the Benjamini-Hochberg false-discovery-rate adjusted `p` over the non-null `p` of one `fit_wt` or `fit_pair` call. Null for reference and failed rows. It depends on the set of comparisons fitted together.

## Engines

**Engine**:
How the regressions are run: `vectorised` (default, batched over many sequences or pairs) or `scalar` (one at a time).

**Oracle engine**:
The scalar engine, the slow, simple reference that the vectorised engine is checked against and that re-fits any result that is not a clean convergence.

## Validation

**Export**:
The full-data tables (counts, samples, sequences, published WT table, pair list and published edge table) that R writes from the Papkou et al. 2023 data into a scratch directory, for the full-data checks.

**Gate**:
A full-data comparison against the published tables and R that applies the tolerance rule (`just full-wt-gate` for the WT fit, `just full-pair-gate` for the pair fit). It fails only on an asserted failure; other findings are reported.

**Tolerance rule**:
The comparison rule as one value: absolute and relative tolerance, a noise floor on `m`, and a clause that treats `p` below 1e-300 on both sides as equal.

**Gate rule**:
The tolerance rule together with the tie fraction, the tie limits on `se` and `p`, and the scalar tolerance that a gate applies.

**Threshold tie**:
A sequence whose iteration count differs between two implementations because its deviance convergence ratio sits at the `1e-8` stopping threshold. Its `se` can differ by up to about 1e-5. A gate counts and bounds such ties rather than holding them to 1e-6.

**Run config**:
The committed file `validation.toml` of result-affecting parameters of the full-data checks (chunk size, tolerances, caps). Machine-specific paths and core counts stay out of it.

**Machine environment**:
The machine-specific values of a run: the data path (`PAPKOU2023RUGGED_DATA`), the scratch directory (`BARSEQFIT_SCRATCH`) and the core count (`BARSEQFIT_CORES`), taken from the shell and a repository-root `.env` file. A shell value, even empty, wins; empty means unset.
