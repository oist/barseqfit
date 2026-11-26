# Validation

barseqfit is validated by reproducing published estimates, not by simulation.
The repository contains no simulation test. The claim is that `fit_wt` and
`fit_pair` reproduce the fitness estimates published with Papkou et al. 2023
(Zenodo record 13837982) within a stated tolerance, and that `fit_sequences`
agrees with `fit_wt` when the reference sequence is the same.

The raw data are not in this repository. The fixture in `tests/data/` is a small
subset of the published tables; the full-data checks read the Papkou2023 store
from the path in `PAPKOU2023RUGGED_DATA` (see `.env.example`).

The `pseudocount` and `dispersion_model` output columns are Python-only: they
echo the arguments passed to the fit. The R reference package does not output
them, and the comparisons against R and the published tables exclude them.

## Release-candidate rerun

PLACEHOLDER: the numbers in this section will be filled in from a rerun of the
full-data gates on the release candidate before release. It will state the
commit, date, machine and tolerances, and cover:

- `fit_wt` and `fit_pair` against the published tables on the full data
  (261,332 non-reference genotypes; full edge table), through R `glm` and
  directly;
- the `fit_sequences` identity with `fit_wt`, and scalar against vectorised
  agreement;
- the known deviation: the 3 genotypes with a zero before-count, where the
  default pseudocount rule matches the published tables.

Until then the dated records below stand. They were produced on 2026-10-14, and
the rerun replaces them.

## Source of the reference data

- Zenodo record 13837982, DOI 10.5281/zenodo.13837982.
- Creator: Andrei Papkou.
- Licence: CC-BY-4.0.
- The data files come from the Zenodo record. The authors' 13 fitness scripts
  do not: they come from the authors' pipeline repository, archived on Zenodo
  as record 8229097 (DOI 10.5281/zenodo.8229097, CC-BY-4.0, file
  `ampseq_nextflow-master.zip`), not from record 13837982. They are not
  redistributed here. The local store holds a copy under `ampseq-nextflow/`.
  The recipe `just validate-papkou2023` runs it from the local store, after
  checking each of the 13 files against the md5 sums in
  `r/validation/authors_scripts.md5` (paths relative to `ampseq-nextflow/`).
  The sums are those of the Zenodo archive. `lib/fit_model_wt.R` differs from an
  earlier local copy only by one trailing space on line 51; the sum is pinned to
  the Zenodo version.
- `scripts/calculate_fitness/calculate_fitness_from_to_chunk.R` sources
  `lib/fit_model_from_to.v2.R` (a hard-coded absolute path on the authors'
  machine), not `lib/fit_model_from_to.R`. The v2 library is therefore the pair
  model validated against. `lib/fit_model_wt.R` is the WT model used by
  `calculate_fitness_wt*.R`.

## Fixture verdict (record of 2026-10-14)

Produced by `just validate-papkou2023` (`r/validation/run_authors_pipeline.R`,
R 4.1.3, `data.table` from the renv library). It sources the authors'
`fit_model_wt.R` and `fit_model_from_to.v2.R` untouched and fits the 20
non-reference genotypes of the fixture (against `GCCGATCTC`) and the 10 pairs in
`tests/data/`. Tolerance 1e-6; the reference row (m = 0, NA se and p) is
excluded. Maximum absolute residuals against `expected_wt.csv` and
`expected_pair.csv`:

| time | pseudocount | WT m | WT se | WT p | pair m | pair se | pair p | pass |
| :-- | :-- | --: | --: | --: | --: | --: | --: | :-: |
| 3.51662 | both, Diff before | 7.3e-08 | 7.5e-09 | 5.0e-19 | 6.4e-08 | 9.2e-09 | 1.4e-14 | yes |
| 3.51662 | both, Diff after | 8.0e-08 | 7.5e-09 | 5.0e-19 | 6.4e-08 | 9.2e-09 | 6.8e-09 | yes |
| 3.51662 | after-only, Diff after | 8.0e-08 | 7.5e-09 | 5.0e-19 | 6.4e-08 | 9.2e-09 | 6.8e-09 | yes |
| 3.51662 | after-only, Diff before | 7.3e-08 | 7.5e-09 | 5.0e-19 | 6.4e-08 | 9.2e-09 | 1.4e-14 | yes |
| 3.51662 | none | 2.0e-01 | 3.7e-02 | 3.1e-17 | 2.0e-01 | 3.2e-02 | 4.3e-02 | no |
| 3.5166197103891137 | both, Diff before | 4.6e-14 | 2.4e-15 | 4.6e-19 | 4.1e-15 | 4.6e-16 | 4.7e-14 | yes |
| 3.5166197103891137 | both, Diff after | 6.8e-09 | 1.4e-12 | 4.6e-19 | 6.8e-09 | 9.7e-16 | 6.8e-09 | yes |
| 3.5166197103891137 | after-only, Diff after | 6.8e-09 | 1.4e-12 | 4.6e-19 | 6.8e-09 | 9.7e-16 | 6.8e-09 | yes |
| 3.5166197103891137 | after-only, Diff before | 4.6e-14 | 2.4e-15 | 4.6e-19 | 4.1e-15 | 4.6e-16 | 4.7e-14 | yes |
| 3.5166197103891137 | none | 2.0e-01 | 3.7e-02 | 3.1e-17 | 2.0e-01 | 3.2e-02 | 4.3e-02 | no |

Definitions:

- "Diff before" means `Diff = Total - Nreads` computed before the pseudocount (as
  in `01.prepare_data_for_fitness.r`). "Diff after" means it is recomputed after
  the pseudocount (as in the older `lib/fit_model_from_to.R`).
- "both" and "after-only" say which timepoints get `Nreads[Nreads == 0] <- 1`.
- `Total` is `sample_total` in `tests/data/samples.csv`, the per-sample
  mapped-read row sum of the full unfiltered matrices (before: all 261,378
  genotypes; after: all 261,333).
- The first time is the rounded 3.51662, the second the recovered value (next
  section). Both are applied explicitly to the after samples.
- The script also recomputes the totals from the store as a cross-check. Those
  rows are identical to the fixture rows above, so the fixture totals equal the
  unfiltered store row sums.

Best candidate (residuals about 1e-14, the exact published numbers): total equal
to the unfiltered mapped-read row sum; `Diff = Total - Nreads` computed before
the pseudocount, then `Nreads[Nreads == 0] <- 1`; time 3.5166197103891137. The
rounded 3.51662 leaves residuals up to 7e-8, still within 1e-6. Every candidate
except "none" passes at 1e-6. Using no pseudocount fails because the after
counts contain zeros.

What the fixture can and cannot distinguish:

- It cannot distinguish the both-timepoint pseudocount from the after-only one:
  no before-count is zero in the subset (only after counts are), so the two rules
  give identical numbers. The both-timepoint choice follows the authors'
  `01.prepare_data_for_fitness.r`, not the fixture. The full-data records below
  do distinguish them.
- It does distinguish the Diff order: "Diff before" reaches about 1e-14 while
  "Diff after" stalls at about 6.8e-9.

## Generation time (record of 2026-10-14)

The after samples use the generation time 3.5166197103891137, recovered so that
the authors' code reproduces the published numbers. The check is the fixture
verdict above: with the rounded 3.51662 the residuals reach 7e-8, with
3.5166197103891137 they fall to about 1e-14 ("both, Diff before"). The recovered value is stored
as `time` in `tests/data/samples.csv` (numeric, in generations) and as the
constant `.papkou_generation_time` in `r/R/papkou.R`, and the fits use it
directly.

## Full-data WT verdict (record of 2026-10-14)

Produced by `just validate-full-wt` (export, R run, Python gate) on 2026-10-14.
Inputs are regenerated from the store into the git-ignored `$BARSEQFIT_SCRATCH`
and never committed. 261,333 genotypes (261,332 non-reference) by 12 samples.
Tolerance: the fixture rule (abs 1e-6, rel 1e-6 when |expected| < 1; p below
1e-300 on both sides counts as equal), with a threshold-tie rule for `se` and
`p`. Verdict: PASS.

| comparison | m max abs | se max abs | p max abs | rows outside (se/p) | note |
| :-- | --: | --: | --: | --: | :-- |
| Python vs published | 5.7e-09 | 4.3e-06 | 4.9e-13 | 8 | all excused (tie or R outlier) |
| Python vs R | 1.4e-08 | 9.8e-06 | 1.2e-13 | 18 | 18 threshold ties; 0 on 261,314 matched rows |
| R vs published | | 9.8e-06 | | 16 | same BLAS-level mechanism |

Scalar against vectorised (1,000 seeded genotypes): 0 entries over 1e-9, maximum
5.6e-14. Status: 261,332 converged, 0 disagreements, no warnings, 0 scalar
fallbacks. Pseudocount (3 genotypes with a zero before-count): the both-timepoint
rule matches the published tables to 1e-14; the after-only rule is off by 2.6e-3
in m.

Timings (one machine, Apple silicon, not a budget):

| step | time |
| :-- | --: |
| export (R) | about 5 s |
| R `fit_wt`, 4 cores (188% CPU) | 1281 s |
| Python vectorised, 261,333 genotypes | about 11 s on an 8-core Apple Silicon machine, single run |
| Python scalar, 1,000 genotypes | 0.6 s |

WT evidence refreshed on 2026-10-14, after the Parquet and run-config migration:
same gate, PASS with a noise floor of 1e-13 on m. R `fit_wt` on 4 cores took
285.9 s (R 4.1.3), Python vectorised 9.0 s, the scalar sample (1,000 genotypes)
0.55 s, with 0 fallbacks. The 1281 s R timing above is from the earlier run.

### Export manifest (record of 2026-10-14)

`export_manifest.json`, R 4.1.3, md5 from the store's `MANIFEST.md5`:

| store file | md5 |
| :-- | :-- |
| `zenodo/in_data/sequence_reads_counts.RDS` | 821a5866d34c540eb528de805624a9e6 |
| `zenodo/in_data/fitness_data_wt.rds` | e636390a6681aeda8ed90dc54fdc17d2 |
| `zenodo/in_data/fitness_diff_wald_test_data.RDS` | 20bfb173c6f0ea8815b08d140b617f7a |

## Edge table (published pair fits)

- File: `zenodo/in_data/fitness_diff_wald_test_data.RDS`, md5
  20bfb173c6f0ea8815b08d140b617f7a (from `export_manifest.json`, taken from the
  store's `MANIFEST.md5`).
- 3,518,696 rows; columns `from`, `to`, `m`, `m_se`, `m_p.value`; one
  orientation per unordered Hamming-1 pair; 261,333 distinct genotypes; no
  status column. `edge_rows_dropped` is 0 in the export manifest.
- The orientation is copied as data. It cannot be derived: string order predicts
  `from` for 44.9% of rows, and WT-table row order for 44.9%.

## Full pair verdict (record of 2026-10-14)

`verdict_full_pair.md`. Verdict: PASS, 0 asserted failures. Noise floor 1e-13 on
m.

- Rows: 3,518,696 fitted, 3,518,696 published and R rows. Set differences
  (fitted, published, R, in all six directions): 0.
- m against published and against R: 0 rows outside tolerance (maximum abs
  3.8e-9, threshold ties with |m| about 0.09).
- Threshold ties (Python and R iteration counts differ): 49 rows, cap 1,759
  (0.05%); maximum se residual 3.4e-6, maximum p residual 2.1e-5 (caps 1e-4).
- R against published: 45 rows with se or p outside tolerance (cap 1,759).
  Python against published: 44 rows outside, all excused by a tie or an R
  outlier. On the 3,518,627 rows that are neither, 0 are outside tolerance. 0
  are outside on the 3,518,647 rows matched against R.
- p underflow (both sides below 1e-300, treated as equal): 1,758 rows against R
  and against published.
- Pseudocount check (66 pairs touch the 3 genotypes with a zero before-count):
  the both-timepoint rule matches published (maximum abs m 5.6e-15, se 8.9e-16,
  p 1.8e-14, 0 rows outside); the after-only rule is off by 2.6e-3 in m with all
  66 rows outside.
- Scalar against vectorised (1,000 seeded pairs): 0 entries over 1e-9, maximum
  abs difference 2.5e-13; 0 status or iteration mismatches; 0 scalar fallbacks.

Pair timings (one machine, Apple silicon, not a budget):

| step | time |
| :-- | --: |
| R `fit_pair`, 4 cores, R 4.1.3 (`timing_pair_r.json`) | 7050 s (about 2 h) |
| Python vectorised, 3,518,696 pairs (`timing_pair_python.json`) | 122.8 s |
| Python scalar, 1,000 sampled pairs | 4.4 s |

The first full R run was abandoned after 3 h because of parent-process memory
pressure; binding the R fit results per chunk in the workers fixed it.
