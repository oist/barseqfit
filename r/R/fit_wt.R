#' Fit wild-type-relative fitness.
#'
#' One binomial GLM per sequence against the reference genotype:
#' `cbind(Nreads, Diff) ~ 1 + type * time`, with `m` the `type:time`
#' coefficient, `se` its standard error and `p` its Wald p-value. The
#' reference genotype itself returns `m = 0` with `NA` `se` and `p`.
#' The fit uses `samples$time` directly (numeric time in `time_unit`), so `m`
#' is per `time_unit`; `elapsed_time` must equal `max(time) - min(time)`
#' within 1e-9 or the call errors. `status` is `"converged"`, `"not_converged"`
#' (estimate kept), `"failed"` (error or non-finite estimate; `NA` results) or
#' `"reference"`. GLM warnings do not fail a fit; their text goes in
#' `diagnostics`, a `key=value; ...` string (converged, iterations, warning).
#' The WT fit is the set of pairs (reference genotype, x) for every other
#' sequence x, plus the reference row.
#'
#' @param counts data.frame: sequence_id, sample_id, count.
#' @param samples data.frame: sample_id, condition, replicate_id, time, sample_total.
#' @param sequences data.frame: sequence_id, sequence.
#' @param reference_id sequence_id of the reference genotype.
#' @param elapsed_time finite positive number.
#' @param time_unit nonempty string.
#' @param cores positive integer number of worker processes (default 1; no
#'   environment variable is read). Above 1, the non-reference sequences are
#'   fitted in chunks of `chunk_size` with `parallel::mclapply` (forked;
#'   `mc.cores > 1` errors on Windows). Output is identical.
#' @param chunk_size positive integer number of sequences per `mclapply` task
#'   (used only when `cores > 1`; does not change the output).
#' @param pseudocount required; the only accepted value is `"all_zeros"`: every
#'   zero count is replaced by 1 at both time points (the published rule). `Diff`
#'   is computed from the raw count. Any other value, or no value, is an error.
#' @param dispersion required; the only accepted value is `"fixed"`: the binomial
#'   dispersion is fixed at 1 and `p` is a Wald z test. Any other value, or no
#'   value, is an error. The `dispersion` column of the result holds the Pearson
#'   chi-square over the residual df for information only (`NA` for the reference
#'   row and failed fits), and `df_residual` the residual df (rows - 4).
#' @return data.frame: sequence_id, reference_id, m, se, p, p_adj, status,
#'   dispersion, df_residual, flags, diagnostics, elapsed_time, time_unit.
#'   `p_adj` is the Benjamini-Hochberg adjusted `p` over the non-`NA` `p` of the
#'   call (`NA` for the reference row and failed fits).
#'   `flags` is a string of sorted, `|`-joined data-quality flags (`censored`:
#'   every later-time raw count is 0), `NA` when none. A
#'   comparison carries the union over both genotypes and the reference row the
#'   reference genotype's own flags. Flags never change `m`, `se`, `p` or `status`.
#' @export
fit_wt <- function(counts, samples, sequences, reference_id, elapsed_time, time_unit,
                   pseudocount, dispersion, cores = 1L, chunk_size = 2000L) {
  .check_options(pseudocount, dispersion)
  .check_time_args(elapsed_time, time_unit)
  cores <- .check_cores(cores)
  .check_chunk_size(chunk_size)
  by_seq <- .assemble_counts(counts, samples, sequences, elapsed_time)
  if (!(length(reference_id) == 1L && reference_id %in% names(by_seq))) {
    stop("reference_id ", format(reference_id), " is absent from sequences", call. = FALSE)
  }
  ids <- unique(as.character(sequences$sequence_id))
  tests <- ids[ids != reference_id]
  fitted <- if (length(tests)) {
    .fit_comparisons(by_seq, rep(reference_id, length(tests)), tests, cores, chunk_size)
  }
  ref_flags <- .join_flags(.genotype_flags(by_seq)[[reference_id]])
  ref_row <- data.frame(from = reference_id, to = reference_id, m = 0, se = NA_real_,
                        p = NA_real_, status = "reference", dispersion = NA_real_,
                        df_residual = NA_integer_, flags = ref_flags,
                        diagnostics = "reference=TRUE",
                        stringsAsFactors = FALSE)
  res <- rbind(ref_row, fitted)
  # Keep the order of `sequences`: the reference row goes back to its position.
  res <- res[match(ids, res$to), , drop = FALSE]
  names(res)[names(res) == "to"] <- "sequence_id"
  names(res)[names(res) == "from"] <- "reference_id"
  res <- res[c("sequence_id", "reference_id", "m", "se", "p", "status", "dispersion",
                "df_residual", "flags", "diagnostics")]
  res <- .add_p_adj(res)
  res$elapsed_time <- as.numeric(elapsed_time)
  res$time_unit <- time_unit
  rownames(res) <- NULL
  stopifnot(nrow(res) == length(ids))
  res
}
