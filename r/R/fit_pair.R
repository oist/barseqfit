#' Fit pairwise (edge) fitness.
#'
#' One binomial GLM per pair, in the given orientation: `from` is the
#' reference level and `to` the test level. `m` is the `type:time`
#' coefficient (fitness of `to` relative to `from`). Reversing a pair negates
#' `m` and leaves `se` and `p` unchanged.
#'
#' @param pairs data.frame: from, to (sequence_ids).
#' @inheritParams fit_wt
#' @return data.frame: from, to, m, se, p, p_adj, status, dispersion, df_residual,
#'   flags, diagnostics, elapsed_time, time_unit.
#' @export
fit_pair <- function(counts, samples, sequences, pairs, elapsed_time, time_unit,
                     pseudocount, dispersion, cores = 1L, chunk_size = 2000L) {
  .check_options(pseudocount, dispersion)
  .check_time_args(elapsed_time, time_unit)
  cores <- .check_cores(cores)
  .check_chunk_size(chunk_size)
  pairs <- .check_frame(pairs, "pairs", c("from", "to"))
  by_seq <- .assemble_counts(counts, samples, sequences, elapsed_time)
  if (any(as.character(pairs$from) == as.character(pairs$to))) {
    stop("pairs must not compare a sequence with itself (from == to)", call. = FALSE)
  }
  if (anyDuplicated(data.frame(from = as.character(pairs$from), to = as.character(pairs$to)))) {
    stop("pairs must not contain duplicate (from, to) rows", call. = FALSE)
  }
  ids <- c(as.character(pairs$from), as.character(pairs$to))
  if (!all(ids %in% names(by_seq))) {
    stop("pairs contains sequence_id values absent from sequences", call. = FALSE)
  }
  res <- .fit_comparisons(by_seq, as.character(pairs$from), as.character(pairs$to),
                          cores, chunk_size)
  res <- .add_p_adj(res)
  res$elapsed_time <- as.numeric(elapsed_time)
  res$time_unit <- time_unit
  res
}
