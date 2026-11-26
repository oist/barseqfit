# Internal helpers shared by fit_wt() and fit_pair(). Base R only.

.check_frame <- function(x, name, columns) {
  if (!is.data.frame(x)) stop(name, " must be a data.frame", call. = FALSE)
  missing <- setdiff(columns, names(x))
  if (length(missing)) {
    stop(name, " is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (nrow(x) == 0L) stop(name, " must be nonempty", call. = FALSE)
  x
}

.check_time_args <- function(elapsed_time, time_unit) {
  if (!is.numeric(elapsed_time) || length(elapsed_time) != 1L ||
      !is.finite(elapsed_time) || elapsed_time <= 0) {
    stop("elapsed_time must be a finite positive number", call. = FALSE)
  }
  if (!is.character(time_unit) || length(time_unit) != 1L || is.na(time_unit) ||
      !nzchar(trimws(time_unit))) {
    stop("time_unit must be a nonempty string", call. = FALSE)
  }
}

# Validate the required fit options shared by fit_wt() and fit_pair(). Each
# accepts exactly one value; there is no default and no partial matching.
.check_options <- function(pseudocount, dispersion) {
  if (missing(pseudocount)) {
    stop("pseudocount is required; the only accepted value is 'all_zeros'", call. = FALSE)
  }
  if (missing(dispersion)) {
    stop("dispersion is required; the only accepted value is 'fixed'", call. = FALSE)
  }
  if (!(is.character(pseudocount) && length(pseudocount) == 1L &&
        identical(pseudocount, "all_zeros"))) {
    stop("pseudocount must be 'all_zeros' (the only accepted value)", call. = FALSE)
  }
  if (!(is.character(dispersion) && length(dispersion) == 1L &&
        identical(dispersion, "fixed"))) {
    stop("dispersion must be 'fixed' (the only accepted value)", call. = FALSE)
  }
  invisible(TRUE)
}

# Validate inputs and assemble a per-sequence list of data frames with
# Nreads, Diff (from the raw count) and time. Diff = sample_total - count is
# computed before the pseudocount: every zero count of Nreads becomes 1 (the
# "all_zeros" rule, at both time points).
.assemble_counts <- function(counts, samples, sequences, elapsed_time) {
  counts <- .check_frame(counts, "counts", c("sequence_id", "sample_id", "count"))
  samples <- .check_frame(samples, "samples", c("sample_id", "condition", "replicate_id", "time", "sample_total"))
  sequences <- .check_frame(sequences, "sequences", c("sequence_id", "sequence"))
  if (anyDuplicated(samples$sample_id)) stop("samples.sample_id must be unique", call. = FALSE)
  if (anyDuplicated(sequences$sequence_id)) stop("sequences.sequence_id must be unique", call. = FALSE)
  if (anyDuplicated(counts[c("sequence_id", "sample_id")])) {
    stop("counts must have one row per sequence_id and sample_id", call. = FALSE)
  }
  if (!all(counts$sample_id %in% samples$sample_id)) {
    stop("counts contains sample_id values absent from samples", call. = FALSE)
  }
  if (!all(counts$sequence_id %in% sequences$sequence_id)) {
    stop("counts contains sequence_id values absent from sequences", call. = FALSE)
  }
  if (length(unique(samples$condition)) != 1L) {
    stop("exactly one condition is currently supported", call. = FALSE)
  }
  if (!is.numeric(samples$time) || !all(is.finite(samples$time))) {
    stop("samples$time must be numeric and finite", call. = FALSE)
  }
  if (length(unique(samples$time)) != 2L) {
    stop("exactly two sample time points are required", call. = FALSE)
  }
  span <- max(samples$time) - min(samples$time)
  if (abs(span - elapsed_time) > 1e-9) {
    stop("elapsed_time (", format(elapsed_time, digits = 17), ") must equal max(time) - min(time) (",
         format(span, digits = 17), ") within 1e-9", call. = FALSE)
  }
  whole <- function(v) is.numeric(v) && all(is.finite(v)) && all(v == floor(v))
  if (!whole(counts$count) || any(counts$count < 0)) {
    stop("count must be a nonnegative integer", call. = FALSE)
  }
  if (!whole(samples$sample_total) || any(samples$sample_total <= 0)) {
    stop("sample_total must be a positive integer", call. = FALSE)
  }
  expected <- nrow(sequences) * nrow(samples)
  if (nrow(counts) != expected) {
    stop("counts must provide complete sequence/sample coverage", call. = FALSE)
  }
  joined <- merge(counts[c("sequence_id", "sample_id", "count")],
                  samples[c("sample_id", "time", "sample_total")],
                  by = "sample_id", sort = FALSE)
  if (any(joined$count > joined$sample_total)) {
    stop("count cannot exceed sample_total", call. = FALSE)
  }
  joined$Diff <- joined$sample_total - joined$count
  joined$time <- as.numeric(joined$time)
  later <- joined$time == max(joined$time)
  joined$Nreads <- ifelse(joined$count == 0, 1, joined$count)
  # droplevels: unused factor levels must not yield empty entries, which would
  # have no raw totals.
  by_seq <- split(joined[c("Nreads", "Diff", "time")],
                  droplevels(as.factor(joined$sequence_id)))
  # Raw later-time count total per genotype (before the pseudocount) for the
  # censored flag; `all_after_zero` is TRUE when every later-time raw count is 0.
  # Kept as an attribute so by_seq stays a plain list of frames.
  seq_id <- as.character(joined$sequence_id)
  after <- tapply(joined$count[later], seq_id[later], sum)
  attr(by_seq, "raw_totals") <- list(all_after_zero = after[names(by_seq)] == 0)
  by_seq
}

# Data-quality flags per genotype from the raw counts attached by
# .assemble_counts (never the pseudocount-patched counts). Returns a named list
# (by sequence_id) of character vectors, character(0) when clean. The only flag
# is "censored": every later-time raw count is 0.
.genotype_flags <- function(by_seq) {
  censored <- as.vector(attr(by_seq, "raw_totals")$all_after_zero)
  out <- lapply(censored, function(x) if (isTRUE(x)) "censored" else character(0))
  names(out) <- names(by_seq)
  out
}

# Union of flag vectors as one sorted "|"-joined string; NA_character_ if none.
.join_flags <- function(...) {
  u <- sort(unique(unlist(list(...))))
  if (length(u)) paste(u, collapse = "|") else NA_character_
}

# One GLM: reference-level rows first, test-level rows second. `time` is the
# numeric sample time in time_unit, so `m` is per time_unit.
# Returns list(m, se, p, status, dispersion, df_residual, diagnostics). status is "converged",
# "not_converged" (estimate kept) or "failed" (error or non-finite estimate;
# m, se, p, dispersion and df_residual are NA). GLM warnings never fail the fit: the estimate is kept and
# the warning text is recorded in diagnostics ("key=value; ..." string).
# `dispersion` in the result is the raw Pearson phi = sum(w * r^2) / df.residual
# (as summary.glm), reported for information; `df_residual` is rows - 4. The
# dispersion is fixed at 1: se and p are the Wald z values. The diagnostics
# string never mentions dispersion.
.fit_one <- function(ref_data, test_data) {
  failed <- function(msg) {
    list(m = NA_real_, se = NA_real_, p = NA_real_, status = "failed",
         dispersion = NA_real_, df_residual = NA_integer_,
         diagnostics = paste0("converged=NA; iterations=NA; error=", msg))
  }
  warns <- character()
  tryCatch({
    d <- rbind(ref_data, test_data)
    d$type <- factor(rep(c("ref", "test"), c(nrow(ref_data), nrow(test_data))),
                     levels = c("ref", "test"))
    fit <- withCallingHandlers(
      stats::glm(cbind(Nreads, Diff) ~ 1 + type * time,
                 family = stats::binomial(link = "logit"), data = d,
                 control = stats::glm.control(epsilon = 1e-8, maxit = 25)),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    wmsg <- if (length(warns)) gsub(";", ",", paste(unique(warns), collapse = " | ")) else "none"
    diag <- paste0("converged=", isTRUE(fit$converged), "; iterations=", fit$iter,
                   "; warning=", wmsg)
    co <- summary(fit)$coefficients["typetest:time", ]
    est <- co[["Estimate"]]
    se <- co[["Std. Error"]]
    p <- co[["Pr(>|z|)"]]
    # Same quantity as summary.glm: final-iteration working weights and residuals.
    phi_raw <- sum(fit$weights * fit$residuals^2) / fit$df.residual
    # No residual df (one replicate per time point): phi is undefined.
    phi_report <- if (is.finite(phi_raw)) phi_raw else NA_real_
    if (!all(is.finite(c(est, se, p)))) {
      return(list(m = NA_real_, se = NA_real_, p = NA_real_, status = "failed",
                  dispersion = NA_real_, df_residual = NA_integer_,
                  diagnostics = paste0(diag, "; error=non-finite estimate")))
    }
    list(m = unname(est), se = unname(se), p = unname(p),
         status = if (isTRUE(fit$converged)) "converged" else "not_converged",
         dispersion = phi_report, df_residual = as.integer(fit$df.residual),
         diagnostics = diag)
  }, error = function(e) failed(gsub(";", ",", conditionMessage(e))))
}

# Bind a list of equally-named scalar result rows into a data.frame, column-wise
# (one vapply per field; no per-row data.frame). Returns NULL for an empty list.
# Column types follow the first row (character, numeric or integer fields).
.bind_results <- function(rows) {
  if (!length(rows)) return(NULL)
  fields <- names(rows[[1]])
  cols <- lapply(fields, function(f) {
    proto <- rows[[1]][[f]]
    vapply(rows, function(r) r[[f]], FUN.VALUE = proto, USE.NAMES = FALSE)
  })
  names(cols) <- fields
  as.data.frame(cols, stringsAsFactors = FALSE, optional = TRUE)
}

# Shared per-pair loop behind fit_wt() and fit_pair(). `from_ids` and `to_ids`
# are parallel character vectors of ids present in `by_seq`; `from` is the
# reference level of each comparison. Returns a data.frame of one row per pair,
# in order: from, to, m, se, p, status, dispersion, df_residual, flags, diagnostics. Callers add the
# elapsed_time and time_unit columns and any renaming. Returns NULL for zero
# pairs.
.fit_comparisons <- function(by_seq, from_ids, to_ids, cores, chunk_size) {
  stopifnot(length(from_ids) == length(to_ids))
  chunks <- .parallel_map(seq_along(from_ids), function(i) {
    from <- from_ids[i]
    to <- to_ids[i]
    c(list(from = from, to = to), .fit_one(by_seq[[from]], by_seq[[to]]))
  }, cores, chunk_size, collect = .bind_results)
  res <- do.call(rbind, chunks)
  if (!is.null(res)) {
    flags <- .genotype_flags(by_seq)
    res$flags <- vapply(seq_along(from_ids), function(i) {
      .join_flags(flags[[from_ids[i]]], flags[[to_ids[i]]])
    }, character(1))
    cols <- names(res)
    res <- res[c(setdiff(cols, c("flags", "diagnostics", "elapsed_time", "time_unit")),
                 "flags", intersect(c("diagnostics", "elapsed_time", "time_unit"), cols))]
  }
  res
}

# Benjamini-Hochberg adjusted p-values over the non-NA entries (NA passes
# through). stats::p.adjust drops NAs before its default `n = length(p)` is
# evaluated, so n is the number of non-NA p, as in barseqfit._adjust.
.bh <- function(p) {
  stats::p.adjust(p, method = "BH")
}

# Insert `p_adj` right after `p` (BH over the non-NA p of the call).
.add_p_adj <- function(res) {
  res$p_adj <- .bh(res$p)
  cols <- names(res)
  res[c(setdiff(cols[seq_len(match("p", cols))], "p_adj"), "p_adj",
        setdiff(cols[-seq_len(match("p", cols))], "p_adj"))]
}
