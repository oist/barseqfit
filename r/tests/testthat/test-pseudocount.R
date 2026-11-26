# Zero counts: Diff = sample_total - raw count (so Diff = total for a zero),
# then Nreads 0 -> 1 at BOTH timepoints (authors' 01.prepare_data_for_fitness.r).
test_that("zero counts get Nreads = 1 at both timepoints with Diff from the raw count", {
  samples <- data.frame(
    sample_id = c("b1", "b2", "a1", "a2"), condition = "c",
    replicate_id = c("r1", "r2", "r1", "r2"), time = c(0, 0, 2, 2),
    sample_total = c(100, 120, 110, 90), stringsAsFactors = FALSE)
  sequences <- data.frame(sequence_id = c("wt", "v"), sequence = c("AAAA", "AAAT"))
  cnt <- list(wt = c(b1 = 10, b2 = 0, a1 = 12, a2 = 9),
              v  = c(b1 = 0, b2 = 5, a1 = 0, a2 = 20))
  counts <- do.call(rbind, lapply(names(cnt), function(id) {
    data.frame(sequence_id = id, sample_id = samples$sample_id,
               count = unname(cnt[[id]][samples$sample_id]), stringsAsFactors = FALSE)
  }))
  res <- fit_wt(counts, samples, sequences, "wt", elapsed_time = 2, time_unit = "g",
                pseudocount = "all_zeros", dispersion = "fixed")
  v <- res[res$sequence_id == "v", ]

  manual <- function(nreads_fun, diff_fun) {
    d <- do.call(rbind, lapply(names(cnt), function(id) {
      raw <- unname(cnt[[id]][samples$sample_id])
      data.frame(Nreads = nreads_fun(raw, samples$time), Diff = diff_fun(raw, samples$sample_total),
                 time = samples$time, type = factor(id, levels = c("wt", "v")))
    }))
    co <- summary(stats::glm(cbind(Nreads, Diff) ~ 1 + type * time, family = binomial, data = d,
                             control = stats::glm.control(epsilon = 1e-8, maxit = 25)))$coefficients
    co["typev:time", c("Estimate", "Std. Error")]
  }
  both <- function(raw, tm) ifelse(raw == 0, 1, raw)
  after_only <- function(raw, tm) ifelse(raw == 0 & tm > 0, 1, raw)
  diff_raw <- function(raw, total) total - raw                       # Diff = total for a zero
  diff_after <- function(raw, total) total - ifelse(raw == 0, 1, raw) # total - 1 for a zero

  ref <- manual(both, diff_raw)
  expect_lt(abs(v$m - ref[["Estimate"]]), 1e-10)
  expect_lt(abs(v$se - ref[["Std. Error"]]), 1e-10)
  # The alternatives are distinguishable on this data, so the match above is informative.
  expect_gt(abs(v$m - manual(both, diff_after)[["Estimate"]]), 1e-6)
  expect_gt(abs(v$m - manual(after_only, diff_raw)[["Estimate"]]), 1e-6)
})

# The only rule is all_zeros: every zero cell of Nreads becomes 1, Diff stays raw.
scope_data <- function() {
  samples <- data.frame(
    sample_id = c("b1", "b2", "a1", "a2"), condition = "c",
    replicate_id = c("r1", "r2", "r1", "r2"), time = c(0, 0, 2, 2),
    sample_total = 100, stringsAsFactors = FALSE)
  cnt <- list(wt = c(10, 12, 11, 13),         # no zeros
              sporadic = c(0, 5, 7, 9),       # one zero before
              afterzero = c(8, 6, 0, 0),      # all-zero after
              beforezero = c(0, 0, 4, 6))     # all-zero before
  counts <- do.call(rbind, lapply(names(cnt), function(id) {
    data.frame(sequence_id = id, sample_id = samples$sample_id, count = cnt[[id]],
               stringsAsFactors = FALSE)
  }))
  list(samples = samples, counts = counts,
       sequences = data.frame(sequence_id = names(cnt), sequence = "AAAA"))
}

test_that("all_zeros patches every zero cell and Diff comes from the raw count", {
  d <- scope_data()
  expected <- list(sporadic = c(1, 5, 7, 9), afterzero = c(8, 6, 1, 1),
                   beforezero = c(1, 1, 4, 6))
  by_seq <- .assemble_counts(d$counts, d$samples, d$sequences, 2)
  for (id in names(expected)) {
    got <- by_seq[[id]]
    exp <- expected[[id]]
    # Rows are in merge order; compare within each time (before = b1, b2).
    expect_equal(sort(got$Nreads[got$time == 0]), sort(exp[1:2]), info = id)
    expect_equal(sort(got$Nreads[got$time == 2]), sort(exp[3:4]), info = id)
    expect_equal(sort(got$Diff), sort(100 - d$counts$count[d$counts$sequence_id == id]),
                 info = id)
  }
})

test_that("pseudocount and dispersion are required and accept exactly one value each", {
  d <- scope_data()
  args <- list(d$counts, d$samples, d$sequences, "wt", elapsed_time = 2, time_unit = "g")
  pargs <- list(d$counts, d$samples, d$sequences, data.frame(from = "wt", to = "sporadic"),
                elapsed_time = 2, time_unit = "g")
  for (fn in list(list(fit_wt, args), list(fit_pair, pargs))) {
    f <- fn[[1]]
    a <- fn[[2]]
    expect_error(do.call(f, c(a, dispersion = "fixed")), "pseudocount is required")
    expect_error(do.call(f, c(a, pseudocount = "all_zeros")), "dispersion is required")
    expect_error(do.call(f, a), "pseudocount is required")
    for (bad in c("after_zeros", "separated_only", "none", "all", "")) {
      expect_error(do.call(f, c(a, pseudocount = bad, dispersion = "fixed")),
                   "pseudocount must be 'all_zeros'", info = bad)
    }
    for (bad in c("quasi", "fix", "")) {
      expect_error(do.call(f, c(a, pseudocount = "all_zeros", dispersion = bad)),
                   "dispersion must be 'fixed'", info = bad)
    }
    expect_error(do.call(f, c(a, list(pseudocount = c("all_zeros", "all_zeros"),
                                   dispersion = "fixed"))), "pseudocount must be")
    expect_error(do.call(f, c(a, pseudocount = NA_character_, dispersion = "fixed")),
                 "pseudocount must be")
  }
  expect_error(.check_options("none", "fixed"), "pseudocount")
  expect_error(.check_options("all_zeros", "quasi"), "dispersion")
  expect_true(.check_options("all_zeros", "fixed"))
})

test_that("the removed read-filter arguments are rejected", {
  d <- scope_data()
  expect_error(fit_wt(d$counts, d$samples, d$sequences, "wt", elapsed_time = 2, time_unit = "g",
                      pseudocount = "all_zeros", dispersion = "fixed", min_before_reads = 10),
               "unused argument")
})
