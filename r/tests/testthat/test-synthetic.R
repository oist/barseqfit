# Deterministic synthetic fixture (no biological data): six replicate samples
# at each of two times, eight reads per sample. WT 2/8 at both times; variant
# 2/8 at time 0 and 4/8 at time 1. Analytic check: m = log(3),
# se = sqrt(5/12).
test_that("synthetic 24-row fixture matches the reference numbers", {
  samples <- do.call(rbind, lapply(c(0, 1), function(tm) {
    data.frame(sample_id = paste0("t", tm, "-r", 0:5), condition = "constant",
               replicate_id = paste0("r", 0:5), time = tm, sample_total = 8, stringsAsFactors = FALSE)
  }))
  mk <- function(id, per_time) {
    data.frame(sequence_id = id, sample_id = samples$sample_id,
               count = per_time[as.character(samples$time)], stringsAsFactors = FALSE)
  }
  counts <- rbind(mk("wt", c("0" = 2, "1" = 2)), mk("variant", c("0" = 2, "1" = 4)))
  sequences <- data.frame(sequence_id = c("wt", "variant"), sequence = c("AAAA", "AAAT"))
  res <- fit_wt(counts, samples, sequences, reference_id = "wt",
                elapsed_time = 1, time_unit = "arbitrary",
                pseudocount = "all_zeros", dispersion = "fixed")
  v <- res[res$sequence_id == "variant", ]
  expect_equal(v$status, "converged")
  expect_equal(v$m, log(3), tolerance = 1e-5)
  expect_equal(v$se, sqrt(5 / 12), tolerance = 1e-5)
  expect_equal(v$p, 2 * pnorm(-abs(log(3) / sqrt(5 / 12))), tolerance = 1e-5)
  expect_equal(v$elapsed_time, 1)
  expect_equal(v$time_unit, "arbitrary")
})
