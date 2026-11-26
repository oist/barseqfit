flag_data <- function() {
  sids <- c("b1", "b2", "a1", "a2")
  tab <- list(wt = c(100, 110, 120, 130), ok = c(90, 95, 200, 210),
              cens = c(50, 60, 0, 0), low = c(2, 3, 4, 5),
              lowafter = c(90, 95, 1, 2), sporadic = c(0, 80, 150, 170))
  samples <- data.frame(sample_id = sids, replicate_id = c("r1", "r2", "r1", "r2"),
                        condition = "c", time = c(0, 0, 2, 2), sample_total = 1000,
                        stringsAsFactors = FALSE)
  counts <- do.call(rbind, lapply(names(tab), function(g) {
    data.frame(sequence_id = g, sample_id = sids, count = tab[[g]], stringsAsFactors = FALSE)
  }))
  sequences <- data.frame(sequence_id = names(tab), sequence = "AAAA", stringsAsFactors = FALSE)
  list(counts = counts, samples = samples, sequences = sequences, tab = tab)
}
wt_fit <- function(d, ...) {
  suppressWarnings(fit_wt(d$counts, d$samples, d$sequences, "wt", elapsed_time = 2,
                          time_unit = "g", pseudocount = "all_zeros", dispersion = "fixed", ...))
}
pair_fit <- function(d, from, to, ...) {
  suppressWarnings(fit_pair(d$counts, d$samples, d$sequences,
                            data.frame(from = from, to = to, stringsAsFactors = FALSE),
                            elapsed_time = 2, time_unit = "g", pseudocount = "all_zeros",
                            dispersion = "fixed", ...))
}
flags_of <- function(res) setNames(res$flags, res$sequence_id)

test_that("only censored genotypes are flagged and the column is in place", {
  res <- wt_fit(flag_data())
  fl <- flags_of(res)
  expect_identical(fl[["cens"]], "censored")
  expect_true(all(is.na(fl[setdiff(names(fl), "cens")])))
  expect_type(res$flags, "character")
  n <- ncol(res)
  expect_equal(names(res)[(n - 4):n],
               c("df_residual", "flags", "diagnostics", "elapsed_time", "time_unit"))
})

test_that("reference row carries the reference genotype's own flags", {
  d <- flag_data()
  expect_true(is.na(wt_fit(d)$flags[wt_fit(d)$sequence_id == "wt"]))
  d$counts$count[d$counts$sequence_id == "wt" & d$counts$sample_id %in% c("a1", "a2")] <- 0
  res <- wt_fit(d)
  expect_identical(res$status[res$sequence_id == "wt"], "reference")
  expect_identical(flags_of(res)[["wt"]], "censored")
  expect_identical(flags_of(res)[["ok"]], "censored")
})

test_that("pair flags are symmetric under reversal", {
  d <- flag_data()
  res <- pair_fit(d, c("ok", "cens", "low", "ok"),
                  c("cens", "ok", "ok", "sporadic"))
  expect_identical(res$flags, c("censored", "censored", NA, NA))
  both <- pair_fit(d, c("cens", "low"), c("low", "cens"))
  expect_identical(both$flags, rep("censored", 2))
})

test_that("censored follows the raw counts, not the pseudocount", {
  d <- flag_data()
  fl <- flags_of(wt_fit(d))
  expect_identical(fl[["cens"]], "censored")
  expect_true(is.na(fl[["sporadic"]]))
  pr <- pair_fit(d, c("ok", "cens"), c("cens", "ok"))
  expect_identical(pr$flags, c("censored", "censored"))
})

test_that("fit_wt and fit_pair agree on flags for (reference, x)", {
  d <- flag_data()
  wt <- wt_fit(d)
  wt <- wt[wt$sequence_id != "wt", ]
  pr <- pair_fit(d, rep("wt", nrow(wt)), wt$sequence_id)
  expect_identical(pr$flags, wt$flags)
})

test_that("fixture flags are NA under defaults and no fixture genotype is censored", {
  res <- fit_wt(td$counts, td$samples, td$sequences, td$reference_id, elapsed_time = td$elapsed_time,
                time_unit = "generations", pseudocount = "all_zeros", dispersion = "fixed")
  expect_true(all(is.na(res$flags)))
  pr <- fit_pair(td$counts, td$samples, td$sequences, td$pairs,
                 elapsed_time = td$elapsed_time, time_unit = "generations",
                 pseudocount = "all_zeros", dispersion = "fixed")
  expect_true(all(is.na(pr$flags)))
  by_seq <- .assemble_counts(td$counts, td$samples, td$sequences, td$elapsed_time)
  expect_false(any(attr(by_seq, "raw_totals")$all_after_zero))
})

test_that("flags are correct with more genotypes than chunk_size", {
  sids <- c("b1", "b2", "a1", "a2")
  tab <- list(wt = c(100, 110, 120, 130))
  expected <- character(0)
  for (i in 0:49) {
    gid <- sprintf("g%02d", i)
    if (i %% 5 == 0) {
      tab[[gid]] <- c(50 + i, 60, 0, 0)
      expected[gid] <- "censored"
    } else if (i %% 7 == 0) {
      tab[[gid]] <- c(2, 3, 4, 5)
      expected[gid] <- NA_character_
    } else {
      tab[[gid]] <- c(90 + i, 95, 200, 210)
      expected[gid] <- NA_character_
    }
  }
  d <- flag_data()
  d$counts <- do.call(rbind, lapply(names(tab), function(g) {
    data.frame(sequence_id = g, sample_id = sids, count = tab[[g]], stringsAsFactors = FALSE)
  }))
  d$sequences <- data.frame(sequence_id = names(tab), sequence = "AAAA", stringsAsFactors = FALSE)
  res <- wt_fit(d, cores = 2, chunk_size = 7)
  got <- flags_of(res)
  expect_true(is.na(got[["wt"]]))
  expect_identical(got[names(expected)], expected)
  pr <- pair_fit(d, rep("wt", length(expected)), names(expected))
  expect_identical(setNames(pr$flags, pr$to), expected)
})

test_that("flag parity literals match the Python test (same hand counts)", {
  d <- flag_data()
  expected <- c(wt = NA, ok = NA, cens = "censored", low = NA, lowafter = NA, sporadic = NA)
  got <- flags_of(wt_fit(d))
  expect_identical(got[names(expected)], expected)
  others <- setdiff(names(expected), "wt")
  pr <- pair_fit(d, rep("wt", length(others)), others)
  expect_identical(setNames(pr$flags, pr$to), expected[others])
})

test_that("unused factor levels in sequence_id do not break the fit", {
  d <- flag_data()
  present <- unique(d$counts$sequence_id)
  d$counts$sequence_id <- factor(d$counts$sequence_id, levels = c(present, "zz"))
  res <- wt_fit(d)
  expect_identical(res$status[res$sequence_id == "ok"], "converged")
  expect_identical(flags_of(res)[["cens"]], "censored")
})
