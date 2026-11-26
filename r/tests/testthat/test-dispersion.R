elapsed <- td$elapsed_time
fit_disp <- function(pairs, dispersion) {
  fit_pair(td$counts, td$samples, td$sequences, pairs,
           elapsed_time = elapsed, time_unit = "generations",
           pseudocount = "all_zeros", dispersion = dispersion)
}

test_that("fixed dispersion leaves the published pair table unchanged", {
  default <- fit_disp(td$pairs, "fixed")
  exp <- td$expected_pair
  expect_lt(max(abs(default$se - exp$se)), 1e-6)
  expect_lt(max(abs(default$p - exp$p)), 1e-6)
  expect_false(anyNA(default$dispersion))
  expect_true(all(default$df_residual == 20L))
  expect_false(grepl("dispersion", paste(default$diagnostics, collapse = " ")))
})

test_that("reference and failed rows have NA dispersion and df_residual", {
  d <- data.frame(Nreads = c(-1, 2, 3, 4), Diff = c(4, 3, 2, 1), time = c(0, 0, 1, 1))
  bad <- .fit_one(d, d)  # negative count: glm() errors
  expect_equal(bad$status, "failed")
  expect_true(is.na(bad$dispersion) && is.na(bad$df_residual))
  ref <- fit_wt(td$counts, td$samples, td$sequences, td$reference_id, elapsed,
                "generations", pseudocount = "all_zeros", dispersion = "fixed")
  r <- ref[ref$status == "reference", ]
  expect_true(is.na(r$dispersion) && is.na(r$df_residual))
})

saturated_data <- function() {
  samples <- data.frame(sample_id = c("b", "a"), replicate_id = "r1", condition = "c",
                        time = c(0, 2), sample_total = 1000, stringsAsFactors = FALSE)
  tab <- list(wt = c(100, 150), x = c(80, 260))
  counts <- do.call(rbind, lapply(names(tab), function(g) {
    data.frame(sequence_id = g, sample_id = c("b", "a"), count = tab[[g]],
               stringsAsFactors = FALSE)
  }))
  sequences <- data.frame(sequence_id = names(tab), sequence = "AAAA", stringsAsFactors = FALSE)
  list(counts = counts, samples = samples, sequences = sequences, tab = tab)
}

test_that("zero residual df keeps the fixed fit", {
  d <- saturated_data()
  lo <- function(n) log(n / (1000 - n))
  m <- ((lo(d$tab$x[2]) - lo(d$tab$x[1])) - (lo(d$tab$wt[2]) - lo(d$tab$wt[1]))) / 2
  se <- sqrt(sum(vapply(unlist(d$tab), function(n) 1 / n + 1 / (1000 - n), 0))) / 2
  fixed <- fit_wt(d$counts, d$samples, d$sequences, "wt", elapsed_time = 2, time_unit = "g",
                  pseudocount = "all_zeros", dispersion = "fixed")
  row <- fixed[fixed$sequence_id == "x", ]
  expect_identical(row$status, "converged")
  expect_equal(row$m, m, tolerance = 1e-8)
  expect_equal(row$se, se, tolerance = 1e-8)
  expect_equal(row$p, 2 * pnorm(-abs(m / se)), tolerance = 1e-8)
  expect_true(is.na(row$dispersion))
  expect_identical(row$df_residual, 0L)
})
