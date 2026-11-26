fit_td <- function(pairs) {
  fit_pair(td$counts, td$samples, td$sequences, pairs,
           elapsed_time = td$elapsed_time, time_unit = "generations",
           pseudocount = "all_zeros", dispersion = "fixed")
}

test_that("fit_pair reproduces the published pair table", {
  res <- fit_td(td$pairs)
  exp <- td$expected_pair
  expect_equal(nrow(res), 10L)
  expect_equal(res$from, exp$from)
  expect_equal(res$to, exp$to)
  expect_true(all(res$status == "converged"))
  expect_lt(max(abs(res$m - exp$m)), 1e-6)
  expect_lt(max(abs(res$se - exp$se)), 1e-6)
  expect_lt(max(abs(res$p - exp$p)), 1e-6)
  expect_equal(names(res), c("from", "to", "m", "se", "p", "p_adj", "status", "dispersion",
                             "df_residual", "flags", "diagnostics", "elapsed_time",
                             "time_unit"))
})

test_that("reversing a pair negates m and keeps se and p", {
  fwd <- fit_td(td$pairs)
  rev <- fit_td(data.frame(from = td$pairs$to, to = td$pairs$from))
  expect_lt(max(abs(rev$m + fwd$m)), 1e-10)
  expect_lt(max(abs(rev$se - fwd$se)), 1e-10)
  expect_lt(max(abs(rev$p - fwd$p)), 1e-10)
})

test_that("fit_pair rejects self-pairs and duplicate pairs", {
  a <- td$pairs[1, ]
  expect_error(fit_td(data.frame(from = a$from, to = a$from)), "from == to")
  expect_error(fit_td(rbind(a, a)), "duplicate")
})

test_that("fit_pair cores = 2 matches cores = 1", {
  skip_on_os("windows")
  a <- list(td$counts, td$samples, td$sequences, td$pairs, td$elapsed_time, "generations",
       pseudocount = "all_zeros", dispersion = "fixed")
  one <- do.call(fit_pair, c(a, cores = 1))
  two <- do.call(fit_pair, c(a, cores = 2, chunk_size = 3L))
  expect_identical(two, one)
  expect_error(do.call(fit_pair, c(a, chunk_size = 0)), "chunk_size")
  expect_error(do.call(fit_pair, c(a, cores = 0)), "cores")
  expect_error(do.call(fit_pair, c(a, cores = 1.5)), "cores")
})
