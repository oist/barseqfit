test_that("fit_wt reproduces the published WT table", {
  res <- fit_wt(td$counts, td$samples, td$sequences, reference_id = td$reference_id,
                elapsed_time = td$elapsed_time, time_unit = "generations",
                pseudocount = "all_zeros", dispersion = "fixed")
  exp <- td$expected_wt
  expect_equal(sort(res$sequence_id), sort(exp$sequence_id))
  res <- res[match(exp$sequence_id, res$sequence_id), ]
  nonref <- exp$sequence_id != td$reference_id
  expect_equal(res$status[nonref], rep("converged", sum(nonref)))
  expect_lt(max(abs(res$m[nonref] - exp$m[nonref])), 1e-6)
  expect_lt(max(abs(res$se[nonref] - exp$se[nonref])), 1e-6)
  expect_lt(max(abs(res$p[nonref] - exp$p[nonref])), 1e-6)
  wt <- res[res$sequence_id == td$reference_id, ]
  expect_equal(wt$status, "reference")
  expect_equal(wt$m, 0)
  expect_true(is.na(wt$se))
  expect_true(is.na(wt$p))
  expect_equal(names(res), c("sequence_id", "reference_id", "m", "se", "p", "p_adj", "status",
                             "dispersion", "df_residual", "flags", "diagnostics", "elapsed_time",
                             "time_unit"))
  expect_true(is.na(wt$dispersion) && is.na(wt$df_residual))
  expect_equal(unique(res$time_unit), "generations")
})

test_that("fit_wt validates arguments", {
  a <- list(td$counts, td$samples, td$sequences)
  expect_error(fit_wt(a[[1]], a[[2]], a[[3]], "nope", td$elapsed_time, "g", pseudocount = "all_zeros", dispersion = "fixed"), "absent")
  expect_error(fit_wt(a[[1]], a[[2]], a[[3]], td$reference_id, -1, "g", pseudocount = "all_zeros", dispersion = "fixed"), "elapsed_time")
  expect_error(fit_wt(a[[1]], a[[2]], a[[3]], td$reference_id, td$elapsed_time, " ", pseudocount = "all_zeros", dispersion = "fixed"), "time_unit")
})

test_that("fit_wt errors when elapsed_time disagrees with samples$time", {
  expect_error(fit_wt(td$counts, td$samples, td$sequences, td$reference_id, 3.5, "generations",
                      pseudocount = "all_zeros", dispersion = "fixed"),
               "max\\(time\\) - min\\(time\\)")
  bad <- td$samples; bad$time <- as.character(bad$time)
  expect_error(fit_wt(td$counts, bad, td$sequences, td$reference_id, td$elapsed_time, "g",
                      pseudocount = "all_zeros", dispersion = "fixed"),
               "numeric and finite")
})

test_that("fit_wt requires replicate_id", {
  expect_error(fit_wt(td$counts, td$samples[setdiff(names(td$samples), "replicate_id")],
                      td$sequences, td$reference_id, td$elapsed_time, "g", pseudocount = "all_zeros", dispersion = "fixed"),
               "replicate_id")
})

test_that("fit_wt cores = 2 matches cores = 1", {
  skip_on_os("windows")
  a <- list(td$counts, td$samples, td$sequences, td$reference_id, td$elapsed_time, "generations",
       pseudocount = "all_zeros", dispersion = "fixed")
  one <- do.call(fit_wt, c(a, cores = 1))
  two <- do.call(fit_wt, c(a, cores = 2, chunk_size = 3L))
  expect_identical(two, one)
  expect_error(do.call(fit_wt, c(a, chunk_size = 0)), "chunk_size")
  expect_error(do.call(fit_wt, c(a, cores = 0)), "cores")
  expect_error(do.call(fit_wt, c(a, cores = 1.5)), "cores")
})

test_that("fit_wt keeps the reference row at its position in sequences", {
  ids <- td$sequences$sequence_id
  ref <- ids[3]
  res <- fit_wt(td$counts, td$samples, td$sequences, reference_id = ref,
                elapsed_time = td$elapsed_time, time_unit = "generations",
                pseudocount = "all_zeros", dispersion = "fixed")
  expect_equal(res$sequence_id, ids)
  expect_equal(res$status[3], "reference")
})
