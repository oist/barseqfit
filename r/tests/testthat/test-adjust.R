# Benjamini-Hochberg adjusted p; mirrors tests/test_adjust.py.

# Generated with R 4.1.3:
#   set.seed(42); p <- round(runif(20)^2, 6); dput(p); dput(p.adjust(p, "BH"))
r_p <- c(0.83687, 0.87811, 0.081876, 0.689643, 0.411837, 0.269461, 0.542562,
         0.018135, 0.431639, 0.497116, 0.209528, 0.517122, 0.873612, 0.065244,
         0.213715, 0.883627, 0.956927, 0.013803, 0.225622, 0.313973)
r_p_adj <- c(0.930133684210526, 0.930133684210526, 0.40938, 0.919524, 0.775088571428571,
             0.6736525, 0.775088571428571, 0.18135, 0.775088571428571, 0.775088571428571,
             0.644634285714286, 0.775088571428571, 0.930133684210526, 0.40938,
             0.644634285714286, 0.930133684210526, 0.956927, 0.18135, 0.644634285714286,
             0.697717777777778)

test_that(".bh matches hand values", {
  expect_equal(.bh(c(0.01, 0.04, 0.03, 0.2)), c(0.04, 0.16 / 3, 0.16 / 3, 0.2),
               tolerance = 1e-15)
})

test_that(".bh passes NA through and counts only non-NA p", {
  got <- .bh(c(0.01, NA, 0.04, 0.03, 0.2))
  expect_true(is.na(got[2]))
  expect_equal(got[-2], c(0.04, 0.16 / 3, 0.16 / 3, 0.2), tolerance = 1e-15)
  expect_length(.bh(numeric(0)), 0)
  expect_true(all(is.na(.bh(c(NA_real_, NA_real_)))))
})

test_that(".bh handles ties and is monotone and at most 1", {
  expect_equal(.bh(c(0.02, 0.02, 0.02, 0.5)), c(rep(0.02 * 4 / 3, 3), 0.5), tolerance = 1e-15)
  set.seed(0)
  p <- runif(500)^0.5
  got <- .bh(p)
  expect_lte(max(got), 1)
  expect_equal(max(got), max(p))
  expect_true(all(got >= p))
  expect_false(is.unsorted(got[order(p)]))
})

test_that(".bh equals the pinned R p.adjust vector", {
  expect_equal(.bh(r_p), r_p_adj, tolerance = 1e-12)
  expect_equal(stats::p.adjust(r_p, "BH"), r_p_adj, tolerance = 1e-12)
})

test_that("p_adj sits right after p and the fixture matches the pinned R output", {
  wt <- fit_wt(td$counts, td$samples, td$sequences, td$reference_id, td$elapsed_time,
               "generations", pseudocount = "all_zeros", dispersion = "fixed")
  expect_equal(names(wt)[which(names(wt) == "p") + 1L], "p_adj")
  expect_true(is.na(wt$p_adj[wt$sequence_id == td$reference_id]))
  nonref <- wt$sequence_id != td$reference_id
  expect_equal(wt$p_adj[nonref], stats::p.adjust(wt$p[nonref], "BH"), tolerance = 1e-15)
  pinned <- read.csv(file.path(test_data_dir(), "expected_wt_padj_r.csv"),
                     na.strings = c("", "NA"), stringsAsFactors = FALSE)
  expect_identical(pinned$sequence_id, wt$sequence_id)
  expect_equal(wt$p_adj, pinned$p_adj, tolerance = 1e-9)
  pr <- fit_pair(td$counts, td$samples, td$sequences, td$pairs, td$elapsed_time,
                 "generations", pseudocount = "all_zeros", dispersion = "fixed")
  pinned <- read.csv(file.path(test_data_dir(), "expected_pair_padj_r.csv"),
                     na.strings = c("", "NA"), stringsAsFactors = FALSE)
  expect_equal(pr$p_adj, pinned$p_adj, tolerance = 1e-9)
})

test_that("p_adj is identical for fit_wt and fit_pair on the same (reference, x) pairs", {
  el <- td$elapsed_time
  wt <- fit_wt(td$counts, td$samples, td$sequences, td$reference_id, el, "generations",
               pseudocount = "all_zeros", dispersion = "fixed")
  wt <- wt[wt$sequence_id != td$reference_id, ]
  pr <- fit_pair(td$counts, td$samples, td$sequences,
                 data.frame(from = td$reference_id, to = wt$sequence_id, stringsAsFactors = FALSE),
                 el, "generations", pseudocount = "all_zeros", dispersion = "fixed")
  expect_identical(pr$p_adj, wt$p_adj)
})

test_that("p_adj is NA where p is NA (failed rows) and sits right after p", {
  res <- data.frame(m = c(1, NA, 2), se = c(1, NA, 1), p = c(0.01, NA, 0.04),
                    status = c("converged", "failed", "converged"))
  out <- .add_p_adj(res)
  expect_equal(names(out), c("m", "se", "p", "p_adj", "status"))
  expect_equal(out$p_adj, c(0.02, NA, 0.04), tolerance = 1e-15)
})
