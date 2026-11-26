test_that(".load_run_config reads validation.toml", {
  cfg <- .load_run_config()
  expect_identical(cfg$fit$chunk_size, 2000L)
  expect_equal(cfg$gate$abs_tol, 1e-6)
  expect_equal(cfg$gate$tie_fraction, 0.0005)
  expect_error(.load_run_config(file.path(tempdir(), "missing.toml")), "not found")
})

# Golden timing files, in the byte format scripts/run_full_wt.R and run_full_pair.R wrote.
timing_cfg <- list(fit = list(chunk_size = 2000L), gate = list(abs_tol = 1e-6))
timing_manifest <- list(
  timestamp = "2026-10-14T17:32:00+0900",
  store_files = data.frame(md5 = c("821a5866d34c540eb528de805624a9e6",
                                   "e636390a6681aeda8ed90dc54fdc17d2",
                                   "20bfb173c6f0ea8815b08d140b617f7a")),
  config = list(fit = list(chunk_size = 1L)))
timing_tail <- paste0("\"export_timestamp\": \"2026-10-14T17:32:00+0900\", \"store_md5\": ",
                      "[\"821a5866d34c540eb528de805624a9e6\", \"e636390a6681aeda8ed90dc54fdc17d2\", ",
                      "\"20bfb173c6f0ea8815b08d140b617f7a\"], ",
                      "\"config\": {\"fit\": {\"chunk_size\": 2000}, \"gate\": {\"abs_tol\": 1e-06}}}\n")

timing_bytes <- function(...) {
  path <- tempfile(fileext = ".json")
  on.exit(unlink(path))
  .write_timing(path, ...)
  rawToChar(readBin(path, "raw", 1e6))
}

test_that(".write_timing reproduces the WT timing bytes", {
  got <- timing_bytes(285.9254, 4L, list(n_genotypes = 261333L), timing_manifest, timing_cfg)
  expect_identical(got, paste0(
    sprintf("{\"wall_seconds\": 285.925, \"cores\": 4, \"r_version\": \"%s\", ", R.version.string),
    "\"n_genotypes\": 261333, ", timing_tail))
})

test_that(".write_timing reproduces the pair timing bytes, with null or a pair limit", {
  prefix <- sprintf("{\"wall_seconds\": 7050.474, \"cores\": 4, \"r_version\": \"%s\", ", R.version.string)
  full <- timing_bytes(7050.4742, 4, list(n_genotypes = 261333L, n_pairs = 3518696L, pair_limit = NA),
                       timing_manifest, timing_cfg)
  expect_identical(full, paste0(prefix, "\"n_genotypes\": 261333, \"n_pairs\": 3518696, ",
                                "\"pair_limit\": null, ", timing_tail))
  smoke <- timing_bytes(7050.4742, 4, list(n_genotypes = 11L, n_pairs = 10L, pair_limit = 10L),
                        timing_manifest, timing_cfg)
  expect_identical(smoke, paste0(prefix, "\"n_genotypes\": 11, \"n_pairs\": 10, ",
                                 "\"pair_limit\": 10, ", timing_tail))
})

test_that(".write_timing writes the passed config, not the manifest's, and one final newline", {
  got <- timing_bytes(1, 1L, list(n_genotypes = 1L), timing_manifest, timing_cfg)
  expect_match(got, "\"chunk_size\": 2000", fixed = TRUE)
  expect_false(grepl("\"chunk_size\": 1}", got, fixed = TRUE))
  expect_match(got, "\\}\n$")
  expect_false(grepl("\n\n", got, fixed = TRUE))
})
