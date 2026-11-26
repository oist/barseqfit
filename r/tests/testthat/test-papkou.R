# Synthetic Papkou2023 store built from the committed fixture: the three RDS
# files plus a MANIFEST.md5 that lists two of them.
make_store <- function(.env = parent.frame()) {
  store <- file.path(withr::local_tempdir(.local_envir = .env), "store")
  in_dir <- file.path(store, "zenodo", "in_data")
  dir.create(in_dir, recursive = TRUE)
  ids <- td$sequences$sequence_id
  filler <- "TTTTTTTTT"
  only <- "CCCCCCCCC"
  reps <- sub("^after_", "", td$samples$sample_id[grepl("^after_", td$samples$sample_id)])
  mat <- function(prefix, genotypes) {
    m <- matrix(0L, length(reps), length(genotypes),
                dimnames = list(paste0("Sample_", reps), genotypes))
    k <- startsWith(td$counts$sample_id, prefix)
    m[cbind(match(sub(paste0("^", prefix, "_"), "", td$counts$sample_id[k]), reps),
            match(td$counts$sequence_id[k], genotypes))] <- as.integer(td$counts$count[k])
    m
  }
  total <- function(prefix) td$samples$sample_total[match(paste0(prefix, "_", reps), td$samples$sample_id)]
  # Filler genotype: sample_total minus the fixture sums, so the unfiltered row
  # sums equal sample_total (the before-only genotype holds 5 reads per sample).
  before <- cbind(mat("before", ids), 0L, 0L)
  colnames(before) <- c(ids, filler, only)
  before[, filler] <- as.integer(total("before") - rowSums(before[, ids]) - 5L)
  before[, only] <- 5L
  after <- cbind(mat("after", ids), 0L)
  colnames(after) <- c(ids, filler)
  after[, filler] <- as.integer(total("after") - rowSums(after[, ids]))
  # 60-character names with the id at 52-60; before columns in another order.
  long <- function(m) {
    colnames(m) <- paste0(strrep("x", 51), colnames(m))
    m
  }
  before <- long(before)[, rev(seq_len(ncol(before)))]
  after <- long(after)
  saveRDS(list(before = before, "M9+TMP" = after), file.path(in_dir, "sequence_reads_counts.RDS"))
  ep <- td$expected_pair
  edge <- data.frame(from = c(ep$from, ids[2]), to = c(ep$to, only),
                     m = c(ep$m, 1), m_se = c(ep$se, 1), m_p.value = c(ep$p, 1),
                     stringsAsFactors = FALSE)
  saveRDS(edge, file.path(in_dir, "fitness_diff_wald_test_data.RDS"))
  ew <- td$expected_wt
  wt <- data.frame(SV = c(rev(ew$sequence_id), filler), m = c(rev(ew$m), NA),
                   m_se = c(rev(ew$se), NA), m_p.value = c(rev(ew$p), NA),
                   stringsAsFactors = FALSE)
  saveRDS(wt, file.path(in_dir, "fitness_data_wt.rds"))
  writeLines(c("0123456789abcdef0123456789abcdef  ./zenodo/in_data/sequence_reads_counts.RDS",
               "fedcba9876543210fedcba9876543210 *./zenodo/in_data/fitness_data_wt.rds",
               "00000000000000000000000000000000  ./README.md"),
             file.path(store, "MANIFEST.md5"))
  store
}

test_that("read_papkou_store builds the Export tables and the manifest", {
  store <- make_store()
  x <- read_papkou_store(paste0(store, "/"))
  expect_identical(x$manifest$store, normalizePath(store))
  expect_equal(x$manifest$edge_rows_dropped, 1)
  sf <- x$manifest$store_files
  in_dir <- file.path(normalizePath(store), "zenodo", "in_data")
  expect_identical(sf$path, file.path(in_dir, c("sequence_reads_counts.RDS", "fitness_data_wt.rds",
                                                "fitness_diff_wald_test_data.RDS")))
  expect_identical(sf$md5_source, c("MANIFEST.md5", "MANIFEST.md5", "computed"))
  expect_identical(sf$md5[1:2], c("0123456789abcdef0123456789abcdef", "fedcba9876543210fedcba9876543210"))
  expect_identical(sf$md5[3], unname(tools::md5sum(sf$path[3])))
  expect_match(sf$mtime, "^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}[+-]\\d{4}$")
  expect_identical(x$reference_id, td$reference_id)
  expect_identical(x$elapsed_time, td$elapsed_time)
  expect_identical(x$time_unit, td$time_unit)
  # Filler genotype and the before-only id.
  expect_identical(x$before_only_ids, data.frame(sequence_id = "CCCCCCCCC"))
  expect_identical(x$sequences$sequence_id, c(td$sequences$sequence_id, "TTTTTTTTT"))
  # expected_wt keeps RDS order, not the fixture order.
  expect_identical(x$expected_wt$sequence_id, c(rev(td$expected_wt$sequence_id), "TTTTTTTTT"))
  expect_identical(x$pairs, td$pairs)
  expect_equal(x$samples, td$samples)
  expect_equal(nrow(x$counts), 12L * length(x$sequences$sequence_id))
})

test_that(".papkou_subset reproduces tests/data from the synthetic store", {
  x <- read_papkou_store(paste0(make_store(), "/"))
  sub <- .papkou_subset(x, seq_len(nrow(x$pairs)))
  for (nm in c("counts", "samples", "sequences", "pairs", "expected_wt", "expected_pair")) {
    expect_equal(sub[[nm]], td[[nm]], info = nm)
  }
})

test_that(".papkou_subset(td, seq_len(10)) reproduces tests/data", {
  sub <- .papkou_subset(td, seq_len(10))
  for (nm in c("counts", "samples", "sequences", "pairs", "expected_wt", "expected_pair")) {
    expect_equal(sub[[nm]], td[[nm]], info = nm)
  }
})

test_that(".papkou_subset keeps the chosen pairs, the reference and the time-block order", {
  sub <- .papkou_subset(td, c(3L, 1L))
  ids <- c(td$reference_id, setdiff(sort(unique(unlist(td$pairs[c(3L, 1L), ]))), td$reference_id))
  expect_equal(sub$pairs, td$pairs[c(3L, 1L), ], ignore_attr = TRUE)
  expect_identical(unique(sub$counts$sequence_id), ids)
  expect_identical(sub$expected_wt$sequence_id, ids)
  expect_identical(sub$sequences$sequence_id, ids)
  expect_true(all(startsWith(sub$counts$sample_id[seq_len(nrow(sub$counts) / 2)], "before_")))
  expect_equal(sub$samples, td$samples)
})

test_that("the Export round trips through write_papkou_export and read_papkou_export", {
  x <- read_papkou_store(paste0(make_store(), "/"))
  dir <- file.path(withr::local_tempdir(), "export")
  cfg <- list(fit = list(chunk_size = 2000L), gate = list(abs_tol = 1e-6))
  .write_papkou_export(x, dir, cfg, "2026-10-14T17:32:00+0900")
  expect_setequal(list.files(dir), c(paste0(.papkou_export_tables, ".parquet"), "export_manifest.json"))
  got <- read_papkou_export(dir, tables = .papkou_export_tables)
  for (nm in .papkou_export_tables) {
    expect_s3_class(got[[nm]], "data.frame")
    expect_false(inherits(got[[nm]], "tbl_df"))
    expect_equal(got[[nm]], x[[nm]], info = nm)
  }
  expect_identical(got$reference_id, td$reference_id)
  expect_identical(got$elapsed_time, td$elapsed_time)
  expect_identical(got$time_unit, td$time_unit)
  m <- got$manifest
  expect_identical(m$timestamp, "2026-10-14T17:32:00+0900")
  expect_identical(m$store, x$manifest$store)
  expect_identical(m$store_files$md5, x$manifest$store_files$md5)
  expect_identical(m$config$fit$chunk_size, 2000L)
  expect_identical(m$row_counts$counts, nrow(x$counts))
  expect_identical(m$row_counts$edge_rows_dropped, 1L)
  # The default tables.
  expect_named(read_papkou_export(dir), c("counts", "samples", "sequences", "pairs",
                                          "reference_id", "elapsed_time", "time_unit", "manifest"))
})

test_that("read_papkou_export refuses missing files and does not compare row counts", {
  x <- .papkou_subset(td, seq_len(10))
  x$manifest <- list(store = "/s", store_files = data.frame(path = "/s/a", mtime = "t", md5 = "m",
                                                           md5_source = "computed"),
                     edge_rows_dropped = 0L)
  x$before_only_ids <- data.frame(sequence_id = character())
  dir <- file.path(withr::local_tempdir(), "export")
  .write_papkou_export(x, dir, list(a = 1L), "ts")
  # A trimmed Parquet file with the full manifest is accepted.
  nanoparquet::write_parquet(x$pairs[1:2, ], file.path(dir, "pairs.parquet"))
  expect_equal(nrow(read_papkou_export(dir)$pairs), 2L)
  file.remove(file.path(dir, "counts.parquet"))
  expect_error(read_papkou_export(dir), "just export-papkou")
  expect_no_error(read_papkou_export(dir, tables = "pairs"))
  file.remove(file.path(dir, "export_manifest.json"))
  expect_error(read_papkou_export(dir, tables = "pairs"), "just export-papkou")
  expect_error(read_papkou_export(dir, tables = "nope"), "unknown")
  expect_error(read_papkou_export(file.path(dir, "absent")), "just export-papkou")
})

test_that("read_papkou_store refuses an unset or missing store", {
  expect_error(read_papkou_store(NA_character_), "not set")
  expect_error(read_papkou_store(file.path(tempdir(), "absent-store")), "does not exist")
  empty <- withr::local_tempdir()
  expect_error(read_papkou_store(empty), "not found")
})

test_that("export_manifest.json bytes match the format of scripts/export_papkou.R", {
  x <- list(
    counts = data.frame(sequence_id = c("A", "B"), sample_id = "s", count = 1:2),
    samples = data.frame(sample_id = "s", condition = "M9+TMP", replicate_id = "F01", time = 0,
                         sample_total = 3),
    sequences = data.frame(sequence_id = "A", sequence = "A"),
    expected_wt = data.frame(sequence_id = c("A", "B", "C"), m = 0, se = 0, p = 0),
    before_only_ids = data.frame(sequence_id = character()),
    pairs = data.frame(from = "A", to = "B"),
    expected_pair = data.frame(from = "A", to = "B", m = 0, se = 0, p = 0),
    manifest = list(
      store = "/data/st\"ore",
      store_files = data.frame(
        path = c("/data/store/a.RDS", "/data/store/b.rds"),
        mtime = c("2023-02-14T23:09:08+0900", "2023-08-14T00:23:23+0900"),
        md5 = c("821a5866d34c540eb528de805624a9e6", "e636390a6681aeda8ed90dc54fdc17d2"),
        md5_source = c("MANIFEST.md5", "computed")),
      edge_rows_dropped = 0L))
  dir <- file.path(withr::local_tempdir(), "export")
  cfg <- list(fit = list(chunk_size = 2000L), gate = list(abs_tol = 1e-6, tie_fraction = 0.0005))
  .write_papkou_export(x, dir, cfg, "2026-10-14T17:32:00+0900")
  golden <- c(
    "{",
    "  \"store\": \"/data/st\\\"ore\",",
    "  \"store_files\": [",
    "    {\"path\": \"/data/store/a.RDS\", \"mtime\": \"2023-02-14T23:09:08+0900\", \"md5\": \"821a5866d34c540eb528de805624a9e6\", \"md5_source\": \"MANIFEST.md5\"},",
    "    {\"path\": \"/data/store/b.rds\", \"mtime\": \"2023-08-14T00:23:23+0900\", \"md5\": \"e636390a6681aeda8ed90dc54fdc17d2\", \"md5_source\": \"computed\"}",
    "  ],",
    sprintf("  \"r_version\": \"%s\",", R.version.string),
    "  \"timestamp\": \"2026-10-14T17:32:00+0900\",",
    "  \"format\": \"parquet\",",
    "  \"config\": {\"fit\": {\"chunk_size\": 2000}, \"gate\": {\"abs_tol\": 1e-06, \"tie_fraction\": 0.0005}},",
    "  \"row_counts\": {",
    "    \"counts\": 2,",
    "    \"samples\": 1,",
    "    \"sequences\": 1,",
    "    \"expected_wt\": 3,",
    "    \"before_only_ids\": 0,",
    "    \"pairs\": 1,",
    "    \"expected_pair\": 1,",
    "    \"edge_rows_dropped\": 0",
    "  }",
    "}")
  got <- readBin(file.path(dir, "export_manifest.json"), "raw", 1e6)
  expect_identical(rawToChar(got), paste0(paste(golden, collapse = "\n"), "\n"))
})
