# Papkou et al. 2023 constants and the Export format. The raw store is
# read only by read_papkou_store(); everything downstream reads the Export.

# Reference genotype and the Generation time recovered for the M9+TMP after samples
# (see docs/validation.md, Generation time).
.papkou_reference_id <- "GCCGATCTC"
.papkou_generation_time <- 3.5166197103891137
.papkou_condition <- "M9+TMP"
.papkou_time_unit <- "generations"

# Store files, relative to the store directory, in manifest order.
.papkou_store_files <- file.path("zenodo", "in_data",
  c(reads = "sequence_reads_counts.RDS", wt = "fitness_data_wt.rds",
    edge = "fitness_diff_wald_test_data.RDS"))

# Export tables in write order: Parquet file names, row_counts keys and the
# reader's `tables` argument.
.papkou_export_tables <- c("counts", "samples", "sequences", "expected_wt",
                           "before_only_ids", "pairs", "expected_pair")

# Build the Export tables from the raw store objects. `raw` is list(reads, wt, edge):
# the reads list (matrices samples x genotypes, named "before" and the condition),
# the WT table and the edge table. Returns list(tables, edge_rows_dropped).
.papkou_tables <- function(raw) {
  before <- raw$reads[["before"]]
  after <- raw$reads[[.papkou_condition]]
  # Per-sample totals: row sum of the FULL UNFILTERED matrices (before genotype subsetting).
  total_before <- rowSums(before)
  total_after <- rowSums(after)
  colnames(before) <- substr(colnames(before), 52, 60)
  colnames(after) <- substr(colnames(after), 52, 60)
  stopifnot(!anyDuplicated(colnames(before)), !anyDuplicated(colnames(after)))
  stopifnot(all(colnames(after) %in% colnames(before)))
  before_only <- setdiff(colnames(before), colnames(after))
  # Mirror the authors' filter: keep genotypes present in the after matrix.
  before <- before[, colnames(after), drop = FALSE]
  stopifnot(identical(rownames(before), rownames(after)))
  reps <- sub("^Sample_", "", rownames(after))
  ids <- colnames(after)

  # Long counts, zeros kept. Matrices are samples x genotypes, so as.integer(mat)
  # runs sample-fastest within each genotype.
  mk <- function(mat, prefix) {
    data.frame(sequence_id = rep(ids, each = nrow(mat)),
               sample_id = rep(paste0(prefix, "_", reps), times = length(ids)),
               count = as.integer(mat), stringsAsFactors = FALSE)
  }
  counts <- rbind(mk(before, "before"), mk(after, "after"))
  samples <- data.frame(
    sample_id = c(paste0("before_", reps), paste0("after_", reps)),
    condition = .papkou_condition,
    replicate_id = rep(reps, 2),
    time = rep(c(0, .papkou_generation_time), each = length(reps)),
    sample_total = c(total_before, total_after),
    stringsAsFactors = FALSE)
  wt <- as.data.frame(raw$wt)
  edge <- as.data.frame(raw$edge)
  keep <- edge$from %in% ids & edge$to %in% ids
  edge <- edge[keep, ]
  from <- as.character(edge$from)
  to <- as.character(edge$to)
  list(
    tables = list(
      counts = counts,
      samples = samples,
      sequences = data.frame(sequence_id = ids, sequence = ids, stringsAsFactors = FALSE),
      expected_wt = data.frame(sequence_id = as.character(wt$SV), m = as.double(wt$m),
                               se = as.double(wt$m_se), p = as.double(wt$m_p.value),
                               stringsAsFactors = FALSE),
      before_only_ids = data.frame(sequence_id = sort(before_only), stringsAsFactors = FALSE),
      pairs = data.frame(from = from, to = to, stringsAsFactors = FALSE),
      expected_pair = data.frame(from = from, to = to, m = as.double(edge$m),
                                 se = as.double(edge$m_se), p = as.double(edge$m_p.value),
                                 stringsAsFactors = FALSE)),
    edge_rows_dropped = sum(!keep))
}

#' Read the Papkou2023 store into the Export tables.
#'
#' The only function that reads the raw RDS files. Returns a named list of the
#' seven Export tables (`counts`, `samples`, `sequences`, `expected_wt`,
#' `before_only_ids`, `pairs`, `expected_pair`) plus `reference_id`,
#' `elapsed_time`, `time_unit` and `manifest` (the store path, the store files
#' with mtime and md5, and `edge_rows_dropped`). The md5 comes from the store's
#' `MANIFEST.md5` when it lists the file, otherwise it is computed.
#' @param store path to the store directory; default `papkou_store()`.
#' @export
read_papkou_store <- function(store = papkou_store()) {
  store <- .papkou_store_path(store, "store")
  paths <- file.path(store, .papkou_store_files)
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("store file(s) not found: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  tabs <- .papkou_tables(list(reads = readRDS(paths[[1]]), wt = readRDS(paths[[2]]),
                              edge = readRDS(paths[[3]])))

  man_path <- file.path(store, "MANIFEST.md5")
  man <- if (file.exists(man_path)) {
    ln <- readLines(man_path, warn = FALSE)
    ln <- ln[nzchar(ln)]
    stats::setNames(sub("\\s.*$", "", ln), sub("^\\S+\\s+\\*?", "", ln))
  } else character()
  key <- paste0("./", .papkou_store_files)
  listed <- key %in% names(man)
  store_files <- data.frame(
    path = paths,
    mtime = format(file.mtime(paths), "%Y-%m-%dT%H:%M:%S%z"),
    md5 = ifelse(listed, unname(man[key]), unname(tools::md5sum(paths))),
    md5_source = ifelse(listed, "MANIFEST.md5", "computed"),
    stringsAsFactors = FALSE)
  rownames(store_files) <- NULL
  c(tabs$tables,
    list(reference_id = .papkou_reference_id, elapsed_time = .papkou_generation_time,
         time_unit = .papkou_time_unit,
         manifest = list(store = store, store_files = store_files,
                         edge_rows_dropped = tabs$edge_rows_dropped)))
}

# Restrict a reader result to the genotypes of pair rows `pair_rows`: keeps those
# pairs, then the reference genotype and every genotype they name. Samples are
# unchanged. Counts keep their order (time block, then genotype in `ids` order).
.papkou_subset <- function(x, pair_rows) {
  pairs <- x$pairs[pair_rows, , drop = FALSE]
  expected_pair <- x$expected_pair[pair_rows, , drop = FALSE]
  ids <- c(x$reference_id,
           setdiff(sort(unique(c(pairs$from, pairs$to))), x$reference_id))
  stopifnot(all(ids %in% x$sequences$sequence_id))
  counts <- x$counts[x$counts$sequence_id %in% ids, , drop = FALSE]
  block <- x$samples$time[match(counts$sample_id, x$samples$sample_id)]
  counts <- counts[order(block, match(counts$sequence_id, ids)), , drop = FALSE]
  sequences <- x$sequences[x$sequences$sequence_id %in% ids, , drop = FALSE]
  sequences <- sequences[order(match(sequences$sequence_id, ids)), , drop = FALSE]
  expected_wt <- x$expected_wt[match(ids, x$expected_wt$sequence_id), , drop = FALSE]
  stopifnot(identical(expected_wt$sequence_id, ids))
  out <- x
  out$counts <- counts
  out$sequences <- sequences
  out$expected_wt <- expected_wt
  out$pairs <- pairs
  out$expected_pair <- expected_pair
  for (nm in c("counts", "sequences", "expected_wt", "pairs", "expected_pair")) {
    rownames(out[[nm]]) <- NULL
  }
  out
}

# Write the Export tables as Parquet (full double precision, NA as null) and a
# hand-built export_manifest.json. `x` is a read_papkou_store() result, `config`
# the Run config, `timestamp` a "%Y-%m-%dT%H:%M:%S%z" string. row_counts are the
# nrow() of the tables written, plus edge_rows_dropped.
.write_papkou_export <- function(x, dir, config, timestamp) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  for (nm in .papkou_export_tables) {
    nanoparquet::write_parquet(x[[nm]], file.path(dir, paste0(nm, ".parquet")))
  }
  js <- function(s) paste0("\"", gsub("([\"\\\\])", "\\\\\\1", s), "\"")
  sf <- x$manifest$store_files
  entries <- paste0("    {\"path\": ", js(sf$path), ", \"mtime\": ", js(sf$mtime),
                    ", \"md5\": ", js(sf$md5), ", \"md5_source\": ", js(sf$md5_source), "}")
  n_rows <- vapply(.papkou_export_tables, function(nm) nrow(x[[nm]]), integer(1))
  counts <- c(paste0("    \"", names(n_rows), "\": ", n_rows, ","),
              paste0("    \"edge_rows_dropped\": ", x$manifest$edge_rows_dropped))
  manifest <- c(
    "{",
    paste0("  \"store\": ", js(x$manifest$store), ","),
    "  \"store_files\": [",
    paste(entries, collapse = ",\n"),
    "  ],",
    paste0("  \"r_version\": ", js(R.version.string), ","),
    paste0("  \"timestamp\": ", js(timestamp), ","),
    "  \"format\": \"parquet\",",
    paste0("  \"config\": ", .config_json(config), ","),
    "  \"row_counts\": {",
    counts,
    "  }",
    "}")
  writeLines(manifest, file.path(dir, "export_manifest.json"))
  invisible(dir)
}

#' Read an Export directory written by `just export-papkou`.
#'
#' Returns the requested tables as plain data.frames plus `reference_id`,
#' `elapsed_time`, `time_unit` and `manifest` (`export_manifest.json`, parsed).
#' Refuses a directory missing a requested file or the manifest. Row counts are
#' not compared with the manifest, because smoke runs trim Parquet files.
#' @param dir the Export directory (the scratch directory).
#' @param tables Export tables to read.
#' @export
read_papkou_export <- function(dir, tables = c("counts", "samples", "sequences", "pairs")) {
  bad <- setdiff(tables, .papkou_export_tables)
  if (length(bad)) {
    stop("unknown Export table(s): ", paste(bad, collapse = ", "), call. = FALSE)
  }
  files <- c(file.path(dir, paste0(tables, ".parquet")), file.path(dir, "export_manifest.json"))
  missing <- files[!file.exists(files)]
  if (length(missing)) {
    stop("Export file(s) missing: ", paste(missing, collapse = ", "),
         ". Run `just export-papkou`.", call. = FALSE)
  }
  out <- stats::setNames(lapply(tables, function(nm) {
    as.data.frame(nanoparquet::read_parquet(file.path(dir, paste0(nm, ".parquet"))))
  }), tables)
  c(out, list(reference_id = .papkou_reference_id, elapsed_time = .papkou_generation_time,
              time_unit = .papkou_time_unit,
              manifest = jsonlite::fromJSON(file.path(dir, "export_manifest.json"))))
}
