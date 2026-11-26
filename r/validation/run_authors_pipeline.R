# Run the Papkou et al. 2023 authors' fitness code on the tests/data fixture
# and compare against the published numbers (expected_*.csv).
# Usage (repo root): Rscript r/validation/run_authors_pipeline.R
# Needs data.table (`just validate-papkou2023` installs it into the renv library).
# The authors' scripts are read from the Papkou2023 store ($PAPKOU2023RUGGED_DATA,
# directory ampseq-nextflow/) and each file is checked against the md5 list in
# r/validation/authors_scripts.md5 before it is sourced; a missing file or a
# mismatch stops the run. See docs/validation.md.
# tests/data/samples.csv carries the unfiltered mapped-read totals and the generation
# time (.papkou_generation_time). The store-derived totals are also run as a
# cross-check that the fixture totals equal the store's.
# The verdict is written to $VERDICT_OUT (default: <tempdir>/verdict.csv).
# The authors' files are sourced untouched. Both libs define get_data_fit_model_get_results
# and get_model_results, so each is sourced into its own environment.
# The package files in r/R are sourced for the Papkou2023 constants and read_papkou_store().

suppressPackageStartupMessages(library(data.table))
tol <- 1e-6

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
root <- if (length(file_arg)) normalizePath(file.path(dirname(file_arg[1]), "..", "..")) else getwd()
for (f in list.files(file.path(root, "r", "R"), pattern = "\\.R$", full.names = TRUE)) source(f)
wt_id <- .papkou_reference_id
fixture_time <- .papkou_generation_time
dat <- file.path(root, "tests", "data")

# ---- authors' scripts: locate in the store and verify md5 ----
store <- .resolve_machine_env(root)$store
if (is.na(store) || !nzchar(store) || !dir.exists(store)) {
  stop("PAPKOU2023RUGGED_DATA is unset or not a directory (set it in the shell or in .env); ",
       "the authors' scripts are read from its ampseq-nextflow/ directory.", call. = FALSE)
}
authors_dir <- file.path(store, "ampseq-nextflow")
md5_lines <- readLines(file.path(root, "r", "validation", "authors_scripts.md5"))
md5_lines <- md5_lines[nzchar(md5_lines)]
want <- setNames(sub("^[0-9a-f]+  ", "", md5_lines), sub("^([0-9a-f]+)  .*$", "\\1", md5_lines))
for (i in seq_along(want)) {
  path <- file.path(authors_dir, want[[i]])
  if (!file.exists(path)) stop("Authors' script missing from the store: ", path, call. = FALSE)
  got <- unname(tools::md5sum(path))
  if (!identical(got, names(want)[i])) {
    stop("md5 mismatch for ", path, ": expected ", names(want)[i], ", found ", got, call. = FALSE)
  }
}
cat("Verified md5 of", length(want), "authors' scripts in", authors_dir, "\n")
lib <- file.path(authors_dir, "lib")

rd <- function(f) read.csv(file.path(dat, f), stringsAsFactors = FALSE)
counts <- rd("counts.csv"); samples <- rd("samples.csv")
pairs <- rd("pairs.csv"); ewt <- rd("expected_wt.csv"); epair <- rd("expected_pair.csv")

env_wt <- new.env(); sys.source(file.path(lib, "fit_model_wt.R"), envir = env_wt)
env_pair <- new.env(); sys.source(file.path(lib, "fit_model_from_to.v2.R"), envir = env_pair)

# ---- totals (per sample_id) ----
# Fixture sample_total = row sum of the full UNFILTERED count matrices (before any
# genotype subsetting), i.e. the validated total.
totals <- list(fixture_sample_total = setNames(samples$sample_total, samples$sample_id))
# Cross-check: recompute the same unfiltered row sums from the store RDS.
st <- read_papkou_store(store)$samples
totals$unfiltered_rowsum_STORE_crosscheck <- setNames(st$sample_total, st$sample_id)

# ---- data frame as the authors' code expects ----
make_df <- function(total, rule, after_time) {
  d <- merge(counts, samples[, c("sample_id", "time")], by = "sample_id")
  d$time[d$time != 0] <- after_time  # override the fixture time with the candidate
  d$Total <- unname(total[d$sample_id])
  d$Nreads <- as.numeric(d$count)
  is_after <- d$time != 0
  zero <- d$Nreads == 0
  # rule = <which timepoints get pseudocount>_<Diff computed before|after pseudocount>
  switch(rule,
    both_diff_before = { d$Diff <- d$Total - d$Nreads; d$Nreads[zero] <- 1 },        # 01.prepare_data_for_fitness.r
    both_diff_after = { d$Nreads[zero] <- 1; d$Diff <- d$Total - d$Nreads },
    after_diff_after = { d$Nreads[zero & is_after] <- 1; d$Diff <- d$Total - d$Nreads }, # lib/fit_model_from_to.R
    after_diff_before = { d$Diff <- d$Total - d$Nreads; d$Nreads[zero & is_after] <- 1 },
    none = { d$Diff <- d$Total - d$Nreads },
    stop("bad rule"))
  dt <- as.data.table(d[, c("sequence_id", "Nreads", "Diff", "time")])
  setnames(dt, "sequence_id", "SV")
  setkey(dt, SV)
  dt
}

rules <- c("both_diff_before", "both_diff_after", "after_diff_after", "after_diff_before", "none")
times <- list(`3.51662` = 3.51662)
times[[format(fixture_time, digits = 17)]] <- fixture_time

maxres <- function(x, y) if (all(is.na(x) & is.na(y))) NA_real_ else max(abs(x - y), na.rm = TRUE)

rows <- list()
for (tn in names(times)) for (tot in names(totals)) for (rule in rules) {
  dt <- make_df(totals[[tot]], rule, times[[tn]])
  # WT: every subset genotype vs GCCGATCTC (WT itself has m = 0, NA se/p: excluded)
  ids <- setdiff(ewt$sequence_id, wt_id)
  w <- do.call(rbind, lapply(ids, function(s) suppressWarnings(env_wt$get_data_fit_model_get_results(dt, s))[, c("SV", "m", "m_se", "m_p.value")]))
  e <- ewt[match(w$SV, ewt$sequence_id), ]
  p <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) suppressWarnings(env_pair$get_data_fit_model_get_results(dt, pairs$from[i], pairs$to[i]))[, c("from", "to", "m", "m_se", "m_p.value")]))
  ep <- epair[match(paste(p$from, p$to), paste(epair$from, epair$to)), ]
  r <- data.frame(time = tn, total = tot, pseudocount = rule,
    wt_n = nrow(w), wt_m = maxres(w$m, e$m), wt_se = maxres(w$m_se, e$se), wt_p = maxres(w$m_p.value, e$p),
    pair_n = nrow(p), pair_m = maxres(p$m, ep$m), pair_se = maxres(p$m_se, ep$se), pair_p = maxres(p$m_p.value, ep$p))
  r$pass <- all(c(r$wt_m, r$wt_se, r$wt_p, r$pair_m, r$pair_se, r$pair_p) < tol)
  rows[[length(rows) + 1]] <- r
}
verdict <- do.call(rbind, rows)

out <- Sys.getenv("VERDICT_OUT", unset = file.path(tempdir(), "verdict.csv"))
write.csv(verdict, out, row.names = FALSE)
cat("Verdict written to", out, "\n")
options(width = 250, digits = 4)
print(verdict, row.names = FALSE)
cat("\nPassing candidates (all max abs residuals <", tol, "):\n")
print(verdict[verdict$pass, c("time", "total", "pseudocount")], row.names = FALSE)
