# Data-store and fixture access helpers. Base R only.

# Repo root: the nearest ancestor directory holding both tests/data and
# r/DESCRIPTION. Walks up from the working directory, then from the installed
# package location (useful for a package loaded from a source checkout).
.is_repo_root <- function(d) {
  dir.exists(file.path(d, "tests", "data")) && file.exists(file.path(d, "r", "DESCRIPTION"))
}

.repo_root <- function() {
  pkg <- tryCatch(find.package("barseqfit", quiet = TRUE), error = function(e) character())
  starts <- c(getwd(), if (length(pkg)) pkg)
  for (start in starts) {
    d <- normalizePath(start, mustWork = FALSE)
    repeat {
      if (.is_repo_root(d)) return(d)
      parent <- dirname(d)
      if (identical(parent, d)) break
      d <- parent
    }
  }
  stop("Cannot locate the repository root (a directory with tests/data and r/DESCRIPTION) ",
       "above the working directory.", call. = FALSE)
}

# Validate a resolved store path (NA or "" when unset) and return it normalised.
# `what` names the source in messages: the variable, or "store" for an explicit path.
.papkou_store_path <- function(path, what = "PAPKOU2023RUGGED_DATA") {
  if (is.na(path) || !nzchar(path)) {
    stop(what, " is not set.",
         if (what == "PAPKOU2023RUGGED_DATA")
           " Export it in your shell or put it in a repo-root .env file (see .env.example).",
         call. = FALSE)
  }
  if (!dir.exists(path)) {
    stop(what, " points to a directory that does not exist: ", path, call. = FALSE)
  }
  normalizePath(path)
}

#' Path to the Papkou2023 data store.
#'
#' Resolves `PAPKOU2023RUGGED_DATA` from the shell, falling back to a repo-root
#' `.env` file (shell wins, `.env` fills gaps), and checks that the directory
#' exists. A relative path resolves against the repo root.
#' @export
papkou_store <- function() {
  root <- tryCatch(.repo_root(), error = function(e) NULL)
  .papkou_store_path(.resolve_machine_env(root)$store)
}

#' Repo-root `tests/data`, the single source of truth for fixtures.
#' @export
test_data_dir <- function() {
  file.path(.repo_root(), "tests", "data")
}

#' Read the six test-data CSVs from `dir` (default: repo-root tests/data).
#'
#' Empty fields become NA; ids are character, numbers numeric. Also carries
#' `reference_id`, `elapsed_time` and `time_unit`, as the Export readers do.
#' @export
read_test_data <- function(dir = test_data_dir()) {
  chr <- "character"; num <- "numeric"
  rd <- function(file, classes) {
    path <- file.path(dir, file)
    got <- utils::read.csv(path, colClasses = classes, na.strings = c("", "NA"),
                           stringsAsFactors = FALSE)
    if (!identical(names(got), names(classes))) {
      stop(file, " has columns ", paste(names(got), collapse = ", "), "; expected ",
           paste(names(classes), collapse = ", "), call. = FALSE)
    }
    got
  }
  list(
    counts = rd("counts.csv", c(sequence_id = chr, sample_id = chr, count = num)),
    samples = rd("samples.csv", c(sample_id = chr, condition = chr, replicate_id = chr,
                                  time = num, sample_total = num)),
    sequences = rd("sequences.csv", c(sequence_id = chr, sequence = chr)),
    pairs = rd("pairs.csv", c(from = chr, to = chr)),
    expected_wt = rd("expected_wt.csv", c(sequence_id = chr, m = num, se = num, p = num)),
    expected_pair = rd("expected_pair.csv", c(from = chr, to = chr, m = num, se = num, p = num)),
    reference_id = .papkou_reference_id,
    elapsed_time = .papkou_generation_time,
    time_unit = .papkou_time_unit
  )
}
