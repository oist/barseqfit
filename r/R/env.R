# Machine environment: store path, scratch directory and core count,
# from the shell and a repo-root `.env`. Base R only. Pure: nothing here calls
# Sys.setenv, so only the scripts and papkou_store() ever look at the environment.

# Project variables that may live in the shell or a repo-root `.env`.
.env_vars <- c("PAPKOU2023RUGGED_DATA", "BARSEQFIT_SCRATCH", "BARSEQFIT_CORES")

# Parse a `.env` file into a named character vector of every key in it (no
# dependency). Grammar, shared with the Python reader:
#   - blank lines and lines whose first non-space character is `#` are skipped;
#   - `KEY=value` or `export KEY=value`, KEY matching [A-Za-z_][A-Za-z0-9_]*,
#     whitespace around key and value trimmed;
#   - one pair of matching quotes around the value is removed, nothing is
#     escaped or interpolated, and `#` inside a value is literal;
#   - lines without `=` or with an invalid key are ignored, a later line for a
#     key overrides an earlier one.
# A missing file (or path "") gives character().
.read_env_file <- function(path) {
  if (!utils::file_test("-f", path)) return(character())
  txt <- readChar(path, file.size(path), useBytes = TRUE)
  lines <- sub("\r$", "", strsplit(txt, "\n", fixed = TRUE)[[1]])
  pat <- "^\\s*(?:export\\s+)?([A-Za-z_][A-Za-z0-9_]*)\\s*=\\s*(.*?)\\s*$"
  hits <- regmatches(lines, regexec(pat, lines, perl = TRUE))
  out <- character()
  for (h in hits[lengths(hits) == 3L]) {
    out[[h[2]]] <- sub("^([\"'])(.*)\\1$", "\\2", h[3], perl = TRUE)
  }
  out
}

# Resolve the machine environment. A variable set in the shell, even to "", is
# used and `.env` is not consulted; otherwise the `.env` value is used. Empty means
# unset. `root` is the repo root: without it `.env` is not read and the scratch
# default is NA. Relative store and scratch paths resolve against `root`; `~` expands.
# Returns list(store, scratch, scratch_source, cores):
#   store          path, or NA when unset (not validated here);
#   scratch        path, defaulting to <root>/.scratch/papkou2023;
#   scratch_source "shell", "env_file" or "default";
#   cores          as.numeric(value), 1 when unset; not validated (see .check_cores).
.resolve_machine_env <- function(root = NULL, shell = Sys.getenv(.env_vars, unset = NA)) {
  env_file <- if (is.null(root)) character() else .read_env_file(file.path(root, ".env"))
  lookup <- function(var) {
    from_shell <- if (var %in% names(shell)) shell[[var]] else NA_character_
    if (!is.na(from_shell)) return(list(value = from_shell, source = "shell"))
    if (var %in% names(env_file)) return(list(value = env_file[[var]], source = "env_file"))
    list(value = "", source = "default")
  }
  as_path <- function(value) {
    path <- path.expand(value)
    if (!is.null(root) && !grepl("^(/|[A-Za-z]:)", path)) file.path(root, path) else path
  }
  store <- lookup("PAPKOU2023RUGGED_DATA")
  scratch <- lookup("BARSEQFIT_SCRATCH")
  cores <- lookup("BARSEQFIT_CORES")
  has_scratch <- nzchar(scratch$value)
  list(
    store = if (nzchar(store$value)) as_path(store$value) else NA_character_,
    scratch = if (has_scratch) as_path(scratch$value)
              else if (!is.null(root)) file.path(root, ".scratch", "papkou2023")
              else NA_character_,
    scratch_source = if (has_scratch) scratch$source else "default",
    cores = if (nzchar(cores$value)) suppressWarnings(as.numeric(cores$value)) else 1)
}

# A subset run (BARSEQFIT_PAIR_LIMIT, read from the shell by the script) must name
# its throwaway scratch directory in the shell, so a partial table cannot overwrite
# the default full-run outputs. `limit` is NA when unset.
.check_pair_limit <- function(limit, menv) {
  if (!is.na(limit) && menv$scratch_source != "shell") {
    stop("BARSEQFIT_PAIR_LIMIT requires an explicit BARSEQFIT_SCRATCH (throwaway directory)",
         call. = FALSE)
  }
  invisible(limit)
}
