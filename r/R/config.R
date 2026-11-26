# Loader for the committed run config validation.toml.
# Only the validation scripts call this; fit_wt and fit_pair take their
# parameters as arguments.

# Parse `path` (default: <repo root>/validation.toml) into a nested list.
.load_run_config <- function(path = file.path(.repo_root(), "validation.toml")) {
  if (!file.exists(path)) stop("run config not found: ", path, call. = FALSE)
  RcppTOML::parseTOML(path, fromFile = TRUE)
}

# Serialise a parsed config (nested named list of numbers and strings) to JSON
# text, for embedding in manifests and timing files. Doubles use %.15g.
.config_json <- function(x) {
  if (is.list(x)) {
    keys <- names(x)
    if (length(x) > 0 && (is.null(keys) || any(!nzchar(keys)))) {
      stop("config lists must be fully named", call. = FALSE)
    }
    body <- vapply(seq_along(x), function(i) {
      sprintf("\"%s\": %s", keys[[i]], .config_json(x[[i]]))
    }, character(1))
    return(paste0("{", paste(body, collapse = ", "), "}"))
  }
  if (length(x) != 1) stop("config values must be length-1 scalars", call. = FALSE)
  if (is.numeric(x) && !is.finite(x)) stop("config numbers must be finite", call. = FALSE)
  if (is.character(x)) return(sprintf("\"%s\"", gsub("([\"\\\\])", "\\\\\\1", x)))
  if (is.logical(x)) return(tolower(as.character(x)))
  if (is.integer(x)) return(as.character(x))
  sprintf("%.15g", x)
}

# Write a timing file (one JSON line plus a final newline). `counts` is a named
# list of the run's sizes, written after r_version in the given order: an integer
# per entry, `null` when NA (the pair run's pair_limit). `export_manifest` (the
# parsed export_manifest.json, or NULL) supplies only its timestamp and the store
# md5s; `config` is the Run config of this run.
.write_timing <- function(path, wall_seconds, cores, counts, export_manifest, config) {
  stamp <- export_manifest$timestamp
  md5s <- export_manifest$store_files$md5
  sizes <- vapply(names(counts), function(nm) {
    v <- counts[[nm]]
    sprintf("\"%s\": %s", nm, if (is.na(v)) "null" else sprintf("%d", as.integer(v)))
  }, character(1))
  parts <- c(
    sprintf("\"wall_seconds\": %.3f", wall_seconds),
    sprintf("\"cores\": %d", as.integer(cores)),
    sprintf("\"r_version\": \"%s\"", R.version.string),
    sizes,
    sprintf("\"export_timestamp\": %s",
            if (length(stamp) == 1L) sprintf("\"%s\"", stamp) else "null"),
    sprintf("\"store_md5\": %s",
            if (length(md5s)) sprintf("[%s]", paste(sprintf("\"%s\"", md5s), collapse = ", ")) else "null"),
    sprintf("\"config\": %s", .config_json(config)))
  writeLines(paste0("{", paste(parts, collapse = ", "), "}\n"), path, sep = "")
  invisible(path)
}
