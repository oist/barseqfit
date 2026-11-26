# Export full-data Papkou et al. 2023 tables from the raw Zenodo store into a scratch dir.
# Usage: just export-papkou   (or PAPKOU2023RUGGED_DATA=/path Rscript scripts/export_papkou.R)
# Output dir: $BARSEQFIT_SCRATCH (default <repo>/.scratch/papkou2023). Gitignored; never commit.
# Tables are written as Parquet via nanoparquet (full double precision, NA as null).
# validation.toml is copied into export_manifest.json.

t0 <- Sys.time()
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
root <- if (length(file_arg)) normalizePath(file.path(dirname(file_arg[1]), "..")) else getwd()
source(file.path(root, "scripts", "_setup.R"))

# Shell wins, .env fills gaps (PAPKOU2023RUGGED_DATA, BARSEQFIT_SCRATCH, BARSEQFIT_CORES).
menv <- .resolve_machine_env(root)
config <- .load_run_config(file.path(root, "validation.toml"))
x <- read_papkou_store(menv$store)
if (x$manifest$edge_rows_dropped != 0L) {
  warning(x$manifest$edge_rows_dropped, " edge rows dropped (id not in the after matrix); expected 0")
}
dir.create(menv$scratch, showWarnings = FALSE, recursive = TRUE)
out_dir <- normalizePath(menv$scratch)
.write_papkou_export(x, out_dir, config, format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))

message("Wrote export to ", out_dir, " (", nrow(x$counts), " count rows, ", nrow(x$sequences),
        " genotypes, ", nrow(x$before_only_ids), " before-only) in ",
        round(as.numeric(Sys.time() - t0, units = "secs"), 1), " s")
