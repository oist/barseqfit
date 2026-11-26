# Fit the full Papkou et al. 2023 WT table with the R reference implementation.
# Usage: BARSEQFIT_CORES=6 Rscript scripts/run_full_wt.R   (needs `just export-papkou` first)
# Reads the Export in $BARSEQFIT_SCRATCH (default .scratch/papkou2023/) and writes
# fit_wt_full_r.parquet and timing_r.json there.
# timing_r.json records the export manifest timestamp and store md5s the fit was made from.
# Variables are read from the shell, with a repo-root .env filling any that are unset.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
root <- if (length(file_arg)) normalizePath(file.path(dirname(file_arg[1]), "..")) else getwd()
source(file.path(root, "scripts", "_setup.R"))

# Shell wins, .env fills gaps (PAPKOU2023RUGGED_DATA, BARSEQFIT_SCRATCH, BARSEQFIT_CORES).
menv <- .resolve_machine_env(root)
# Result-affecting parameters; recorded in timing_r.json.
config <- .load_run_config(file.path(root, "validation.toml"))
x <- read_papkou_export(menv$scratch, c("counts", "samples", "sequences"))
n_ids <- nrow(x$sequences)
stopifnot(x$reference_id %in% x$sequences$sequence_id)

cores <- .check_cores(menv$cores)
message("Fitting ", n_ids, " genotypes with cores = ", cores)
t0 <- proc.time()[["elapsed"]]
res <- fit_wt(x$counts, x$samples, x$sequences, reference_id = x$reference_id,
              elapsed_time = x$elapsed_time, time_unit = x$time_unit,
              pseudocount = "all_zeros", dispersion = "fixed", cores = cores,
              chunk_size = as.integer(config$fit$chunk_size))
wall <- proc.time()[["elapsed"]] - t0
stopifnot(nrow(res) == n_ids)  # reference row included

nanoparquet::write_parquet(res, file.path(menv$scratch, "fit_wt_full_r.parquet"))
# Tie this fit to the export it was made from: the gate refuses a table whose
# timestamp or store md5s differ from the current export_manifest.json.
.write_timing(file.path(menv$scratch, "timing_r.json"), wall, cores,
              list(n_genotypes = nrow(res)), x$manifest, config)
message("Wrote ", menv$scratch, " (", nrow(res), " rows, ", sprintf("%.1f", wall), " s)")
print(table(res$status))
