# Fit the full Papkou et al. 2023 edge table (3.5M pairs) with the R reference implementation.
# Usage: BARSEQFIT_CORES=6 Rscript scripts/run_full_pair.R   (needs `just export-papkou` first)
# Reads the Export in $BARSEQFIT_SCRATCH (default .scratch/papkou2023/) and writes
# fit_pair_full_r.parquet and timing_pair_r.json there.
# BARSEQFIT_PAIR_LIMIT=N fits only the first N pairs (smoke runs; refused unless BARSEQFIT_SCRATCH is
# set explicitly in the shell, so a partial table cannot overwrite the default full-run outputs).
# Progress is logged to stderr as each chunk finishes (chunks complete out of order under mclapply). timing_pair_r.json records the export manifest timestamp and store md5s the fit was made from.
# Variables are read from the shell, with a repo-root .env filling any that are unset.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
root <- if (length(file_arg)) normalizePath(file.path(dirname(file_arg[1]), "..")) else getwd()
source(file.path(root, "scripts", "_setup.R"))

# Shell wins, .env fills gaps (PAPKOU2023RUGGED_DATA, BARSEQFIT_SCRATCH, BARSEQFIT_CORES).
menv <- .resolve_machine_env(root)
# A subset run must name its throwaway directory in the shell, not in .env.
limit <- suppressWarnings(as.integer(Sys.getenv("BARSEQFIT_PAIR_LIMIT", unset = "")))
.check_pair_limit(limit, menv)
# Result-affecting parameters; recorded in timing_pair_r.json.
config <- .load_run_config(file.path(root, "validation.toml"))
x <- read_papkou_export(menv$scratch)
n_ids <- nrow(x$sequences)
stopifnot(x$reference_id %in% x$sequences$sequence_id)

pairs <- x$pairs
if (!is.na(limit)) pairs <- pairs[seq_len(min(limit, nrow(pairs))), , drop = FALSE]
n_pairs <- nrow(pairs)

cores <- .check_cores(menv$cores)
chunk_size <- as.integer(config$fit$chunk_size)
ts <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
t0 <- proc.time()[["elapsed"]]
options(barseqfit.progress = function(k, n_chunks) {
  el <- proc.time()[["elapsed"]] - t0
  message(sprintf("[%s] chunk %d finished of %d (elapsed %.0f s)", ts(), k, n_chunks, el))
})
message("[", ts(), "] Fitting ", n_pairs, " pairs (", n_ids, " genotypes) with cores = ", cores,
        ", chunk_size = ", chunk_size, " (~", chunk_size, " pairs/chunk)")
res <- fit_pair(x$counts, x$samples, x$sequences, pairs,
                elapsed_time = x$elapsed_time, time_unit = x$time_unit,
                pseudocount = "all_zeros", dispersion = "fixed", cores = cores,
                chunk_size = chunk_size)
wall <- proc.time()[["elapsed"]] - t0
stopifnot(nrow(res) == n_pairs)

nanoparquet::write_parquet(res, file.path(menv$scratch, "fit_pair_full_r.parquet"))
# Tie this fit to the export it was made from: the gate refuses a table whose
# timestamp or store md5s differ from the current export_manifest.json.
.write_timing(file.path(menv$scratch, "timing_pair_r.json"), wall, cores,
              list(n_genotypes = n_ids, n_pairs = n_pairs, pair_limit = limit),
              x$manifest, config)
message("Wrote ", menv$scratch, " (", nrow(res), " rows, ", sprintf("%.1f", wall), " s)")
print(table(res$status))
