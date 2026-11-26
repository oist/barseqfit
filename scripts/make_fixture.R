# Extract a small, seeded Papkou et al. 2023 subset into tests/data/.
# Usage: just make-fixture   (reads the store from $PAPKOU2023RUGGED_DATA or .env)
# No filtering on zero counts or fitness values.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
root <- if (length(file_arg)) normalizePath(file.path(dirname(file_arg[1]), "..")) else getwd()
source(file.path(root, "scripts", "_setup.R"))

menv <- .resolve_machine_env(root)
out_dir <- file.path(root, "tests", "data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

x <- read_papkou_store(menv$store)
set.seed(20261014)
fx <- .papkou_subset(x, sort(sample.int(nrow(x$pairs), 10L)))

num <- function(v) ifelse(is.na(v), "", sprintf("%.17g", v))
wr <- function(df, name) {
  write.table(df, file.path(out_dir, name), sep = ",", row.names = FALSE,
              quote = FALSE, na = "", eol = "\n")
}
fx$samples$time <- num(fx$samples$time)
wr(fx$counts, "counts.csv")
wr(fx$samples, "samples.csv")
wr(fx$sequences, "sequences.csv")
wr(fx$pairs, "pairs.csv")
wr(data.frame(sequence_id = fx$expected_wt$sequence_id, m = num(fx$expected_wt$m),
              se = num(fx$expected_wt$se), p = num(fx$expected_wt$p)), "expected_wt.csv")
wr(data.frame(from = fx$expected_pair$from, to = fx$expected_pair$to, m = num(fx$expected_pair$m),
              se = num(fx$expected_pair$se), p = num(fx$expected_pair$p)), "expected_pair.csv")

message("Wrote fixture to ", out_dir)
print(fx$pairs)
