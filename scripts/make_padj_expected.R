# Generate tests/data/expected_{wt,pair}_padj_r.csv: R's Benjamini-Hochberg p_adj
# for every fixture comparison under the default (fixed) dispersion. The R
# package is sourced from r/R, so the numbers come from fit_wt/fit_pair.
# Run: just make-padj-expected
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
root <- if (length(file_arg)) normalizePath(file.path(dirname(file_arg[1]), "..")) else getwd()
source(file.path(root, "scripts", "_setup.R"))
data_dir <- file.path(root, "tests", "data")
td <- read_test_data(data_dir)
out <- function(df, name) {
  utils::write.table(df, file.path(data_dir, name), sep = ",", row.names = FALSE,
                     quote = FALSE, na = "",
                     col.names = TRUE, qmethod = "double", dec = ".")
}
num17 <- function(df, cols) {
  for (col in cols) df[[col]] <- ifelse(is.na(df[[col]]), NA, format(df[[col]], digits = 17, trim = TRUE))
  df
}
wt <- fit_wt(td$counts, td$samples, td$sequences, td$reference_id, td$elapsed_time, td$time_unit,
              pseudocount = "all_zeros", dispersion = "fixed")
out(num17(wt[c("sequence_id", "p_adj")], "p_adj"), "expected_wt_padj_r.csv")
pair <- fit_pair(td$counts, td$samples, td$sequences, td$pairs, td$elapsed_time, td$time_unit,
                pseudocount = "all_zeros", dispersion = "fixed")
out(num17(pair[c("from", "to", "p_adj")], "p_adj"), "expected_pair_padj_r.csv")
