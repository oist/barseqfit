# When run via test_dir(load_package = "source") the package R/ files are
# loaded by testthat; otherwise (plain test_dir) source them here.
if (!exists("fit_wt", mode = "function")) {
  for (f in list.files(file.path("..", "..", "R"), full.names = TRUE, pattern = "\\.R$")) source(f)
}
td <- read_test_data()
