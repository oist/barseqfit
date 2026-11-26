set shell := ["bash", "-eu", "-o", "pipefail", "-c"]

# Rscript to use; R 4.1.3 (the Papkou2023 reference version) is required. Override with RSCRIPT.
rscript := env_var_or_default("RSCRIPT", "Rscript")

# Install the default development environment.
sync:
    uv sync

# Run the Python test suite.
test:
    uv run pytest tests

# Restore the renv library (RcppTOML, nanoparquet, ...) under R 4.1.3.
renv-restore:
    cd r && {{rscript}} -e 'renv::restore(prompt = FALSE)'

# Run the R package testthat suite under R 4.1.3 (restores the renv library first).
test-r: renv-restore
    cd r && {{rscript}} -e 'testthat::test_dir("tests/testthat", load_package = "source", stop_on_failure = TRUE)'

# Run the Python and R test suites.
test-all: test test-r

# Alias for test-all.
check: test-all

# Regenerate tests/data from the Papkou2023 store ($PAPKOU2023RUGGED_DATA).
make-fixture: renv-restore
    {{rscript}} scripts/make_fixture.R

# Regenerate tests/data/expected_*_padj_r.csv (R's BH-adjusted p) from the R package.
make-padj-expected: renv-restore
    {{rscript}} scripts/make_padj_expected.R

# Export full-data Papkou2023 tables to $BARSEQFIT_SCRATCH (default .scratch/papkou2023/) from the store.
export-papkou: renv-restore
    {{rscript}} scripts/export_papkou.R

# Run the authors' Papkou2023 fitness code from the store ($PAPKOU2023RUGGED_DATA, md5-checked) on tests/data; installs data.table into the renv library if missing.
validate-papkou2023:
    cd r && {{rscript}} -e 'if (!requireNamespace("data.table", quietly=TRUE)) install.packages("data.table", repos="https://cloud.r-project.org")'
    cd r && {{rscript}} validation/run_authors_pipeline.R

# Fit the full Papkou2023 WT table with R (set BARSEQFIT_CORES for parallelism; writes to $BARSEQFIT_SCRATCH).
full-wt-r: renv-restore
    {{rscript}} scripts/run_full_wt.R

# Full-data WT gate: Python vectorised vs scalar sample vs published vs R, plus the fit_sequences identity (needs `just export-papkou` first).
full-wt-gate:
    uv run python scripts/full_wt_gate.py

# Export, fit in R, then run the Python gate.
validate-full-wt: export-papkou full-wt-r full-wt-gate

# Fit the full Papkou2023 edge table with R (set BARSEQFIT_CORES; ~2 h on 4 cores; writes to $BARSEQFIT_SCRATCH).
full-pair-r: renv-restore
    {{rscript}} scripts/run_full_pair.R

# Full-data pair gate: Python vectorised vs scalar sample vs published vs R (needs `just export-papkou` first).
full-pair-gate:
    uv run python scripts/full_pair_gate.py

# Export, fit pairs in R, then run the Python gate.
validate-full-pair: export-papkou full-pair-r full-pair-gate
