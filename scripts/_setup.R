# Shared preamble for the R scripts. The caller sets `root` (the repo root) and
# sources this file: it activates the r/ renv library (RcppTOML, nanoparquet,
# jsonlite) and sources every r/R/*.R file.
local({
  old <- setwd(file.path(root, "r"))
  on.exit(setwd(old))
  source("renv/activate.R")
})
for (f in list.files(file.path(root, "r", "R"), pattern = "\\.R$", full.names = TRUE)) source(f)
