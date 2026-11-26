new_root <- function(env_lines = NULL, .env = parent.frame()) {
  root <- withr::local_tempdir(.local_envir = .env)
  if (!is.null(env_lines)) writeLines(env_lines, file.path(root, ".env"))
  root
}

shell_vars <- function(store = NA_character_, scratch = NA_character_, cores = NA_character_) {
  c(PAPKOU2023RUGGED_DATA = store, BARSEQFIT_SCRATCH = scratch, BARSEQFIT_CORES = cores)
}

test_that(".read_env_file parses the shared .env grammar", {
  path <- file.path(withr::local_tempdir(), ".env")
  lines <- c("# comment line", "   # indented comment", "", "PLAIN=value",
             "export EXPORTED=exp", "  SPACED  =  spaced value  ",
             "DQ=\"double quoted\"", "SQ='single quoted'", "QEMPTY=\"\"",
             "QPAD=\"  padded  \"", "MISMATCH=\"mismatched'", "HASH=a # not a comment",
             "EMPTY=", "DUP=first", "DUP=second", "INTERP=${HOME}/x",
             "no equals sign here", "1BAD=invalid key")
  writeBin(charToRaw(paste0(paste0(lines, "\n", collapse = ""), "CRLF=crlf\r\n")), path)
  expected <- list(PLAIN = "value", EXPORTED = "exp", SPACED = "spaced value",
                   DQ = "double quoted", SQ = "single quoted", QEMPTY = "",
                   QPAD = "  padded  ", MISMATCH = "\"mismatched'",
                   HASH = "a # not a comment", EMPTY = "", DUP = "second",
                   INTERP = "${HOME}/x", CRLF = "crlf")
  got <- as.list(.read_env_file(path))
  expect_equal(got[order(names(got))], expected[order(names(expected))])
})

test_that(".read_env_file handles a missing file, an empty file and no final newline", {
  dir <- withr::local_tempdir()
  expect_identical(.read_env_file(file.path(dir, "missing")), character())
  expect_identical(.read_env_file(""), character())
  writeBin(raw(), file.path(dir, "empty"))
  expect_identical(.read_env_file(file.path(dir, "empty")), character())
  writeBin(charToRaw("A=1\nB=2"), file.path(dir, "unterminated"))
  expect_equal(.read_env_file(file.path(dir, "unterminated")), c(A = "1", B = "2"))
})

test_that("the shell wins over .env, even when empty", {
  root <- new_root(c("PAPKOU2023RUGGED_DATA=/from/env"))
  expect_equal(.resolve_machine_env(root, shell_vars(store = "/from/shell"))$store, "/from/shell")
  expect_equal(.resolve_machine_env(root, shell_vars())$store, "/from/env")
  expect_true(is.na(.resolve_machine_env(root, shell_vars(store = ""))$store))
})

test_that("an unset store is NA and is not validated", {
  root <- new_root()
  expect_true(is.na(.resolve_machine_env(root, shell_vars())$store))
  expect_equal(.resolve_machine_env(root, shell_vars(store = "/no/such/dir"))$store, "/no/such/dir")
})

test_that("scratch comes from the shell, then .env, then the default", {
  root <- new_root(c("BARSEQFIT_SCRATCH=/from/env"))
  default <- file.path(root, ".scratch", "papkou2023")
  got <- .resolve_machine_env(root, shell_vars(scratch = "/from/shell"))
  expect_equal(got$scratch, "/from/shell")
  expect_identical(got$scratch_source, "shell")
  got <- .resolve_machine_env(root, shell_vars())
  expect_equal(got$scratch, "/from/env")
  expect_identical(got$scratch_source, "env_file")
  got <- .resolve_machine_env(root, shell_vars(scratch = ""))
  expect_equal(got$scratch, default)
  expect_identical(got$scratch_source, "default")
})

test_that("an empty .env scratch gives the default", {
  root <- new_root(c("BARSEQFIT_SCRATCH="))
  got <- .resolve_machine_env(root, shell_vars())
  expect_equal(got$scratch, file.path(root, ".scratch", "papkou2023"))
  expect_identical(got$scratch_source, "default")
})

test_that("with neither shell nor .env, scratch is the default", {
  root <- new_root()
  got <- .resolve_machine_env(root, shell_vars())
  expect_equal(got$scratch, file.path(root, ".scratch", "papkou2023"))
  expect_identical(got$scratch_source, "default")
})

test_that("store and scratch expand ~ and resolve relative paths against the root", {
  root <- new_root(c("BARSEQFIT_SCRATCH=rel/dir"))
  expect_equal(.resolve_machine_env(root, shell_vars())$scratch, file.path(root, "rel/dir"))
  got <- .resolve_machine_env(root, shell_vars(scratch = "~/x"))
  expect_equal(got$scratch, path.expand("~/x"))
  expect_identical(got$scratch_source, "shell")
  expect_equal(.resolve_machine_env(root, shell_vars(store = "~/store"))$store, path.expand("~/store"))
  expect_equal(.resolve_machine_env(root, shell_vars(store = "data/store"))$store,
               file.path(root, "data/store"))
})

test_that("cores is the numeric value, 1 when unset or empty, NA when not a number", {
  root <- new_root()
  cores <- function(value) .resolve_machine_env(root, shell_vars(cores = value))$cores
  expect_equal(cores(NA_character_), 1)
  expect_equal(cores(""), 1)
  expect_equal(cores("4"), 4)
  expect_true(is.na(cores("abc")))
  expect_equal(.resolve_machine_env(new_root("BARSEQFIT_CORES=3"), shell_vars())$cores, 3)
  expect_equal(.resolve_machine_env(new_root("BARSEQFIT_CORES=3"), shell_vars(cores = ""))$cores, 1)
})

test_that("a root without a .env file is fine", {
  root <- new_root()
  expect_false(file.exists(file.path(root, ".env")))
  got <- .resolve_machine_env(root, shell_vars(store = "/s", cores = "2"))
  expect_equal(got$store, "/s")
  expect_equal(got$cores, 2)
})

test_that("without a root .env is not read and scratch is NA", {
  root <- new_root(c("PAPKOU2023RUGGED_DATA=/from/env", "BARSEQFIT_SCRATCH=/from/env",
                     "BARSEQFIT_CORES=8"))
  withr::local_dir(root)
  got <- .resolve_machine_env(NULL, shell_vars())
  expect_true(is.na(got$store))
  expect_true(is.na(got$scratch))
  expect_identical(got$scratch_source, "default")
  expect_equal(got$cores, 1)
  got <- .resolve_machine_env(NULL, shell_vars(store = "~/s", scratch = "~/x", cores = "2"))
  expect_equal(got$store, path.expand("~/s"))
  expect_equal(got$scratch, path.expand("~/x"))
  expect_identical(got$scratch_source, "shell")
  expect_equal(got$cores, 2)
  expect_equal(.resolve_machine_env(NULL, shell_vars(scratch = "rel/dir"))$scratch, "rel/dir")
  expect_identical(.resolve_machine_env(NULL, shell_vars(scratch = ""))$scratch_source, "default")
})

test_that("resolving does not change the process environment", {
  root <- new_root(c("PAPKOU2023RUGGED_DATA=/from/env", "BARSEQFIT_SCRATCH=/from/env",
                     "BARSEQFIT_CORES=8"))
  before <- Sys.getenv(.env_vars, unset = NA)
  .resolve_machine_env(root, shell_vars())
  .resolve_machine_env(root)
  expect_identical(Sys.getenv(.env_vars, unset = NA), before)
})

test_that(".check_pair_limit refuses a limit unless the scratch came from the shell", {
  root <- new_root(c("BARSEQFIT_SCRATCH=/from/env"))
  from_env <- .resolve_machine_env(root, shell_vars())
  by_default <- .resolve_machine_env(new_root(), shell_vars())
  from_shell <- .resolve_machine_env(root, shell_vars(scratch = "/from/shell"))
  msg <- "BARSEQFIT_PAIR_LIMIT requires an explicit BARSEQFIT_SCRATCH \\(throwaway directory\\)"
  expect_error(.check_pair_limit(10L, from_env), msg)
  expect_error(.check_pair_limit(10L, by_default), msg)
  expect_identical(.check_pair_limit(10L, from_shell), 10L)
  for (menv in list(from_env, by_default, from_shell)) {
    expect_identical(.check_pair_limit(NA_integer_, menv), NA_integer_)
  }
})

test_that("papkou_store reports an unset variable", {
  withr::local_envvar(PAPKOU2023RUGGED_DATA = "")
  expect_error(papkou_store(), paste0("PAPKOU2023RUGGED_DATA is not set\\. Export it in your ",
                                      "shell or put it in a repo-root \\.env file ",
                                      "\\(see \\.env\\.example\\)\\."))
})

test_that("papkou_store reports a missing directory", {
  missing <- file.path(withr::local_tempdir(), "no_such_store")
  withr::local_envvar(PAPKOU2023RUGGED_DATA = missing)
  expect_error(papkou_store(), paste0("PAPKOU2023RUGGED_DATA points to a directory that does ",
                                      "not exist: ", missing), fixed = TRUE)
})

test_that("papkou_store returns the normalised store path", {
  store <- withr::local_tempdir()
  withr::local_envvar(PAPKOU2023RUGGED_DATA = store)
  expect_identical(papkou_store(), normalizePath(store))
})

test_that(".papkou_store_path rejects NA and empty paths", {
  expect_error(.papkou_store_path(NA_character_), "is not set")
  expect_error(.papkou_store_path(""), "is not set")
})

test_that(".papkou_store_path names the store, not the variable, for an explicit path", {
  missing <- file.path(tempdir(), "absent-store")
  expect_error(.papkou_store_path(missing, "store"),
               paste0("store points to a directory that does not exist: ", missing), fixed = TRUE)
  expect_error(.papkou_store_path(NA_character_, "store"), "^store is not set\\.$")
})

test_that(".check_cores validates and returns an integer", {
  expect_identical(.check_cores(4), 4L)
  expect_identical(.check_cores(1L), 1L)
  expect_error(.check_cores(4.5), "cores must be a positive integer")
  expect_error(.check_cores(NA), "cores must be a positive integer")
  expect_error(.check_cores(NA_real_), "cores must be a positive integer")
  expect_error(.check_cores(0), "cores must be a positive integer")
  expect_error(.check_cores(NULL), "cores must be a positive integer")
})

test_that("the fits do not read BARSEQFIT_CORES", {
  elapsed <- max(td$samples$time) - min(td$samples$time)
  fit_all <- function() {
    list(fit_wt(td$counts, td$samples, td$sequences, td$reference_id, elapsed, "generations",
                pseudocount = "all_zeros", dispersion = "fixed"),
         fit_pair(td$counts, td$samples, td$sequences, td$pairs, elapsed, "generations",
                 pseudocount = "all_zeros", dispersion = "fixed"))
  }
  base <- fit_all()
  withr::local_envvar(BARSEQFIT_CORES = "0")
  expect_identical(fit_all(), base)
})
