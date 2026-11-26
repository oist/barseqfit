# Shared chunked parallel map for fit_wt and fit_pair.

.check_cores <- function(cores) {
  if (!(is.numeric(cores) && length(cores) == 1L && is.finite(cores) &&
        cores >= 1 && cores == round(cores))) {
    stop("cores must be a positive integer", call. = FALSE)
  }
  as.integer(cores)
}

.check_chunk_size <- function(chunk_size) {
  if (!(is.numeric(chunk_size) && length(chunk_size) == 1L && is.finite(chunk_size) &&
        chunk_size >= 1 && chunk_size == round(chunk_size))) {
    stop("chunk_size must be a positive integer", call. = FALSE)
  }
  invisible(as.integer(chunk_size))
}

# Apply `fn` to each element of `items` and return the flat list of results, in
# order. With cores > 1, chunks of `chunk_size` items run under
# parallel::mclapply (forked); the output does not depend on cores or chunk_size.
# If option `barseqfit.progress` is a function, it is called as f(chunk, n_chunks)
# after each chunk finishes (in the worker, so it should only log).
# If `collect` is a function it is applied to each chunk's list of results inside
# the worker, and the list of per-chunk values (in order) is returned instead of
# the flat list; this keeps the parent from holding millions of small lists.
.parallel_map <- function(items, fn, cores, chunk_size, collect = NULL) {
  progress <- getOption("barseqfit.progress")
  n_items <- length(items)
  if (cores <= 1L && !is.function(progress)) {
    res <- lapply(items, fn)
    return(if (is.function(collect)) list(collect(res)) else res)
  }
  chunks <- split(seq_len(n_items), ceiling(seq_len(n_items) / chunk_size))
  n_chunks <- length(chunks)
  run_chunk <- function(k) {
    out <- lapply(items[chunks[[k]]], fn)
    if (is.function(collect)) out <- collect(out)
    if (is.function(progress)) progress(k, n_chunks)
    out
  }
  out <- if (cores > 1L) {
    parallel::mclapply(seq_len(n_chunks), run_chunk, mc.cores = cores, mc.preschedule = FALSE)
  } else {
    lapply(seq_len(n_chunks), run_chunk)
  }
  # mclapply yields NULL for a child that died and a try-error for one that errored.
  bad <- vapply(out, function(x) is.null(x) || inherits(x, "try-error"), logical(1))
  if (any(bad)) {
    first <- out[[which(bad)[1]]]
    stop("parallel fit failed: ",
         if (is.null(first)) "a worker process died (NULL result)" else as.character(first),
         call. = FALSE)
  }
  if (is.function(collect)) return(out)
  unlist(out, recursive = FALSE, use.names = FALSE)
}
