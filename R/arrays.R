array_drop <- function(x, i, call = NULL) {
  dx <- dim2(x)
  rank <- length(dx)
  if (rank < max(i)) {
    cli::cli_abort(
      "Can't update dimension {max(i)}, array only has {rank} dimension{?s}",
      call = call)
  }
  if (!all(dx[i] == 1)) {
    err <- i[dx[i] != 1]
    ## I can't make cli's pluralisation work here easily, doing that
    ## manually:
    if (length(err) == 1L) {
      cli::cli_abort(
        "Can't drop dimension {err} as it is length {dx[err]}, not 1",
        call = call)
    } else {
      cli::cli_abort(
        "Can't drop dimensions {err} as they are length {dx[err]}, not 1",
        call = call)
    }
  }
  dn <- dimnames(x)
  dim(x) <- dim(x)[-i]
  if (!is.null(dn)) {
    if (rank == 2) {
      x <- c(x)
      names(x) <- dn[[-i]]
    } else {
      dimnames(x) <- dn[-i]
    }
  }
  x
}
