date_to_int <- function(date, origin = "1970-01-01") {
  int <- as.integer(as.Date(date) - as.Date(origin))
  
  d <- dim(date)
  if (!is.null(d)) {
    int <- array(int, d)
  }
  
  int
}


int_to_date <- function(int, origin = "1970-01-01") {
  int + as.Date(origin)
}


data_frame_to_array <- function(df) {
  array(unlist(df), dim(df))
}


`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


squote <- function(x) {
  sprintf("'%s'", x)
}


visFALSE <- function(x) {
  !is.na(x) & x == FALSE
}

#' Clean group names from delay map
#' @param raw_group The group element from a delay_map row
#' @noRd
clean_group_name <- function(raw_group) {
  raw_group <- as.character(raw_group)
  
  if (length(raw_group) == 1 && grepl("^c\\(", raw_group)) {
    raw_group <- gsub("^c\\(|\\)$", "", raw_group)
    raw_group <- gsub("[\"']", "", raw_group)
    raw_group <- trimws(strsplit(raw_group, ",")[[1]])
  }
  
  clean_group <- paste(raw_group, collapse = ", ")
  clean_group <- gsub("[-_]", " ", clean_group)
  
  tools::toTitleCase(clean_group)
}

#' Clean event names from delay map
#' @param raw_event The from/to element from a delay_map row
#' @noRd
clean_event_name <- function(raw_event) {
  clean_name <- gsub("[-_]", " ", as.character(raw_event))
  tools::toTitleCase(clean_name)
}

vnapply <- function(...) {
  vapply(..., FUN.VALUE = 1)
}

suggested_package_status <- function(pkg) {
  status <- ifelse(pkg %in% loadedNamespaces(), "loaded", "missing")
  i <- status == "missing"
  if (any(i)) {
    status[i] <- ifelse(pkg[i] %in% .packages(TRUE), "installed", "missing")
  }
  set_names(status, pkg)
}

set_names <- function(x, nms) {
  if (length(nms) == 1 && length(x) != 1) {
    nms <- rep_len(nms, length(x))
  }-
  names(x) <- nms
  x
}


array_reshape <- function(x, i, d, call = NULL) {
  dx <- dim2(x)
  if (length(dx) < i) {
    cli::cli_abort(
      "array only has {length(dx)} dimension{?s}, can't update dimension {i}",
      call = call)
  }
  if (dx[[i]] != prod(d)) {
    dim_str <- paste(d, collapse = ", ")
    cli::cli_abort(
      paste("New dimensions ({dim_str}) imply dimension {i} has length",
            "{prod(d)} but found {dx[[i]]}"),
      call = call)
  }
  
  dn <- dimnames2(x)
  
  ## The actual reshape is easy:
  dim(x) <- append(dx[-i], d, i - 1L)
  
  ## Can't preserve dimension names on modified dimensions
  if (!is.null(dn)) {
    dimnames(x) <- append(dn[-i], rep(list(NULL), length(d)), i - 1L)
  }
  x
}


dim2 <- function(x) {
  dim(x) %||% length(x)
}


dimnames2 <- function(x) {
  if (!is.null(dim(x))) {
    dimnames(x)
  } else {
    nms <- names(x)
    if (is.null(nms)) NULL else list(nms)
  }
}
