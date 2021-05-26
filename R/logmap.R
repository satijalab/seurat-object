#' @include zzz.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' A Logical Map
#'
#' A simple container for storing mappings of values using logical matrices.
#' Keeps track of which values (rows) are present in which observations
#' (columns). \code{LogMap} objects can be created with \code{LogMap}; queries
#' can be performed with \code{[[} and observations can be added or removed
#' with \code{[[<-}
#'
#' @slot .Data A logical matrix with at least one row
#'
setClass(
  Class = 'LogMap',
  contains = 'matrix'
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param y A character vector
#'
#' @return \code{LogMap}: A new \code{LogMap} object with zero columns and
#' \code{length(x = x)} rows; rownames are set to \code{x}
#'
#' @rdname LogMap-class
#'
#' @export
#'
LogMap <- function(y) {
  if (!is.character(x = y)) {
    y <- as.character(x = y)
  }
  return(new(
    Class = 'LogMap',
    .Data = matrix(nrow = length(x = y), ncol = 0, dimnames = list(y, NULL))
  ))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname LogMap-class
#'
#' @param x A \code{LogMap} object
#' @param i A character vector of length 1, or (for \code{[[}) \code{NULL}
#'
#' @return \code{[[}: if \code{i} is a character vector, the rownames that are
#' mapped to \code{i}; otherwise the rownames of \code{x}
#'
#' @export
#'
setMethod(
  f = '[[',
  signature = c(x = 'LogMap', i = 'character', j = 'missing'),
  definition = function(x, i, ...) {
    i <- i[1]
    i <- match.arg(arg = i, choices = colnames(x = x))
    return(rownames(x = x)[x[, i, drop = TRUE]])
  }
)

#' @rdname LogMap-class
#'
#' @export
#'
setMethod(
  f = '[[',
  signature = c(x = 'LogMap', i = 'NULL', j = 'missing'),
  definition = function(x, i, ...) {
    return(rownames(x = x))
  }
)

#' @rdname LogMap-class
#'
#' @param value A character or integer vector of values to record in the map
#' for \code{i}, or \code{NULL} to remove the record for \code{i}
#'
#' @return \code{[[<-}: If \code{value} is \code{NULL}, then \code{x} without
#' the observations for \code{i}; otherwise, \code{x} with a new column for
#' \code{i} recording a \code{TRUE} for all values present in \code{value}
#'
#' @export
#'
setMethod(
  f = '[[<-',
  signature = c(
    x = 'LogMap',
    i = 'character',
    j = 'missing',
    value = 'character'
  ),
  definition = function(x, i, ..., value) {
    value <- MatchCells(new = rownames(x = x), orig = value)
    x[[i]] <- value
    return(x)
  }
)

#' @rdname LogMap-class
#'
#' @export
#'
setMethod(
  f = '[[<-',
  signature = c(
    x = 'LogMap',
    i = 'character',
    j = 'missing',
    value = 'integer'
  ),
  definition = function(x, i, ..., value) {
    if (i %in% colnames(x = x)) {
      x[[i]] <- NULL
    }
    value <- MatchCells(new = value, rownames(x = x))
    mat <- slot(object = x, name = '.Data')
    if (length(x = value)) {
      mat <- cbind(
        mat,
        matrix(data = FALSE, nrow = nrow(x = x), dimnames = list(NULL, i))
      )
      mat[value, i] <- TRUE
    }
    slot(object = x, name = '.Data') <- mat
    validObject(object = x)
    return(x)
  }
)

#' @rdname LogMap-class
#'
#' @export
#'
setMethod(
  f = '[[<-',
  signature = c(x = 'LogMap', i = 'character', j = 'missing', value = 'NULL'),
  definition = function(x, i, ..., value) {
    mat <- slot(object = x, name = '.Data')
    for (name in i) {
      idx <- which(x = colnames(x = mat) == name)
      if (length(x = idx)) {
        mat <- mat[, -idx, drop = FALSE]
      }
    }
    slot(object = x, name = '.Data') <- mat
    validObject(object = x)
    return(x)
  }
)

#' @rdname LogMap-class
#'
#' @export
#'
setMethod(
  f = '[[<-',
  signature = c(
    x = 'LogMap',
    i = 'character',
    j = 'missing',
    value = 'numeric'
  ),
  definition = function(x, i, ..., value) {
    value <- as.integer(x = value)
    x[[i]] <- value
    return(x)
  }
)

setMethod(
  f = 'show',
  signature = 'LogMap',
  definition = function(object) {
    cat(
      "A logical map for",
      nrow(x = object),
      "values across",
      ncol(x = object),
      "observations"
    )
    return(invisible(x = NULL))
  }
)

setValidity(
  Class = 'LogMap',
  method = function(object) {
    valid <- NULL
    # Ensure we have a logical matrix
    if (!is.logical(x = object)) {
      valid <- c(valid, "The map must be a logical matrix")
    }
    # Check rownames
    if (is.null(x = rownames(x = object))) {
      valid <- c(valid, "Rownames must be supplied")
    } else if (anyDuplicated(x = rownames(x = object))) {
      valid <- c(valid, "Duplicate rownames not allowed")
    }
    # Check colnames
    if (!is.null(x = colnames(x = object))) {
      if (!all(nchar(x = colnames(x = object)))) {
        valid <- c(valid, "Columnn names cannot be empty strings")
      }
      if (anyDuplicated(x = colnames(x = object))) {
        valid <- c(valid, "Duplicate colnames not allowed")
      }
    }
    return(valid %||% TRUE)
  }
)
