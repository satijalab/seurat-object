#' @include zzz.R
#' @importFrom methods as callNextMethod
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setGeneric(
  name = '.GetLayerData',
  def = function(x, ...) {
    standardGeneric(f = '.GetLayerData')
  },
  signature = c('x')
)

setGeneric(
  name = '.PrepLayerData',
  def = function(x, target, transpose = NULL, ...) {
    standardGeneric(f = '.PrepLayerData')
  },
  signature = c('x', 'target')
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.CheckDimnames <- function(dnames, dims) {
  if (is.null(x = dnames)) {
    return(NULL)
  }
  if (!is.numeric(x = dims)) {
    stop("'dims' must be a numeric vector", call. = FALSE)
  } else if (!is.list(x = dnames) || length(x = dnames) != length(x = dims)) {
    stop("'dnames' must be a list of length ", length(x = dims), call. = FALSE)
  }
  dims <- dims %/% 1L
  didx <- match(
    x = vapply(
      X = dnames,
      FUN = length,
      FUN.VALUE = numeric(length = 1L),
      USE.NAMES = FALSE
    ),
    table = dims
  )
  didx[duplicated(x = didx)] <- NA
  didx[is.na(x = didx)] <- setdiff(
    x = seq_along(along.with = dims),
    y = didx
  )
  return(dnames[didx])
}

.CheckFmargin <- function(fmargin) {
  fmargin <- fmargin %/% 1L
  if (!fmargin %in% seq.int(from = 1L, to = 2L)) {
    stop("'fmargin' must be either 1 or 2", call. = FALSE)
  }
  return(fmargin)
}

.Cmargin <- function(fmargin) {
  fmargin <- .CheckFmargin(fmargin = fmargin)
  return(setdiff(x = seq.int(from = 1L, to = 2L), y = fmargin))
}

#' Get and Prepare Layer Data
#'
#' Assemble layer data pulled from and prepare layer data for addition to
#' v5 assays; v5 assays allow layers to be in multiple formats and support
#' both regular and transposed orientations
#'
#' When transposition is required, \code{.GetLayerData2} and
#' \code{.PrepLayerData2} will attempt to
#' \link[.GetMethod]{determine the optimal method} of \code{\link[base:t]{t()}}
#' to use; if no optimal method is found, \code{base::t.default} will be used
#' for transposition, which may be slow
#'
#' @param x A matrix-like object
#' @param dnames An optional list with feature and cell names
#' (in order for \code{.GetLayerData2})
#' @param fmargin Margin for features (1 or 2); for \code{.GetLayerData2}, if
#' \code{fmargin} is 2, \code{x} will be transposed (see details)
#'
#' @return \code{.GetLayerData2}: \code{x}, potentially transposed and
#' potentially with \code{dnames} set as the \code{\link{dimnames}}
#'
#' @keywords internal
#'
#' @noRd
#'
.GetLayerData2 <- function(x, dnames = NULL, fmargin = 1L) {
  # Check dimnames
  if (!is.null(x = dnames)) {
    ndims <- length(x = dim(x = x))
    if (!is_bare_list(x = dnames, n = ndims)) {
      abort(message = paste("'dnames' must be a list of length", ndims))
    }
    didx <- match(
      x = sapply(X = dnames, FUN = length, USE.NAMES = FALSE),
      table = dim(x = x)
    )
    didx[duplicated(x = didx)] <- NA
    didx[is.na(x = didx)] <- setdiff(x = seq_len(length.out = ndims), y = didx)
    dnames <- dnames[didx]
  }
  # Check fmargin
  fmargin <- fmargin %/% 1L
  if (!fmargin %in% c(1L, 2L)) {
    abort(message = "'fmargin' must be either 1 or 2")
  }
  # Do we transpose
  if (fmargin == 2L) {
    tf <- .GetMethod(fxn = 't', cls = class(x = x))
    x <- tf(x)
    dnames <- rev(x = dnames)
  }
  suppressWarnings(expr = suppressMessages(expr = dimnames(x = x) <- dnames))
  return(x)
}

#' @param target An optional two-length integer vector with dimensions of
#' the v5 assay that \code{x} will be added to; used only if
#' \code{transpose} is \code{NULL}
#' @param transpose Transpose \code{x} before returning it; if \code{NULL} and
#' \code{target} is provided, will attempt to determine if transposition is
#' necessary (see details)
#'
#' @return \code{.PrepLayerData2}: \code{x} with \code{\link{dimnames}} removed
#' and \code{dnames} added as attributes \dQuote{\code{features}} and
#' \dQuote{\code{cells}} and potentially transposed
#'
#' @rdname dot-GetLayerData2
#'
#' @keywords internal
#'
#' @noRd
#'
.PrepLayerData2 <- function(
  x,
  target = NULL,
  transpose = FALSE,
  dnames = NULL,
  fmargin = 1L
) {
  if (is.null(x = x)) {
    return(x)
  }
  # Check fmargin
  fmargin <- fmargin %/% 1L
  if (!fmargin %in% c(1L, 2L)) {
    abort(message = "'fmargin' must be either 1 or 2")
  }
  cmargin <- c(2L, 1L)[fmargin]
  # Auto-check transposition
  if (!is.null(x = target) && is.null(x = transpose)) {
    if (!is_bare_integerish(x = target, n = 2L, finite = TRUE)) {
      abort(message = "'target' must be a two-length integer vector")
    }
    xdim <- dim(x)[c(fmargin, cmargin)]
    if (all(rev(xdim) == target)) {
      transpose <- TRUE
    }
  }
  # Check dimnames
  if (is.null(x = dnames)) {
    dnames <- dimnames(x = x)
  } else if (!is_bare_list(x = dnames, n = 2L)) {
    abort(message = "'dnames' must be a two-length list")
  }
  # Handle transposition
  if (isTRUE(x = transpose)) {
    tf <- .GetMethod(fxn = 't', cls = class(x = x))
    x <- tf(x)
    dnames <- rev(x = dnames)
  }
  x <- suppressMessages(expr = unname(x))
  attr(x = x, which = 'features') <- dnames[[fmargin]]
  attr(x = x, which = 'cells') <- dnames[[cmargin]]
  return(x)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setMethod(
  f = '.GetLayerData',
  signature = c(x = 'ANY'),
  definition = function(x, dnames = NULL, fmargin = 1L, ...) {
    dnames <- .CheckDimnames(dnames = dnames, dims = dim(x = x))
    fmargin <- .CheckFmargin(fmargin = fmargin)
    if (fmargin == 2L) {
      x <- base::t(x = x)
      dnames <- rev(x = dnames)
    }
    dimnames(x = x) <- dnames
    return(x)
  }
)

setMethod(
  f = '.GetLayerData',
  signature = c(x = 'data.frame'),
  definition = function(x, dnames = NULL, fmargin = 1L, ...) {
    x <- callNextMethod()
    return(as.data.frame(x = x))
  }
)

#' @importFrom Matrix t
#' @importClassesFrom Matrix Matrix
#'
setMethod(
  f = '.GetLayerData',
  signature = c(x = 'Matrix'),
  definition = function(x, dnames = NULL, fmargin = 1L, ...) {
    dnames <- .CheckDimnames(dnames = dnames, dims = dim(x = x))
    fmargin <- .CheckFmargin(fmargin = fmargin)
    if (fmargin == 2L) {
      x <- Matrix::t(x = x)
      dnames <- rev(x = dnames)
    }
    dimnames(x = x) <- dnames
    return(x)
  }
)

#' @importFrom spam t
#' @importClassesFrom spam spam
#'
setMethod(
  f = '.GetLayerData',
  signature = c(x = 'spam'),
  definition = function(x, fmargin = 1L, ...) {
    fmargin <- .CheckFmargin(fmargin = fmargin)
    if (fmargin == 2L) {
      x <- spam::t(x = x)
    }
    return(x)
  }
)

setMethod(
  f = '.PrepLayerData',
  signature = c(x = 'ANY', target = 'numeric'),
  definition = function(
    x,
    target,
    transpose = NULL,
    dnames = NULL,
    fmargin = 1L,
    ...
  ) {
    # Check the value of target
    # target should be c(nfeatures, ncells)
    if (length(x = target) != 2L) {
      stop("'target' must be a two-length numeric vector", call. = FALSE)
    }
    # If transpose is NULL, try to determine if we should transpose
    if (is.null(x = transpose)) {
      fmargin <- .CheckFmargin(fmargin = fmargin)
      cmargin <- .Cmargin(fmargin = fmargin)
      xdim <- dim(x = x)[c(fmargin, cmargin)]
      if (all(rev(x = xdim) == target)) {
        transpose <- TRUE
      }
    }
    return(.PrepLayerData(
      x = x,
      target = NULL,
      transpose = transpose,
      dnames = dnames,
      fmargin = fmargin,
      ...
    ))
  }
)

setMethod(
  f = '.PrepLayerData',
  signature = c(x = 'ANY', target = 'NULL'),
  definition = function(
    x,
    target,
    transpose = NULL,
    dnames = NULL,
    fmargin = 1L,
    tf = base::t,
    ...
  ) {
    fmargin <- .CheckFmargin(fmargin = fmargin)
    cmargin <- .Cmargin(fmargin = fmargin)
    if (isTRUE(x = transpose)) {
      x <- tf(x)
      dnames <- rev(x = dnames)
    }
    x <- suppressMessages(expr = unname(x))
    attr(x = x, which = 'features') <- dnames[[fmargin]]
    attr(x = x, which = 'cells') <- dnames[[cmargin]]
    return(x)
  }
)

#' @importFrom Matrix t
#' @importClassesFrom Matrix Matrix
#'
setMethod(
  f = '.PrepLayerData',
  signature = c(x = 'Matrix', target = 'NULL'),
  definition = function(x, target, transpose = NULL, ...) {
    return(callNextMethod(
      x = x,
      target = target,
      transpose = transpose,
      tf = Matrix::t,
      ...
    ))
  }
)

setMethod(
  f = '.PrepLayerData',
  signature = c(x = 'NULL', target = 'ANY'),
  definition = function(x, target, transpose = NULL, ...) {
    return(x)
  }
)

#' @importFrom spam t
#' @importClassesFrom spam spam
#'
setMethod(
  f = '.PrepLayerData',
  signature = c(x = 'spam', target = 'NULL'),
  definition = function(x, target, transpose = NULL, ...) {
    return(callNextMethod(
      x = x,
      target = target,
      transpose = transpose,
      tf = spam::t,
      ...
    ))
  }
)
