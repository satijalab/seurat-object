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
    attr(x = x, which = 'features') <- dnames[[fmargin]]
    attr(x = x, which = 'cells') <- dnames[[cmargin]]
    x <- suppressMessages(expr = unname(x))
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
