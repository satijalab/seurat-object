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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
      fmargin <- fmargin %/% 1L
      if (!fmargin %in% seq.int(from = 1L, to = 2L)) {
        stop("'fmargin' must be either 1 or 2", call. = FALSE)
      }
      cmargin <- setdiff(x = seq.int(from = 1L, to = 2L), y = fmargin)
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
    fmargin <- fmargin %/% 1L
    if (!fmargin %in% seq.int(from = 1L, to = 2L)) {
      stop("'fmargin' must be either 1 or 2", call. = FALSE)
    }
    cmargin <- setdiff(x = seq.int(from = 1L, to = 2L), y = fmargin)
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
