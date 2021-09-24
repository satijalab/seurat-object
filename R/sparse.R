#' @include zzz.R
#' @include generics.R
#' @importClassesFrom spam spam
#' @importClassesFrom Matrix CsparseMatrix RsparseMatrix
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @export
#'
.SparseSlots <- function(x, type = c('pointers', 'indices', 'entries')) {
  UseMethod(generic = '.SparseSlots', object = x)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @export
#'
IsSparse <- function(x) {
  if (!isS4(x)) {
    return(FALSE)
  }
  classkey <- unlist(x = strsplit(
    x = ClassKey(class = class(x = x)),
    split = ':'
  ))
  cls <- classkey[[2L]]
  pkg <- classkey[[1L]]
  sparse <- cls %in% sparse.classes[[pkg]]
  if (!sparse) {
    sparse <- any(sparse.classes[[pkg]] %in% .Contains(object = x))
  }
  return(sparse)
}

#' Register Sparse Matrix Classes
#'
#' @inheritParams ClassKey
#'
#' @return Invisibly returns \code{NULL}
#'
#' @export
#'
RegisterSparseMatrix <- function(class, package = NULL) {
  classkey <- unlist(x = strsplit(
    x = ClassKey(class = class, package = package),
    split = ':'
  ))
  sparse.classes[[classkey[[1L]]]] <- unique(c(
    sparse.classes[[classkey[[1L]]]],
    classkey[[2L]]
  ))
  return(invisible(x = NULL))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @method .MARGIN CsparseMatrix
#' @export
#'
.MARGIN.CsparseMatrix <- function(x, ...) {
  return(2L)
}

#' @method .MARGIN RsparseMatrix
#' @export
#'
.MARGIN.RsparseMatrix <- function(x, ...) {
  return(1L)
}

#' @method .MARGIN spam
#' @export
#'
.MARGIN.spam <- .MARGIN.RsparseMatrix

#' @method .SparseSlots CsparseMatrix
#' @export
#'
.SparseSlots.CsparseMatrix <- function(
  x,
  type = c('pointers', 'entries', 'indices')
) {
  type <- type[1L]
  type <- match.arg(arg = type)
  return(switch(
    EXPR = type,
    'pointers' = 'p',
    'indices' = 'i',
    'entries' = 'x'
  ))
}

#' @method .SparseSlots RsparseMatrix
#' @export
#'
.SparseSlots.RsparseMatrix <- function(
  x,
  type = c('pointers', 'indices', 'entries')
) {
  type <- type[1L]
  type <- match.arg(arg = type)
  return(switch(
    EXPR = type,
    'pointers' = 'p',
    'indices' = 'j',
    'entries' = 'x'
  ))
}

#' @method .SparseSlots spam
#' @export
#'
.SparseSlots.spam <- function(x, type = c('pointers', 'indices', 'entries')) {
  type <- type[1L]
  type <- match.arg(arg = type)
  return(switch(
    EXPR = type,
    'pointers' = 'rowpointers',
    'indices' = 'colindices',
    'entries' = 'entries'
  ))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 Methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
