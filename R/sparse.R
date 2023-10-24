#' @include zzz.R
#' @include generics.R
#' @importClassesFrom spam spam
#' @importClassesFrom Matrix CsparseMatrix RsparseMatrix
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Identify Sparse Slots
#'
#' @param x A sparse matrix
#' @param type ...
#'
#' @return ...
#'
#' @keywords internal
#'
#' @export
#'
#' @family sparse
#'
.SparseSlots <- function(x, type = c('pointers', 'indices', 'entries')) {
  UseMethod(generic = '.SparseSlots', object = x)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Is a Matrix Sparse
#'
#' @param x A matrix
#'
#' @return ...
#'
#' @keywords internal
#'
#' @export
#'
#' @family sparse
#'
#' @examples
#' IsSparse(matrix())
#' IsSparse(LayerData(pbmc_small, "counts"))
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
#' @keywords internal
#'
#' @export
#'
#' @family sparse
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

#' @rdname dot-SparseSlots
#' @method .SparseSlots CsparseMatrix
#' @export
#'
.SparseSlots.CsparseMatrix <- function(
  x,
  type = c('pointers', 'entries', 'indices')
) {
  type <- arg_match(arg = type)
  return(switch(
    EXPR = type,
    'pointers' = 'p',
    'indices' = 'i',
    'entries' = 'x'
  ))
}

#' @rdname dot-SparseSlots
#' @method .SparseSlots RsparseMatrix
#' @export
#'
.SparseSlots.RsparseMatrix <- function(
  x,
  type = c('pointers', 'indices', 'entries')
) {
  type <- arg_match(arg = type)
  return(switch(
    EXPR = type,
    'pointers' = 'p',
    'indices' = 'j',
    'entries' = 'x'
  ))
}

#' @rdname dot-SparseSlots
#' @method .SparseSlots spam
#' @export
#'
.SparseSlots.spam <- function(x, type = c('pointers', 'indices', 'entries')) {
  check_installed(pkg = 'spam')
  type <- arg_match(arg = type)
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
