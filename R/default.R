#' @include generics.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @method .AssayClass StdAssay
#' @export
#'
.AssayClass.StdAssay <- function(object) {
  return('Assay (v5)')
}

#' @method .MARGIN default
#' @export
#'
.MARGIN.default <- function(object, type = c('features', 'cells')) {
  type <- type[1]
  type <- match.arg(arg = type)
  return(unname(obj = c(features = 1L, cells = 2L)[type]))
}

#' @rdname Cells
#' @method Cells default
#' @export
#'
Cells.default <- function(x, ...) {
  return(colnames(x = x))
}

#' @rdname IsGlobal
#' @method IsGlobal default
#' @export
#'
IsGlobal.default <- function(object, ...) {
  return(FALSE)
}

#' @method Key character
#' @export
#'
Key.character <- function(object, quiet = FALSE, ...) {
  f <- ifelse(test = isTRUE(x = quiet), yes = suppressWarnings, no = identity)
  return(f(UpdateKey(key = object)))
}
