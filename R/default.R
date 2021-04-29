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
