#' @include generics.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname Cells
#' @export
#'
Cells.default <- function(x) {
  return(colnames(x = x))
}

#' @rdname IsGlobal
#' @export
#' @method IsGlobal default
#'
IsGlobal.default <- function(object, ...) {
  return(FALSE)
}
