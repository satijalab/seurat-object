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

#' @importFrom stats na.omit
#'
#' @rdname MatchCells
#' @method MatchCells character
#' @export
#'
MatchCells.character <- function(new, orig, ordered = FALSE) {
  cmatch <- as.vector(x = na.omit(object = match(x = orig, table = new)))
  if (!isTRUE(x = ordered)) {
    cmatch <- sort(x = cmatch)
  }
  # return(as.vector(x = na.omit(object = )))
  return(cmatch)
}

#' @rdname MatchCells
#' @method MatchCells numeric
#' @export
#'
MatchCells.numeric <- function(new, orig, ordered = FALSE) {
  new <- unique(x = new[new <= length(x = orig)])
  if (isTRUE(x = ordered)) {
    new <- sort(x = new)
  }
  return(new)
}
