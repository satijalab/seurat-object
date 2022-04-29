#' @include generics.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname Cells
#' @export
#'
Cells.default <- function(x, ...) {
  return(colnames(x = x))
}

#' @rdname IsGlobal
#' @export
#' @method IsGlobal default
#'
IsGlobal.default <- function(object, ...) {
  return(FALSE)
}

#' @importFrom stats na.omit
#'
#' @rdname MatchCells
#' @method MatchCells character
#' @export
#'
MatchCells.character <- function(new, orig, ordered = FALSE) {
  cmatch <- as.vector(x = na.omit(object = match(x = orig, table = new)))
  if (!length(x = cmatch)) {
    return(NULL)
  }
  if (!isTRUE(x = ordered)) {
    cmatch <- sort(x = cmatch)
  }
  return(cmatch)
}

#' @rdname MatchCells
#' @method MatchCells NULL
#' @export
#'
MatchCells.NULL <- function(new, orig, ordered = FALSE) {
  return(NULL)
}

#' @rdname MatchCells
#' @method MatchCells numeric
#' @export
#'
MatchCells.numeric <- function(new, orig, ordered = FALSE) {
  new <- unique(x = new[new <= length(x = orig)])
  if (!length(x = new)) {
    return(NULL)
  }
  if (isTRUE(x = ordered)) {
    new <- sort(x = new)
  }
  return(new)
}
