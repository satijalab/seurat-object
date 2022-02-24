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
.MARGIN.default <- function(object, type = c('features', 'cells'), ...) {
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Head and Tail Object Metadata
#'
#' Internal \code{\link[utils]{head}} and \code{\link[utils]{tail}} definitions
#'
#' @param x An object
#' @param n Number of rows to return
#' @inheritDotParams utils::head
#'
#' @return The first or last \code{n} rows of object metadata
#'
#' @keywords internal
#'
#' @noRd
#'
.head <- function(x, n = 10L, ...) {
  return(head(x = x[[]], n = n, ...))
}

#' Miscellaneous Data
#'
#' Internal functions for getting and setting miscellaneous data
#'
#' @param object An object
#' @param slot Name of miscellaneous data to get or set
#' @param ... Arguments passed to other methods
#'
#' @return \code{.Misc}: If \code{slot} is \code{NULL}, all miscellaneous
#' data, otherwise the miscellaneous data for \code{slot}
#'
#' @keywords internal
#'
#' @noRd
#'
.Misc <- function(object, slot = NULL, ...) {
  CheckDots(...)
  if (is.null(x = slot)) {
    return(slot(object = object, name = 'misc'))
  }
  return(slot(object = object, name = 'misc')[[slot]])
}

#' @param value Data to add
#'
#' @return \code{.Misc<-}: \code{object} with \code{value} added to the
#' miscellaneous data slot \code{slot}
#'
#' @rdname dot-Misc
#'
#' @noRd
#'
".Misc<-" <- function(object, slot, ..., value) {
  CheckDots(...)
  if (slot %in% names(x = Misc(object = object))) {
    warning(
      "Overwriting miscellanous data for ",
      slot,
      call. = FALSE,
      immediate. = TRUE
    )
  }
  if (is.list(x = value)) {
    slot(object = object, name = 'misc')[[slot]] <- c(value)
  } else {
    slot(object = object, name = 'misc')[[slot]] <- value
  }
  return(object)
}

#' @rdname dot-head
#'
#' @noRd
#'
.tail <- function(x, n = 10L, ...) {
  return(tail(x = x[[]], n = n, ...))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Add Object Metadata
#'
#' Internal \code{\link{AddMetaData}} definition
#'
#' @param object An object
#' @param metadata A vector, list, or data.frame with metadata to add
#' @param col.name A name for meta data if not a named list or data.frame
#'
#' @return object with metadata added
#'
#' @keywords internal
#'
#' @noRd
#'
.AddMetaData <- function(object, metadata, col.name = NULL) {
  if (is.null(x = col.name) && is.atomic(x = metadata)) {
    stop("'col.name' must be provided for atomic metadata types (eg. vectors)")
  }
  if (inherits(x = metadata, what = c('matrix', 'Matrix'))) {
    metadata <- as.data.frame(x = metadata)
  }
  col.name <- col.name %||% names(x = metadata) %||% colnames(x = metadata)
  if (is.null(x = col.name)) {
    stop("No metadata name provided and could not infer it from metadata object")
  }
  object[[col.name]] <- metadata
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
