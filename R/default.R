#' @include zzz.R
#' @include generics.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname dot-AssayClass
#' @method .AssayClass StdAssay
#' @export
#'
.AssayClass.StdAssay <- function(object) {
  cls <- gsub(
    pattern = '5$|v5$',
    replacement = '',
    x = class(x = object)[1L],
    ignore.case = TRUE
  )
  return(paste(cls, '(v5)'))
}

#' @method .MARGIN default
#' @export
#'
.MARGIN.default <- function(x, type = c('features', 'cells'), ...) {
  type <- arg_match(arg = type)
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
  if (is.null(x = col.name) && (is.atomic(x = metadata) && !is.matrix(x = metadata))) {
    abort(message = "'col.name' must be provided for atomic meta data")
  }
  if (inherits(x = metadata, what = c('matrix', 'Matrix'))) {
    metadata <- as.data.frame(x = metadata)
  }
  col.name <- col.name %||% names(x = metadata) %||% colnames(x = metadata)
  if (is.null(x = col.name)) {
    abort(message = "No metadata name provided and could not infer it from metadata object")
  }
  object[[col.name]] <- metadata
  return(object)
}

#' Internal Cropping Function
#'
#' @inheritParams Crop
#'
#' @return ...
#'
#' @keywords internal
#'
#' @noRd
#'
.Crop <- function(object, x = NULL, y = NULL, coords = c('plot','tissue'), ...) {
  if (is.null(x = x) && is.null(x = y)) {
    return(object)
  }
  compact <- .hasSlot(object = object, name = 'compact') && slot(object = object, name = 'compact')
  coords <- coords[1L]
  coords <- match.arg(arg = coords)
  switch(
    EXPR = coords,
    'plot' = {
      cx <- 'y'
      cy <- 'x'
    },
    'tissue' = {
      cx <- 'x'
      cy <- 'y'
    }
  )
  x <- range(x %||% bbox(obj = object)[cx, , drop = TRUE])
  y <- range(y %||% bbox(obj = object)[cy, , drop = TRUE])
  idx <- c(max = 1L, min = 2L)[[getOption(
    x = 'Seurat.coords.short_range',
    default = Seurat.options$Seurat.coords.short_range
  )]]
  if (x[1L] == x[2L]) {
    x[idx] <- bbox(obj = object)[cx, idx]
  }
  if (y[1L] == y[2L]) {
    y[idx] <- bbox(obj = object)[cy, idx]
  }
  args <- list(x, y)
  names(x = args) <- switch(
    EXPR = coords,
    'plot' = c('y', 'x'),
    'tissue' = c('x', 'y')
  )
  args <- args[c('x', 'y')]
  df <- do.call(what = expand.grid, args = args)
  df <- df[c(1, 3, 4, 2), ]
  df$cell <- 'cell'
  return(Overlay(x = object, y = CreateSegmentation(coords = df, compact = compact)))
}

#' Test Finiteness of Centroids
#'
#' Determines if a \code{\link{Centroids}} object should be finite; for
#' \code{Centroids}, this means if their \code{nsides} slot is an integer >= 3
#'
#' @param x A \code{\link{Centroids}} object
#'
#' @return \code{TRUE} if the \code{Centroids} are finite; otherwise
#' \code{FALSE}
#'
#' @keywords internal
#'
#' @noRd
#'
.FiniteCentroids <- function(x) {
  return(as.logical(x = length(x = x)))
}

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

.tail <- function(x, n = 10L, ...) {
  return(tail(x = x[[]], n = n, ...))
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

.OverBbox <- function(x, y, invert = FALSE, ...) {
  df <- .BboxDF(x = bbox(obj = y))
  df$cell <- 'cell'
  return(Overlay(
    x = x,
    y = CreateSegmentation(coords = df),
    invert = invert,
    ...
  ))
}

#' Internal Overlay Method
#'
#' @param x Query spatial object
#' @param y Target spatial object
#' @param ... Ignored
#'
#' @return \code{x} with only the components that fall within
#' the bounds of \code{y}
#'
#' @keywords internal
#'
#' @noRd
#'
.Overlay <- function(x, y, ...) {
  idx <- over(x = x, y = y)
  idx <- idx[!is.na(x = idx)]
  names(x = idx) <- vapply(
    X = strsplit(x = names(x = idx), split = '\\.'),
    FUN = '[[',
    FUN.VALUE = character(length = 1L),
    1L,
    USE.NAMES = FALSE
  )
  return(x[names(x = idx)])
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
