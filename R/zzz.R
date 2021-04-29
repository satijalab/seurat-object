#' @include utils.R
#' @importFrom methods new setClass setClassUnion setGeneric setOldClass
#' setValidity slot slot<-
#' @importClassesFrom Matrix dgCMatrix
#'
NULL

#' @docType package
#' @name SeuratObject-package
#' @rdname SeuratObject-package
#'
"_PACKAGE"

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Reexports
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom Matrix colMeans
#' @export
#'
Matrix::colMeans

#' @importFrom Matrix colSums
#' @export
#'
Matrix::colSums

#' @importFrom Matrix rowMeans
#' @export
#'
Matrix::rowMeans

#' @importFrom Matrix rowSums
#' @export
#'
Matrix::rowSums

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setClassUnion(name = 'AnyMatrix', members = c("matrix", "dgCMatrix"))
setClassUnion(name = 'OptionalCharacter', members = c('NULL', 'character'))
setClassUnion(name = 'OptionalList', members = c('NULL', 'list'))

setOldClass(Classes = 'package_version')

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

#' Check List Names
#'
#' Check to see if a list has names; also check to enforce that all names are
#' present and unique
#'
#' @param x A list
#' @param all.unique Require that all names are unique from one another
#' @param allow.empty Allow empty (\code{nchar = 0}) names
#'
#' @return \code{TRUE} if ..., otherwise \code{FALSE}
#'
#' @importFrom rlang is_bare_list
#'
#' @keywords internal
#'
#' @noRd
#'
IsNamedList <- function(x, all.unique = TRUE, allow.empty = FALSE) {
  if (!is_bare_list(x = x)) {
    stop("'x' must be a list")
  }
  n <- names(x = x)
  named <- !is.null(x = n)
  if (!isTRUE(x = allow.empty)) {
    named <- named && all(vapply(
      X = n,
      FUN = nchar,
      FUN.VALUE = integer(length = 1L)
    ))
  }
  if (isTRUE(x = all.unique)) {
    named <- named && (length(x = n) == length(x = unique(x = n)))
  }
  return(named)
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
