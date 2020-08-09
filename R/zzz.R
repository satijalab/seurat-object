#' @include utils.R
#' @importFrom methods setOldClass setClassUnion slot slot<-
#'
NULL

#' @docType package
#' @name SeuratObject-package
#' @rdname SeuratObject-package
#'
"_PACKAGE"

setOldClass(Classes = 'package_version')
setClassUnion(name = 'AnyMatrix', members = c("matrix", "dgCMatrix"))
setClassUnion(name = 'OptionalCharacter', members = c('NULL', 'character'))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Object Metadata
#'
#' Internal \code{\link{AddMetaData}} defintion
#'
#' @param object An object
#' @param metadata A vector, list, or data.frame with metadata to add
#' @param col.name A name for meta data if not a named list or data.frame
#'
#' @return object with metadata added
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

#' Check whether an assay has been processed by sctransform
#'
#' @param assay assay to check
#'
#' @return Boolean
#'
#' @keywords internal
#'
IsSCT <- function(assay) {
  if (is.list(x = assay)) {
    sct.check <- lapply(X = assay, FUN = function(x) {
      return(!is.null(x = Misc(object = x, slot = 'vst.out')) | !is.null(x = Misc(object = x, slot = 'vst.set')))
    })
    return(unlist(x = sct.check))
  }
  return(!is.null(x = Misc(object = assay, slot = 'vst.out')) | !is.null(x = Misc(object = assay, slot = 'vst.set')))
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
