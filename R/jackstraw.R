#' @include zzz.R
#' @include generics.R
#' @importFrom methods slot slot<- slotNames
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The JackStrawData Class
#'
#' The JackStrawData is used to store the results of a JackStraw computation.
#'
#' @slot empirical.p.values Empirical p-values
#' @slot fake.reduction.scores Fake reduction scores
#' @slot empirical.p.values.full Empirical p-values on full
#' @slot overall.p.values Overall p-values from ScoreJackStraw
#'
#' @name JackStrawData-class
#' @rdname JackStrawData-class
#' @exportClass JackStrawData
#'
JackStrawData <- setClass(
  Class = "JackStrawData",
  slots = list(
    empirical.p.values = "matrix",
    fake.reduction.scores = "matrix",
    empirical.p.values.full = "matrix",
    overall.p.values = "matrix"
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname JS
#' @export
#' @method JS JackStrawData
#'
JS.JackStrawData <- function(object, slot, ...) {
  CheckDots(...)
  slot <- switch(
    EXPR = slot,
    'empirical' = 'empirical.p.values',
    'fake' = 'fake.reduction.scores',
    'full' = 'empirical.p.values.full',
    'overall' = 'overall.p.values',
    slot
  )
  return(slot(object = object, name = slot))
}

#' @rdname JS
#' @export
#' @method JS<- JackStrawData
#'
"JS<-.JackStrawData" <- function(object, slot, ..., value) {
  CheckDots(...)
  slot <- switch(
    EXPR = slot,
    'empirical' = 'empirical.p.values',
    'fake' = 'fake.reduction.scores',
    'full' = 'empirical.p.values.full',
    'overall' = 'overall.p.values',
    slot
  )
  slot(object = object, name = slot) <- value
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' \code{JackStrawData} Methods
#'
#' Methods for \code{\link{JackStrawData}} objects for generics defined in
#' other packages
#'
#' @param x,object A \code{\link{JackStrawData}} object
#'
#' @name JackStrawData-methods
#' @rdname JackStrawData-methods
#'
#' @concept jackstraw
#'
NULL

#' @describeIn JackStrawData-methods Autocompletion for \code{$} access on a
#' \code{JackStrawData} object
#'
#' @inheritParams utils::.DollarNames
#'
#' @importFrom utils .DollarNames
#' @export
#' @method .DollarNames JackStrawData
#'
".DollarNames.JackStrawData" <- function(x, pattern = '') {
  slotnames <- as.list(x = slotNames(x = x))
  names(x = slotnames) <- unlist(x = slotnames)
  return(.DollarNames(x = slotnames, pattern = pattern))
}

#' @describeIn JackStrawData-methods Access data from a \code{JackStrawData}
#' object
#'
#' @param i A \code{JackStrawData} slot name
#'
#' @return \code{$}: Slot \code{i} from \code{x}
#' @export
#'
"$.JackStrawData" <- function(x, i, ...) {
  return(slot(object = x, name = i))
}

#' @describeIn JackStrawData-methods Have empirical p-values for a
#' \code{JackStrawData} object been calculated
#'
#' @return \code{as.logical}: \code{TRUE} if empirical p-values have been
#' calculated otherwise \code{FALSE}
#'
#' @export
#' @method as.logical JackStrawData
#'
as.logical.JackStrawData <- function(x, ...) {
  CheckDots(...)
  empP <- JS(object = x, slot = 'empirical')
  return(!(all(dim(x = empP) == 0) || all(is.na(x = empP))))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @describeIn JackStrawData-methods Overview of a \code{JackStrawData} object
#'
#' @return \code{show}: Prints summary to \code{\link[base]{stdout}} and
#' invisibly returns \code{NULL}
#'
#' @importFrom utils head
#' @importFrom methods show
#'
#' @export
#'
setMethod(
  f = 'show',
  signature = 'JackStrawData',
  definition = function(object) {
    empp <- object$empirical.p.values
    scored <- object$overall.p.values
    cat(
      "A JackStrawData object simulated on",
      nrow(x = empp),
      "features for",
      ncol(x = empp),
      "dimensions.\n",
      "Scored for:",
      nrow(x = scored),
      "dimensions.\n"
    )
    return(invisible(x = NULL))
  }
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
