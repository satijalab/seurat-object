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

#' @importFrom utils .DollarNames
#' @export
#' @method .DollarNames JackStrawData
#'
".DollarNames.JackStrawData" <- function(x, pattern = '') {
  slotnames <- as.list(x = slotNames(x = x))
  names(x = slotnames) <- unlist(x = slotnames)
  return(.DollarNames(x = slotnames, pattern = pattern))
}

#' @export
#'
"$.JackStrawData" <- function(x, i, ...) {
  return(slot(object = x, name = i))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setMethod(
  f = 'show',
  signature = 'JackStrawData',
  definition = function(object) {
    # empp <- GetJS(object = object, slot = "empirical.p.values")
    empp <- object$empirical.p.values
    # scored <- GetJS(object = object, slot = "overall.p.values")
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
  }
)
