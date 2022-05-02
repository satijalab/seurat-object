#' @include zzz.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' A Mixin for Keyed objects
#'
#' A mixin (virtual class) for enabling keyed objects; provides consistent
#' behavior for getting, setting, and validating keys
#'
#' @slot key A one-length character vector with the object's key; keys must
#' be one or more alphanumeric characters followed by an underscore
#' \dQuote{\code{_}} (regex pattern \dQuote{\code{^[[:alnum:]]+_$}})
#'
#' @keywords internal
#'
#' @exportClass KeyMixin
#'
#' @aliases KeyMixin
#'
setClass(
  Class = 'KeyMixin',
  contains = 'VIRTUAL',
  slots = list(key = 'character')
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param object An object
#' @param quiet Suppress warnings when updating characters to keys
#' @param ... Ignored
#' @param value A key to set
#'
#' @details \code{Key.character}: Update a character to a key
#'
#' @return \code{Key.character}: \code{object} but as a syntactically-valid key
#'
#' @rdname KeyMixin-class
#' @method Key character
#' @export
#'
Key.character <- function(object, quiet = FALSE, ...) {
  f <- ifelse(test = isTRUE(x = quiet), yes = suppressWarnings, no = identity)
  return(f(UpdateKey(key = object)))
}

#' @details \code{Key.KeyMixin}: Get the key of a keyed object
#'
#' @return \code{Key.KeyMixin}: The key from \code{object}
#'
#' @rdname KeyMixin-class
#' @method Key KeyMixin
#' @export
#'
Key.KeyMixin <- function(object, ...) {
  return(slot(object = object, name = 'key'))
}

#' @details \code{Key<-}: Set the key of a keyed object
#'
#' @return \code{Key<-}: \code{object} with the key set to \code{value}
#'
#' @rdname KeyMixin-class
#' @method Key<- KeyMixin
#' @export
#'
"Key<-.KeyMixin" <- function(object, ..., value) {
  slot(object = object, name = 'key') <- Key(object = value, ...)
  validObject(object = object)
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setValidity(
  Class = 'KeyMixin',
  method = function(object) {
    valid <- NULL
    pattern <- '^[[:alnum:]]+_$'
    key <- Key(object = object)
    # Ensure key has length of 1
    if (length(x = key) != 1L) {
      valid <- c(valid, "Keys must be a one-length character vector")
    } else if (nchar(x = key) && !grepl(pattern = pattern, x = key)) {
      # Ensure proper key composition
      valid <- c(valid, paste0("Keys must match the pattern '", pattern, "'"))
    }
    return(valid %||% TRUE)
  }
)
