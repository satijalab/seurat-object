#' @include zzz.R
#' @include generics.R
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

#' Generate a Random Key
#'
#' @inheritParams RandomName
#'
#' @return ...
#'
#' @examples
#' set.seed(42L)
#' RandomKey()
#'
#' @export
#'
RandomKey <- function(length = 7L, ...) {
  return(Key(
    object = RandomName(
      length = length,
      chars = c(letters, LETTERS, seq.int(from = 0L, to = 9L)),
      ...
    ),
    quiet = TRUE
  ))
}

#' Update a Key
#'
#' @param key A character to become a Seurat Key
#'
#' @return An updated Key that's valid for Seurat
#'
#' @keywords internal
#'
#' @noRd
#'
UpdateKey <- function(key) {
  key.msg <- 'Keys should be one or more alphanumeric characters followed by an underscore'
  if (isTRUE(x = grepl(pattern = '^[[:alnum:]]+_$', x = key))) {
    return(key)
  }
  new.key <- regmatches(
    x = key,
    m = gregexpr(pattern = '[[:alnum:]]+', text = key)
  )
  new.key <- paste0(paste(unlist(x = new.key), collapse = ''), '_')
  if (new.key == '_') {
    new.key <- paste0(RandomName(length = 3), '_')
  }
  warning(
    key.msg,
    ", setting key from ",
    key,
    " to ",
    new.key,
    call. = FALSE,
    immediate. = TRUE
  )
  return(new.key)
}

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
  key <- slot(object = object, name = 'key')
  if (!length(x = key)) {
    key <- NULL
  }
  return(key)
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
    if (!is.null(x = key)) {
      if (length(x = key) != 1L) {
        valid <- c(valid, "Keys must be a one-length character vector")
      } else if (nzchar(x = key) && !grepl(pattern = pattern, x = key)) {
        # Ensure proper key composition
        valid <- c(valid, paste0("Keys must match the pattern '", pattern, "'"))
      }
    }
    return(valid %||% TRUE)
  }
)
