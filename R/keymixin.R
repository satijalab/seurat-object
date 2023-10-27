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
#' @template slot-key
#'
#' @keywords internal
#'
#' @exportClass KeyMixin
#'
#' @aliases KeyMixin
#'
#' @family key
#'
setClass(
  Class = 'KeyMixin',
  contains = 'VIRTUAL',
  slots = list(key = 'character')
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Regex Pattern for Keys
#'
#' @return Returns the regex pattern for keys
#' (\dQuote{\Sexpr[stage=build]{SeuratObject:::.KeyPattern()}})
#'
#' @keywords internal
#'
#' @export
#'
#' @family key
#'
.KeyPattern <- function() {
  return('^[a-zA-Z][a-zA-Z0-9]*_$')
}

#' Generate a Random Key
#'
#' @inheritParams RandomName
#'
#' @return Returns a valid key
#'
#' @keywords internal
#'
#' @export
#'
#' @family key
#'
#' @examples
#' set.seed(42L)
#' .RandomKey()
#'
.RandomKey <- function(length = 7L, ...) {
  return(Key(
    object = RandomName(
      length = length,
      chars = c(letters, LETTERS, seq.int(from = 0L, to = 9L)),
      ...
    ),
    quiet = TRUE
  ))
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
#' @return \code{Key.KeyMixin}: The key from \code{object}; if no key set,
#' returns \code{NULL}
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

#' @method Key NULL
#' @export
#'
Key.NULL <- function(object, ...) {
    return(NULL)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Internal Key Methods
#'
#' Internal key methods for classes that inherit from \code{\link{KeyMixin}};
#' these functions are designed to be used as the body for methods for
#' \code{Key()} and \code{Key<-()} when an immediate public method is required.
#' Generally speaking, classes that inherit from \code{KeyMixin} should use the
#' \code{KeyMixin} methods for \code{Key()} and \code{Key<-()}
#'
#' @inheritParams Key
#'
#' @inherit Key return
#'
#' @keywords internal
#'
#' @noRd
#'
.Key <- function(object, ...) {
  CheckDots(...)
  return(NextMethod())
}

#' @rdname dot-Key
#'
#' @noRd
#'
".Key<-" <- function(object, ..., value) {
  CheckDots(...)
  object <- UpdateSlots(object = object)
  object <- NextMethod()
  return(object)
}

#' Update a Key
#'
#' @param key A character to become a Seurat Key
#'
#' @return An updated Key that's valid for Seurat
#'
#' @keywords internal
#'
#' @family key
#'
#' @noRd
#'
UpdateKey <- function(key) {
  key.msg <- 'Keys should be one or more alphanumeric characters followed by an underscore'
  if (isTRUE(x = grepl(pattern = .KeyPattern(), x = key))) {
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
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Key Validity
#'
#' Validation of \code{\link{KeyMixin}} objects is handled by
#' \code{\link[methods]{validObject}}
#'
#' @section Key Validation:
#' Keys must be a one-length character vector; a key must be composed of one
#' of the following:
#' \itemize{
#'  \item An empty string (eg. \dQuote{\code{''}}) where \code{nzchar() == 0}
#'  \item An string composed of one or more alphanumeric values
#'  (both lower- and upper-case) that ends with an underscore
#'  (\dQuote{\code{_}}); the first character must be a letter
#' }
#' Keys that are not empty strings are validated with the regex
#' \dQuote{\code{\Sexpr[stage=build]{SeuratObject:::.KeyPattern()}}}
#'
#' @importFrom rlang is_scalar_character
#'
#' @keywords internal
#'
#' @name Key-validity
#'
#' @family key
#'
setValidity(
  Class = 'KeyMixin',
  method = function(object) {
    if (isFALSE(x = getOption(x = "Seurat.object.validate", default = TRUE))) {
      warn(
        message = paste("Not validating", class(x = object)[1L], "objects"),
        class = 'validationWarning'
      )
      return(TRUE)
    }
    valid <- NULL
    key <- Key(object = object)
    # Ensure key has length of 1
    if (!is.null(x = key) && .GetSeuratCompat() >= '5.0.0') {
      if (!is_scalar_character(x = key)) {
        valid <- c(valid, "Keys must be a one-length character vector")
      } else if (is_na(x = key)) {
        valid <- c(valid, "Keys may not be 'NA'")
      } else if (nzchar(x = key) && !grepl(pattern = .KeyPattern(), x = key)) {
        # Ensure proper key composition
        valid <- c(
          valid,
          paste0("Keys must match the pattern '", .KeyPattern(), "'")
        )
      }
    }
    return(valid %||% TRUE)
  }
)

.CheckKey <- function(key, existing = NULL, name = NULL) {
  if (rlang::is_missing(x = key) || !length(x = key) || !nzchar(x = key)) {
    key <- Key(object = tolower(name) %||% RandomName(), quiet = TRUE)
  }
  if (!is.null(x = names(x = existing)) && !is.null(x = name)) {
    existing <- existing[setdiff(x = names(x = existing), y = name)]
  }
  if (key %in% existing) {
    old <- key
    key <- Key(object = tolower(x = name %||% RandomName()), quiet = TRUE)
    i <- 1L
    n <- 5L
    while (key %in% existing) {
      key <- Key(object = RandomName(length = n), quiet = TRUE)
      i <- i + 1L
      if (!i %% 7L) {
        n <- n + 2L
      }
    }
    warn(
      message = paste(
        "Key",
        sQuote(x = old),
        "taken, using",
        sQuote(x = key),
        "instead"
      )
    )
  }
  return(key)
}
