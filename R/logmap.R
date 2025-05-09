#' @include zzz.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' A Logical Map
#'
#' A simple container for storing mappings of values using logical matrices.
#' Keeps track of which values (rows) are present in which observations
#' (columns). \code{LogMap} objects can be created with \code{LogMap()};
#' queries can be performed with \code{[[} and observations can be added
#' or removed with \code{[[<-}
#'
#' @slot .Data A logical matrix with at least one row
#'
#' @exportClass LogMap
#'
#' @family logmap
#'
setClass(
  Class = 'LogMap',
  contains = 'matrix'
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname LogMap-class
#'
#' @param y A character vector
#'
#' @return \code{LogMap}: A new \code{LogMap} object with zero columns and
#' \code{length(x = x)} rows; rownames are set to \code{x}
#'
#' @export
#'
#' @order 1
#'
#' @examples
#' # Create a LogMap
#' map <- LogMap(letters[1:10])
#' map
#'
#' # Get the names of values in the LogMap
#' map[[NULL]]
#' rownames(map)
#'
#' # Add an observation to the LogMap
#' map[['obs']] <- c(1, 3, 7)
#' map[['entry']] <- c(2, 7, 10)
#' map
#'
#' # Get the names of observations in the LogMap
#' colnames(map)
#'
#' # Fetch an observation from the LogMap
#' map[['obs']]
#'
#' # Get the full logical matrix
#' map[[]]
#'
#' # Remove an observation from the LogMap
#' map[['obs']] <- NULL
#' map[['entry']] <- NULL
#' map
#'
LogMap <- function(y) {
  if (!is.character(x = y)) {
    y <- as.character(x = y)
  }
  return(new(
    Class = 'LogMap',
    .Data = matrix(nrow = length(x = y), ncol = 0, dimnames = list(y, NULL))
  ))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Coerce Logical Maps to Matrices
#'
#' Coerce a logical map to a matrix; this removes all
#' \link[=LogMap]{logical map} class capabilities from
#' the object and returns a base-R matrix object
#'
#' @param x A \code{\link{LogMap}} object
#' @template param-dots-ignored
#'
#' @return A base-R matrix created from \code{x}
#'
#' @method as.matrix LogMap
#' @export
#'
#' @family logmap
#'
#' @examples
#' map <- LogMap(letters[1:10])
#' map[['obs']] <- c(1, 3, 7)
#' mat <- as.matrix(map)
#' mat
#' class(mat)
#'
as.matrix.LogMap <- function(x, ...) {
  return(as(object = x, Class = 'matrix'))
}

#' Drop Unused Logical Map Values
#'
#' Remove any unused values from a \link[=LogMap]{logical map}
#'
#' @template param-dots-ignored
#' @param x A \code{LogMap} object
#'
#' @return \code{x} with values not present in any
#' observation removed
#'
#' @method droplevels LogMap
#' @export
#'
#' @family logmap
#'
#' @examples
#' map <- LogMap(letters[1:10])
#' map[['obs']] <- c(1, 3, 7)
#' map[['entry']] <- c(2, 7, 10)
#'
#' # Remove unused values
#' map <- droplevels(map)
#' map
#' map[[]]
#'
droplevels.LogMap <- function(x, ...) {
  fidx <- which(x = apply(
    X = x,
    MARGIN = 1L,
    FUN = function(row) {
      return(all(vapply(
        X = row,
        FUN = isFALSE,
        FUN.VALUE = logical(length = 1L)
      )))
    }
  ))
  if (length(x = fidx)) {
    x <- as(object = x[-fidx, , drop = FALSE], Class = 'LogMap')
  }
  validObject(object = x)
  return(x)
}

#' Find Common Logical Map Values
#'
#' Identify values in a \link[=LogMap]{logical map} that are
#' common to every observation
#'
#' @inheritParams droplevels.LogMap
#' @param y Ignored
#'
#' @return The values of \code{x} that are present in \strong{every} observation
#'
#' @method intersect LogMap
#' @export
#'
#' @family logmap
#'
#' @examples
#' map <- LogMap(letters[1:10])
#' map[['obs']] <- c(1, 3, 7)
#' map[['entry']] <- c(2, 7, 10)
#'
#' # Identify values that are present in every observation
#' intersect(map)
#'
intersect.LogMap <- function(x, y = missing_arg(), ...) {
  if (!is_missing(x = y)) {
    abort(message = "'y' must not be provided")
  }
  idx <- which(x = apply(X = x, MARGIN = 1L, FUN = all))
  return(rownames(x = x)[idx])
}

#' Find Observations by Value
#'
#' Identify the observations that contain a specific
#' value in a \link[=LogMap]{logical map}
#'
#' @template param-dots-ignored
#' @param object A \code{\link{LogMap}} object
#' @param values A vector of values to find observations for
#' @param select Observation selection method; choose from:
#' \itemize{
#'  \item \dQuote{\code{first}}: the first observation the value is found in
#'  \item \dQuote{\code{last}}: the last observation the value is found in
#'  \item \dQuote{\code{common}}: the first most-common observation the value
#'  is found in; most-common is determined by the observation that contains
#'  the most of the values requested
#'  \item \dQuote{\code{all}}: all observations the value is found in
#' }
#' @param simplify Simplify the resulting list to a vector
#'
#' @return \code{labels}: A list, or vector if \code{simplify} is \code{TRUE},
#' of all values and the observations they're found in, according
#' to the value of \code{select}
#'
#' @method labels LogMap
#' @export
#'
#' @family logmap
#'
#' @examples
#' map <- LogMap(letters[1:10])
#' map[['obs']] <- c(1, 3, 7)
#' map[['entry']] <- c(2, 7, 10)
#'
#' # Find observations for a set of values
#' labels(map, c('a', 'b', 'g'))
#'
labels.LogMap <- function(
  object,
  values,
  select = c('first', 'last', 'common', 'all'),
  simplify = TRUE,
  ...
) {
  select <- select[1L]
  select <- match.arg(arg = select)
  values <- intersect(x = values, y = rownames(x = object))
  mat <- as.matrix(object)
  cols <- colnames(object)
  idx <- match(values, rownames(object))
  obs <- vector("list", length(values))
  names(obs) <- values
  for (i in seq_along(idx)) {
    id <- idx[i]
    vals <- cols[mat[id, , drop = FALSE]]
    if (length(vals) > 0) {
      obs[[i]] <- vals
    }
  }
  obs <- Filter(f = length, x = obs)
  obs <- switch(
    EXPR = select,
    'first' = lapply(X = obs, FUN = '[[', 1L),
    'last' = lapply(
      X = obs,
      FUN = function(x) {
        return(x[[length(x = x)]])
      }
    ),
    common = {
      counts <- table(unlist(x = obs))
      tmp <- obs
      obs <- vector(mode = 'character', length = length(x = tmp))
      names(x = obs) <- names(x = tmp)
      for (i in seq_along(along.with = obs)) {
        obs[i] <- names(x = which.max(
          x = counts[names(x = counts) %in% tmp[[i]]]
        ))
      }
      obs
    },
    obs
  )
  if (isTRUE(x = simplify)) {
    tmp <- obs
    obs <- unlist(x = tmp)
    names(x = obs) <- make.unique(names = rep.int(
      x = names(x = tmp),
      times = vapply(X = tmp, FUN = length, FUN.VALUE = numeric(length = 1L))
    ))
  }
  return(obs)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Matrix-like Subsetting for \link[=LogMap]{Logical Maps}
#'
#' @inheritParams base::`[`
#' @inheritParams LogMap-class
#' @param i,j Vectors of values (\code{i}) and observations (\code{j}) to pull
#' from \code{x}
#' @template param-dots-method
#'
#' @note \code{i} is not reordable; passing a different order for \code{i}
#' will return a subset with rows in the same order as \code{x}
#'
#' @name [,LogMap
#' @rdname sub-LogMap-method
#'
#' @keywords internal
#'
#' @examples
#' map <- LogMap(letters[1:10])
#' map[['obs']] <- c(1, 3, 7)
#' map[['entry']] <- c(2, 7, 10)
#'
#' map[]
#' map[1:5, 2L]
#' map[c("b", "c", "f"), "obs"]
#'
#' # Pass `drop = TRUE` to cast to `matrix`
#' map[1:3, , drop = TRUE]
#'
#' # Note that `i` is non-reordable
#' rownames(map)[1:3]
#' map[c("b", "c", "a"), , drop = TRUE]
#'
NULL

#' @rdname sub-LogMap-method
#'
setMethod(
  f = '[',
  signature = c(x = 'LogMap', i = 'missing', j = 'missing'),
  definition = function(x, i, j, ..., drop = FALSE) {
    return(x)
  }
)

#' @rdname sub-LogMap-method
#'
setMethod(
  f = '[',
  signature = c(x = 'LogMap', i = 'character', j = 'character'),
  definition = function(x, i, j, ..., drop = FALSE) {
    i <- i[MatchCells(new = i, orig = rownames(x = x), ordered = TRUE)]
    x <- as.matrix(x = x)[i, j, drop = drop]
    if (!isTRUE(x = drop)) {
      x <- as(object = x, Class = 'LogMap')
    }
    return(x)
  }
)

#' @rdname sub-LogMap-method
#'
setMethod(
  f = '[',
  signature = c(x = 'LogMap', i = 'character', j = 'missing'),
  definition = function(x, i, j, ..., drop = FALSE) {
    i <- i[MatchCells(new = i, orig = rownames(x = x), ordered = TRUE)]
    x <- as.matrix(x = x)[i, , drop = drop]
    if (!isTRUE(x = drop)) {
      x <- as(object = x, Class = 'LogMap')
    }
    return(x)
  }
)

#' @rdname sub-LogMap-method
#'
setMethod(
  f = '[',
  signature = c(x = 'LogMap', i = 'missing', j = 'character'),
  definition = function(x, i, j, ..., drop = FALSE) {
    x <- as.matrix(x = x)[, j, drop = drop]
    if (!isTRUE(x = drop)) {
      x <- as(object = x, Class = 'LogMap')
    }
    return(x)
  }
)

#' @rdname sub-LogMap-method
#'
setMethod(
  f = '[',
  signature = c(x = 'LogMap', i = 'numeric', j = 'missing'),
  definition = function(x, i, j, ..., drop = FALSE) {
    stopifnot(is_bare_integerish(x = i))
    if (is.unsorted(x = i)) {
      i <- sort(x = i)
    }
    i <- rownames(x = x)[i]
    return(callNextMethod(x, i, ..., drop = drop))
  }
)

#' @rdname sub-LogMap-method
#'
setMethod(
  f = '[',
  signature = c(x = 'LogMap', i = 'missing', j = 'numeric'),
  definition = function(x, i, j, ..., drop = FALSE) {
    stopifnot(is_bare_integerish(x = j))
    j <- colnames(x = x)[j]
    return(callNextMethod(x, , j, ..., drop = drop))
  }
)

#' @rdname sub-LogMap-method
#'
setMethod(
  f = '[',
  signature = c(x = 'LogMap', i = 'numeric', j = 'numeric'),
  definition = function(x, i, j, ..., drop = FALSE) {
    stopifnot(is_bare_integerish(x = i), is_bare_integerish(x = j))
    if (is.unsorted(x = i)) {
      i <- sort(x = i)
    }
    i <- rownames(x = x)[i]
    j <- colnames(x = x)[j]
    return(callNextMethod(x, i, j, ..., drop = drop))
  }
)

#' @rdname LogMap-class
#'
#' @param x A \code{LogMap} object
#' @param i A character vector of length 1, or \code{NULL}
#' @param j Not used
#' @param ... Ignored
#'
#' @return \code{[[}: if \code{i} is a character vector, the rownames that are
#' mapped to \code{i}; otherwise the rownames of \code{x}
#'
#' @export
#'
#' @order 2
#'
setMethod(
  f = '[[',
  signature = c(x = 'LogMap', i = 'character', j = 'missing'),
  definition = function(x, i, ...) {
    i <- i[1]
    i <- match.arg(arg = i, choices = colnames(x = x))
    return(rownames(x = x)[x[, i, drop = TRUE]])
  }
)

#' \code{\link{LogMap}} Interaction Methods
#'
#' Additional methods for using \code{[[} with \code{\link{LogMap}} objects
#'
#' @inheritParams LogMap
#' @param i An integer or numeric vector of length 1
#'
#' @return The rownames that are mapped to \code{i}
#'
#' @rdname sub-sub-LogMap-internal-method
#'
#' @keywords internal
#'
setMethod(
  f = '[[',
  signature = c(x = 'LogMap', i = 'integer', j = 'missing'),
  definition = function(x, i, ...) {
    return(x[[colnames(x = x)[i]]])
  }
)

#' @rdname LogMap-class
#'
#' @export
#'
#' @order 3
#'
setMethod(
  f = '[[',
  signature = c(x = 'LogMap', i = 'missing', j = 'missing'),
  definition = function(x, ...) {
    return(slot(object = x, name = '.Data'))
  }
)

#' @rdname sub-sub-LogMap-internal-method
#'
setMethod(
  f = '[[',
  signature = c(x = 'LogMap', i = 'numeric', j = 'missing'),
  definition = function(x, i, ...) {
    return(x[[as.integer(x = i)]])
  }
)

#' @rdname LogMap-class
#'
#' @export
#'
#' @order 4
#'
setMethod(
  f = '[[',
  signature = c(x = 'LogMap', i = 'NULL', j = 'missing'),
  definition = function(x, i, ...) {
    return(rownames(x = x))
  }
)

#' @rdname LogMap-class
#'
#' @param value A character or integer vector of values to record in the map
#' for \code{i}, or \code{NULL} to remove the record for \code{i}
#'
#' @return \code{[[<-}: If \code{value} is \code{NULL}, then \code{x} without
#' the observations for \code{i}; otherwise, \code{x} with a new column for
#' \code{i} recording a \code{TRUE} for all values present in \code{value}
#'
#' @export
#'
#' @order 5
#'
setMethod(
  f = '[[<-',
  signature = c(
    x = 'LogMap',
    i = 'character',
    j = 'missing',
    value = 'character'
  ),
  definition = function(x, i, ..., value) {
    value <- MatchCells(new = rownames(x = x), orig = value)
    x[[i]] <- value
    return(x)
  }
)

#' @rdname LogMap-class
#'
#' @export
#'
#' @order 6
#'
setMethod(
  f = '[[<-',
  signature = c(
    x = 'LogMap',
    i = 'character',
    j = 'missing',
    value = 'integer'
  ),
  definition = function(x, i, ..., value) {
    if (i %in% colnames(x = x)) {
      x[[i]] <- NULL
    }
    value <- MatchCells(new = value, rownames(x = x))
    mat <- slot(object = x, name = '.Data')
    if (length(x = value)) {
      mat <- cbind(
        mat,
        matrix(data = FALSE, nrow = nrow(x = x), dimnames = list(NULL, i))
      )
      mat[value, i] <- TRUE
    }
    slot(object = x, name = '.Data') <- mat
    validObject(object = x)
    return(x)
  }
)

#' @rdname LogMap-class
#'
#' @export
#'
#' @order 7
#'
setMethod(
  f = '[[<-',
  signature = c(x = 'LogMap', i = 'character', j = 'missing', value = 'NULL'),
  definition = function(x, i, ..., value) {
    mat <- slot(object = x, name = '.Data')
    for (name in i) {
      idx <- which(x = colnames(x = mat) == name)
      if (length(x = idx)) {
        mat <- mat[, -idx, drop = FALSE]
      }
    }
    slot(object = x, name = '.Data') <- mat
    validObject(object = x)
    return(x)
  }
)

#' @rdname LogMap-class
#'
#' @export
#'
#' @order 8
#'
setMethod(
  f = '[[<-',
  signature = c(
    x = 'LogMap',
    i = 'character',
    j = 'missing',
    value = 'numeric'
  ),
  definition = function(x, i, ..., value) {
    value <- as.integer(x = value)
    x[[i]] <- value
    return(x)
  }
)

#' \code{\link{LogMap}} Object Overview
#'
#' Overview of a \code{\link{LogMap}} object
#'
#' @param object A \code{\link{LogMap}} object
#'
#' @template return-show
#'
#' @concept logmap
#'
setMethod(
  f = 'show',
  signature = 'LogMap',
  definition = function(object) {
    cat(
      "A logical map for",
      nrow(x = object),
      "values across",
      ncol(x = object),
      "observations"
    )
    return(invisible(x = NULL))
  }
)

#' Logical Map Validity
#'
#' @templateVar cls LogMap
#' @template desc-validity
#'
#' @section Data Validation:
#' Logical maps must be a logical matrix containing only TRUE or FALSE values
#'
#' @section Value Validation:
#' All values must be named within the rownames of the object. Duplicate or
#' empty (\code{""}) values are not allowed
#'
#' @section Observation Validation:
#' All observations must be named within the column names of the object.
#' Duplicate or empty (\code{""}) observations are not allowed
#'
#' @name LogMap-validity
#'
#' @family logmap
#' @seealso \code{\link[methods]{validObject}}
#'
#' @examples
#' map <- LogMap(letters[1:10])
#' map[['obs']] <- c(1, 3, 7)
#' map[['entry']] <- c(2, 7, 10)
#' validObject(map)
#'
setValidity(
  Class = 'LogMap',
  method = function(object) {
    if (isFALSE(x = getOption(x = "Seurat.object.validate", default = TRUE))) {
      warn(
        message = paste("Not validating", class(x = object)[1L], "objects"),
        class = 'validationWarning'
      )
      return(TRUE)
    }
    valid <- NULL
    # Ensure we have a logical matrix
    if (!is.logical(x = object)) {
      valid <- c(valid, "The map must be a logical matrix")
    }
    if (any(is.na(x = object))) {
      valid <- c(valid, "The may may not contain NAs")
    }
    # Check rownames
    if (is.null(x = rownames(x = object))) {
      valid <- c(valid, "Rownames must be supplied")
    }
    if (any(!nzchar(x = rownames(x = object)))) {
      valid <- c(valid, "Rownames cannot be empty strings")
    }
    if (anyDuplicated(x = rownames(x = object))) {
      valid <- c(valid, "Duplicate rownames not allowed")
    }
    # Check colnames
    if (!is.null(x = colnames(x = object))) {
      if (any(!nzchar(x = colnames(x = object)))) {
        valid <- c(valid, "Columnn names cannot be empty strings")
      }
      if (anyDuplicated(x = colnames(x = object))) {
        valid <- c(valid, "Duplicate colnames not allowed")
      }
    } else if (ncol(x = object)) {
      valid <- c(valid, "Colnames must be supplied")
    }
    return(valid %||% TRUE)
  }
)
