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
#' (columns). \code{LogMap} objects can be created with \code{LogMap()}; queries
#' can be performed with \code{[[} and observations can be added or removed
#' with \code{[[<-}
#'
#' @slot .Data A logical matrix with at least one row
#'
#' @exportClass LogMap
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
#' # Find observations for a set of values
#' map[['entry']] <- c(2, 7, 10)
#' labels(map, c('a', 'b', 'g'))
#'
#' # Remove unused values
#' map <- droplevels(map)
#' map
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

#' @method as.matrix LogMap
#' @export
#'
as.matrix.LogMap <- function(x, ...) {
  return(as(object = x, Class = 'matrix'))
}

#' @rdname LogMap-class
#'
#' @param x,object A \code{LogMap} object
#'
#' @return \code{droplevels}: \code{x} with values not present in any
#' observation removed
#'
#' @method droplevels LogMap
#' @export
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

#' @rdname LogMap-class
#'
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
  obs <- pbapply::pbsapply(
    X = values,
    FUN = function(x) {
      return(colnames(x = object)[which(x = object[x, , drop = TRUE])])
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
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

setMethod(
  f = '[',
  signature = c(x = 'LogMap', i = 'missing', j = 'missing'),
  definition = function(x, i, j, ..., drop = FALSE) {
    return(x)
  }
)

setMethod(
  f = '[',
  signature = c(x = 'LogMap', i = 'character', j = 'character'),
  definition = function(x, i, j, ..., drop = FALSE) {
    x <- as.matrix(x = x)[i, j, drop = drop]
    if (!isTRUE(x = drop)) {
      x <- as(object = x, Class = 'LogMap')
    }
    return(x)
  }
)

setMethod(
  f = '[',
  signature = c(x = 'LogMap', i = 'character', j = 'missing'),
  definition = function(x, i, j, ..., drop = FALSE) {
    x <- as.matrix(x = x)[i, , drop = drop]
    if (!isTRUE(x = drop)) {
      x <- as(object = x, Class = 'LogMap')
    }
    return(x)
  }
)

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

setMethod(
  f = '[',
  signature = c(x = 'LogMap', i = 'numeric', j = 'missing'),
  definition = function(x, i, j, ..., drop = FALSE) {
    i <- rownames(x = x)[i]
    return(callNextMethod(x, i, ..., drop = drop))
  }
)

setMethod(
  f = '[',
  signature = c(x = 'LogMap', i = 'missing', j = 'numeric'),
  definition = function(x, i, j, ..., drop = FALSE) {
    j <- colnames(x = x)[j]
    return(callNextMethod(x, , j, ..., drop = drop))
  }
)


setMethod(
  f = '[',
  signature = c(x = 'LogMap', i = 'numeric', j = 'numeric'),
  definition = function(x, i, j, ..., drop = FALSE) {
    i <- rownames(x = x)[i]
    j <- colnames(x = x)[j]
    return(callNextMethod(x, i, j, ..., drop = drop))
  }
)


#' @rdname LogMap-class
#'
#' @param x,object A \code{LogMap} object
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
    # Check rownames
    if (is.null(x = rownames(x = object))) {
      valid <- c(valid, "Rownames must be supplied")
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
