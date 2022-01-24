#' @importFrom utils head tail
#' @importFrom sp bbox over
#' @importFrom methods new setClass setClassUnion setMethod setOldClass
#' setValidity slot slot<- validObject
#' @importClassesFrom Matrix dgCMatrix
#' @useDynLib SeuratObject
#'
NULL

#' @docType package
#' @name SeuratObject-package
#' @rdname SeuratObject-package
#'
"_PACKAGE"

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Options
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Seurat.options <- list(
  Seurat.input.sparse_ratio = 0.4,
  Seurat.coords.short_range = 'max',
  progressr.clear = FALSE
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Reexports
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom future plan
#' @export
#'
future::plan

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

#' @importFrom progressr handlers
#' @export
#'
progressr::handlers

#' @importFrom progressr with_progress
#' @export
#'
progressr::with_progress

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

#' @keywords internal
#'
#' @noRd
#'
.AutoRadius <- function(coords) {
  return(0.01 * mean(x = apply(
    X = apply(X = coords, MARGIN = 2L, FUN = range),
    MARGIN = 2L,
    FUN = diff
  )))
}

#' @keywords internal
#'
#' @noRd
#'
.BboxDF <- function(x) {
  df <- expand.grid(x = x['x', ], y = x['y', ])
  df <- df[c(1, 3, 4, 2), ]
  return(df)
}

#' Test Intersections of Bounding Boxes
#'
#' @param i,j \link[sp::bbox]{Bounding boxes}
#' @param constraint Type of intersection to perform; choose from:
#' \itemize{
#'  \item \dQuote{\code{intersect}}: \code{i} must fall at least
#'   partially within the bounds of \code{j} for the dimensions
#'   specified by \code{MARGIN}
#'  \item \dQuote{\code{contained}}: \code{i} must fall completely
#'   within the bounds of \code{j} for the dimensions specified
#'   by \code{MARIGN}
#'  \item \dQuote{\code{overlap}}: \code{i} must fall at least partially
#'   within \code{j}, or \code{j} must fall at least partially within
#'   \code{i}, for the dimensions specified by \code{MARGIN}; essentially
#'   \code{.BboxIntersect(i, j, 'intersect', MARGIN) || .BboxIntersect(j, i, 'intersect', MARGIN)}
#' }
#' @param MARGIN Direction of intersection; choose from:
#' \itemize{
#'  \item \code{1L}: intersect along the x-dimension
#'  \item \code{2L}: intersect along the y-dimension
#'  \item \code{3L}: intersect along both the x- and y-dimensions
#' }
#'
#' @return \code{TRUE} if \code{i} intersects with \code{j};
#' otherwise \code{FALSE}
#'
#' @keywords internal
#'
#' @noRd
#'
.BboxIntersect <- function(
  i,
  j,
  constraint = c('intersect', 'contained', 'overlap'),
  MARGIN = 3L
) {
  constraint <- constraint[1L]
  constraint <- match.arg(arg = constraint)
  if (!MARGIN %in% seq.int(from = 1L, to = 3L)) {
    stop(".MARGIN must be 1, 2, or 3")
  } else if (MARGIN == 3L) {
    MARGIN <- seq.int(from = 1L, to = 2L)
  }
  check <- vector(mode = 'logical', length = length(x = MARGIN))
  names(x = check) <- c('x', 'y')[MARGIN]
  for (x in names(x = check)) {
    check[[x]] <- switch(
      EXPR = constraint,
      'intersect' = {
        (i[x, 'min'] >= j[x, 'min'] && i[x, 'min'] <= j[x, 'max']) ||
          (i[x, 'max'] >= j[x, 'min'] && i[x, 'max'] <= j[x, 'max'])
      },
      'contained' = i[x, 'min'] >= j[x, 'min'] && i[x, 'max'] <= j[x, 'max'],
      'overlap' = {
        .margin <- c(x = 1L, y = 2L)[x]
        .BboxIntersect(i = i, j = j, constraint = 'i', MARGIN = .margin) ||
          .BboxIntersect(i = j, j = i, constraint = 'i', MARGIN = .margin)
      },
      stop("Constraint '", constraint, "' not yet implemented")
    )
  }
  return(all(check))
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
.Crop <- function(object, x = NULL, y = NULL, coords = c('plot',' tissue'), ...) {
  if (is.null(x = x) && is.null(x = y)) {
    return(object)
  }
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
  return(Overlay(x = object, y = CreateSegmentation(coords = df)))
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

.PruneLogMap <- function(x) {
  fidx <- which(x = apply(
    X = x,
    MARGIN = 1L,
    FUN = function(row) {
      return(all(vapply(X = row, FUN = isFALSE, FUN.VALUE = logical(length = 1L))))
    }
  ))
  if (length(x = fidx)) {
    x <- as(object = x[-fidx, , drop = FALSE], Class = 'LogMap')
  }
  validObject(object = x)
  return(x)
}

#' Indexes from Run Length Encodings
#'
#' Generate an index for subsetting from a \link[base:rle]{run length encoding}
#'
#' @inheritParams base::lengths
#' @param x An \code{\link[base:rle]{rle}} object
#'
#' @return A list where each entry is the indices a particular value
#'
#' @keywords internal
#'
#' @noRd
#'
.RleIndex <- function(x, use.names = TRUE) {
  idx <- lapply(
    X = seq_len(length.out = length(x = x$values)),
    FUN = function(i) {
      from <- (x$lengths[i] * i) - (x$lengths[i] - 1L)
      return(seq.int(from = from, to = from + x$lengths[i] - 1L))
    }
  )
  if (isTRUE(x = use.names)) {
    names(x = idx) <- x$values
  }
  return(idx)
}

#' Index of Names
#'
#' Get the index of row- or column-names
#'
#' @param x A two-dimensional object
#' @param names A vector of names to index
#' @param MARGIN Either \code{1L} for row-names or \code{2L} for column-names
#'
#' @return A named integer vector of length \code{length(names)}; the names are
#' \code{names} and the values are the index of \code{names} in the row- or
#' column-names. If no name is found, uses the lowest available index
#'
#' @importFrom stats na.omit
#'
#' @keywords internal
#'
#' @noRd
#'
NameIndex <- function(x, names, MARGIN) {
  if (!MARGIN %in% c(1L, 2L)) {
    stop("MARGIN must be either 1 or 2", call. = FALSE)
  }
  if (!length(x = dim(x = x)) == 2L) {
    stop("'x' must be a two-dimensional object", call. = FALSE)
  }
  nfunc <- list(rownames, colnames)[[MARGIN]]
  xnames <- nfunc(x = x)
  if (length(x = names) > length(x = xnames)) {
    stop(
      "Too many names requested (",
      length(x = names),
      " requested, ",
      length(x = xnames),
      " provided)",
      call. = FALSE
    )
  }
  idx <- vector(mode = 'integer', length = length(x = names))
  names(x = idx) <- names
  for (i in names) {
    idx[[i]] <- ifelse(
      test = i %in% xnames,
      yes = which(x = xnames == i)[1],
      no = NA_integer_
    )
  }
  idx.na <- which(x = is.na(x = idx))
  xind <- setdiff(
    x = seq_len(length.out = ncol(x = x)),
    y = na.omit(object = idx)
  )
  for (i in idx.na) {
    idx[[i]] <- xind[[i]]
  }
  return(idx)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Hooks
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.onLoad <- function(libname, pkgname) {
  toset <- setdiff(x = names(x = Seurat.options), y = names(x = options()))
  if (length(x = toset)) {
    options(Seurat.options[toset])
  }
  return(invisible(x = NULL))
}
