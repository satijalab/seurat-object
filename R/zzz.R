#' @importFrom sp bbox over
#' @importFrom Rcpp evalCpp
#' @importFrom utils head tail
#' @importFrom rlang abort arg_match check_installed inform is_bare_character
#' is_bare_integerish is_bare_list is_bare_numeric is_na warn
#' @importFrom lifecycle deprecated deprecate_soft deprecate_stop
#' deprecate_warn is_present
#' @importFrom methods new setClass setClassUnion setGeneric setMethod
#' setOldClass setValidity show slot slot<- validObject
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
# Package options
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

default.options <- list()

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

# #' @importFrom Matrix colMeans
# #' @export
# #'
# Matrix::colMeans

#' @importFrom progressr handlers
#' @export
#'
progressr::handlers

#' @importFrom progressr with_progress
#' @export
#'
progressr::with_progress

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Environments
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sparse.classes <- new.env()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setClassUnion(name = 'OptionalCharacter', members = c('NULL', 'character'))
setClassUnion(name = 'OptionalList', members = c('NULL', 'list'))

setOldClass(Classes = 'package_version')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

#' Test Future Compatibility with \pkg{Seurat}
#'
#' Check to see if \pkg{SeuratObject} and/or \pkg{Seurat} are at least a
#' specific version or if they're configured to act as if they're a
#' specific version (see details below). This allows testing compatibility with
#' future requirements for both \pkg{SeuratObject} and \pkg{Seurat}
#'
#' Blah blah blah
#'
#' @inheritParams utils::packageVersion
#' @param version A version string or object of class
#' \code{\link{package_version}}
#'
#' @return \code{TRUE} if \pkg{SeuratObject} and/or \pkg{Seurat}
#'
#' @importFrom utils packageVersion
#'
#' @keywords internal
#'
#' @export
#'
#' @aliases IsFutureSeurat
#'
.IsFutureSeurat <- function(version, lib.loc = NULL) {
  version <- package_version(x = version)
  opt <- paste0(
    'Seurat.future.v',
    gsub(pattern = '\\.', replacement = '_', x = as.character(x = version))
  )
  future <- isTRUE(x = getOption(x = opt, default = FALSE)) ||
    packageVersion(pkg = 'SeuratObject', lib.loc = lib.loc) >= version
  if (requireNamespace('Seurat', quietly = TRUE)) {
    future <- future ||
      packageVersion(pkg = 'Seurat', lib.loc = lib.loc) >= version
  }
  return(future)
}

#' Create a List with a Serial Comma
#'
#' @param ... A character vector to join
#' @param cnj Conjunction to use for final entry
#' @param quote Quote the entries of \code{...}; choose from:
#' \itemize{
#'  \item \dQuote{\code{single}}: regular single quotes
#'  \item \dQuote{\code{fancysingle}}: fancy single quotes
#'  \item \dQuote{\code{double}}: regular double quotes
#'  \item \dQuote{\code{fancydouble}}: fancy double quotes
#'  \item \dQuote{\code{none}}: no extra quoting
#' }
#'
#' @return \code{...} arranged into an English list with a serial comma
#' when needed
#'
#' @keywords internal
#'
#' @seealso \code{\link[base]{sQuote}}
#'
#' @examples
#' .Oxford('cell')
#' .Oxford('cell', 'ident')
#' .Oxford('cell', 'ident', 'gene')
#'
#' @noRd
#'
.Oxford <- function(
  ...,
  cnj = c('or', 'and'),
  quote = c('single', 'fancysingle', 'double', 'fancydouble', 'none')
) {
  x <- as.character(x = c(...))
  cnj <- arg_match(arg = cnj)
  quote <- arg_match(arg = quote)
  x <- switch(
    EXPR = quote,
    single = sQuote(x = x, q = FALSE),
    fancysingle = sQuote(x = x, q = TRUE),
    double = dQuote(x = x, q = FALSE),
    fancydouble = dQuote(x = x, q = TRUE),
    x
  )
  if (length(x = x) <= 1L) {
    return(x)
  } else if (length(x = x) == 2L) {
    return(paste(x, collapse = paste0(' ', cnj, ' ')))
  }
  return(paste(
    paste0(paste(x[1:(length(x = x) - 1L)], collapse = ', '), ','),
    cnj,
    x[length(x = x)]
  ))
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

.onAttach <- function(libname, pkgname) {
  AttachDeps(deps = 'sp')
  return(invisible(x = NULL))
}

.onLoad <- function(libname, pkgname) {
  toset <- setdiff(x = names(x = Seurat.options), y = names(x = options()))
  if (length(x = toset)) {
    options(Seurat.options[toset])
  }
  return(invisible(x = NULL))
}
