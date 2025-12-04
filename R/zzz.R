#' @importFrom sp bbox over
#' @importFrom Rcpp evalCpp
#' @importFrom Matrix nnzero
#' @importFrom progressr progressor
#' @importFrom lifecycle deprecated is_present
#' @importFrom utils head packageVersion tail upgrade
#' @importFrom methods new setClass setClassUnion setGeneric setMethod
#' setOldClass setValidity show slot slot<- validObject
#' @importFrom rlang abort arg_match arg_match0 caller_env check_installed
#' enquo eval_tidy have_name inform is_bare_character is_bare_integerish
#' is_bare_list is_bare_numeric is_missing is_na is_named is_quosure
#' missing_arg warn
#' @importClassesFrom Matrix dgCMatrix CsparseMatrix dsparseMatrix generalMatrix
#' dMatrix sparseMatrix compMatrix Matrix
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

#' \pkg{Seurat} Options
#'
#' Various options used in \pkg{Seurat}
#'
#' @section Package Options:
#' \subsection{Seurat.coords.short_range}{
#' Defaults to
#' \dQuote{\Sexpr[stage=build]{SeuratObject:::Seurat.options$Seurat.coords.short_range}}\cr
#' Currently set to \dQuote{\Sexpr[stage=render]{getOption("Seurat.coords.short_range")}}
#' }
#' \subsection{Seurat.input.sparse_ratio}{
#' Defaults to
#' \dQuote{\Sexpr[stage=build]{SeuratObject:::Seurat.options$Seurat.input.sparse_ratio}}\cr
#' Currently set to \dQuote{\Sexpr[stage=render]{getOption("Seurat.input.sparse_ratio")}}
#' }
#' \subsection{Seurat.io.rds.strict}{
#' Defaults to
#' \dQuote{\Sexpr[stage=build]{SeuratObject:::Seurat.options$Seurat.io.rds.strict}}\cr
#' Currently set to \dQuote{\Sexpr[stage=render]{getOption("Seurat.io.rds.strict")}}
#' }
#' \subsection{Seurat.object.assay.calcn}{
#' Run \code{CalcN} when adding assay data to a \code{Seurat} object\cr
#' Defaults to
#' \dQuote{\Sexpr[stage=build]{SeuratObject:::Seurat.options$Seurat.object.assay.calcn}}\cr
#' Currently set to \dQuote{\Sexpr[stage=render]{getOption("Seurat.object.assay.calcn")}}
#' }
#' \subsection{Seurat.object.assay.version}{
#' Defaults to
#' \dQuote{\Sexpr[stage=build]{SeuratObject:::Seurat.options$Seurat.object.assay.version}}\cr
#' Currently set to \dQuote{\Sexpr[stage=render]{getOption("Seurat.object.assay.version")}}
#' }
#' \subsection{Seurat.object.assay.v3.missing_layer}{
#' Defaults to
#' \dQuote{\Sexpr[stage=build]{SeuratObject:::Seurat.options$Seurat.object.assay.v3.missing_layer}}\cr
#' Currently set to \dQuote{\Sexpr[stage=render]{getOption("Seurat.object.assay.v3.missing_layer")}}
#' }
#' \subsection{Seurat.object.project}{
#' Default project for new \code{\link{Seurat}} objects\cr
#' Defaults to
#' \dQuote{\Sexpr[stage=build]{SeuratObject:::Seurat.options$Seurat.object.project}}\cr
#' Currently set to \dQuote{\Sexpr[stage=render]{getOption("Seurat.object.project")}}
#' }
#'
#' @name SeuratObject-options
#'
#' @keywords internal
#'
NULL

Seurat.options <- list(
  Seurat.coords.short_range = 'max',
  Seurat.input.sparse_ratio = 0.4,
  Seurat.io.rds.strict = FALSE,
  Seurat.object.assay.brackets = 'v5',
  Seurat.object.assay.calcn = NULL,
  Seurat.object.assay.version = 'v5',
  Seurat.object.assay.v3.missing_layer = 'matrix',
  Seurat.object.project = 'SeuratProject',
  progressr.clear = FALSE
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Built With
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.BuiltWith <- c(
  R = format(x = getRversion()),
  Matrix = format(x = packageVersion(pkg = "Matrix"))
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Reexports
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom future plan
#' @export
#'
future::plan

#' @importFrom generics intersect
#' @export
#'
generics::intersect

# #' @importFrom Matrix colMeans
# #' @export
# #'
# Matrix::colMeans

#' @importFrom Matrix t
#' @export
#'
#' @noRd
#'
Matrix::t

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
setClassUnion(name = 'OptionalLogical', members = c('NULL', 'logical'))

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
#' @param i,j \link[sp:bbox]{Bounding boxes}
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

.IsDataFrame <- function(x) {
  return(length(x = class(x = x)) == 1L && inherits(x = x, what = 'data.frame'))
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
  # If `Seurat` is installed, check against it's version.
  # This is a (probably bad/hacky) workaround to avoid calling `requireNamespace`
  # and forcing us to list `Seurat` as a dependency.
  if (nzchar(system.file(package = "Seurat"))) {
    future <- future ||
      packageVersion("Seurat", lib.loc = lib.loc) >= version
  }
  return(future)
}

.IsNull <- function(x) {
  return(vapply(X = x, FUN = is.null, FUN.VALUE = logical(length = 1L)))
}

.MakeNames <- function(x, strict = FALSE, type = c('layers')) {
  if (isTRUE(x = strict)) {
    return(make.names(names = x, unique = TRUE))
  }
  type <- type[[1L]]
  type <- match.arg(arg = type)
  x <- switch(
    EXPR = type,
    layers = {
      # Remove white spaces
      x <- gsub(pattern = '[[:space:]]+', replacement = '_', x = x)
      # Remove illegal characters
      x <- gsub(
        pattern = '[\\;\\:\\!\\@\\#\\$\\%\\^\\&\\*\\(\\)\\{\\}\\[]',
        replacement = '',
        x = x
      )
      x <- gsub(pattern = '\\]', replacement = '', x = x)
      x
    }
  )
  return(x)
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

#' Round Version Information
#'
#' @param current A package version
#'
#' @return ...
#'
#' @keywords internal
#'
#' @noRd
#'
.RoundVersion <- function(current) {
  current <- as.character(x = numeric_version(x = current, strict = TRUE))
  current <- unlist(x = strsplit(x = current, split = '\\.'))
  if (length(x = current) > 4L) {
    if (length(x = current) > 4L) {
      current[4L] <- paste(
        current[seq.int(from = 4L, to = length(x = current))],
        collapse = '.'
      )
      current <- current[1:4]
    }
  }
  names(x = current) <- c('major', 'minor', 'patch', 'devel')[seq_along(along.with = current)]
  if (!is_na(x = current['devel'])) {
    if (all(current[c('minor', 'patch')] == '9')) {
      current['major'] <- as.character(x = as.integer(x = current['major']) + 1L)
      current[c('minor', 'patch')] <- '0'
    } else if (current['patch'] == '0') {
      current['minor'] <- as.character(x = as.integer(x = current['minor']) + 1L)
      current['patch'] <- '0'
    } else {
      current['patch'] <- as.character(x = as.integer(x = current['patch']) + 1L)
    }
    current <- current[c('major', 'minor', 'patch')]
  }
  current <- vapply(
    X = current,
    FUN = as.integer,
    FUN.VALUE = integer(length = 1L),
    USE.NAMES = TRUE
  )
  if (!is_na(x = current['devel'])) {
    if (all(current[c('minor', 'patch')] == '9')) {
      current['major'] <- as.character(x = as.integer(x = current['major']) + 1L)
      current[c('minor', 'patch')] <- '0'
    } else if (current['patch'] == '0') {
      current['minor'] <- as.character(x = as.integer(x = current['minor']) + 1L)
      current['patch'] <- '0'
    } else {
      current['patch'] <- as.character(x = as.integer(x = current['patch']) + 1L)
    }
    current <- current[c('major', 'minor', 'patch')]
  }
  return(current)
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
  for (i in names(x = .BuiltWith)) {
    current <- switch(EXPR = i, R = getRversion(), packageVersion(pkg = i))
    if (current > .BuiltWith[i]) {
      msg <- paste(
        sQuote(x = pkgname),
        "was built",
        switch(
          EXPR = i,
          R = "under R",
          paste("with package", sQuote(x = i))
        ),
        .BuiltWith[i],
        "but the current version is",
        paste0(current, ';'),
        "it is recomended that you reinstall ",
        sQuote(x = pkgname),
        " as the ABI for",
        switch(EXPR = i, R = i, sQuote(x = i)),
        "may have changed"
      )
      packageStartupMessage(paste(strwrap(x = msg), collapse = '\n'))
    }
  }
}

.onLoad <- function(libname, pkgname) {
  toset <- setdiff(x = names(x = Seurat.options), y = names(x = options()))
  if (length(x = toset)) {
    options(Seurat.options[toset])
  }
  setHook(
    hookName = packageEvent(pkgname = 'Seurat', event = 'onLoad'),
    value = .SetSeuratCompat
  )
  setHook(
    hookName = packageEvent(pkgname = 'Signac', event = 'onLoad'),
    value = .SetSeuratCompat
  )
  setHook(
    hookName = packageEvent(pkgname = 'Seurat', event = 'attach'),
    value = .SeuratCompatMessage
  )
  setHook(
    hookName = packageEvent(pkgname = 'Signac', event = 'attach'),
    value = .SeuratCompatMessage
  )
  return(invisible(x = NULL))
}
