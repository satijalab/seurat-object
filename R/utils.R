#' @include zzz.R
#' @include generics.R
#' @importFrom Rcpp evalCpp
#' @useDynLib SeuratObject
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Set a default value depending on if an object is \code{NULL}
#'
#' @param x An object to test
#' @param y A default value
#'
#' @return For \code{\%||\%}: \code{y} if \code{x} is \code{NULL} otherwise
#' \code{x}
#'
#' @importFrom rlang %||%
#'
#' @name set-if-null
#' @rdname set-if-null
#'
#' @export
#'
#' @concept utils
#'
#' @examples
#' 1 %||% 2
#' NULL %||% 2
#'
rlang::`%||%`

#' @rdname set-if-null
#'
#' @return For \code{\%iff\%}: \code{y} if \code{x} is \strong{not}
#' \code{NULL}; otherwise \code{x}
#'
#' @importFrom rlang is_null
#'
#' @export
#'
#' @examples
#' 1 %iff% 2
#' NULL %iff% 2
#'
`%iff%` <- function(x, y) {
  if (!is_null(x = x)) {
    return(y)
  }
  return(x)
}

#' Attach Required Packages
#'
#' Helper function to attach required packages. Detects if a package is already
#' attached and if so, skips it. Should be called in \code{\link[base]{.onAttach}}
#'
#' @param deps A character vector of packages to attach
#'
#' @return Invisibly returns \code{NULL}
#'
#' @export
#'
#' @concept utils
#'
#' @examples
#' # Use in your .onAttach hook
#' if (FALSE) {
#'   .onAttach <- function(libname, pkgname) {
#'     AttachDeps(c("SeuratObject", "rlang"))
#'   }
#' }
#'
AttachDeps <- function(deps) {
  for (d in deps) {
    if (!paste0('package:', d) %in% search()) {
      packageStartupMessage("Attaching ", d)
      attachNamespace(ns = d)
    }
  }
  return(invisible(x = NULL))
}

#' Conditional Garbage Collection
#'
#' Call \code{gc} only when desired
#'
#' @param option ...
#'
#' @return Invisibly returns \code{NULL}
#'
#' @export
#'
#' @concept utils
#'
CheckGC <- function(option = 'SeuratObject.memsafe') {
  if (isTRUE(x = getOption(x = option, default = FALSE))) {
    gc(verbose = FALSE)
  }
  return(invisible(x = NULL))
}

#' Find the default \code{\link{DimReduc}}
#'
#' Searches for \code{\link{DimReduc}s} matching \dQuote{umap}, \dQuote{tsne},
#' or \dQuote{pca}, case-insensitive, and in that order. Priority given to
#' \code{\link{DimReduc}s} matching the \code{DefaultAssay} or assay specified
#' (eg. \dQuote{pca} for the default assay weights higher than \dQuote{umap}
#' for a non-default assay)
#'
#' @param object A \code{\link{Seurat}} object
#' @param assay Name of assay to use; defaults to the default assay of the object
#'
#' @return The default \code{\link{DimReduc}}, if possible
#'
#' @export
#'
#' @concept utils
#'
#' @examples
#' DefaultDimReduc(pbmc_small)
#'
DefaultDimReduc <- function(object, assay = NULL) {
  object <- UpdateSlots(object = object)
  assay <- assay %||% DefaultAssay(object = object)
  drs.use <- c('umap', 'tsne', 'pca')
  dim.reducs <- FilterObjects(object = object, classes.keep = 'DimReduc')
  drs.assay <- Filter(
    f = function(x) {
      return(DefaultAssay(object = object[[x]]) == assay)
    },
    x = dim.reducs
  )
  if (length(x = drs.assay) > 0) {
    index <- lapply(
      X = drs.use,
      FUN = grep,
      x = drs.assay,
      ignore.case = TRUE
    )
    index <- Filter(f = length, x = index)
    if (length(x = index) > 0) {
      return(drs.assay[min(index[[1]])])
    }
  }
  index <- lapply(
    X = drs.use,
    FUN = grep,
    x = dim.reducs,
    ignore.case = TRUE
  )
  index <- Filter(f = length, x = index)
  if (length(x = index) < 1) {
    stop(
      "Unable to find a DimReduc matching one of '",
      paste(drs.use[1:(length(x = drs.use) - 1)], collapse = "', '"),
      "', or '",
      drs.use[length(x = drs.use)],
      "', please specify a dimensional reduction to use",
      call. = FALSE
    )
  }
  return(dim.reducs[min(index[[1]])])
}

#' Generate a Class Key
#'
#' Generate class keys for S4 classes. A class key follows the following
#' structure: \dQuote{\code{package:class}}
#'
#' @param class Class name
#' @param package Optional name of package; by default, will search namespaces
#' of loaded packages to determine the providing package
#'
#' @return The class key (\dQuote{\code{package:class}})
#'
#' @importFrom methods getClass slot
#'
#' @export
#'
#' @examples
#' ClassKey("Seurat")
#'
ClassKey <- function(class, package = NULL) {
  class <- class[1L]
  package <- package %||% slot(
    object = getClass(Class = class),
    name = 'package'
  )
  return(paste(package, class, sep = ':'))
}


#' @name s4list
#' @rdname s4list
#'
#' @return \code{IsS4List}: \code{TRUE} if \code{x} is a list with an S4 class
#' definition attribute
#'
#' @export
#'
IsS4List <- function(x) {
  return(
    inherits(x = x, what = 'list') &&
      isTRUE(x = grepl(
        pattern = '^[[:alnum:]]+:[[:alnum:]]+$',
        x = attr(x = x, which = 'classDef')
      ))
  )
}

#' @name s4list
#' @rdname s4list
#'
#' @return \code{ListToS4}: An S4 object as defined by the S4 class definition
#' attribute
#'
#' @importFrom methods getClassDef new
#'
#' @export
#'
ListToS4 <- function(x) {
  if (!inherits(x = x, what = 'list')) {
    return(x)
  }
  for (i in seq_along(along.with = x)) {
    if (!is.null(x = x[[i]])) {
      x[[i]] <- ListToS4(x = x[[i]])
    }
  }
  classdef <- attr(x = x, which = 'classDef')
  x <- Filter(f = Negate(f = is.function), x = x)
  attr(x = x, which = 'classDef') <- classdef
  if (!IsS4List(x = x)) {
    return(x)
  }
  classdef <- unlist(x = strsplit(
    x = attr(x = x, which = 'classDef'),
    split = ':'
  ))
  pkg <- classdef[1]
  cls <- classdef[2]
  formal <- getClassDef(Class = cls, package = pkg, inherits = FALSE)
  return(do.call(what = new, args = c(list(Class = formal), x)))
}

#' Check the existence of a package
#'
#' @param ... Package names
#' @param error If true, throw an error if the package doesn't exist
#'
#' @return Invisibly returns boolean denoting if the package is installed
#'
#' @export
#'
#' @concept utils
#'
#' @examples
#' PackageCheck("SeuratObject", error = FALSE)
#'
PackageCheck <- function(..., error = TRUE) {
  pkgs <- unlist(x = c(...), use.names = FALSE)
  package.installed <- vapply(
    X = pkgs,
    FUN = requireNamespace,
    FUN.VALUE = logical(length = 1L),
    quietly = TRUE
  )
  if (error && any(!package.installed)) {
    stop(
      "Cannot find the following packages: ",
      paste(pkgs[!package.installed], collapse = ', '),
      ". Please install"
    )
  }
  invisible(x = package.installed)
}

#' Generate a random name
#'
#' Make a name from randomly sampled characters, pasted together with no spaces
#'
#' @param length How long should the name be
#' @param chars A vector of 1-length characters to use to generate the name
#' @param ... Extra parameters passed to \code{\link[base]{sample}}
#'
#' @return A character with \code{nchar == length} of randomly sampled letters
#'
#' @seealso \code{\link[base]{sample}}
#'
#' @export
#'
#' @concept utils
#'
#' @examples
#' set.seed(42L)
#' RandomName()
#' RandomName(7L, replace = TRUE)
#'
RandomName <- function(length = 5L, chars = letters, ...) {
  CheckDots(..., fxns = 'sample')
  chars <- unique(x = unlist(x = strsplit(
    x = as.character(x = chars),
    split = ''
  )))
  return(paste(sample(x = chars, size = length, ...), collapse = ''))
}

#' Merge Sparse Matrices by Row
#'
#' Merge two or more sparse matrices by rowname.
#'
#' @details
#' Shared matrix rows (with the same row name) will be merged, and unshared
#' rows (with different names) will be filled with zeros in the matrix not
#' containing the row.
#'
#' @param mat1 First matrix
#' @param mat2 Second matrix or list of matrices
#'
#' @return Returns a sparse matrix
#'
#' @importFrom methods as
#
#' @export
#'
#' @concept utils
#'
RowMergeSparseMatrices <- function(mat1, mat2) {
  all.mat <- c(list(mat1), mat2)
  all.colnames <- all.rownames <- vector(
    mode = 'list',
    length = length(x = all.mat)
  )
  for (i in seq_along(along.with = all.mat)) {
    if (is.data.frame(x = all.mat[[1]])) {
      all.mat[[i]] <- as.matrix(x = all.mat[[i]])
    }
    all.rownames[[i]] <- rownames(x = all.mat[[i]])
    all.colnames[[i]] <- colnames(x = all.mat[[i]])
  }
  use.cbind <- all(duplicated(x = all.rownames)[2:length(x = all.rownames)])
  if (isTRUE(x = use.cbind)) {
    new.mat <- do.call(what = cbind, args = all.mat)
  } else {
    all.mat <- lapply(X = all.mat, FUN = as, Class = "RsparseMatrix")
    all.names <- unique(x = unlist(x = all.rownames))
    new.mat <- RowMergeMatricesList(
      mat_list = all.mat,
      mat_rownames = all.rownames,
      all_rownames = all.names
    )
    rownames(x = new.mat) <- make.unique(names = all.names)
  }
  colnames(x = new.mat) <- make.unique(names = unlist(x = all.colnames))
  return(new.mat)
}

#' Update slots in an object
#'
#' @param object An object to update
#'
#' @return \code{object} with the latest slot definitions
#'
#' @importFrom methods slotNames slot
#'
#' @export
#'
#' @concept utils
#'
UpdateSlots <- function(object) {
  object.list <- sapply(
    X = slotNames(x = object),
    FUN = function(x) {
      return(tryCatch(
        expr = slot(object = object, name = x),
        error = function(...) {
          return(NULL)
        }
      ))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  object.list <- Filter(f = Negate(f = is.null), x = object.list)
  object.list <- c('Class' = class(x = object)[1], object.list)
  object <- do.call(what = 'new', args = object.list)
  for (x in setdiff(x = slotNames(x = object), y = names(x = object.list))) {
    xobj <- slot(object = object, name = x)
    if (is.vector(x = xobj) && !is.list(x = xobj) && length(x = xobj) == 0) {
      slot(object = object, name = x) <- vector(
        mode = class(x = xobj),
        length = 1L
      )
    }
  }
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param row.names \code{NULL} or a character vector giving the row names for
#' the data; missing values are not allowed
#'
#' @rdname as.sparse
#' @export
#' @method as.sparse data.frame
#'
as.sparse.data.frame <- function(x, row.names = NULL, ...) {
  CheckDots(...)
  dnames <- list(row.names %||% rownames(x = x), colnames(x = x))
  if (length(x = dnames[[1]]) != nrow(x = x)) {
    stop("Differing numbers of rownames and rows", call. = FALSE)
  }
  x <- as.data.frame(x = x)
  dimnames(x = x) <- dnames
  return(as.sparse(x = as.matrix(x = x)))
}

#' @importFrom methods as
#'
#' @rdname as.sparse
#' @export
#' @method as.sparse Matrix
#'
as.sparse.Matrix <- function(x, ...) {
  CheckDots(...)
  return(as(object = x, Class = 'dgCMatrix'))
}

#' @rdname as.sparse
#' @export
#' @method as.sparse matrix
#'
as.sparse.matrix <- as.sparse.Matrix

#' @rdname CheckMatrix
#' @method CheckMatrix default
#' @export
#'
CheckMatrix.default <- function(object, checks, ...) {
  return(invisible(x = NULL))
}

#' @rdname CheckMatrix
#' @method CheckMatrix dMatrix
#' @export
#'
CheckMatrix.dMatrix <- function(
  object,
  checks = c('infinite', 'logical', 'integer', 'na'),
  ...
) {
  checks <- match.arg(arg = checks, several.ok = TRUE)
  x <- slot(object = object, name = 'x')
  for (i in checks) {
    switch(
      EXPR = i,
      'infinite' = if (any(is.infinite(x = x))) {
        warning("Input matrix contains infinite values")
      },
      'logical' = if (any(is.logical(x = x))) {
        warning("Input matrix contains logical values")
      },
      'integer' = if (!all(round(x = x) == x, na.rm = TRUE)) {
        warning("Input matrix contains non-integer values")
      },
      'na' = if (anyNA(x = x)) {
        warning("Input matrix contains NA/NaN values")
      },
    )
  }
  return(invisible(x = NULL))
}

#' @rdname CheckMatrix
#' @method CheckMatrix lMatrix
#' @export
#'
CheckMatrix.lMatrix <- function(
  object,
  checks = c('infinite', 'logical', 'integer', 'na'),
  ...
) {
  warning("Input matrix contains logical values")
  return(invisible(x = NULL))
}

#' @rdname IsMatrixEmpty
#' @export
#' @method IsMatrixEmpty default
#'
IsMatrixEmpty.default <- function(x) {
  matrix.dims <- dim(x = x)
  matrix.na <- all(matrix.dims == 1) && all(is.na(x = x))
  return(all(matrix.dims == 0) || matrix.na)
}

#' @importFrom methods slotNames
#'
#' @rdname s4list
#' @export
#' @method S4ToList default
#'
S4ToList.default <- function(object) {
  obj.list <- sapply(
    X = slotNames(x = object),
    FUN = function(x) {
      return(S4ToList(object = slot(object = object, name = x)))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  attr(x = obj.list, which = 'classDef') <- paste(
    c(
      attr(x = class(x = object), which = 'package'),
      class(x = object)
    ),
    collapse = ':'
  )
  return(obj.list)
}

#' @rdname s4list
#' @export
#' @method S4ToList list
#'
S4ToList.list <- function(object) {
  if (length(x = object)) {
    for (i in seq_along(along.with = object)) {
      if (!is.null(x = object[[i]])) {
        object[[i]] <- S4ToList(object = object[[i]])
      }
    }
  }
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.CheckNames <- function(x, n) {
  stopifnot(length(x = x) == length(x = n))
  if (is.null(x = names(x = x))) {
    names(x = x) <- n
  }
  if (any(!nzchar(x = names(x = x)))) {
    idx <- which(x = !nzchar(x = names(x = x)))
    n2 <- setdiff(x = n, y = names(x = x))
    if (length(x = idx) != length(x = n2)) {
      stop("Not all provided names fit with the values provided", call. = FALSE)
    }
    names(x = x)[idx] <- n2
  }
  return(x)
}

#' Identify Object Collections
#'
#' Find all collection (named lists) slots in an S4 object
#'
#' @param object An S4 object
#' @param exclude A character vector of slot names to exclude
#' @inheritDotParams IsNamedList
#'
#' @return A character vector of names of collection slots
#'
#' @importFrom methods slotNames
#'
#' @keywords internal
#'
#' @noRd
#'
.Collections <- function(object, exclude = character(length = 0L), ...) {
  if (!isS4(object)) {
    stop("Not an S4 object")
  }
  collections <- slotNames(x = object)
  collections <- Filter(
    f = function(s) {
      return(IsNamedList(x = slot(object = object, name = s), ...))
    },
    x = collections
  )
  if (is.character(x = exclude) && length(x = exclude)) {
    collections <- setdiff(x = collections, y = exclude)
  }
  return(collections)
}

#' Get Parent S4 Classes
#'
#' @param object An \link[methods:Classes_Details]{S4} object
#'
#' @return A vector of class names that \code{object} inherits from
#'
#' @importFrom methods getClass
#'
#' @keywords internal
#'
#' @export
#'
#' @examples
#' .Contains(pbmc_small)
#'
.Contains <- function(object) {
  if (!isS4(object)) {
    stop("'object' not an S4 object")
  }
  return(names(x = slot(
    object = getClass(Class = class(x = object)),
    name = 'contains'
  )))
}

#' Find A Subobject
#'
#' Determine the slot that a subobject is contained in
#'
#' @inheritParams .Collections
#' @param name Name of subobject to find
#'
#' @return The name of the slot that contains \code{name}
#'
#' @keywords internal
#'
#' @noRd
#'
.FindObject <- function(object, name, exclude = c('misc', 'tools'), ...) {
  collections <- .Collections(object = object, exclude = exclude)
  object.names <- sapply(
    X = collections,
    FUN = function(x) {
      return(names(x = slot(object = object, name = x)))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  object.names <- Filter(f = Negate(f = is.null), x = object.names)
  for (i in names(x = object.names)) {
    if (name %in% names(x = slot(object = object, name = i))) {
      return(i)
    }
  }
  return(NULL)
}

#' Find Subobjects Of A Certain Class
#'
#' @inheritParams .Contains
#' @param classes.keep A vector of classes to keep
#'
#' @return A vector of object names that are of class \code{classes.keep}
#'
#' @keywords internal
#'
#' @noRd
#'
.FilterObjects <- function(object, classes.keep = c('StdAssay', 'DimReduc')) {
  collections <- .Collections(object = object, exclude = c('misc', 'tools'))
  subobjects <- unlist(x = lapply(
    X = collections,
    FUN = function(x) {
      return(Filter(
        f = function(i) {
          return(inherits(
            x = slot(object = object, name = x)[[i]],
            what = classes.keep
          ))
        },
        x = names(x = slot(object = object, name = x))
      ))
    }
  ))
  if (!length(x = subobjects)) {
    subobjects <- NULL
  }
  return(subobjects)
}

#' Get An Option
#'
#' @inheritParams base::getOption
#' @param choices A named list of default options; has higher priority
#' than \code{default}
#'
#' @return ...
#'
#' @keywords internal
#'
#' @noRd
#'
.Opt <- function(x, choices = default.options, default = NULL) {
  return(getOption(x = x, default = choices[[x]] %||% default))
}

#' @param pkg Name of package
#' @param external Include packages imported, but not defined, by \code{pkg}
#' @param old Includes S3 classes registered by
#' \code{\link[methods]{setOldClass}}
#' @param unions Include class unions
#'
#' @importFrom methods getClass getClasses isClassUnion isXS3Class
#'
#' @noRd
#'
.PkgClasses <- function(
  pkg = 'SeuratObject',
  external = FALSE,
  old = FALSE,
  unions = FALSE,
  virtual = NA,
  collapse = TRUE,
  include = NULL,
  exclude = NULL
) {
  classes <- getClasses(where = getNamespace(name = pkg))
  include <- intersect(x = include, y = classes)
  # Filter out classes imported, but not defined by pkg
  if (!isTRUE(x = external)) {
    classes <- Filter(
      f = function(x) {
        return(slot(object = getClass(Class = x), name = 'package') == pkg)
      },
      x = classes
    )
  }
  # Filter out S3 classes
  if (!isTRUE(x = old)) {
    classes <- Filter(
      f = function(x) {
        return(!isXS3Class(classDef = getClass(Class = x)))
      },
      x = classes
    )
  }
  # Filter out class unions
  if (!isTRUE(x = unions)) {
    classes <- Filter(f = Negate(f = isClassUnion), x = classes)
  }
  # TODO: Remove virtual classes
  if (isFALSE(x = virtual)) {
    ''
  }
  # TODO: Collapse classes
  if (isTRUE(x = collapse)) {
    ''
  }
  # Add classes back
  classes <- union(x = classes, y = include)
  # Remove excluded classes
  classes <- setdiff(x = classes, y = exclude)
  return(classes)
}

.Vowels <- function() {
  return(c('a', 'e', 'i', 'o', 'u'))
}

#' Check the Use of Dots
#'
#' Function to check the use of unused arguments passed to \code{...}; this
#' function is designed to be called from another function to see if an
#' argument passed to \code{...} remains unused and alert the user if so. Also
#' accepts a vector of function or function names to see if \code{...} can be
#' used in a downstream function
#'
#' Behavior of \code{CheckDots} can be controlled by the following option(s):
#' \describe{
#'  \item{\dQuote{\code{Seurat.checkdots}}}{Control how to alert the presence
#'  of unused arguments in \code{...}; choose from
#'  \itemize{
#'   \item \dQuote{\code{warn}}: emit a warning (default)
#'   \item \dQuote{\code{error}}: throw an error
#'   \item \dQuote{\code{silent}}: no not alert the presence of unused
#'   arguments in \code{...}
#'  }
#'  }
#' }
#'
#' @param ... Arguments passed to a function that fall under \code{...}
#' @param fxns A list/vector of functions or function names
#'
#' @return Emits either an error or warning if an argument passed is unused;
#' invisibly returns \code{NULL}
#'
#' @importFrom utils isS3stdGeneric methods argsAnywhere isS3method
#'
#' @keywords internal
#'
#' @export
#'
#' @examples
#' \dontrun{
#' f <- function(x, ...) {
#'   CheckDots(...)
#'   return(x ^ 2)
#' }
#' f(x = 3, y = 9)
#' }
#'
CheckDots <- function(..., fxns = NULL) {
  args.names <- names(x = list(...))
  if (length(x = list(...)) == 0) {
    return(invisible(x = NULL))
  }
  if (is.null(x = args.names)) {
    stop("No named arguments passed")
  }
  if (length(x = fxns) == 1) {
    fxns <- list(fxns)
  }
  for (f in fxns) {
    if (!(is.character(x = f) || is.function(x = f))) {
      stop("CheckDots only works on characters or functions, not ", class(x = f))
    }
  }
  fxn.args <- suppressWarnings(expr = sapply(
    X = fxns,
    FUN = function(x) {
      x <- tryCatch(
        expr = if (isS3stdGeneric(f = x)) {
          as.character(x = methods(generic.function = x))
        } else {
          x
        },
        error = function(...) {
          return(x)
        }
      )
      x <- if (is.character(x = x)) {
        sapply(X = x, FUN = argsAnywhere, simplify = FALSE, USE.NAMES = TRUE)
      } else if (length(x = x) <= 1) {
        list(x)
      }
      return(sapply(
        X = x,
        FUN = function(f) {
          return(names(x = formals(fun = f)))
        },
        simplify = FALSE,
        USE.NAMES = TRUE
      ))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  ))
  fxn.args <- unlist(x = fxn.args, recursive = FALSE)
  fxn.null <- vapply(
    X = fxn.args,
    FUN = is.null,
    FUN.VALUE = logical(length = 1L)
  )
  if (all(fxn.null) && !is.null(x = fxns)) {
    stop("None of the functions passed could be found", call. = FALSE)
  } else if (any(fxn.null)) {
    warning(
      "The following functions passed could not be found: ",
      paste(names(x = which(x = fxn.null)), collapse = ', '),
      call. = FALSE,
      immediate. = TRUE
    )
    fxn.args <- Filter(f = Negate(f = is.null), x = fxn.args)
  }
  dfxns <- vector(mode = 'logical', length = length(x = fxn.args))
  names(x = dfxns) <- names(x = fxn.args)
  for (i in 1:length(x = fxn.args)) {
    dfxns[i] <- any(grepl(pattern = '...', x = fxn.args[[i]], fixed = TRUE))
  }
  if (any(dfxns)) {
    dfxns <- names(x = which(x = dfxns))
    if (any(nchar(x = dfxns) > 0)) {
      fx <- vapply(
        X = Filter(f = nchar, x = dfxns),
        FUN = function(x) {
          if (isS3method(method = x)) {
            x <- unlist(x = strsplit(x = x, split = '\\.'))
            x <- x[length(x = x) - 1L]
          }
          return(x)
        },
        FUN.VALUE = character(length = 1L)
      )
      message(
        "The following functions and any applicable methods accept the dots: ",
        paste(unique(x = fx), collapse = ', ')
      )
      if (any(nchar(x = dfxns) < 1)) {
        message(
          "In addition, there is/are ",
          length(x = Filter(f = Negate(f = nchar), x = dfxns)),
          " other function(s) that accept(s) the dots"
        )
      }
    } else {
      message("There is/are ", length(x = dfxns), 'function(s) that accept(s) the dots')
    }
  } else {
    unused <- Filter(
      f = function(x) {
        return(!x %in% unlist(x = fxn.args))
      },
      x = args.names
    )
    if (length(x = unused) > 0) {
      msg <- paste0(
        "The following arguments are not used: ",
        paste(unused, collapse = ', ')
      )
      switch(
        EXPR = getOption(x = "Seurat.checkdots", default = 'warn'),
        "warn" = warning(msg, call. = FALSE, immediate. = TRUE),
        "stop" = stop(msg),
        "silent" = NULL,
        stop("Invalid Seurat.checkdots option. Please choose one of warn, stop, silent")
      )
      # unused.hints <- sapply(X = unused, FUN = OldParamHints)
      # names(x = unused.hints) <- unused
      # unused.hints <- na.omit(object = unused.hints)
      # if (length(x = unused.hints) > 0) {
      #   message(
      #     "Suggested parameter: ",
      #     paste(unused.hints, "instead of", names(x = unused.hints), collapse = '; '),
      #     "\n"
      #   )
      # }
    }
  }
  return(invisible(x = NULL))
}

#' Check a list of objects for duplicate cell names
#'
#' @param object.list List of Seurat objects
#' @param verbose Print message about renaming
#' @param stop Error out if any duplicate names exist
#'
#' @return Returns list of objects with duplicate cells renamed to be unique
#'
#' @keywords internal
#'
#' @noRd
#'
CheckDuplicateCellNames <- function(object.list, verbose = TRUE, stop = FALSE) {
  cell.names <- unlist(x = lapply(X = object.list, FUN = colnames))
  if (any(duplicated(x = cell.names))) {
    if (stop) {
      stop("Duplicate cell names present across objects provided.")
    }
    if (verbose) {
      warning("Some cell names are duplicated across objects provided. Renaming to enforce unique cell names.")
    }
    object.list <- lapply(
      X = 1:length(x = object.list),
      FUN = function(x) {
        return(RenameCells(
          object = object.list[[x]],
          new.names = paste0(Cells(x = object.list[[x]]), "_", x)
        ))
      }
    )
  }
  return(object.list)
}

#' Empty Data Frames
#'
#' Create an empty \link[base:data.frame]{data frame} with no row names and
#' zero columns
#'
#' @param n Number of rows for the data frame
#'
#' @return A \link[base:data.frame]{data frame} with \code{n} rows and
#' zero columns
#'
#' @keywords internal
#'
#' @noRd
#'
EmptyDF <- function(n) {
  return(as.data.frame(x = matrix(nrow = n, ncol = 0L)))
}

#' Extract delimiter information from a string.
#'
#' Parses a string (usually a cell name) and extracts fields based
#'  on a delimiter
#'
#' @param string String to parse.
#' @param field Integer(s) indicating which field(s) to extract. Can be a
#' vector multiple numbers.
#' @param delim Delimiter to use, set to underscore by default.
#'
#' @return A new string, that parses out the requested fields, and
#' (if multiple), rejoins them with the same delimiter
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' \donttest{
#' SeuratObject:::ExtractField('Hello World', field = 1, delim = '_')
#' }
#'
ExtractField <- function(string, field = 1, delim = "_") {
  fields <- as.numeric(x = unlist(x = strsplit(
    x = as.character(x = field),
    split = ","
  )))
  if (length(x = fields) == 1) {
    return(strsplit(x = string, split = delim)[[1]][field])
  }
  return(paste(
    strsplit(x = string, split = delim)[[1]][fields],
    collapse = delim
  ))
}

#' Check List Names
#'
#' Check to see if a list has names; also check to enforce that all names are
#' present and unique
#'
#' @param x A list
#' @param all.unique Require that all names are unique from one another
#' @param allow.empty Allow empty (\code{nchar = 0}) names
#' @param pass.zero Pass on zero-length lists
#'
#' @return \code{TRUE} if ..., otherwise \code{FALSE}
#'
#' @importFrom rlang is_bare_list
#'
#' @keywords internal
#'
#' @noRd
#'
IsNamedList <- function(
  x,
  all.unique = TRUE,
  allow.empty = FALSE,
  pass.zero = FALSE
) {
  if (!is_bare_list(x = x)) {
    return(FALSE)
  }
  if (isTRUE(x = pass.zero) && !length(x = x)) {
    return(TRUE)
  }
  n <- names(x = x)
  named <- !is.null(x = n)
  if (!isTRUE(x = allow.empty)) {
    named <- named && all(vapply(
      X = n,
      FUN = nchar,
      FUN.VALUE = integer(length = 1L)
    ))
  }
  if (isTRUE(x = all.unique)) {
    named <- named && (length(x = n) == length(x = unique(x = n)))
  }
  return(named)
}

#' Test Null Pointers
#'
#' Check to see if a C++ pointer is a null pointer on the compiled side
#'
#' @param x An \link[methods:externalptr-class]{external pointer} object
#'
#' @return \code{TRUE} if \code{x} is a null pointer, otherwise \code{FALSE}
#'
#' @importFrom methods is
#'
#' @references \url{https://stackoverflow.com/questions/26666614/how-do-i-check-if-an-externalptr-is-null-from-within-r}
#'
#' @keywords internal
#'
#' @noRd
#'
IsNullPtr <- function(x) {
  stopifnot(is(object = x, class2 = 'externalptr'))
  return(.Call('isnull', x))
}

#' Test Empty Characters
#'
#' Check to see if a \code{\link[base]{character}} vector is empty. A character
#' is empty if it has no length or an \code{nchar == 0}
#'
#' @param x A \code{\link[base]{character}} vector
#' @param mode Stringency of emptiness test:
#' \describe{
#'  \item{\dQuote{each}}{Return a single value for each member of \code{x}}
#'  \item{\dQuote{any}}{Return \code{TRUE} if any member of \code{x} is empty}
#'  \item{\dQuote{all}}{Return \code{TRUE} if \emph{every} member of \code{x} is
#'  empty}
#' }
#' @param na Control how \code{\link[base]{NA}} values are treated:
#' \describe{
#'  \item{\dQuote{empty}}{Treat \code{NA}s as empty values}
#'  \item{\dQuote{keep}}{Keep \code{NA} values and treat them as \code{NA}}
#'  \item{\dQuote{remove}}{Remove \code{NA} values before testing emptiness}
#' }
#'
#' @return If \code{mode} is \dQuote{each}, a vector of logical values denoting
#' the emptiness of of each member of \code{x}; otherwise, a singular
#' \code{\link[base]{logical}} denoting the overall emptiness of \code{x}
#'
#' @keywords internal
#'
#' @noRd
#'
IsCharEmpty <- function(
  x,
  mode = c('each', 'any', 'all'),
  na = c('empty', 'keep', 'remove')
) {
  if (!is.character(x = x)) {
    return(FALSE)
  }
  mode <- mode[1]
  mode <- match.arg(arg = mode)
  na <- na[1]
  na <- match.arg(arg = na)
  switch(
    EXPR = na,
    'empty' = {
      x[is.na(x = x)] <- ''
    },
    'remove' = {
      x <- x[!is.na(x = x)]
    }
  )
  if (!length(x = x)) {
    return(TRUE)
  }
  empty <- vapply(
    X = x,
    FUN = Negate(f = nchar),
    FUN.VALUE = logical(length = 1L),
    USE.NAMES = FALSE
  )
  empty <- switch(
    EXPR = mode,
    'any' = any(empty),
    'all' = all(empty),
    empty
  )
  return(empty)
}

#' Update a Class's Package
#'
#' Swap packages for an object's class definition. As classes move between
#' packages, these functions rescope the namespace of the S4 class. This allows
#' objects to depend only on the new package for class definitions rather than
#' both the new and old packages
#'
#' @inheritParams s4list
#' @param from A vector of one or more packages to limit conversion from
#' @param to A character naming the package to search for new class definitions;
#' defaults to the package of the function calling this function
#'
#' @return \code{SwapClassPkg}: \code{x} with an updated S4 class
#' definition attribute
#'
#' @inheritSection s4list S4 Class Definition Attributes
#'
#' @name classpkg
#' @rdname classpkg
#'
#' @keywords internal
#'
#' @seealso \code{\link{s4list}}
#'
#' @noRd
#'
SwapClassPkg <- function(x, from = NULL, to = NULL) {
  if (!inherits(x = x, what = 'list')) {
    return(x)
  }
  to <- to[1] %||% environmentName(env = environment(
    fun = sys.function(which = 1L)
  ))
  if (!nchar(x = to) || !paste0('package:', to) %in% search()) {
    to <- environmentName(env = environment(fun = sys.function(which = 0L)))
  }
  for (i in seq_along(along.with = x)) {
    if (!is.null(x = x[[i]])) {
      x[[i]] <- SwapClassPkg(x = x[[i]], from = from, to = to)
    }
  }
  if (!IsS4List(x = x)) {
    return(x)
  }
  classdef <- unlist(x = strsplit(
    x = attr(x = x, which = 'classDef'),
    split = ':'
  ))
  pkg <- classdef[1]
  cls <- classdef[2]
  if (is.null(x = from) || pkg %in% from) {
    pkg <- ifelse(
      test = is.null(x = getClassDef(
        Class = cls,
        package = to,
        inherits = FALSE
      )),
      yes = pkg,
      no = to
    )
  }
  attr(x = x, which = 'classDef') <- paste(pkg, cls, sep = ':')
  return(x)
}

#' Get the top
#'
#' @param data Data to pull the top from
#' @param num Pull top \code{num}
#' @param balanced Pull even amounts of from positive and negative values
#'
#' @return The top \code{num}
#'
#' @importFrom utils head tail
#'
#' @keywords internal
#'
#' @noRd
#'
Top <- function(data, num = 20, balanced = FALSE) {
  nr <- nrow(x = data)
  if (num > nr) {
    warning(
      "Requested number is larger than the number of available items (",
      nr,
      "). Setting to ",
      nr ,
      ".",
      call. = FALSE
    )
    num <- nr
  }
  balanced <- ifelse(test = nr == 1, yes = FALSE, no = balanced)
  top <- if (isTRUE(x = balanced)) {
    num <- round(x = num / 2)
    data <- data[order(data, decreasing = TRUE), , drop = FALSE]
    positive <- head(x = rownames(x = data), n = num)
    negative <- rev(x = tail(x = rownames(x = data), n = num))
    # remove duplicates
    if (positive[num] == negative[num]) {
      negative <- negative[-num]
    }
    list(positive = positive, negative = negative)
  } else {
    data <- data[rev(x = order(abs(x = data))), , drop = FALSE]
    top <- head(x = rownames(x = data), n = num)
    top[order(data[top, ])]
  }
  return(top)
}

#' @name classpkg
#' @rdname classpkg
#'
#' @return \code{UpdateClassPkg}: \code{object} with the updated
#' class definition
#'
#' @keywords internal
#'
#' @noRd
#'
UpdateClassPkg <- function(object, from = NULL, to = NULL) {
  if (!isS4(object)) {
    return(object)
  }
  obj.list <- S4ToList(object = object)
  obj.list <- SwapClassPkg(x = obj.list, from = from, to = to)
  # browser()
  return(ListToS4(x = obj.list))
}

#' Update slots in an object
#'
#' @param object An object to update
#'
#' @return \code{object} with the latest slot definitions
#'
#' @importFrom methods slotNames slot
#'
#' @keywords internal
#'
#' @noRd
#'
UpdateSlots <- function(object) {
  object.list <- sapply(
    X = slotNames(x = object),
    FUN = function(x) {
      return(tryCatch(
        expr = slot(object = object, name = x),
        error = function(...) {
          return(NULL)
        }
      ))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  object.list <- Filter(f = Negate(f = is.null), x = object.list)
  object.list <- c('Class' = class(x = object)[1], object.list)
  object <- do.call(what = 'new', args = object.list)
  for (x in setdiff(x = slotNames(x = object), y = names(x = object.list))) {
    xobj <- slot(object = object, name = x)
    if (is.vector(x = xobj) && !is.list(x = xobj) && length(x = xobj) == 0) {
      slot(object = object, name = x) <- vector(
        mode = class(x = xobj),
        length = 1L
      )
    }
  }
  return(object)
}
