#' @importFrom Rcpp evalCpp
#' @useDynLib SeuratObject
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Cast to Sparse
#'
#' Convert dense objects to sparse representations
#'
#' @param x An object
#' @param ... Arguments passed to other methods
#'
#' @return A sparse representation of the input data
#'
#' @rdname as.sparse
#' @export as.sparse
#'
#' @concept utils
#'
as.sparse <- function(x, ...) {
  UseMethod(generic = 'as.sparse', object = x)
}

#' Check if a matrix is empty
#'
#' Takes a matrix and asks if it's empty (either 0x0 or 1x1 with a value of NA)
#'
#' @param x A matrix
#'
#' @return Whether or not \code{x} is empty
#'
#' @rdname IsMatrixEmpty
#' @export IsMatrixEmpty
#'
#' @concept utils
#'
#' @examples
#' IsMatrixEmpty(new("matrix"))
#' IsMatrixEmpty(matrix())
#' IsMatrixEmpty(matrix(1:3))
#'
IsMatrixEmpty <- function(x) {
  UseMethod(generic = 'IsMatrixEmpty', object = x)
}

#' S4/List Conversion
#'
#' Convert S4 objects to lists and vice versa. Useful for declassing an S4
#' object while keeping track of it's class using attributes (see section
#' \strong{S4 Class Definition Attributes} below for more details). Both
#' \code{ListToS4} and \code{S4ToList} are recursive functions, affecting all
#' lists/S4 objects contained as sub-lists/sub-objects.
#'
#' @param x A list with an S4 class definition attribute
#' @param object An S4 object
#'
#' @return \code{S4ToList}: A list with an S4 class definition attribute
#'
#' @section S4 Class Definition Attributes:
#' S4 classes are scoped to the package and class name. In order to properly
#' track which class a list is generated from in order to build a new one,
#' these function use an \code{\link[base:attr]{attribute}} to denote the
#' class name and package of origin. This attribute is stored as
#' \dQuote{classDef} and takes the form of \dQuote{\code{package:class}}.
#'
#' @name s4list
#' @rdname s4list
#'
#' @concept utils
#'
#' @export
#'
S4ToList <- function(object) {
  if (!(isS4(object) || inherits(x = object, what = 'list'))) {
    return(object)
  }
  UseMethod(generic = 'S4ToList', object = object)
}

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

#' Check the use of dots
#'
#' @param ... Arguments passed to a function that fall under \code{...}
#' @param fxns A list/vector of functions or function names
#'
#' @return ...
#'
#' @importFrom utils isS3stdGeneric methods argsAnywhere isS3method
#'
#' @keywords internal
#'
#' @noRd
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
        EXPR = getOption(x = "Seurat.checkdots"),
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

#' Generate a random name
#'
#' Make a name from randomly sampled lowercase letters, pasted together with no
#' spaces or other characters
#'
#' @param length How long should the name be
#' @param ... Extra parameters passed to \code{\link[base]{sample}}
#'
#' @return A character with \code{nchar == length} of randomly sampled letters
#'
#' @seealso \code{\link[base]{sample}}
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#' set.seed(42L)
#' SeuratObject:::RandomName()
#' SeuratObject:::RandomName(7L, replace = TRUE)
#' }
#'
RandomName <- function(length = 5L, ...) {
  CheckDots(..., fxns = 'sample')
  return(paste(sample(x = letters, size = length, ...), collapse = ''))
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

#' Update a Key
#'
#' @param key A character to become a Seurat Key
#'
#' @return An updated Key that's valid for Seurat
#'
#' @section \code{Seurat} Object Keys:
#' blah
#'
#' @keywords internal
#'
#' @noRd
#'
UpdateKey <- function(key) {
  if (grepl(pattern = '^[[:alnum:]]+_$', x = key)) {
    return(key)
  } else {
    new.key <- regmatches(
      x = key,
      m = gregexpr(pattern = '[[:alnum:]]+', text = key)
    )
    new.key <- paste0(paste(unlist(x = new.key), collapse = ''), '_')
    if (new.key == '_') {
      new.key <- paste0(RandomName(length = 3), '_')
    }
    warning(
      "Keys should be one or more alphanumeric characters followed by an underscore, setting key from ",
      key,
      " to ",
      new.key,
      call. = FALSE,
      immediate. = TRUE
    )
    return(new.key)
  }
}
