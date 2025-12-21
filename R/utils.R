#' @include zzz.R
#' @include generics.R
#' @include centroids.R
#' @include segmentation.R
#' @importFrom Rcpp evalCpp
#' @importFrom methods as setAs
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Set If or If Not \code{NULL}
#'
#' Set a default value depending on if an object is \code{NULL}
#'
#' @usage x \%||\% y
#'
#' @param x An object to test
#' @param y A default value
#'
#' @return For \code{\%||\%}: \code{y} if \code{x} is \code{NULL};
#' otherwise \code{x}
#'
#' @name set-if-null
#' @rdname set-if-null
#'
#' @author For \code{\%||\%}: \pkg{rlang} developers
#'
#' @seealso \code{\link[rlang:op-null-default]{rlang::\%||\%}}
#'
#' @aliases %||%
#'
#' @concept utils
#'
#' @examples
#' # Set if NULL
#' 1 %||% 2
#' NULL %||% 2
#'
NULL

#' @importFrom rlang %||%
#' @export
#'
#' @noRd
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
#' # Set if *not* NULL
#' 1 %iff% 2
#' NULL %iff% 2
#'
`%iff%` <- function(x, y) {
  if (!is_null(x = x)) {
    return(y)
  }
  return(x)
}

#' Set If or If Not \code{NA}
#'
#' Set a default value depending on if an object is \code{\link[base]{NA}}
#'
#' @inheritParams set-if-null
#'
#' @return For \code{\%NA\%}: \code{y} if \code{x} is \code{\link[base]{NA}};
#' otherwise \code{x}
#'
#' @name set-if-na
#' @rdname set-if-na
#'
#' @importFrom rlang is_na
#'
#' @keywords internal
#'
#' @export
#'
#' @concept utils
#'
#' @examples
#' # Set if NA
#' 1 %NA% 2
#' NA %NA% 2
#'
`%NA%` <- function(x, y) {
  if (is_na(x = x)) {
    return(y)
  }
  return(x)
}

#' @rdname set-if-na
#'
#' @export
#'
`%na%` <- `%NA%`

#' @return For \code{\%!NA\%}: \code{y} if \code{x} is \strong{not}
#' \code{\link[base]{NA}}; otherwise \code{x}
#'
#' @rdname set-if-na
#'
#' @importFrom rlang is_na
#'
#' @export
#'
#' @examples
#' # Set if *not* NA
#' 1 %!NA% 2
#' NA %!NA% 2
#'
`%!NA%` <- function(x, y) {
  if (is_na(x = x)) {
    return(x)
  }
  return(y)
}

#' @rdname set-if-na
#'
#' @export
#'
`%!na%` <- `%!NA%`

#' \pkg{BPCells} Matrix Mode
#'
#' Get the mode (on-disk, in-memory) of an \code{IterableMatrix} object
#' from \pkg{BPCells}
#'
#' @param object An \code{IterableMatrix}
#' @param simplify Return \dQuote{\code{disk}} for on-disk matrices
#'
#' @return One of the following, depending on the mode of \code{object}:
#' \itemize{
#'  \item \dQuote{\code{memory}}
#'  \item \dQuote{\code{file}}
#'  \item \dQuote{\code{directory}}
#' }
#' If \code{simplify} is \code{TRUE}, returns \dQuote{\code{disk}} instead of
#' \dQuote{\code{file}} or \dQuote{\code{directory}}
#'
#' @keywords internal
#'
#' @export
#'
.BPMatrixMode <- function(object, simplify = FALSE) {
  check_installed(pkg = 'BPCells', reason = 'for working with BPCells')
  if (!inherits(x = object, what = 'IterableMatrix')) {
    return(NULL)
  }
  stopifnot(rlang::is_bare_logical(x = simplify, n = 1L))
  # Get a vector of all the slots in all sub-matrices
  slots <- Reduce(
    f = union,
    x = lapply(
      X = BPCells::all_matrix_inputs(object),
      FUN = \(x) methods::slotNames(x = methods::getClass(Class = class(x = x)))
    )
  )
  # Figure out if any sub-matrix points to a directory or a file path
  type <- c(path = FALSE, dir = FALSE)
  for (s in slots) {
    if (s %in% names(x = type)) {
      type[s] <- TRUE
    }
  }
  # If no matrix points to a directory or file, it's an in-memory one
  if (!any(type)) {
    return('memory')
  }
  # If any matrix points to a directory or file, it's an on-disk matrix
  if (isTRUE(x = simplify) && any(type)) {
    return("disk")
  }
  # Get the exact type; there should only be one
  return(c(path = 'file', dir = 'directory')[[names(x = type)[type]]])
}

#' Identify Object Collections
#'
#' Find all collection (named lists) slots in an S4 object
#'
#' @inheritParams .Contains
#' @param exclude A character vector of slot names to exclude
#' @param ... Arguments passed to \code{\link{IsNamedList}}
#'
#' @return A character vector of names of collection slots
#'
#' @importFrom methods slotNames
#'
#' @keywords internal
#'
#' @export
#'
#' @family subobjects
#' @concept utils
#'
#' @examples
#' .Collections(pbmc_small)
#'
.Collections <- function(object, exclude = character(length = 0L), ...) {
  if (!isS4(object)) {
    abort(message = "'object' is not an S4 object")
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
#' @importFrom methods getClass slot
#'
#' @keywords internal
#'
#' @export
#'
#' @concept utils
#'
#' @examples
#' .Contains(pbmc_small)
#'
.Contains <- function(object) {
  if (!isS4(object)) {
    abort(message = "'object' not an S4 object")
  }
  return(names(x = slot(
    object = getClass(Class = class(x = object)),
    name = 'contains'
  )))
}

#' Find the Default FOV
#'
#' Attempts to find the \dQuote{default} FOV using the revamped
#' spatial framework
#'
#' @param object A \code{{Seurat}} object
#'
#' @return ...
#'
#' @keywords internal
#'
#' @export
#'
#' @concept utils
#'
.DefaultFOV <- function(object, assay = NULL) {
  images <- .FilterObjects(object = object, classes.keep = 'FOV')
  if (!is.null(x = assay)) {
    assays <- c(assay, DefaultAssay(object = object[[assay]]))
    images <- Filter(
      f = function(x) {
        return(DefaultAssay(object = object[[x]]) %in% assays)
      },
      x = images
    )
  }
  if (!length(x = images)) {
    return(NULL)
  }
  return(images)
}

#' Deprecate Functions and Arguments
#'
#' Provides automatic deprecation and defunctation of functions and arguments;
#'
#' @inheritParams lifecycle::deprecate_soft
#' @inheritDotParams lifecycle::deprecate_soft
#' @param pkg Name of package to use for comparison
#' @param env,user_env Managed internally by \code{.Deprecate()}
#'
#' @return Run for its side effect and invisibly returns \code{NULL}
#'
#' @importFrom rlang ns_env_name
#' @importFrom utils packageVersion
#' @importFrom lifecycle deprecate_soft deprecate_stop deprecate_warn
#'
#' @keywords internal
#'
#' @export
#'
#' @seealso \code{\link[lifecycle:deprecate_soft]{lifecycle::deprecate_soft}()}
#' \code{\link[lifecycle:deprecate_warn]{lifecycle::deprecate_warn}()}
#' \code{\link[lifecycle:deprecate_stop]{lifecycle::deprecate_stop}()}
#'
.Deprecate <- function(
  when,
  what,
  with = NULL,
  ...,
  pkg = NULL,
  env = missing_arg(),
  user_env = missing_arg()
) {
  # Figure out current version, rounding up development versions
  caller <- caller_env()
  current <- .RoundVersion(current = packageVersion(
    pkg = ns_env_name(x = caller)
  ))
  cv <- paste(current, collapse = '.')
  # Ensure our 'when' is a valid version
  wv <- when <- as.character(x = numeric_version(x = when, strict = TRUE))
  # If we haven't reached deprecation, exit out silently
  if (cv < wv) {
    return(invisible(x = NULL))
  }
  # Figure out if this is a soft deprecation, a warning deprecation, or a defunct
  when <- unlist(x = strsplit(x = when, split = '\\.'))
  if (length(x = when) > 4L) {
    when[4L] <- paste(
      when[seq.int(from = 4L, to = length(x = when))],
      collapse = '.'
    )
    when <- when[1:4]
  }
  names(x = when) <- c('major', 'minor', 'patch', 'devel')[seq_along(along.with = when)]
  when <- vapply(
    X = when,
    FUN = as.integer,
    FUN.VALUE = integer(length = 1L),
    USE.NAMES = TRUE
  )
  diffs <- abs(current - when)
  if (diffs['major'] >= 1L || diffs['minor'] >= 3L) {
    deprecate_stop(
      when = wv,
      what = what,
      with = with,
      env = caller,
      ...
    )
  }
  fn <- if (diffs['minor'] >= 1L) {
    deprecate_warn
  } else {
    deprecate_soft
  }
  fn(
    when = wv,
    what = what,
    with = with,
    env = caller,
    user_env = caller_env(n = 2L),
    ...
  )
  return(invisible(x = NULL))
}

#' Find Subobjects Of A Certain Class
#'
#' @inheritParams .Collections
#' @param classes.keep A vector of classes to keep
#'
#' @return A vector of object names that are of class \code{classes.keep}
#'
#' @keywords internal
#'
#' @export
#'
#' @family subobjects
#' @concept utils
#'
#' @examples
#' .FilterObjects(pbmc_small)
#' .FilterObjects(pbmc_small, "Graph")
#'
.FilterObjects <- function(
  object,
  classes.keep = c('Assay', 'StdAssay', 'DimReduc')
) {
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

#' Find A Subobject
#'
#' Determine the slot that a subobject is contained in
#'
#' @inheritParams .Collections
#' @param name Name of subobject to find
#'
#' @return The name of the slot that contains \code{name}; returns \code{NULL}
#' if a subobject named \code{name} cannot be found
#'
#' @keywords internal
#'
#' @export
#'
#' @family subobjects
#' @concept utils
#'
#' @examples
#' .FindObject(pbmc_small, "tsne")
#'
.FindObject <- function(object, name, exclude = c('misc', 'tools')) {
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

#' Get a Method
#'
#' @param fxn Name of a function as a character
#' @param cls The class to find a method of \code{fxn} for
#'
#' @return The method of \code{fxn} for class \code{cls}; if no method found,
#' returns the default method. If no default method found; returns \code{NULL}
#'
#' @importFrom utils getS3method isS3stdGeneric
#' @importFrom methods isClass isGeneric selectMethod
#'
#' @keywords internal
#'
#' @export
#'
#' @concept utils
#'
#' @examples
#' .GetMethod('t', 'Matrix')
#' .GetMethod('t', 'data.frame')
#'
.GetMethod <- function(fxn, cls) {
  if (is.function(x = fxn)) {
    fxn <- as.character(x = substitute(expr = fxn))
  }
  if (!(isS3stdGeneric(f = fxn) || isGeneric(f = fxn))) {
    abort(message = paste0("'", fxn, "' is not a generic function"))
  }
  default <- NULL
  if (isGeneric(f = fxn) && isClass(Class = cls[1L])) {
    method <- selectMethod(f = fxn, signature = cls)
    if (!inherits(x = method, what = 'derivedDefaultMethod')) {
      return(slot(object = method, name = '.Data'))
    }
    default <- slot(object = method, name = '.Data')
  }
  method <- NULL
  for (i in c(cls, 'default')) {
    method <- getS3method(f = fxn, class = i, optional = TRUE)
    if (!is.null(x = method)) {
      break
    }
  }
  method <- method %||% default
  if (is.null(x = method)) {
    abort(message = paste0(
      "Unable to find a method for '",
      fxn,
      "' for '",
      cls[1L],
      "' objects"
    ))
  }
  return(method)
}

#' Propagate a List
#'
#' @param x A list or character vector
#' @param names A vector of names to keep from \code{x}
#' @param default A default value for unassigned values of \code{x}
#'
#' @return A named list where the names are present in both \code{x} and
#' \code{names} and the values are either the values from \code{x} or
#' \code{default}
#'
#' @keywords internal
#'
#' @export
#'
#' @concept utils
#'
#' @examples
#' .PropagateList("counts", c("RNA", "ADT", "SCT"))
#' .PropagateList(c("counts", "data"), c("RNA", "ADT", "SCT"))
#' .PropagateList("ADT", c("RNA", "ADT", "SCT"))
#' .PropagateList(c("RNA", "SCT"), c("RNA", "ADT", "SCT"))
#' .PropagateList(c("RNA", ADT = "counts"), c("RNA", "ADT", "SCT"))
#' .PropagateList(list(SCT = c("counts", "data"), ADT = "counts"), c("RNA", "ADT", "SCT"))
#' .PropagateList(list(SCT = c("counts", "data"), "ADT"), c("RNA", "ADT", "SCT"))
#'
.PropagateList <- function(x, names, default = NA) {
  # `names` must be a character vector
  if (!is_bare_character(x = names)) {
    abort(message = "'names' must be a character vector")
  }
  # `x` must be a list or character vector
  if (!(is_bare_list(x = x) || is_bare_character(x = x))) {
    abort(message = "'x' must be either a list or character vector")
  }
  # `x` cannot be empty
  if (!length(x = x)) {
    abort(message = "'x' cannot be empty")
  }
  # `x` is a character vector
  if (is_bare_character(x = x)) {
    if (!all(nzchar(x = x))) {
      abort(message = "'x' cannot be empty")
    }
    # Handle cases where `x` is unnamed
    if (!any(have_name(x = x))) {
      # `x` is a vector with values in `names`
      # Return a list for every value in `x` that's present in `names`
      # with a value of `default`
      if (any(x %in% names)) {
        x <- intersect(x = x, y = names)
        ret <- vector(mode = 'list', length = length(x = x))
        names(x = ret) <- x
        for (i in seq_along(along.with = ret)) {
          ret[[i]] <- default
        }
        return(ret)
      }
      # `x` is a vector of default values
      # Return a list for every value in `names` with a value of `x`
      ret <- vector(mode = 'list', length = length(x = names))
      names(x = ret) <- names
      for (i in seq_along(along.with = ret)) {
        ret[[i]] <- unique(x = x)
      }
      return(ret)
    }
    # `x` is named
    # Turn `x` into a list and continue on
    x <- as.list(x = x)
  }
  # `x` is a list
  # Find entries of `x` that correspond to a value in `names`
  # Assign new value of `default`
  for (i in seq_along(along.with = x)) {
    if (is_scalar_character(x = x[[i]]) && x[[i]] %in% names) {
      names(x = x)[i] <- x[[i]]
      x[[i]] <- default
    }
  }
  # Identify values of `x` in `names`
  x.use <- intersect(x = names(x = x), y = names)
  if (!length(x = x.use) && is_named(x = x)) {
    abort(message = "None of the values of 'x' match with 'names")
  }
  #`Return only values of `x` that are in `names``
  return(x[x.use])
}

#' Get the Subobject Names
#'
#' @inheritParams .Collections
#' @param collapse Collapse the list into a vector
#'
#' @return If \code{collapse = TRUE}, then a vector with the names of all
#' subobjects; otherwise, a named list where the names are the names of the
#' collections and the values are the names of subobjects within the collection
#'
#' @keywords internal
#'
#' @export
#'
#' @family subobjects
#' @keywords utils
#'
#' @examples
#' .Subobjects(pbmc_small)
#'
.Subobjects <- function(
  object,
  exclude = c('misc', 'tools'),
  collapse = TRUE,
  ...
) {
  subobjects <- sapply(
    X = .Collections(object = object, exclude = exclude, ...),
    FUN = function(x) {
      return(names(x = slot(object = object, name = x)))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  if (isTRUE(x = collapse)) {
    subobjects <- unlist(x = subobjects, use.names = FALSE)
  }
  return(subobjects)
}

#' Attach Required Packages
#'
#' Helper function to attach required packages. Detects if a package is already
#' attached and if so, skips it. Should be called in \code{\link[base]{.onAttach}}
#'
#' @param deps A character vector of packages to attach
#'
#' @template return-null
#'
#' @export
#'
#' @concept utils
#'
#' @template lifecycle-superseded
#' @section Lifecycle:
#' \code{AttachDeps} has been superseded as of \pkg{SeuratObject} v5.0.0;
#' as an alternative, list dependencies in the \code{Depends} section of
#' \code{DESCRIPTION}
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
#' @concept utils
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
    abort(message = "No named arguments passed")
  }
  if (length(x = fxns) == 1) {
    fxns <- list(fxns)
  }
  for (f in fxns) {
    if (!(is.character(x = f) || is.function(x = f))) {
      abort(message = paste(
        "CheckDots only works on characters or functions, not",
        class(x = f)[1L]
      ))
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

#' Check features names format
#'
#' @param data a matrix input, rownames(data) are feature names
#'
#' @return \code{data} with update feature names
#'
#' @keywords internal
#'
#' @export
#'
CheckFeaturesNames <- function(data) {
  if (any(grepl(pattern = "_", x = rownames(x = data)))) {
    warning(
      "Feature names cannot have underscores ('_'), replacing with dashes ('-')",
      call. = FALSE,
      immediate. = TRUE
    )
    rownames(x = data) <- gsub(
      pattern = "_",
      replacement = "-",
      x = rownames(x = data)
    )
  }
  if (any(grepl(pattern = "|", x = rownames(x = data), fixed = TRUE))) {
    warning(
      "Feature names cannot have pipe characters ('|'), replacing with dashes ('-')",
      call. = FALSE,
      immediate. = TRUE
    )
    rownames(x = data) <- gsub(
      pattern = "|",
      replacement = "-",
      x = rownames(x = data),
      fixed = TRUE
    )
  }
  return(data)
}

#' Conditional Garbage Collection
#'
#' Call \code{gc} only when desired
#'
#' @param option ...
#'
#' @template return-null
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

#' Check layers names for the input list
#'
#'
#' @param matrix.list A list of matrices
#' @param layers.type layers type, such as counts or data
#'
#'
#' @export
#'
#' @concept utils
#'
CheckLayersName <- function(
  matrix.list,
  layers.type = c('counts', 'data')
) {
  layers.type <- match.arg(arg = layers.type)
  if (is.null(x = matrix.list)) {
    return(matrix.list)
  }
  if (!inherits(x = matrix.list, what = 'list')) {
    matrix.list <- list(matrix.list)
  }
  if (length(x = matrix.list) == 1) {
    names(x = matrix.list) <- layers.type
  } else {
    endings <- seq_along(along.with = matrix.list)
    for (i in 1:length(x = matrix.list)) {
      name <- names(x = matrix.list)[i]
      if (!is.null(name) && nzchar(x = name)) {
        if (grepl(pattern = paste0('^', layers.type, '[._\\0-9-]+'), x = name)) {
          name <- gsub(
            pattern = paste0(layers.type, '[._\\0-9-]+'),
            replacement = "",
            x = name
          )
          # If replacement leaves empty string
          if (!nzchar(x = name)) {
            name <- i
          }
        }
        endings[i] <- name
      }
    }
    names(x = matrix.list) <- paste0(paste0(layers.type, '.'), endings)
    names(x = matrix.list) <- make.unique(names = names(x = matrix.list), sep = '')
  }
  return(matrix.list)
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
#' @keywords internal
#'
#' @export
#'
#' @concept utils
#' @family s4list
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

  # check if stored defaults and if not use original behavior
  defaults <- Tool(object, slot = "`DefaultDimReduc<-.Seurat`")
  if (!is.null(x = defaults) && assay %in% names(x = defaults)) {
    return(defaults[[assay]])
  } else {
    drs.use <- c('umap', 'tsne', 'pca')
    dim.reducs <- .FilterObjects(object = object, classes.keep = 'DimReduc')
    drs.assay <- Filter(
      f = function(x) {
        return(DefaultAssay(object = object[[x]]) == assay)
      },
      x = dim.reducs
    )
    if (length(x = drs.assay)) {
      index <- lapply(
        X = drs.use,
        FUN = grep,
        x = drs.assay,
        ignore.case = TRUE
      )
      index <- Filter(f = length, x = index)
      if (length(x = index)) {
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
    if (!length(x = index)) {
      abort(message = paste0(
        "Unable to find a DimReduc matching one of ",
        .Oxford(drs.use),
        "; please specify a dimensional reduction to use"
      ))
    }
    return(dim.reducs[min(index[[1]])])
  }
}

#' @rdname DefaultDimReduc
#' @export
#' @method DefaultDimReduc<-  Seurat
#'
#' @concept utils
#'
#' @examples
#' \dontrun{
#' # Set UMAP as default for RNA assay
#' DefaultDimReduc(seurat_obj) <- "umap"
#'
#' # Clear the set default
#' DefaultDimReduc(seurat_obj) <- NULL
#' }

"DefaultDimReduc<-.Seurat" <- function(object, ..., value) {
  assay <- DefaultAssay(object = object)
  # Get existing defaults or create empty named vector
  defaults <- Tool(object, slot = "`DefaultDimReduc<-.Seurat`")
  if (is.null(x = defaults)) {
    defaults <- character(0)
  }

  if (is.null(x = value)) {
    # Clear default for current assay
    new_defaults <- defaults[names(x = defaults) != assay]
    message(paste0('Removing the set default DimReduc for "', assay, '" assay.'))
  } else {
    # Validate that the dim reduc exists
    if (!value %in% Reductions(object = object)) {
      stop(paste0("DimReduc ", value, " not present in object."))
    }

    # Check reduction is associated with current assay
    reduc_assay <- DefaultAssay(object[[value]])
    if (assay != reduc_assay) {
      stop(paste0('The reduction "', value, '" is not associated with current assay "', assay, '". No change made to default DimReduc.'))
    }

    # Set default for current assay
    message(paste0('Setting "', value, '" as default DimReduc for "', assay, '" assay.'))
    new_defaults <- c(defaults[names(x = defaults) != assay], setNames(value, assay))
  }

  # Store back using Tool<-
  Tool(object) <- new_defaults

  return(object)
}


#' Radian/Degree Conversions
#'
#' Convert degrees to radians and vice versa
#'
#' @param rad Angle in radians
#'
#' @return \code{Degrees}: \code{rad} in degrees
#'
#' @name Angles
#' @rdname angles
#'
#' @keywords internal
#'
#' @export
#'
#' @concept utils
#' @family angles
#'
#' @examples
#' Degrees(pi)
#'
Degrees <- function(rad) {
  return(rad * (180 / pi))
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
#' @export
#'
#' @concept utils
#'
#' @examples
#' EmptyDF(4L)
#'
EmptyDF <- function(n) {
  return(as.data.frame(x = matrix(nrow = n, ncol = 0L)))
}

#' Empty Matrices
#'
#' Create empty 0x0 matrices of varying types
#'
#' @param repr Representation of empty matrix; choose from:
#' \itemize{
#'  \item \dQuote{\code{C}} for a
#'   \code{\link[Matrix:CsparseMatrix-class]{CsparseMatrix}}
#'  \item \dQuote{\code{T}} for a
#'   \code{\link[Matrix:TsparseMatrix-class]{TsparseMatrix}}
#'  \item \dQuote{\code{R}} for an
#'   \code{\link[Matrix:RsparseMatrix-class]{RsparseMatrix}}
#'  \item \dQuote{\code{e}} for an
#'   \code{\link[Matrix:unpackedMatrix-class]{unpackedMatrix}}
#'  \item \dQuote{\code{d}} for a dense S3 \code{\link[base]{matrix}}
#'  \item \dQuote{\code{spam}} for a \code{\link[spam]{spam}} matrix
#' }
#' @param type Type of resulting matrix to return, choose from:
#' \itemize{
#'  \item \dQuote{\code{d}} for numeric matrices
#'  \item \dQuote{\code{l}} for logical matrices
#'  \item \dQuote{\code{n}} for pattern matrices
#' }
#' Note, when \code{repr} is \dQuote{\code{spam}}, \code{type} must be
#' \dQuote{\code{d}}; when \code{repr} is \dQuote{\code{d}}, setting \code{type}
#' to \dQuote{\code{n}} returns a logical matrix
#'
#' @return A 0x0 matrix of the specified representation and type
#'
#' @export
#'
#' @concept utils
#'
#' @seealso \code{\link{IsMatrixEmpty}()}
#'
#' @examples
#' EmptyMatrix()
#' EmptyMatrix("spam")
#'
EmptyMatrix <- function(repr = 'C', type = 'd' ) {
  repr <- arg_match(arg = repr, values = c('C', 'T', 'R', 'e', 'd', 'spam'))
  type <- arg_match(
    arg = type,
    values = switch(
      EXPR = repr,
      spam = 'd',
      c('d', 'l', 'n')
    )
  )
  return(switch(
    EXPR = repr,
    spam = spam::spam(x = 0L, nrow = 0L, ncol = 0L),
    d = matrix(
      data = vector(
        mode = switch(EXPR = type, d = 'numeric', 'logical'),
        length = 0L
      ),
      nrow = 0L,
      ncol = 0L
    ),
    new(Class = paste0(type, 'g', repr, 'Matrix'))
  ))
}

#' Extract delimiter information from a string.
#'
#' Parses a string (usually a cell name) and extracts fields based
#' on a delimiter
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
#' @export
#'
#' @concept utils
#'
#' @examples
#' ExtractField('Hello World', field = 1, delim = '_')
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
#' @export
#'
#' @concept utils
#'
#' @examples
#' IsNamedList(list())
#' IsNamedList(list(), pass.zero = TRUE)
#' IsNamedList(list(1, 2, 3))
#' IsNamedList(list(a = 1, b = 2, c = 3))
#' IsNamedList(list(a = 1, 2, c = 3))
#' IsNamedList(list(a = 1, 2, c = 3), allow.empty = TRUE)
#' IsNamedList(list(a = 1, a = 2, a = 3))
#' IsNamedList(list(a = 1, a = 2, a = 3), all.unique = FALSE)
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

#' @name s4list
#' @rdname s4list
#'
#' @return \code{IsS4List}: \code{TRUE} if \code{x} is a list with an S4 class
#' definition attribute
#'
#' @export
#'
#' @examples
#' IsS4List(pbmc.list)
#'
IsS4List <- function(x) {
  return(
    is_bare_list(x = x) &&
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
#' @examples
#' pbmc2 <- ListToS4(pbmc.list)
#' pbmc2
#' class(pbmc2)
#' Reductions(pbmc2)
#' validObject(pbmc2)
#'
ListToS4 <- function(x) {
  if (!is_bare_list(x = x)) {
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
  pkg <- classdef[1L]
  cls <- classdef[2L]
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
#' @section Lifecycle:
#'
#' \Sexpr[stage=build,results=rd]{lifecycle::badge("deprecated")}
#'
#' \code{PackageCheck} was deprecated in version 5.0.0; please use
#' \code{\link[rlang:check_installed]{rlang::check_installed}()} instead
#'
#'
PackageCheck <- function(..., error = TRUE) {
  .Deprecate(
    when = '5.0.0',
    what = 'PackageCheck()',
    with = 'rlang::check_installed()'
  )
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

#' Polygon Vertices
#'
#' Calculate the vertices of a regular polygon given the number of sides and
#' its radius (distance from center to vertex). Also permits transforming the
#' resulting coordinates by moving the origin and altering the initial angle
#'
#' @param n Number of sides of the polygon
#' @param r Radius of the polygon
#' @param xc,yc X/Y coordinates for the center of the polygon
#' @param t1 Angle of the first vertex in degrees
#'
#' @return A \code{\link[base]{data.frame}} with \code{n} rows and two columns:
#' \describe{
#'  \item{\code{x}}{X positions of each coordinate}
#'  \item{\code{y}}{Y positions of each coordinate}
#' }
#'
#' @keywords internal
#'
#' @export
#'
#' @concept utils
#' @family angles
#'
#' @references \url{https://stackoverflow.com/questions/3436453/calculate-coordinates-of-a-regular-polygons-vertices}
#'
#' @examples
#' (coords <- PolyVtx(5, t1 = 90))
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   ggplot2::ggplot(coords, ggplot2::aes(x = x, y = y)) + ggplot2::geom_polygon()
#' }
#'
PolyVtx <- function(n, r = 1L, xc = 0L, yc = 0L, t1 = 0) {
  if (!is_bare_integerish(x = n, n = 1L, finite = TRUE)) {
    abort(message = "'n' must be a single integer")
  } else if (n < 3L) {
    abort(message = "'n' must be greater than or equal to 3")
  }
  stopifnot(
    "'r' must be a single, finite number" = is_bare_numeric(x = r, n = 1L) &&
      is.finite(x = r),
    "'xc' must be a single, finite number" = is_bare_numeric(x = xc, n = 1L) &&
      is.finite(x = xc),
    "'yc' must be a single, finite number" = is_bare_numeric(x = yc, n = 1L) &&
      is.finite(x = yc),
    "'t1' must be a single, finite number" = is_bare_numeric(x = t1, n = 1L) &&
      is.finite(x = t1)
  )
  t1 <- Radians(deg = t1)
  coords <- matrix(data = 0, nrow = n, ncol = 2)
  colnames(x = coords) <- c('x', 'y')
  for (i in seq_len(length.out = n)) {
    theta <- 2 * pi * (i - 1) / n + t1
    coords[i, ] <- c(
      xc + r * cos(x = theta),
      yc + r * sin(x = theta)
    )
  }
  return(as.data.frame(x = coords))
}

#' @param deg Angle in degrees
#'
#' @return \code{Radians}: \code{deg} in radians
#'
#' @rdname angles
#'
#' @keywords internal
#'
#' @export
#'
#' @examples
#' Radians(180)
#'
Radians <- function(deg) {
  return(deg * (pi / 180))
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

#' Improve S4 validity error messages
#'
#' Catch errors from validObject to allow for more informative error messages.
#'
#' (R's internal validation checking--including making sure the object has all slots in the class definition--
#' occurs before the custom validity method for a class is run. Errors originating from an internal check may be
#' confusing to users, hence they can be modified here to provide more helpful messages.)
#'
#' @keywords internal
#' @noRd
#'
safeValidityCheck <- function(object) {
  tryCatch(
    expr = {
      validObject(object = object)
    },
    error = function(e) {
      if (grepl(pattern = "slots in class definition but not in object", x = e$message)) {
        e <- simpleError(message = paste0("Consider running UpdateSeuratObject; ", conditionMessage(e)),
                        call = conditionCall(e))
      }
      stop(e)
    }
  )
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname dot-AssayClass
#' @method .AssayClass default
#' @export
#'
.AssayClass.default <- function(object) {
  return(class(x = object)[1L])
}

#' @importFrom methods getClass
#'
#' @rdname dot-ClassPkg
#'
#' @method .ClassPkg default
#' @export
#'
.ClassPkg.default <- function(object) {
  if (!isS4(object)) {
    return(NA_character_)
  }
  return(slot(object = getClass(Class = class(x = object)), name = 'package'))
}

#' @rdname dot-ClassPkg
#' @method .ClassPkg DelayedArray
#' @export
#'
.ClassPkg.DelayedArray <- function(object) {
  check_installed(
    pkg = 'DelayedArray',
    reason = 'for working with delayed arrays'
  )
  return(.ClassPkg(object = DelayedArray::seed(x = object)))
}

#' @rdname dot-ClassPkg
#' @method .ClassPkg R6
#' @export
#'
.ClassPkg.R6 <- function(object) {
  for (cls in class(x = object)) {
    x <- eval(expr = as.symbol(x = cls))
    if (inherits(x = x, what = 'R6ClassGenerator')) {
      return(.ClassPkg(object = x))
    }
  }
  warn(message = "No r6")
  return('R6')
}

#' @rdname dot-ClassPkg
#' @method .ClassPkg R6ClassGenerator
#' @export
#'
.ClassPkg.R6ClassGenerator <- function(object) {
  return(environmentName(env = object$parent_env))
}

#' @rdname dot-DiskLoad
#' @method .DiskLoad default
#' @export
#'
.DiskLoad.default <- function(x) {
  return(NULL)
}

#' @rdname dot-DiskLoad
#' @method .DiskLoad 10xMatrixH5
#' @export
#'
.DiskLoad.10xMatrixH5 <- function(x) {
  abort(message = "Unable to determine the feature type of 10x-based BPCells matrices")
  check_installed(
    pkg = 'BPCells',
    reason = 'for working with BPCells matrices'
  )
  f <- paste(
    'function(x)',
    'BPCells::open_matrix_10x_hdf5(path = x, feature_type =',
    sQuote(x = '', q = FALSE),
    ')'
  )
  return(f)
}

#' @rdname dot-DiskLoad
#' @method .DiskLoad AnnDataMatrixH5
#' @export
#'
.DiskLoad.AnnDataMatrixH5 <- function(x) {
  check_installed(
    pkg = 'BPCells',
    reason = 'for working with BPCells matrices'
  )
  f <- paste(
    'function(x)',
    'BPCells::open_matrix_anndata_hdf5(path = x, group =',
    sQuote(x = slot(object = x, name = 'group'), q = FALSE),
    ')'
  )
  return(f)
}

#' @rdname dot-DiskLoad
#' @method .DiskLoad DelayedMatrix
#' @export
#'
.DiskLoad.DelayedMatrix <- function(x) {
  check_installed(
    pkg = 'DelayedArray',
    reason = 'for working with delayed matrices'
  )
  seed <- DelayedArray::seed(x = x)
  return(.DiskLoad(x = DelayedArray::DelayedArray(seed = seed)))
}

#' @rdname dot-DiskLoad
#' @method .DiskLoad H5ADMatrix
#' @export
#'
.DiskLoad.H5ADMatrix <- function(x) {
  check_installed(
    pkg = 'HDF5Array',
    reason = 'for working with H5AD matrices'
  )
  sparse <- DelayedArray::is_sparse(x = x)
  layer <- if (isTRUE(x = sparse)) {
    slot(object = DelayedArray::seed(x = x), name = 'group')
  } else {
    slot(object = DelayedArray::seed(x = x), name = 'name')
  }
  layer <- if (layer == '/X') {
    NULL
  } else {
    basename(path = layer)
  }
  f <- paste(
    "function(x)",
    "HDF5Array::H5ADMatrix(filepath = x",
    if (!is.null(x = layer)) {
      paste(", layer =", sQuote(x = layer, q = FALSE))
    },
    ")"
  )
  return(f)
}

#' @rdname dot-DiskLoad
#' @method .DiskLoad HDF5Matrix
#' @export
#'
.DiskLoad.HDF5Matrix <- function(x) {
  check_installed(
    pkg = 'HDF5Array',
    reason = 'for working with HDF5 matrices'
  )
  sparse <- DelayedArray::is_sparse(x = x)
  name <- slot(object = DelayedArray::seed(x = x), name = 'name')
  f <- paste(
    "function(x)",
    "HDF5Array::HDF5Array(filepath = x, name =",
    sQuote(x = name, q = FALSE),
    ", as.sparse =",
    sparse,
    ")"
  )
  return(f)
}

#' @rdname dot-DiskLoad
#' @method .DiskLoad IterableMatrix
#' @export
#'
.DiskLoad.IterableMatrix <- function(x) {
  check_installed(
    pkg = 'BPCells',
    reason = 'for working with BPCells matrices'
  )
  fxns <- lapply(
    X = BPCells::all_matrix_inputs(x = x),
    FUN = .DiskLoad
  )
  fxns <- Filter(f = Negate(f = is.null), x = fxns)
  if (!length(x = fxns)) {
    return(NULL)
  }
  fn <- if (length(x = fxns) > 1L) {
    # fxns <- paste('list(', paste(sQuote(x = fxns, q = FALSE), collapse = ', '), ')')
    fn <- paste(
      "function(x) {",
      "paths <- unlist(x = strsplit(x = x, split = ','));",
      "fxns <- list(", paste(sQuote(x = fxns, q = FALSE), collapse = ', '), ");",
      "mats <- vector(mode = 'list', length = length(x = paths));",
      "for (i in seq_along(paths)) {",
      "fn <- eval(str2lang(fxns[[i]]));",
      "mats[[i]] <- fn(paths[i]);",
      "};",
      "return(Reduce(cbind, mats));",
      "}"
    )
    fn
    # abort(message = "too many matrices")
  } else {
    fxns[[1L]]
  }
  return(fn)
}

#' @rdname dot-DiskLoad
#' @method .DiskLoad MatrixDir
#' @export
#'
.DiskLoad.MatrixDir <- function(x) {
  check_installed(
    pkg = 'BPCells',
    reason = 'for working with BPCells matrices'
  )
  f <- paste(
    'function(x)',
    'BPCells::open_matrix_dir(dir = x)'
  )
  return(f)
}

#' @rdname dot-DiskLoad
#' @method .DiskLoad MatrixH5
#' @export
#'
.DiskLoad.MatrixH5 <- function(x) {
  check_installed(
    pkg = 'BPCells',
    reason = 'for working with BPCells matrices'
  )
  f <- paste(
    'function(x)',
    'BPCells::open_matrix_hdf5(path = x, group =',
    sQuote(x = slot(object = x, name = 'group'), q = FALSE),
    ')'
  )
  return(f)
}

#' @rdname dot-DiskLoad
#' @method .DiskLoad TileDBMatrix
#' @export
#'
.DiskLoad.TileDBMatrix <- function(x) {
  check_installed(
    pkg = 'TileDBArray',
    reason = 'for working with TileDB matrices'
  )
  tdb.attr <- slot(object = DelayedArray::seed(x = x), name = 'attr')
  f <- paste(
    'function(x)',
    'TileDBArray::TileDBArray(x = x, attr =',
    sQuote(x = tdb.attr, q = FALSE),
    ')'
  )
  return(f)
}

#' @rdname dot-FilePath
#' @method .FilePath default
#' @export
#'
.FilePath.default <- function(x) {
  return(NULL)
}

#' @rdname dot-FilePath
#' @method .FilePath DelayedMatrix
#' @export
#'
.FilePath.DelayedMatrix <- function(x) {
  check_installed(
    pkg = 'DelayedArray',
    reason = 'for working with delayed matrices'
  )
  path <- tryCatch(
    expr = normalizePath(path = DelayedArray::path(object = x)),
    error = \(...) NULL
  )
  if (is.null(x = path)) {
    warn(message = "The matrix provided does not exist on-disk")
  }
  return(path)
}

#' @rdname dot-FilePath
#' @method .FilePath IterableMatrix
#' @export
#'
.FilePath.IterableMatrix <- function(x) {
  check_installed(pkg = "BPCells", reason = "for working with BPCells matrices")
  matrices <- BPCells::all_matrix_inputs(x = x)
  paths <- vector(mode = 'character', length = length(x = matrices))
  for (i in seq_along(along.with = matrices)) {
    mode <- .BPMatrixMode(object = matrices[[i]])
    paths[i] <- switch(
      EXPR = mode,
      memory = '',
      file = slot(object = matrices[[i]], name = "path"),
      directory = slot(object = matrices[[i]], name = 'dir'),
      abort(message = paste("Unknown BPCells matrix mode:", sQuote(x = mode)))
    )
  }
  if (length(paths) > 1){
    paths <- paste(paths, collapse = ",")
  }
  return(paths)
}

#' @rdname as.Centroids
#' @method as.Centroids Segmentation
#' @export
#'
as.Centroids.Segmentation <- function(
  x,
  nsides = NULL,
  radius = NULL,
  theta = NULL,
  ...
) {
  coords <- as(object = x, Class = 'Centroids')
  if (!is.null(x = nsides)) {
    slot(object = coords, name = 'nsides') <- nsides
  }
  if (!is.null(x = theta)) {
    slot(object = coords, name = 'theta') <- theta
  }
  if (is.null(x = radius)) {
    radius <- vapply(
      X = Cells(x = x),
      FUN = function(i) {
        area <- slot(
          object = slot(object = x, name = 'polygons')[[i]],
          name = 'area'
        )
        return(sqrt(x = area / pi))
      },
      FUN.VALUE = numeric(length = 1L),
      USE.NAMES = FALSE
    )
  }
  slot(object = coords, name = 'radius') <- radius
  validObject(object = coords)
  return(coords)
  # x <- c()
  # y <- c()
  # radius <- c()
  # nsides <- 0
  # for (cell in Cells(x)) {
  #   a <- x@polygons[[cell]]@area
  #   radius <- c(radius, sqrt(a / pi))
  #   x <- c(x, x@polygons[[cell]]@labpt[1])
  #   y <- c(y, x@polygons[[cell]]@labpt[2])
  # }
  # coords <- data.frame(x, y)
  # rownames(x = coords) = Cells(x)
  # return(
  #   CreateCentroids(
  #     coords,
  #     radius = radius,
  #     theta = rep(0, length(radius)),
  #     nsides = rep(0, length(radius))
  #   )
  # )
}

#' @rdname as.Centroids
#' @method as.Segmentation Centroids
#' @export
#'
as.Segmentation.Centroids <- function(x, ...) {
  return(as(object = x, Class = 'Segmentation'))
}

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
  return(as(object = as(object = as(object = x, Class = "dMatrix"), Class = "generalMatrix"), Class = "CsparseMatrix"))
}

#' @rdname as.sparse
#' @export
#' @method as.sparse matrix
#'
as.sparse.matrix <- function(x, ...) {
  if (is.character(x = x)) {
    dnames <- dimnames(x = x)
    nc <- ncol(x = x)
    x <- matrix(data = as.numeric(x = x), ncol = nc)
    dimnames(x = x) <- dnames
  }
  x <- as(object = x, Class = "Matrix")
  return(as.sparse.Matrix(x, ...))
}

#' @rdname as.sparse
#' @export
#' @method as.sparse ngCMatrix
#'
as.sparse.ngCMatrix <- function(x, ...) {
  return(as(object = x, Class = "dMatrix"))
}

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
  checks <- arg_match(arg = checks, multiple = TRUE)
  x <- slot(object = object, name = 'x')
  for (i in checks) {
    switch(
      EXPR = i,
      'infinite' = if (any(is.infinite(x = x))) {
        warn(message = "Input matrix contains infinite values")
      },
      'logical' = if (any(is.logical(x = x))) {
        warn(message = "Input matrix contains logical values")
      },
      'integer' = if (!all(round(x = x) == x, na.rm = TRUE)) {
        warn(message = "Input matrix contains non-integer values")
      },
      'na' = if (anyNA(x = x)) {
        warn(message = "Input matrix contains NA/NaN values")
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
  warn(message = "Input matrix contains logical values")
  return(invisible(x = NULL))
}

#' @rdname IsMatrixEmpty
#' @export
#' @method IsMatrixEmpty default
#'
IsMatrixEmpty.default <- function(x) {
  matrix.dims <- dim(x = x)
  if (is.null(x = matrix.dims)) {
    return(FALSE)
  }
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

#' Simplify segmentations by reducing the number of vertices
#'
#' @param coords A `Segmentation` object
#' @param tol Numerical tolerance value to be used by the Douglas-Peuker algorithm
#' @param topologyPreserve Logical determining if the algorithm should attempt to preserve the topology of the original geometry
#'
#' @return A `Segmentation` object with simplified segmentation vertices
#'
#' @rdname Simplify
#' @method Simplify Spatial
#' @export
#'
Simplify.Spatial <- function(coords, tol, topologyPreserve = TRUE) {
  check_installed(pkg = 'sf', reason = 'to simplify spatial data')
  class.orig <- class(x = coords)
  coords.orig <- coords
  dest <- ifelse(
    test = grepl(pattern = "^Spatial", x = class.orig),
    yes = class.orig,
    no = grep(pattern = "^Spatial", x = .Contains(object = coords), value = TRUE)[1L]
  )
  x <- sf::st_as_sfc(as(object = coords, Class = dest))
  coords <- sf::st_simplify(
    x = x,
    dTolerance = as.numeric(x = tol),
    preserveTopology = isTRUE(x = topologyPreserve))
  coords <- sf::st_sf(geometry = coords)
  coords <- as(coords, Class = "Spatial")
  coords <- as(coords, Class = "Segmentation")
  slot(object = coords, name = "polygons") <- mapply(
    FUN = function(x, y) {
      slot(object = x, name = "ID") <- y
      return(x)
    },
    slot(object = coords, name = "polygons"),
    Cells(coords.orig))
  return(coords)
}

#' Generate empty dgC sparse matrix
#'
#' @param ncol,nrow Number of columns and rows in matrix
#' @param rownames,colnames Optional row- and column names for the matrix
#'
#' @keywords internal
#'
#' @export
#'
SparseEmptyMatrix <- function(nrow, ncol, rownames = NULL, colnames = NULL) {
  return(new(
    Class = 'dgCMatrix',
    p = integer(length = ncol + 1L),
    Dim = c(as.integer(x = nrow), as.integer(x = ncol)),
    Dimnames = list(rownames, colnames)
  ))
}

#' @method StitchMatrix default
#' @export
#'
StitchMatrix.default <- function(x, y, rowmap, colmap, ...) {
  abort(message = paste(
    "Stitching matrices of class",
    dQuote(x = class(x = x)[1L]),
    "is not yet supported"
  ))
}

#' @method StitchMatrix dgCMatrix
#' @export
#'
StitchMatrix.dgCMatrix <- function(x, y, rowmap, colmap, ...) {
  on.exit(expr = CheckGC())
  if (!is_bare_list(x = y)) {
    y <- list(y)
  }
  rowmap <- droplevels(x = rowmap)
  colmap <- droplevels(x = colmap)
  stopifnot(ncol(rowmap) == length(y) + 1L)
  stopifnot(ncol(colmap) == length(y) + 1L)
  stopifnot(identical(x = colnames(x = rowmap), y = colnames(x = colmap)))
  dimnames(x = x) <- list(rowmap[[1L]], colmap[[1L]])
  for (i in seq_along(along.with = y)) {
    j <- i + 1L
    y[[i]] <- as(object = y[[i]], Class = 'dgCMatrix')
    dimnames(x = y[[i]]) <- list(rowmap[[j]], colmap[[j]])
  }
  return(RowMergeSparseMatrices(mat1 = x, mat2 = y))
}

#' @method StitchMatrix IterableMatrix
#' @export
#'
StitchMatrix.IterableMatrix <- function(x, y,  rowmap, colmap, ...) {
  on.exit(expr = CheckGC())
  if (!is_bare_list(x = y)) {
    y <- list(y)
  }
  rowmap <- droplevels(x = rowmap)
  colmap <- droplevels(x = colmap)
  stopifnot(ncol(rowmap) == length(y) + 1L)
  stopifnot(ncol(colmap) == length(y) + 1L)
  stopifnot(identical(x = colnames(x = rowmap), y = colnames(x = colmap)))
  y <- c(x, y)
  for (i in seq_along(along.with = y)) {
    #expand matrix to the same size
    missing_row <- setdiff(x = rownames(x = rowmap), y = rowmap[[i]])
    if (length(x = missing_row) > 0) {
      zero_i <- SparseEmptyMatrix(
        nrow = length(x = missing_row),
        ncol = ncol(x = y[[i]]),
        colnames = colmap[[i]],
        rownames = missing_row
      )
      zero_i <- as(object = zero_i, Class = 'IterableMatrix')
      y[[i]] <- rbind(y[[i]], zero_i)[rownames(rowmap),]
    }
  }
  m <- Reduce(f = cbind, x = y)
  return(m)
}


#' @method StitchMatrix matrix
#' @export
#'
StitchMatrix.matrix <- function(x, y, rowmap, colmap, ...) {
  on.exit(expr = CheckGC())
  if (!is_bare_list(x = y)) {
    y <- list(y)
  }
  rowmap <- droplevels(x = rowmap)
  colmap <- droplevels(x = colmap)
  stopifnot(ncol(rowmap) == length(y) + 1L)
  stopifnot(ncol(colmap) == length(y) + 1L)
  stopifnot(identical(x = colnames(x = rowmap), y = colnames(x = colmap)))
  m <- matrix(
    data = 0,
    nrow = nrow(x = rowmap),
    ncol = nrow(x = colmap),
    dimnames = list(rownames(x = rowmap), rownames(x = colmap))
  )
  m[rowmap[[1L]], colmap[[1L]]] <- x
  for (i in seq_along(along.with = y)) {
    j <- i + 1L
    m[rowmap[[j]], colmap[[j]]] <- as.matrix(x = y[[i]])
  }
  return(m)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @method t spam
#' @export
#'
t.spam <- spam::t

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

#' @importFrom stats median
#'
.FeatureRank <- function(features, flist, ranks = FALSE) {
  franks <- vapply(
    X = features,
    FUN = function(x) {
      return(median(x = unlist(x = lapply(
        X = flist,
        FUN = function(fl) {
          if (x %in% fl) {
            return(which(x = x == fl))
          }
          return(NULL)
        }
      ))))
    },
    FUN.VALUE = numeric(length = 1L)
  )
  franks <- sort(x = franks)
  if (!isTRUE(x = ranks)) {
    franks <- names(x = franks)
  }
  return(franks)
}

#' Move Files and Directories
#'
#' Move files and directories with \pkg{fs}; includes a handler for when
#' \code{path} is a directory on a different filesystem than \code{new_path}
#' by explicitly copying and deleting \code{path}
#'
#' @inherit fs::file_move params return
#' @inheritParams rlang::caller_env
#'
#' @keywords internal
#'
#' @export
#'
#' @templateVar pkg fs
#' @template note-reqdpkg
#'
#' @seealso \code{\link[fs:file_move]{fs::file_move}()}
#'
.FileMove <- function(path, new_path, overwrite = FALSE, n = 1L) {
  check_installed(pkg = "fs", reason = "for moving on-disk files")
  stopifnot(
    is_scalar_character(x = path),
    is_scalar_character(x = new_path),
    rlang::is_bare_logical(x = overwrite, n = 1L),
    is_bare_integerish(x = n, n = 1L, finite = TRUE) && n > 0
  )
  eexist <- function(err) {
    warn(
      message = paste(
        strwrap(x = paste(
          "Trying to move",
          sQuote(x = path),
          "to itself, skipping"
        )),
        collapse = '\n'
      ),
      class = c('WEXIST', 'EEXIST')
    )
    return(fs::as_fs_path(x = path))
  }
  hndlr <- function(err) {
    abort(
      message = err$message,
      class = class(x = err),
      call = caller_env(n = 4L + n)
    )
  }
  if (fs::is_dir(path = path)) {
    path <- fs::path_expand(path = path)
    new_path <- fs::path_expand(path = new_path)
    new_path <- fs::dir_create(path = new_path)
    dest <- tryCatch(
      expr = fs::dir_copy(
        path = path,
        new_path = new_path,
        overwrite = overwrite
      ),
      EEXIST = eexist,
      error = hndlr
    )
  } else if (fs::is_file(path = path)) {
    dest <- tryCatch(
      expr = fs::file_copy(
        path = path,
        new_path = new_path,
        overwrite = overwrite
      ),
      EEXIST = eexist,
      error = hndlr
    )
  } else {
    abort(
      message = paste(
        strwrap(x = paste0(
          "Can't find path: ",
          sQuote(x = path),
          "; if path is relative, change working directory"
        )),
        sep = '\n'
      ),
      call = caller_env(n = 1L + n)
    )
  }
  return(invisible(x = dest))
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

#' Get English Vowels
#'
#' @return A vector with English vowels in lower case
#'
#' @keywords internal
#'
#' @examples
#' .Vowels()
#'
#' @noRd
#'
.Vowels <- function() {
  return(c('a', 'e', 'i', 'o', 'u'))
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
  if (anyDuplicated(x = cell.names)) {
    if (isTRUE(x = stop)) {
      stop("Duplicate cell names present across objects provided.", call. = FALSE)
    }
    if (verbose) {
      warning(
        "Some cell names are duplicated across objects provided. Renaming to enforce unique cell names.",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    for (i in seq_along(along.with = object.list)) {
      object.list[[i]] <- RenameCells(
        object = object.list[[i]],
        new.names = paste(
        colnames(x = object.list[[i]]),
        i,
        sep = '_'
      ))
    }
  }
  return(object.list)
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
#' is empty if it has no length or an \code{nzchar == FALSE}
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
  mode <- arg_match(arg = mode)
  na <- arg_match(arg = na)
  x <- switch(
    EXPR = na,
    empty = x[is.na(x = x)] <- '',
    remove = x <- x[!is.na(x = x)],
    x
  )
  if (!length(x = x)) {
    return(TRUE)
  }
  empty <- vapply(
    X = x,
    FUN = Negate(f = nzchar),
    FUN.VALUE = logical(length = 1L),
    USE.NAMES = FALSE
  )
  empty <- switch(
    EXPR = mode,
    any = any(empty),
    all = all(empty),
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
  if (!is_bare_list(x = x)) {
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
#' @concept utils
#'
#' @export
#'
UpdateSlots <- function(object) {
  if (!isS4(object)) {
    return(object)
  }
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
  op <- options(Seurat.object.validate = FALSE)
  on.exit(expr = options(op), add = TRUE)
  object <- suppressWarnings(expr = do.call(what = 'new', args = object.list))
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
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setAs(
  from = 'Centroids',
  to = 'Segmentation',
  def = function(from) {
    if (is.infinite(x = from)) {
      stop("Cannot convert shapeless Centroids", call. = FALSE)
    }
    return(CreateSegmentation(coords = GetTissueCoordinates(
      object = from,
      full = TRUE
    )))
  }
)

setAs(
  from = 'Segmentation',
  to = 'Centroids',
  def = function(from) {
    return(CreateCentroids(coords = GetTissueCoordinates(
      object = from,
      full = FALSE
    )))
  }
)
