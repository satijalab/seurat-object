#' @include zzz.R
#' @include assay5.R
#' @include logmap.R
#'
#' @importFrom methods setGeneric
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setClass(
  Class = 'Seurat5',
  slots = c(
    assays = 'list',
    reductions = 'list',
    cells = 'LogMap',
    meta.data = 'data.frame',
    version = 'package_version'
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setGeneric(
  name = '.AddObject',
  def = function(object, name, value) {
    standardGeneric(f = '.AddObject')
  },
  signature = c('object', 'value')
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @method Assays Seurat5
#' @export
#'
Assays.Seurat5 <- function(object, ...) {
  return(names(x = slot(object = object, name = 'assays')))
}

#' @method Cells Seurat5
#' @export
#'
Cells.Seurat5 <- function(x, assay = NA, ...) {
  assay <- assay[1L] %||% DefaultAssay(object = x)
  if (is.na(x = assay)) {
    return(rownames(x = slot(object = x, name = 'cells')))
  }
  assay <- match.arg(arg = assay, choices = Assays(object = x))
  return(slot(object = x, name = 'cells')[[assay]])
}

#' @method CreateSeurat5Object default
#' @export
#'
CreateSeurat5Object.default <- function(
  counts,
  assay = 'RNA',
  names.field = 1L,
  names.delim = '_',
  meta.data = NULL,
  min.cells = 0,
  min.features = 0,
  ...
) {
  assay.data <- CreateAssay5Object(
    counts = counts,
    min.cells = min.cells,
    min.features = min.features,
   ...
  )
  return(CreateSeurat5Object(
    counts = assay.data,
    assay = assay,
    names.field = names.field,
    names.delim = names.delim,
    meta.data = meta.data
  ))
}

#' @method CreateSeurat5Object StdAssay
#' @export
#'
CreateSeurat5Object.StdAssay <- function(
  counts,
  assay = 'RNA',
  names.field = 1L,
  names.delim = '_',
  meta.data = NULL,
  ...
) {
  cells <- LogMap(y = Cells(x = counts))
  cells[[assay]] <- Cells(x = counts)
  if (IsCharEmpty(x = Key(object = counts))) {
    Key(object = counts) <- Key(object = tolower(x = assay), quiet = TRUE)
  }
  assay.list <- list(counts)
  names(x = assay.list) <- assay
  object <- new(
    Class = 'Seurat5',
    assays = assay.list,
    reductions = list(),
    cells = cells,
    meta.data = EmptyDF(n = nrow(x = cells)),
    version = packageVersion(pkg = 'SeuratObject')
  )
  # TODO: Calculate nCount and nFeature
  n.calc <- CalcN5(object = counts)
  if (!is.null(x = n.calc)) {
    names(x = n.calc) <- paste(names(x = n.calc), assay, sep = '_')
  }
  return(object)
}

#' @method DefaultAssay Seurat5
#' @export
#'
DefaultAssay.Seurat5 <- function(object, ...) {
  return(Assays(object = object)[1L])
}

#' @method DefaultAssay<- Seurat5
#' @export
#'
"DefaultAssay<-.Seurat5" <- function(object, ..., value) {
  value <- value[1L]
  assays <- Assays(object = object)
  value <- match.arg(arg = value, choices = assays)
  idx <- which(x = assays == value)
  slot(object = object, name = 'assays') <- c(
    slot(object = object, name = 'assays')[idx],
    slot(object = object, name = 'assays')[-idx]
  )
  return(object)
}

#' @method Features Seurat5
#' @export
#'
Features.Seurat5 <- function(x, assay = NULL, ...) {
  assay <- assay[1L] %||% DefaultAssay(object = x)
  assay <- match.arg(arg = assay, choices = Assays(object = x))
  return(Features(x = x[[assay]], ...))
}

#' @method Key Seurat5
#' @export
#'
Key.Seurat5 <- function(object, ...) {
  return(sapply(
    X = .FilterObjects(
      object = object,
      classes.keep = c('KeyMixin')
    ),
    FUN = function(x) {
      return(Key(object = object[[x]]))
    }
  ))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @method [[ Seurat5
#' @export
#'
"[[.Seurat5" <- function(x, i, ..., drop = FALSE) {
  if (missing(x = i)) {
    return(slot(object = x, name = 'meta.data'))
  }
  slot.use <- .FindObject(object = x, name = i)
  if (is.null(x = slot.use)) {
    stop("Cannot find '", i, "' in this Seurat object", call. = FALSE)
  }
  data.return <- slot(object = x, name = slot.use)[[i]]
  return(data.return)
}

#' @method dim Seurat5
#' @export
#'
dim.Seurat5 <- function(x) {
  return(c(
    tryCatch(
      expr = nrow(x = x[[DefaultAssay(object = x)]]),
      error = function(...) {
        return(0L)
      }
    ),
    nrow(x = slot(object = x, name = 'cells'))
  ))
}

#' @method dimnames Seurat5
#' @export
#'
dimnames.Seurat5 <- function(x) {
  .NotYetImplemented()
}

#' @method names Seurat5
#' @export
#'
names.Seurat5 <- function(x) {
  return(.FilterObjects(
    object = x,
    classes.keep = c('StdAssay','DimReduc')
  ))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.AddObjectAssay <- function(object, name, value) {
  # Overwrite existing assay
  if (name %in% names(x = object)) {
    if (!inherits(x = object[[name]], what = 'StdAssay')) {
      stop("duplicate names")
    }
    if (!identical(x = dim(x = value), y = dim(x = object[[name]]))) {
      stop("different cells/features from existing")
    }
    if (!all(Cells(x = value) %in% Cells(x = object, assay = name))) {
      stop("different cells")
    }
  } else {
    # Add new assay
    if (!all(Cells(x = value) %in% Cells(x = object))) {
      stop("new cells")
    }
    slot(object = object, name = 'cells')[[name]] <- Cells(x = value)
  }
  if (IsCharEmpty(x = Key(object = value))) {
    Key(object = value) <- Key(object = tolower(x = name))
  }
  if (Key(object = value) %in% Key(object = object) && !name %in% names(x = object)) {
    old <- Key(object = value)
    Key(object = value) <- Key(object = tolower(x = name), quiet = TRUE)
    i <- 1L
    n <- 5L
    while (Key(object = value) %in% Key(object = object)) {
      Key(object = value) <- Key(object = RandomName(length = n), quiet = TRUE)
      i <- i + 1L
      if (!i %% 7L) {
        n <- n + 2L
        i <- 1L
      }
    }
    warning(
      "Key '",
      old,
      "' taken, using '",
      Key(object = value),
      "' instead",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  cmatch <- MatchCells(
    new = colnames(x = object),
    orig = colnames(x = value),
    ordered = TRUE
  )
  value <- subset(x = value, cells = cmatch)
  slot(object = object, name = 'assays')[[name]] <- value
  validObject(object = object)
  return(object)
}

#' Check for Duplicate Keys
#'
#' Check the uniqueness and validity of a key. If the provided key is not
#' unique, \code{.CheckKey} will attempt to create a new key based on the name
#' associated with the key, if provided. If not provided, or still not unique,
#' \code{.CheckKey} will generate a series of random keys until a unique one
#' is found
#'
#' @param key A key to check
#' @param existing A vector of existing keys; can optionally be named with
#' the names of objects corresponding to the keys
#' @param name Name for object that \code{key} will be associated with
#'
#' @return A valid, unique key
#'
#' @keywords internal
#'
#' @noRd
#'
.CheckKey <- function(key, existing = NULL, name = NULL) {
  key <- Key(object = key, quiet = !is.null(x = existing))
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
    warning(
      "Key '",
      old,
      "' taken, using '",
      key,
      "' instead",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  return(key)
}

.DuplicateError <- function(name, cls, error = TRUE) {
  letter <- tolower(x = substr(x = cls[1L], start = 1L, stop = 1L))
  article <- ifelse(test = letter %in% .Vowels(), yes = 'an', no = 'a')
  msg <- paste0("'", name[1L], "' already taken for ", paste(article, cls))
  if (isTRUE(x = error)) {
    stop(msg, call. = FALSE)
  }
  return(msg)
}

#' @importFrom methods slotNames
#' @importFrom rlang is_bare_list
#'
.ObjectNames <- function(object) {
  slots <- setdiff(x = slotNames(x = object), y = c('misc', 'tools'))
  slots <- Filter(
    f = function(x) {
      return(is_bare_list(x = slot(object = object, name = x)))
    },
    x = slots
  )
  objects <- lapply(
    X = slots,
    FUN = function(x) {
      obj.names <- names(x = slot(object = object, name = x))
      cls <- switch(
        EXPR = x,
        'assays' = 'StdAssay',
        'reductions' = 'DimReduc',
        'commands' = 'SeuratCommand',
        x
      )
      objs <- rep_len(x = cls, length.out = length(x = obj.names))
      names(x = objs) <- obj.names
      return(objs)
    }
  )
  return(unlist(x = objects))
}

.RemoveObject <- function(object, name, value) {
  slot.use <- .FindObject(object = object, name = name) %||% 'meta.data'
  switch(
    EXPR = slot.use,
    'meta.data' = {
      .NotYetImplemented()
    },
    'assays' = {
      if (isTRUE(x = name == DefaultAssay(object = object))) {
        stop("Cannot delete default assay", call. = FALSE)
      }
      slot(object = object, name = slot.use)[[name]] <- value
      cmat <- slot(object = object, name = 'cells')
      cmat <- cmat[, -which(x = colnames(x = cmat) == name), drop = FALSE]
      slot(object = object, name = 'cells') <- cmat
    },
    slot(object = object, name = slot.use)[[name]] <- value
  )
  validObject(object = object)
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setMethod(
  f = '[[<-',
  signature = c(x = 'Seurat5', i = 'character', j = 'missing'),
  definition = function(x, i, ..., value) {
    x <- .AddObject(object = x, name = i, value = value)
    return(x)
  }
)

setMethod(
  f = '[[<-',
  signature = c(
    x = 'Seurat5',
    i = 'character',
    j = 'missing',
    value = 'StdAssay'
  ),
  definition = function(x, i, ..., value) {
    i <- make.names(names = i)
    # Checks for if the assay or name already exists
    if (i %in% names(x = x)) {
      # Check for duplicate names
      if (!inherits(x = x[[i]], what = 'StdAssay')) {
        .DuplicateError(name = i, cls = class(x = x[[i]]))
      }
      if (!identical(x = dim(x = value), y = dim(x = x[[i]]))) {
        stop("different cells/features from existing")
      }
      if (!all(Cells(x = value, layer = NA) == Cells(x = x, assay = i))) {
        stop("different cells")
      }
    } else {
      if (!all(Cells(x = value) %in% Cells(x = x))) {
        stop("new cells")
      }
      slot(object = x, name = 'cells')[[i]] <- Cells(x = value, layer = NA)
    }
    # Check keys
    Key(object = value) <- .CheckKey(
      key = Key(object = value),
      existing = Key(object = x),
      name = i
    )
    # Add the assay
    slot(object = x, name = 'assays')[[i]] <- value
    slot(object = x, name = 'assays') <- Filter(
      f = Negate(f = is.null),
      x = slot(object = x, name = 'assays')
    )
    # TODO: Update the cells LogMap
    # Return the Seurat object
    validObject(object = x)
    return(x)
  }
)

setMethod(
  f = '[[<-',
  signature = c(
    x = 'Seurat5',
    i = 'character',
    j = 'missing',
    value = 'DimReduc'
  ),
  definition = function(x, i, ..., value) {
    i <- make.names(names = i)
    # Check for duplicate names
    if (i %in% names(x = x)) {
      # Checks for overwriting DimReducs
      if (inherits(x = x[[i]], what = 'DimReduc')) {
        ''
      } else {
        .DuplicateError(name = i, cls = class(x = x[[i]]))
      }
    }
    # Check keys
    Key(object = value) <- .CheckKey(
      key = Key(object = value),
      existing = Key(object = x),
      name = i
    )
    # TODO: Check cells
    # Add the DimReduc
    slot(object = x, name = 'reductions')[[i]] <- value
    slot(object = x, name = 'reductions') <- Filter(
      f = Negate(f = is.null),
      x = slot(object = x, name = 'reductions')
    )
    return(x)
  }
)

setMethod(
  f = '[[<-',
  signature = c(x = 'Seurat5', i = 'character', j = 'missing', value = 'vector'),
  definition = function(x, i, ..., value) {
    if (length(x = value) > length(x = i)) {
      warning(
        "More values provided than names, using only first ",
        length(x = i),
        " values",
        call. = FALSE,
        immediate. = TRUE
      )
    } else if (!length(x = value)) {
      stop("At least one value must be supplied", call. = FALSE)
    }
    value <- rep_len(x = value, length.out = length(x = i))
    for (idx in seq_along(along.with = i)) {
      x <- .AddObject(object = x, name = i[[idx]], value = value[[idx]])
    }
    return(x)
  }
)

setMethod(
  f = '[[<-',
  signature = c(x = 'Seurat5', i = 'character', j = 'missing', value = 'NULL'),
  definition = function(x, i, ..., value) {
    # return(.AddObject(object = x, name = i, value = value))
    for (name in i) {
      slot.use <- .FindObject(object = x, name = name) %||% 'meta.data'
      switch(
        EXPR = slot.use,
        'meta.data' = {
          .NotYetImplemented()
        },
        'assays' = {
          if (isTRUE(x = name == DefaultAssay(object = x))) {
            stop("Cannot delete default assay", call. = FALSE)
          }
          slot(object = x, name = slot.use)[[x]] <- value
          cmat <- slot(object = x, name = 'cells')
          cmat <- cmat[, -which(x = colnames(x = cmat) == name), drop = FALSE]
          slot(object = x, name = 'cells') <- cmat
        },
        slot(object = x, name = slot.use)[[name]] <- value
      )
    }
    validObject(object = x)
    return(x)
  }
)

setMethod(
  f = '.AddObject',
  signature = c(object = 'Seurat5', value = 'StdAssay'),
  definition = function(object, name, value) {
    return(.AddObjectAssay(object = object, name = name, value = value))
  }
)

setMethod(
  f = '.AddObject',
  signature = c(object = 'Seurat5', value = 'NULL'),
  definition = function(object, name, value) {
    for (i in name) {
      object <- .RemoveObject(object = object, name = i, value = value)
    }
    return(object)
  }
)

setMethod(
  f = 'show',
  signature = c(object = 'Seurat5'),
  definition = function(object) {
    nfeatures <- sum(vapply(
      X = Assays(object = object),
      FUN = function(x) {
        return(nrow(x = object[[x]]))
      },
      FUN.VALUE = numeric(length = 1L),
      USE.NAMES = FALSE
    ))
    nassays <- length(x = Assays(object = object))
    cat(
      "A Seurat (v5) object\n"
    )
    cat(
      nfeatures,
      'features across',
      ncol(x = object),
      'samples within',
      nassays,
      ifelse(test = nassays == 1, yes = 'assay', no = 'assays'),
      '\n'
    )
    cat(
      'Active assay:',
      DefaultAssay(object = object),
      paste0('(', nrow(x = object), ' features)')
      # paste0(
      #   "(",
      #   nrow(x = object),
      #   ' features, ',
      #   length(x = VariableFeatures(object = object)),
      #   ' variable features)'
      # )
    )
    nother <- nassays - 1L
    if (as.logical(x = nother)) {
      cat(
        '\n',
        nother,
        'other',
        ifelse(test = nother == 1, yes = 'assay', no = 'assays'),
        'present:',
        strwrap(x = paste(
          setdiff(
            x = Assays(object = object),
            y = DefaultAssay(object = object)
          ),
          collapse = ', '
        ))
      )
    }
    return(invisible(x = NULL))
  }
)

setValidity(
  Class = 'Seurat5',
  method = function(object) {
    valid <- NULL
    # TODO: Check assays
    if (!IsNamedList(x = slot(object = object, name = 'assays'))) {
      valid <- c(valid, "'assays' must be a named list")
    } else {
      for (assay in Assays(object = object)) {
        if (!inherits(x = object[[assay]], what = 'StdAssay')) {
          valid <- c(
            valid,
            "Objects in the 'assays' list must be 'Assay' objects"
          )
          break
        }
      }
    }
    # TODO: Check reductions
    # TODO: Check metadata
    # TODO: Check version
    if (length(x = slot(object = object, name = 'version')) > 1) {
      valid <- c(valid, "Only one version is allowed")
    }
    return(valid %||% TRUE)
  }
)

