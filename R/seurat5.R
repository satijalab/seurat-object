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
Cells.Seurat5 <- function(x, assay = NULL, ...) {
  if (is.null(x = assay)) {
    return(rownames(x = slot(object = x, name = 'cells')))
  }
  assay <- assay[1L]
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
  return(Features(x = x[[assay]]))
}

#' @method Key Seurat5
#' @export
#'
Key.Seurat5 <- function(object, ...) {
  return(sapply(
    X = .FilterObjects(
      object = object,
      classes.keep = c('StdAssay', 'DimReduc', 'SpatialImage')
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
    if (!all(Cells(x = value) == Cells(x = object, assay = name))) {
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
      Key(object = value) <- Key(object = RandomName(), quiet = TRUE)
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
  slot(object = object, name = 'assays')[[name]] <- value
  validObject(object = object)
  return(object)
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
    return(.AddObject(object = x, name = i, value = value))
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
