#' @include zzz.R
#' @include assay.R
#' @include layers.R
#' @include logmap.R
#' @include keymixin.R
#' @importFrom methods callNextMethod setAs
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Core Assay Infrastructure
#'
#' The \code{StdAssay} class is a virtual class that provides core
#' infrastructure for assay data in \pkg{Seurat}. Assays contain expression
#' data (layers) and associated feature-level meta data. Derived classes
#' (eg. \link[=Assay5]{the v5 Assay}) may optionally
#' define additional functionality
#'
#' @template slot-stdassay
#' @template slot-misc
#' @template slot-key
#'
#' @keywords internal
#'
#' @exportClass StdAssay
#'
#' @aliases StdAssay
#'
#' @family stdassay
#'
#' @seealso \code{\link{Assay5-class}} \code{\link{Assay5T-class}}
#'
setClass(
  Class = 'StdAssay',
  contains = c('VIRTUAL', 'KeyMixin'),
  slots = c(
    layers = 'list',
    cells = 'LogMap',
    features = 'LogMap',
    default = 'integer',
    assay.orig = 'character',
    meta.data = 'data.frame',
    misc = 'list'
  )
)

#' The v5 \code{Assay} Object
#'
#' The v5 \code{Assay} is the typical \code{Assay} class used in \pkg{Seurat}
#' v5; ...
#'
#' @template slot-stdassay
#' @template slot-misc
#' @template slot-key
#'
#' @exportClass Assay5
#'
#' @aliases Assay5
#'
#' @family assay5
#'
setClass(
  Class = 'Assay5',
  contains = 'StdAssay'
)

#' The Transposed v5 \code{Assay} Object
#'
#' @template slot-stdassay
#' @template slot-misc
#' @template slot-key
#'
#' @template lifecycle-experimental
#'
#' @keywords internal
#'
#' @aliases Assay5T
#'
setClass(
  Class = 'Assay5T',
  contains = 'StdAssay'
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Create a v5 Assay object
#'
#' Create an \code{\link{Assay5}} object from a feature expression matrix;
#' the expected format of the matrix is features x cells
#'
#' @inheritParams .CreateStdAssay
#' @param data Optional prenormalized data matrix
#' @template param-dots-method
# @param transpose Create a transposed assay
# @param ... Extra parameters passed to \code{\link{.CreateStdAssay}}
#'
#' @return An \code{\link{Assay5}} object
#'
#' @export
#'
#' @concept assay
#'
CreateAssay5Object <- function(
  counts = NULL,
  data = NULL,
  min.cells = 0,
  min.features = 0,
  csum = NULL,
  fsum = NULL,
  ...
) {
  transpose <- FALSE
  colsums <- Matrix::colSums
  rowsums <- Matrix::rowSums
  type <- 'Assay5'
  csum <- csum %||% colsums
  fsum <- fsum %||% rowsums
  counts <- CheckLayersName(matrix.list = counts, layers.type = 'counts')
  data <- CheckLayersName(matrix.list = data, layers.type = 'data')
  if (!is.null(x = counts) & !is.null(data)) {
    counts.cells <- unlist(
      x = lapply(
        X = counts,
        FUN = function(x) colnames(x = x)
      )
    )
    data.cells <- unlist(
      x = lapply(
        X = data,
        FUN = function(x) colnames(x)
      )
    )
    if (!all(counts.cells == data.cells)) {
      abort(message = 'counts and data input should have the same cells')
    }
  }
  counts <- c(counts, data)
  data <- NULL
  CheckGC()
  return(.CreateStdAssay(
    counts = counts,
    min.cells = min.cells,
    min.features = min.features,
    transpose = transpose,
    type = type,
    csum = csum,
    fsum = fsum,
    ...
  ))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @method .AssayClass Assay5T
#' @export
#'
.AssayClass.Assay5T <- function(object) {
  return('Transposed Assay (v5)')
}

#' @method .CalcN StdAssay
#' @export
#'
.CalcN.StdAssay <- function(object, layer = 'counts', simplify = TRUE, ...) {
  layer <- tryCatch(
    expr = Layers(object = object, search = layer),
    error = \(...) NULL
  ) # %||% DefaultLayer(object = object)
  if (is.null(x = layer)) {
    warn(
      message = "Cannot find the layer(s) specified",
      class = 'missingLayerWarning'
    )
    return(NULL)
  }
  calcn <- vector(mode = 'list', length = length(x = layer))
  names(x = calcn) <- layer
  for (lyr in layer) {
    ldat <- LayerData(object = object, layer = lyr)
    if (IsMatrixEmpty(x = ldat)) {
      next
    }
    calcn[[lyr]] <- .CalcN(object = ldat)
  }
  calcn <- Filter(f = length, x = calcn)
  # If every layer is empty, return `NULL`
  if (!length(x = calcn)) {
    return(NULL)
  } else if (isFALSE(x = simplify)) {
    # If we're not simplifying, return the list as-is
    return(calcn)
  } else if (length(x = calcn) == 1L) {
    # If we're only calculating N for one layer, return those results
    return(calcn[[1L]])
  }
  # Simplify the calcn list for all cells
  all.cells <- Cells(x = object, layer = layer, simplify = TRUE)
  ncells <- length(x = all.cells)
  ncalc <- list(
    nCount = vector(mode = 'numeric', length = ncells),
    nFeature = vector(mode = 'numeric', length = ncells)
  )
  names(x = ncalc$nCount) <- names(x = ncalc$nFeature) <- all.cells
  # For every layer, add the nCount and nFeature counts to existing cells
  for (i in seq_along(along.with = calcn)) {
    lcells <- names(x = calcn[[i]][['nCount']])
    ncalc[['nCount']][lcells] <- calcn[[i]][['nCount']] + ncalc[['nCount']][lcells]
    ncalc[['nFeature']][lcells] <- calcn[[i]][['nFeature']] + ncalc[['nFeature']][lcells]
  }
  return(ncalc)
}

#' @method .CalcN default
#' @export
#'
.CalcN.default <- function(object, ...) {
   return(list(
     nCount = Matrix::colSums(x = object),
     nFeature = Matrix::colSums(x = object > 0)
   ))
}

#' @param layer Name of layer to store \code{counts} as
#'
#' @rdname dot-CreateStdAssay
#' @method .CreateStdAssay default
#' @export
#'
.CreateStdAssay.default <- function(
  counts,
  min.cells = 0,
  min.features = 0,
  cells = NULL,
  features = NULL,
  transpose = FALSE,
  type = 'Assay5',
  layer = 'counts',
  ...
) {
  if (!is_bare_integerish(x = dim(x = counts), n = 2L, finite = TRUE)) {
    abort(message = "'counts' must be a two-dimensional object")
  }
  dnames <- dimnames(x = counts)
  cls <- class(x = counts)

  if (isTRUE(x = transpose)) {
    csum <- .GetMethod(fxn = 'rowSums', cls = cls)
    cells <- cells %||% dnames[[1L]]
    fsum <- .GetMethod(fxn = 'colSums', cls = cls)
    features <- features %||% dnames[[2L]]
  } else {
    csum <- .GetMethod(fxn = 'colSums', cls = cls)
    cells <- cells %||% dnames[[2L]]
    fsum <- .GetMethod(fxn = 'rowSums', cls = cls)
    features <- features %||% dnames[[1L]]
  }
  counts <- list(counts)
  names(x = counts) <- layer
  return(.CreateStdAssay(
    counts = counts,
    min.cells = min.cells,
    layer = layer,
    min.features = min.features,
    cells = cells,
    features = features,
    transpose = transpose,
    type = type,
    fsum = fsum,
    csum = csum,
    ...
  ))
}

#' @param csum Function for calculating cell sums
#' @param fsum Function for calculating feature sums
#'
#' @importFrom methods getClass
#' @importFrom utils getS3method methods
#'
#' @rdname dot-CreateStdAssay
#' @method .CreateStdAssay list
#' @export
#'
.CreateStdAssay.list <- function(
  counts,
  min.cells = 0,
  min.features = 0,
  cells = NULL,
  features = NULL,
  transpose = FALSE,
  type = 'Assay5',
  csum = Matrix::colSums,
  fsum = Matrix::rowSums,
  ...
) {
  # Figure out feature/cell MARGINs
  cdef <- getClass(Class = type)
  contains <- names(x = slot(object = cdef, name = 'contains'))
  if (!'StdAssay' %in% contains) {
    stop("Class '", type, "' does not inherit from StdAssay")
  }
  for (i in c(type, contains, 'default')) {
    fmargin <- getS3method(f = '.MARGIN', class = i, optional = TRUE)
    if (is.function(x = fmargin)) {
      break
    }
  }
  cdim <- fmargin(object = type, type = 'cells')
  fdim <- fmargin(object = type, type = 'features')
  counts <- lapply(X = counts, FUN = function(x) {
    x <- CheckFeaturesNames(data = x)
    return(x)
  })
  # Check cell/feature names for all layers
  if (is.atomic(x = cells) || is.null(x = cells)) {
    cells <- rep_len(x = list(cells), length.out = length(x = counts))
  }
  if (!is_bare_list(x = cells) || length(x = cells) != length(x = counts)) {
    stop("Not enough cells for the counts matrices provided", call. = FALSE)
  }
  cells <- .CheckNames(x = cells, n = names(x = counts))
  if (is.atomic(x = features) || is.null(x = features)) {
    features <- rep_len(x = list(features), length.out = length(x = counts))
  }
  if (!is_bare_list(x = features) || length(x = features) != length(x = counts)) {
    stop("Not enough features for the counts matrices provided", call. = FALSE)
  }
  features <- .CheckNames(x = features, n = names(x = counts))
  for (layer in names(x = counts)) {
    cells[[layer]] <- cells[[layer]] %||%
      dimnames(x = counts[[layer]])[[cdim]] %||%
      paste0('Cell_', seq_len(length.out = dim(x = counts[[layer]])[cdim]))
    features[[layer]] <- features[[layer]] %||%
      dimnames(x = counts[[layer]])[[fdim]] %||%
      paste0('Feature', seq_len(length.out = dim(x = counts[[layer]])[fdim]))
  }
  # Filter based on min.features
  if (min.features > 0) {
    for (layer in names(x = counts)) {
      if (inherits(x = counts[[layer]], what = "IterableMatrix")) {
        check_installed(pkg = 'BPCells', reason = 'for working with BPCells')
        col_stat <- BPCells::matrix_stats(matrix = counts[[layer]], col_stats = 'nonzero')$col_stats
        cells.use <- which(x = col_stat >= min.features)
      } else {
        cells.use <- which(x = csum(counts[[layer]] > 0) >= min.features)
      }
      counts[[layer]] <- if (cdim == 1L) {
        counts[[layer]][cells.use, ]
      } else {
        counts[[layer]][, cells.use]
      }
      cells[[layer]] <- cells[[layer]][cells.use]
    }
  }
  # For now, coerce to dgCMatrix if not dgCMatrix, IterableMatrix, or DelayedArray
  if (!inherits(x = counts[[layer]], what = c('dgCMatrix', 'IterableMatrix', 'DelayedArray'))) {
    warning('Data is of class ', class(counts[[layer]])[1], ". Coercing to dgCMatrix.",
            call. = FALSE, immediate. = TRUE)
    if (inherits(x = counts[[layer]], what = "data.frame")) {
      counts[[layer]] <- as.sparse(x = counts[[layer]], ...)
    } else {
      counts[[layer]] <- as.sparse(x = counts[[layer]])
    }
  }
  # Filter based on min.cells
  if (min.cells > 0) {
    for (layer in names(x = counts)) {
      if (inherits(x = counts[[layer]], what = "IterableMatrix")) {
        check_installed(pkg = 'BPCells', reason = 'for working with BPCells')
        row_stat <- BPCells::matrix_stats(matrix = counts[[layer]], row_stats = 'nonzero')$row_stats
        features.use <- which(x = row_stat >= min.cells)
      } else {
        features.use <- which(x = fsum(counts[[layer]] > 0) >= min.cells)
      }
      counts[[layer]] <- if (fdim == 1L) {
        counts[[layer]][features.use, ]
      } else {
        counts[[layer]][, features.use]
      }
      features[[layer]] <- features[[layer]][features.use]
    }
  }
  features.all <- Reduce(f = union, x = features)
  cells.all <- Reduce(f = union, x = cells)
  calcN_option <- getOption(
    x = 'Seurat.object.assay.calcn',
    default =  Seurat.options$Seurat.object.assay.calcn
  )
  # Create the object
  object <- new(
    Class = type,
    layers = list(),
    default = 0L,
    features = LogMap(y = features.all),
    cells = LogMap(y = cells.all),
    meta.data = EmptyDF(n = length(x = features.all)),
    misc = list(calcN = calcN_option %||% TRUE),
    ...
  )
  for (layer in names(x = counts)) {
    LayerData(
      object = object,
      layer = layer,
      features = features[[layer]],
      cells = cells[[layer]],
      transpose = transpose
    ) <- counts[[layer]]
  }
  DefaultLayer(object = object) <- Layers(object = object)[1L]
  validObject(object = object)
  return(object)
}

#' @rdname dot-CreateStdAssay
#' @method .CreateStdAssay Matrix
#' @export
#'
.CreateStdAssay.Matrix <- function(
  counts,
  min.cells = 0,
  min.features = 0,
  cells = NULL,
  features = NULL,
  transpose = FALSE,
  type = 'Assay5',
  layer = 'counts',
  ...
) {
  counts <- list(counts)
  names(x = counts) <- layer
  if (isTRUE(x = transpose)) {
    csum <- Matrix::rowSums
    fsum <- Matrix::colSums
  } else {
    csum <- Matrix::colSums
    fsum <- Matrix::rowSums
  }
  return(.CreateStdAssay(
    counts = counts,
    layer = layer,
    min.cells = min.cells,
    min.features = min.features,
    cells = cells,
    features = features,
    transpose = transpose,
    type = type,
    csum = csum,
    fsum = fsum,
    ...
  ))
}

#' @rdname dot-CreateStdAssay
#' @method .CreateStdAssay matrix
#' @export
#'
.CreateStdAssay.matrix <- .CreateStdAssay.Matrix

#' @method .MARGIN Assay5T
#' @export
#'
.MARGIN.Assay5T <- function(x, type = c('features', 'cells'), ...) {
  type <- type[1]
  type <- match.arg(arg = type)
  return(unname(obj = c(features = 2L, cells = 1L)[type]))
}

#' @method .SelectFeatures StdAssay
#' @export
#'
.SelectFeatures.StdAssay <- function(object, ...) {
  .NotYetImplemented()
}

#' @templateVar fxn AddMetaData
#' @template method-stdassay
#'
#' @method AddMetaData StdAssay
#' @export
#'
AddMetaData.StdAssay <- AddMetaData.Assay

#' @rdname AddMetaData
#' @method AddMetaData Assay5
#' @export
#'
AddMetaData.Assay5 <- AddMetaData.StdAssay

#' @templateVar fxn CastAssay
#' @template method-stdassay
#'
#' @importFrom methods as
#' @importFrom rlang quo_get_env quo_get_expr
#'
#' @method CastAssay StdAssay
#' @export
#'
CastAssay.StdAssay <- function(object, to, layers = NA, verbose = TRUE, ...) {
  layers <- Layers(object = object, search = layers)
  if (is_quosure(x = to)) {
    to <- eval(
      expr = quo_get_expr(quo = to),
      envir = quo_get_env(quo = to)
    )
  }
  stopifnot(is.character(x = to) || is.function(x = to))
  for (lyr in layers) {
    if (isTRUE(x = verbose)) {
      msg <- paste("Attempting to cast layer", lyr)
      if (is.character(x = to)) {
        msg <- paste(msg, "to", to)
      }
      message(msg)
    }
    clyr <- Cells(x = object, layer = lyr)
    flyr <- Features(x = object, layer = lyr)
    w <- function(e) {
      warn(message = paste0(
        "Unable to cast layer ",
        sQuote(x = lyr),
        ": ",
        e$message
      ))
      return(invisible(x = NULL))
    }
    if (is.function(x = to)) {
      tryCatch(
        expr = LayerData(
          object = object,
          layer = lyr,
          cells = clyr,
          features = flyr
        ) <- to(LayerData(object = object, layer = lyr, fast = TRUE), ...),
        error = w
      )
    } else {
      check <- is(
        object = LayerData(object = object, layer = lyr, fast = TRUE),
        class2 = to
      )
      if (isTRUE(x = check)) {
        next
      }
      tryCatch(
        expr = LayerData(
          object = object,
          layer = lyr,
          cells = clyr,
          features = flyr
        ) <- as(
          object = LayerData(object = object, layer = lyr, fast = TRUE),
          Class = to
        ),
        error = w
      )
    }
  }
  return(object)
}

#' @template param-verbose
#' @param layers A vector of layers to cast; defaults to all layers
#'
#' @rdname CastAssay
#' @method CastAssay Assay5
#' @export
#'
CastAssay.Assay5 <- CastAssay.StdAssay

#' @templateVar fxn Cells
#' @template method-stdassay
#'
#' @method Cells StdAssay
#' @export
#'
Cells.StdAssay <- function(x, layer = NULL, simplify = TRUE, ...) {
  if (any(is.na(x = layer)) || is.null(x = layer)) {
    return(rownames(x = slot(object = x, name = 'cells')))
  }
  layer <- Layers(object = x, search = layer)
  cells <- sapply(
    X = layer,
    FUN = function(lyr) {
      return(slot(object = x, name = 'cells')[[lyr]])
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  if (isFALSE(x = simplify)) {
    return(cells)
  }
  return(Reduce(f = union, x = cells))
}

#' @param layer Layer to pull cells/features for; defaults to default layer;
#' if \code{NA}, returns all cells for the assay
#' @param simplify Simplify the cell/feature names into a single vector; if
#' \code{FALSE}, separates each cell/feature names by layer
#'
#' @rdname Cells
#' @method Cells Assay5
#' @export
#'
Cells.Assay5 <- Cells.StdAssay

#' @templateVar fxn DefaultAssay
#' @template method-stdassay
#'
#' @method DefaultAssay StdAssay
#' @export
#'
DefaultAssay.StdAssay <- function(object, ...) {
  return(slot(object = object, name = 'assay.orig'))
}

#' @rdname DefaultAssay
#' @method DefaultAssay Assay5
#' @export
#'
DefaultAssay.Assay5 <- DefaultAssay.StdAssay

#' @rdname DefaultAssay-StdAssay
#' @method DefaultAssay<- StdAssay
#' @export
#'
"DefaultAssay<-.StdAssay" <- function(object, ..., value) {
  slot(object = object, name = 'assay.orig') <- value
  return(object)
}

#' @rdname DefaultAssay
#' @method DefaultAssay<- Assay5
#' @export
#'
"DefaultAssay<-.Assay5" <- `DefaultAssay<-.StdAssay`

#' @templateVar fxn DefaultLayer
#' @template method-stdassay
#'
#' @method DefaultLayer StdAssay
#' @export
#'
DefaultLayer.StdAssay <- function(object, ...) {
  idx <- slot(object = object, name = 'default')
  if (!length(x = idx) || idx == 0L) {
    idx <- 1L
  }
  return(Layers(object = object)[seq_len(length.out = idx)])
}

#' @rdname DefaultLayer
#' @method DefaultLayer Assay5
#' @export
#'
DefaultLayer.Assay5 <- DefaultLayer.StdAssay

#' @rdname DefaultLayer-StdAssay
#' @method DefaultLayer<- StdAssay
#' @export
#'
"DefaultLayer<-.StdAssay" <- function(object, ..., value) {
  layers <- Layers(object = object)
  value <- Layers(object = object, search = value)
  idx <- MatchCells(new = layers, orig = value, ordered = TRUE)
  slot(object = object, name = 'layers') <- c(
    slot(object = object, name = 'layers')[idx],
    slot(object = object, name = 'layers')[-idx]
  )
  slot(object = object, name = 'default') <- length(x = value)
  validObject(object = object)
  return(object)
}

#' @rdname DefaultLayer
#' @method DefaultLayer<- Assay5
#' @export
#'
"DefaultLayer<-.Assay5" <- `DefaultLayer<-.StdAssay`

#' @rdname Cells-StdAssay
#' @method Features StdAssay
#' @export
#'
Features.StdAssay <- function(x, layer = NULL, simplify = TRUE, ...) {
  if (any(is.na(x = layer)) || is.null(x = layer)) {
    return(rownames(x = slot(object = x, name = 'features')))
  }
  layer <- Layers(object = x, search = layer)
  features <- sapply(
    X = layer,
    FUN = function(lyr) {
      return(slot(object = x, name = 'features')[[lyr]])
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  if (isFALSE(x = simplify)) {
    return(features)
  }
  return(Reduce(f = union, x = features))
}

#' @rdname Cells
#' @method Features Assay5
#' @export
#'
Features.Assay5 <- Features.StdAssay

#' @method FetchData StdAssay
#' @export
#'
FetchData.StdAssay <- function(
  object,
  vars,
  cells = NULL,
  layer = NULL,
  clean = TRUE,
  ...
) {
  # Identify layer(s) to use
  layer.set <- rev(x = Layers(
    object = object,
    search = layer %||% 'data'
  ))
  if (is.null(layer) && length(layer.set) == 1 && layer.set == 'scale.data'){
    warning('Default search for "data" layer yielded no results; utilizing "scale.data" layer instead.')
  }
  if (is.null(layer.set) & is.null(layer) ) {
    warning('data layer is not found and counts layer is used')
    layer.set <- rev(x = Layers(
      object = object,
      search = 'counts'
    ))
  }
  if (is.null(layer.set)) {
  stop('layer "', layer,'" is not found in the object')
  } else {
    layer <- layer.set
  }

  # Identify cells to use
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  cells <- intersect(x = cells, y = colnames(x = object))
  if (!length(x = cells)) {
    abort(message = "None of the cells requested found in this assay")
  }
  # Check vars
  orig <- vars
  vars <- gsub(
    pattern = paste0('^', Key(object = object)),
    replacement = '',
    x = vars
  )
  # Pull expression information
  features <- sapply(
    X = layer,
    FUN = Features,
    x = object,
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  vars <- intersect(x = vars, y = Reduce(f = union, x = features))
  data.fetched <- as.data.frame(x = matrix(
    data = NA_real_,
    nrow = length(x = cells),
    ncol = length(x = vars),
    dimnames = list(cells, vars)
  ))
  for (lyr in layer) {
    lcells <- intersect(x = cells, y = Cells(x = object, layer = lyr))
    lvars <- intersect(x = vars, y = Features(x = object, layer = lyr))
    if (!length(x = lcells) || !length(x = lvars)) {
      next
    }
    data.fetched[lcells, lvars] <- as(t(x = LayerData(
      object = object,
      layer = lyr,
      cells = lcells,
      features = lvars
    )[lvars, lcells, drop = FALSE]),
    "matrix")
  }
  # Clean out missing cells from the expression matrix
  if (isTRUE(x = clean)) {
    no.data <- which(x = apply(
      X = data.fetched,
      MARGIN = 1L,
      FUN = function(x) {
        return(all(is.na(x = x)))
      }
    ))
    if (length(x = no.data)) {
      warn(message = paste(
        "Removing",
        length(x = no.data),
        "cells missing data for features requested"
      ))
      data.fetched <- data.fetched[-no.data, , drop = FALSE]
    }
  }
  # Add keys to keyed vars
  keyed.features <- paste0(Key(object = object), colnames(x = data.fetched))
  keyed.idx <- which(x = keyed.features %in% orig)
  if (length(x = keyed.idx)) {
    colnames(x = data.fetched)[keyed.idx] <- keyed.features[keyed.idx]
  }
  # Check final list of features
  fetched <- setdiff(x = unlist(x = dimnames(x = data.fetched)), y = cells)
  missing <- setdiff(x = orig, y = fetched)
  if (length(x = missing) == length(x = orig)) {
    abort(message = "None of the requested variables found", class = 'varsNotFoundError')
  } else if (length(x = missing)) {
    warn(message = paste(
      "The following variables could not be found:",
      paste(missing, collapse = ', ')
    ))
  }
  return(data.fetched)
  # # Pull feature-level metadata
  # meta.fetch <- c(
  #   grep(pattern = '^md_', x = vars, value = TRUE),
  #   vars[vars %in% colnames(x = object[[]])]
  # )
  # meta.fetch <- setdiff(x = meta.fetch, y = colnames(x = data.fetched))
  # meta.keyed <- which(x = grepl(pattern = '^md', x = meta.fetch))
  # meta.fetch <- gsub(pattern = '^md_', replacement = '', x = meta.fetch)
  # meta.data <- lapply(
  #   X = meta.fetch,
  #   FUN = function(x, f) {
  #     df <- as.data.frame(x = matrix(
  #       data = NA,
  #       nrow = 1L,
  #       ncol = length(x = f),
  #       dimnames = list(x, f)
  #     ))
  #     df[x, ] <- object[[x]][f, , drop = TRUE]
  #     return(df)
  #   },
  #   f = colnames(x = data.fetched)
  # )
  # meta.data <- do.call(what = 'rbind', args = meta.data)
  # if (length(x = meta.keyed)) {
  #   rownames(x = meta.data)[meta.keyed] <- paste0(
  #     'md_',
  #     rownames(x = meta.data)[meta.keyed]
  #   )
  # }
  # keyed.meta <- paste0(Key(object = object), rownames(x = meta.data))
  # keyed.meta.idx <- which(x = keyed.meta %in% orig)
  # if (length(x = keyed.meta.idx)) {
  #   rownames(x = meta.data)[keyed.meta.idx] <- keyed.meta[keyed.meta.idx]
  # }
  # if (nrow(x = data.fetched) && (nrow(x = meta.data) %||% 0)) {
  #   warning(
  #     "Returning both expression and meta data; data types might be different than expected",
  #     call. = FALSE,
  #     immediate. = TRUE
  #   )
  # }
  # data.fetched <- rbind(data.fetched, meta.data)
  # # Add keys to keyed vars
  # keyed.features <- paste0(Key(object = object), colnames(x = data.fetched))
  # keyed.idx <- which(x = keyed.features %in% orig)
  # if (length(x = keyed.idx)) {
  #   colnames(x = data.fetched)[keyed.idx] <- keyed.features[keyed.idx]
  # }
  # # Check final list of features
  # fetched <- setdiff(x = unlist(x = dimnames(x = data.fetched)), y = cells)
  # missing <- setdiff(x = orig, y = fetched)
  # if (length(x = missing) == length(x = orig)) {
  #   stop("None of the requested variables found", call. = FALSE)
  # } else if (length(x = missing)) {
  #   warning(
  #     "The following variables could not be found: ",
  #     paste(missing, collapse = ', '),
  #     call. = FALSE,
  #     immediate. = TRUE
  #   )
  # }
  # return(data.fetched)
}

#' @method FetchData Assay5
#' @export
#'
FetchData.Assay5 <- FetchData.StdAssay

#' @templateVar fxn AssayData
#' @template method-stdassay
#'
#' @method GetAssayData StdAssay
#' @export
#'
GetAssayData.StdAssay <- function(
  object,
  layer = NULL,
  slot = deprecated(),
  ...
) {
  CheckDots(..., fxns = LayerData)
  if (is_present(arg = slot)) {
    .Deprecate(
      when = '5.0.0',
      what = 'GetAssayData(slot = )',
      with = 'GetAssayData(layer = )'
    )
    layer <- slot
  }
  layer_name <- layer[1L] %||% DefaultLayer(object = object)[1L]
  layer.set <- suppressWarnings(expr = Layers(
    object = object,
    search = layer %||% 'data'
  ))
  if (is.null(layer.set) & is.null(layer)) {
    warning('data layer is not found and counts layer is used')
    layer <- rev(x = Layers(
      object = object,
      search = 'counts'
    ))
  } else {
    layer <- layer.set
  }
  if (length(x = layer) > 1) {
    abort("GetAssayData doesn't work for multiple layers in v5 assay.",
         " You can run 'object <- JoinLayers(object = object, layers = layer)'.")
  }
  if (is.null(x = layer)) {
    msg <- paste("Layer", sQuote(x = layer_name), "is empty")
    opt <- getOption(x = "Seurat.object.assay.v3.missing_layer",
                     default = Seurat.options$Seurat.object.assay.v3.missing_layer)
    opt <- tryCatch(expr = arg_match0(
      arg = opt,
      values = c("matrix","null", "error")),
      error = function(...) {
        return(Seurat.options$Seurat.object.assay.v3.missing_layer)
      }
    )
    if (opt == "error") {
      abort(message = msg)
    }
    warn(message = msg)
    return(switch(
      EXPR = opt,
      matrix = switch(
        EXPR = layer_name,
        scale.data = new(Class = "matrix"), new(Class = "dgCMatrix")
      ),
      NULL))
  }
  return(LayerData(object = object, layer = layer, ...))
}

#' @templateVar fxn VariableFeatures
#' @template method-stdassay
#'
#' @importFrom utils adist
#'
#' @method HVFInfo StdAssay
#' @export
#'
HVFInfo.StdAssay <- function(
  object,
  method = NULL,
  status = FALSE,
  layer = NULL,
  strip = TRUE,
  ...
) {
  # Find available HVF methods and layers
  vf.methods.layers <- .VFMethodsLayers(object = object, type = 'hvf')
  #vf.methods <- .VFMethods(object = object, type = 'hvf')
  #vf.layers <- .VFLayers(object = object, type = 'hvf')
  # Determine which method and layer to use
  method <- method[length(methods)] %||% names(vf.methods.layers[length(vf.methods.layers)])
  method <- switch(
    EXPR = tolower(x = method),
    mean.var.plot = 'mvp',
    dispersion = 'disp',
    method
  )
  method <- tryCatch(
    expr = match.arg(arg = method, choices = names(vf.methods.layers)),
    error = function(...) {
      return(NULL)
    }
  )
  # If no methods found, return NULL
  if (is.null(x = method)) {
    return(method)
  }
  vf.methods.layers <- unlist(vf.methods.layers, use.names = FALSE)
  layer <- Layers(object = object, search = layer)
  layer <- vf.methods.layers[which.min(x = adist(x = layer, y = vf.methods.layers))]
  # Find the columns for the specified method and layer
  cols <- grep(
    pattern = paste0(paste('^vf', method, layer, sep = '_'), '_'),
    x = colnames(x = object[[]]),
    value = TRUE
  )
  if (!isTRUE(x = status)) {
    cols <- setdiff(
      x = cols,
      y = paste('vf', method, layer, c('variable', 'rank'), sep = '_')
    )
  }
  hvf.info <- object[[cols]]
  colnames(x = hvf.info) <- gsub(
    pattern = '^vf_',
    replacement = '',
    x = colnames(x = hvf.info)
  )
  if (isTRUE(x = strip)) {
    colnames(x = hvf.info) <- gsub(
      pattern = paste0(paste(method, layer, sep = '_'), '_'),
      replacement = '',
      x = colnames(x = hvf.info)
    )
  }
  return(hvf.info)
}

#' @param layer Layer to pull variable features for
#' @param strip Remove method/layer identifiers from highly variable data frame
#'
#' @rdname VariableFeatures
#' @method HVFInfo Assay5
#' @export
#'
HVFInfo.Assay5 <- HVFInfo.StdAssay

#' @method JoinLayers StdAssay
#' @export
#'
JoinLayers.StdAssay <- function(
  object,
  layers = NULL,
  new = NULL,
  ...
) {
  layers <- layers %||% c('counts', 'data', 'scale.data')
  new <- new %||% layers
  if (length(x = layers) != length(x = new)) {
    stop('Number of layers and new should be the same')
  }
  var.features <- VariableFeatures(object = object)
  for (i in seq_along(layers)) {
    num.layers <- suppressWarnings(
      expr = length(x = Layers(object = object, search = layers[i]))
      )
    if (num.layers > 0L) {
      object <- JoinSingleLayers(
        object = object,
        layers = layers[i],
        new = new[i],
        default = TRUE,
        ...
      )
    }
  }
  VariableFeatures(object = object) <- var.features
 return(object)
}

#' @param layers Names of layers to split or join
#' @param new Name of new layers
#'
#' @rdname SplitLayers
#'
#' @method JoinLayers Assay5
#' @export
#'
JoinLayers.Assay5 <- JoinLayers.StdAssay

#' @rdname Key
#' @method Key Assay5
#' @export
#'
Key.Assay5 <- .Key

#' @rdname Key
#' @method Key<- Assay5
#' @export
#'
"Key<-.Assay5" <- `.Key<-`

#' @templateVar fxn Layers
#' @template method-stdassay
#'
#' @method LayerData StdAssay
#' @export
#'
LayerData.StdAssay <- function(
  object,
  layer = NULL,
  cells = NULL,
  features = NULL,
  fast = FALSE,
  slot = deprecated(),
  ...
) {
  if (is_present(arg = slot)) {
    deprecate_stop(
      when = '5.0.0',
      what = 'LayerData(slot = )',
      with = 'LayerData(layer = )"'
    )
  }
  layer_name <- layer[1L] %||% DefaultLayer(object = object)[1L]
  # Identify layer(s) to use
  layer.set <- suppressWarnings(expr = Layers(
    object = object,
    search = layer %||% 'data'
  ))
  # If layer.set doesnt return anything and layer is not defined
  if (is.null(layer.set) & is.null(layer) ) {
    warning(
      'data layer is not found and counts layer is used',
      call. = F,
      immediate. = T
    )
    layer <- Layers(
      object = object,
      search = 'counts'
    )
  } else {
    layer <- layer.set
  }

  if (length(x = layer) > 1) {
    warning("multiple layers are identified by ",
            paste0(layer, collapse = ' '),
            "\n only the first layer is used")
    layer <- layer[1L]
  }
  # layer <- match.arg(arg = layer, choices = Layers(object = object))
  if (is.null(x = layer) || any(is.na(x = layer))) {
    msg <- paste("Layer", sQuote(x = layer_name), "is empty")
    opt <- getOption(x = "Seurat.object.assay.v3.missing_layer",
                     default = Seurat.options$Seurat.object.assay.v3.missing_layer)
    opt <- tryCatch(expr = arg_match0(
      arg = opt,
      values = c("matrix","null", "error")),
      error = function(...) {
        return(Seurat.options$Seurat.object.assay.v3.missing_layer)
        }
      )
    if (opt == "error") {
      abort(message = msg)
    }
    warn(message = msg)
    return(switch(
      EXPR = opt,
      matrix = switch(
        EXPR = layer_name,
        scale.data = new(Class = "matrix"), new(Class = "dgCMatrix")
        ),
      NULL))
  }
  # Allow cell/feature subsets
  dnames <- list(
    Features(x = object, layer = layer),
    Cells(x = object, layer = layer)
  )
  cells <- cells %||% dnames[[2L]]
  if (is.numeric(x = cells)) {
    cells <- dnames[[2L]][cells]
  }
  cells <- sort(x = MatchCells(
    new = dnames[[2L]],
    orig = cells,
    ordered = TRUE
  ))
  dnames[[2L]] <- dnames[[2L]][cells]
  features <- features %||% dnames[[1L]]
  if (is.numeric(x = features)) {
    features <- dnames[[1L]][features]
  }
  features <- sort(x = MatchCells(
    new = dnames[[1L]],
    orig = features,
    ordered = TRUE
  ))
  dnames[[1L]] <- dnames[[1L]][features]
  if(length(x = dnames[[1L]]) == 0) {
    stop('features are not found')
  }
  # Pull the layer data
  ldat <- if (.MARGIN(x = object) == 1L) {
    methods::slot(object = object, name = 'layers')[[layer]][features, cells, drop = FALSE]
  } else {
    methods::slot(object = object, name = 'layers')[[layer]][cells, features, drop = FALSE]
  }
  # Add dimnames and transpose if requested
  ldat <- if (isTRUE(x = fast)) {
    ldat
  } else if (is_na(x = fast)) {
    .GetLayerData2(x = ldat, dnames = dnames, fmargin = 1L)
    # .GetLayerData(
    #   x = ldat,
    #   dnames = dnames,
    #   fmargin = 1L,
    #   ...
    # )
  } else {
    .GetLayerData2(
      x = ldat,
      dnames = dnames,
      fmargin = .MARGIN(x = object, type = 'features')
    )
    # .GetLayerData(
    #   x = ldat,
    #   dnames = dnames,
    #   fmargin = .MARGIN(x = object, type = 'features'),
    #   ...
    # )
  }
  return(ldat)
}

#' @param features,cells Vectors of features/cells to include
#' @param fast Determine how to return the layer data; choose from:
#' \describe{
#'  \item{\code{FALSE}}{Apply any transpositions and attempt to add
#'   feature/cell names (if supported) back to the layer data}
#'  \item{\code{NA}}{Attempt to add feature/cell names back to the layer data,
#'   skip any transpositions}
#'  \item{\code{TRUE}}{Do not apply any transpositions or add feature/cell
#'   names to the layer data}
#' }
#'
#' @rdname Layers
#' @method LayerData Assay5
#' @export
#'
LayerData.Assay5 <- LayerData.StdAssay

#'
#' @rdname Layers-StdAssay
#' @method LayerData<- StdAssay
#' @export
#'
"LayerData<-.StdAssay" <- function(
  object,
  layer,
  features = NULL,
  cells = NULL,
  ...,
  value
) {
  if (!is_scalar_character(x = layer) || !nzchar(x = layer)) {
    abort(message = "'layer' must be a single non-empty string")
  }
  # Remove a layer
  if (is.null(x = value)) {
    if (length(x = Layers(object = object)) == 1L) {
      stop("Cannot remove only layer")
    } else if (layer %in% DefaultLayer(object = object)) {
      msg <- 'Removing default layer'
      if (length(x = DefaultLayer(object = object)) == 1L) {
        DefaultLayer(object = object) <- Layers(object = object)[2]
        msg <- paste0(
          msg,
          ', setting default to ', DefaultLayer(object = object)
        )
      } else {
        didx <- slot(object = object, name = 'default') - 1L
        slot(object = object, name = 'default') <- didx
      }
      warning(msg, call. = FALSE, immediate. = TRUE)
    }
    slot(object = object, name = 'layers')[[layer]] <- NULL
    if (slot(object = object, name = 'default') > length(x = Layers(object = object)) ||
        !length(x = slot(object = object, name = 'default'))) {
      slot(object = object, name = 'default') <- length(x = Layers(object = object))
    }
    maps <- c(
      'cells',
      'features'
    )
    for (i in maps) {
      slot(object = object, name = i)[[layer]] <- NULL
    }
    validObject(object = object)
    return(object)
  }
  # Add a layer
  fdim <- .MARGIN(x = object, type = 'features')
  cdim <- .MARGIN(x = object, type = 'cells')
  # Assume input matrix is features x cells
  dnames <- list(
    features %||% dimnames(x = value)[[1L]],
    cells %||% dimnames(x = value)[[2L]]
  )
  if (length(x = unique(x = dim(x = value))) > 1L) {
    didx <- match(
      x = vapply(
        X = dnames,
        FUN = length,
        FUN.VALUE = numeric(length = 1L),
        USE.NAMES = FALSE
      ),
      table = dim(x = value)
    )
    dnames <- dnames[didx]
  }
  value <- .PrepLayerData2(
    x = value,
    target = dim(x = object),
    dnames = dnames,
    fmargin = fdim,
    ...
  )
  # value <- .PrepLayerData(
  #   x = value,
  #   target = dim(x = object),
  #   dnames = dnames,
  #   fmargin = fdim,
  #   ...
  # )
  # Check features and cells
  features <- attr(x = value, which = 'features') %||% seq_len(length.out = dim(x = value)[fdim])
  cells <- attr(x = value, which = 'cells') %||% seq_len(length.out = dim(x = value)[cdim])
  fmatch <- MatchCells(
    new = features,
    orig = rownames(x = slot(object = object, name = 'features')),
    ordered = TRUE
  )
  cmatch <- MatchCells(
    new = cells,
    orig = rownames(x = slot(object = object, name = 'cells')),
    ordered = TRUE
  )
  if (is.null(x = fmatch)) {
    stop(
      "No feature overlap between existing object and new layer data",
      call. = FALSE
    )
  } else if (is.null(x = cmatch)) {
    stop(
      "No cell overlap between existing object and new layer data",
      call. = FALSE
    )
  }
  features <- features[fmatch]
  cells <- cells[cmatch]
  # Check for existing layer data
  if (layer %in% Layers(object = object)) {
    fcheck <- if (is.numeric(x = features)) {
      Features(x = object, layer = layer)[features]
    } else {
      features
    }
    ccheck <- if (is.numeric(x = cells)) {
      Cells(x = object, layer = layer)[cells]
    } else {
      cells
    }
    if (!identical(x = fcheck, y = Features(x = object, layer = layer))) {
      warning(
        "Different features in new layer data than already exists for ",
        layer,
        call. = FALSE,
        immediate. = TRUE
      )
    }
    if (!identical(x = ccheck, y = Cells(x = object, layer = layer))) {
      warning(
        "Different cells in new layer data than already exists for ",
        layer,
        call. = FALSE,
        immediate. = TRUE
      )
    }
  }
  # Reorder the layer data
  value <- if (fdim == 1L) {
    value[fmatch, cmatch]
  } else {
    value[cmatch, fmatch]
  }
  # Add the layer
  slot(object = object, name = 'layers')[[layer]] <- value
  # Update the maps
  slot(object = object, name = 'features')[[layer]] <- features
  slot(object = object, name = 'cells')[[layer]] <- cells
  validObject(object = object)
  return(object)
}

#' @rdname Layers
#' @method LayerData<- Assay5
#' @export
#'
"LayerData<-.Assay5" <- `LayerData<-.StdAssay`

#' @rdname Layers-StdAssay
#' @method Layers StdAssay
#' @export
#'
Layers.StdAssay <- function(object, search = NA, ...) {
  if (is.null(x = search)) {
    return(DefaultLayer(object = object))
  }
  layers <- names(x = slot(object = object, name = 'layers'))
  if (!is_na(x = search)) {
    layers <- unique(x = unlist(x = lapply(
      X = search,
      FUN = function(lyr) {
        if (lyr %in% layers) {
          return(lyr)
        }
        patterns <- c(paste0('^', lyr), paste0(lyr, '$'), lyr)
        res <- vector(mode = 'character')
        for (p in patterns) {
          res <- grep(pattern = p, x = layers, value = TRUE, ...)
          if (length(x = res)) {
            break
          }
        }
        return(res)
      }
    )))
    if (!length(x = layers)) {
      warning(message = "No layers found matching search pattern provided",
              call. = FALSE,
              immediate. = TRUE)
      return(NULL)
    }
  }
  return(layers)
}

#' @param search A pattern to search layer names for; pass one of:
#' \itemize{
#'  \item \dQuote{\code{NA}} to pull all layers
#'  \item \dQuote{\code{NULL}} to pull the default layer(s)
#'  \item a \link[base:grep]{regular expression} that matches layer names
#' }
#'
#' @rdname Layers
#' @method Layers Assay5
#' @export
#'
Layers.Assay5 <- Layers.StdAssay

#' @templateVar fxn Misc
#' @template method-stdassay
#'
#' @method Misc StdAssay
#' @export
#'
Misc.StdAssay <- .Misc

#' @rdname Misc
#' @method Misc Assay5
#' @export
#'
Misc.Assay5 <- .Misc

#' @templateVar fxn Misc
#' @template method-stdassay
#'
#' @method Misc<- StdAssay
#' @export
#'
"Misc<-.StdAssay" <- `.Misc<-`

#' @rdname Misc
#' @method Misc<- Assay5
#' @export
#'
"Misc<-.Assay5" <- `.Misc<-`

#' @templateVar fxn RenameCells
#' @template method-stdassay
#'
#' @method RenameCells StdAssay
#' @export
#'
RenameCells.StdAssay <- function(object, new.names = NULL, ...) {
  CheckDots(...)
  colnames(object) <- new.names[colnames(object)]
  return(object)
}

#' @rdname RenameCells
#' @method RenameCells Assay5
#' @export
#'
RenameCells.Assay5 <- RenameCells.StdAssay

#' @rdname AssayData-StdAssay
#' @method SetAssayData StdAssay
#' @export
#'
SetAssayData.StdAssay <- function(
  object,
  layer,
  new.data,
  slot = deprecated(),
  ...
) {
  if (is_present(arg = slot)) {
    .Deprecate(
      when = '5.0.0',
      what = 'SetAssayData(slot = )',
      with = 'SetAssayData(layer = )'
    )
    layer <- slot
  }
  LayerData(object = object, layer = layer) <- new.data
  return(object)
}

#' @rdname VariableFeatures-StdAssay
#' @method VariableFeatures StdAssay
#' @export
#'
VariableFeatures.StdAssay <- function(
  object,
  method = NULL,
  layer = NA,
  simplify = TRUE,
  nfeatures = Inf,
  selection.method = deprecated(),
  ...
) {
  if (is_present(arg = selection.method)) {
    .Deprecate(
      when = '5.0.0',
      what = 'VariableFeatures(selection.method = )',
      with = 'VariableFeatures(method = )'
    )
    method <- selection.method
  }
  nfeatures <- nfeatures %||% Inf
  if ("var.features" %in% colnames(object[[]])) {
    if ("var.features.rank" %in% colnames(object[[]])) {
      var.features <- row.names(x = object[[]])[which(!is.na(object[[]]$var.features.rank))]
      var.features <- var.features[order(object[[]][["var.features.rank"]][which(!is.na(object[[]]$var.features))])]
    }
    else {
      var.features <- as.vector(object["var.features", drop = TRUE])
      var.features <- var.features[!is.na(var.features)]
    }
    if (isTRUE(x = simplify) & (is.null(x = layer) || any(is.na(x = layer))) &
        (is.infinite(x = nfeatures) || length(x = var.features) ==
         nfeatures)) {
      return(var.features)
    }
  }
  msg <- 'No variable features found'
  layer.orig <- layer
  methods <- .VFMethodsLayers(object = object, type = 'hvf', layers = layer)
  layer <- Layers(object = object, search = layer)
  method <- method %||% names(x = methods)[length(x = methods)]
  method <- match.arg(arg = method, choices = names(x = methods))
  if (is_na(x = layer.orig) || is.null(x = layer.orig)) {
    layer <- unlist(methods[method], use.names = FALSE)
  }
  vf <- sapply(
    X = layer,
    FUN = function(lyr) {
      hvf.info <- HVFInfo(
        object = object,
        method = method,
        layer = lyr,
        status = TRUE,
        strip = TRUE
      )
      if (is.null(x = hvf.info)) {
        return(NULL)
      } else if (!'variable' %in% names(x = hvf.info)) {
        return(NA)
      }
      vf <- row.names(x = hvf.info)[which(x = hvf.info$variable)]
      if ('rank' %in% names(x = hvf.info)) {
        vf <- vf[order(hvf.info$rank[which(x = hvf.info$variable)])]
      } else {
        warn(message = paste0(
          "No variable feature rank found for ",
          sQuote(x = lyr),
          ", returning features in assay order"
        ))
      }
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  if (is.null(x = unlist(x = vf))) {
    return(NULL)
  } else if (all(is.na(x = unlist(x = vf)))) {
    abort(message = msg)
  }
  if (isTRUE(x = simplify)) {
    vf <- .SelectFeatures(
      object = vf,
      all.features = intersect(
        x = slot(object = object, name = 'features')[,layer]
      ),
      nfeatures = nfeatures
    )
  }
  return(vf)
  # hvf.info <- HVFInfo(
  #   object = object,
  #   method = method,
  #   layer = layer,
  #   status = TRUE,
  #   strip = TRUE
  # )
  # if (is.null(x = hvf.info)) {
  #   warning(msg, call. = FALSE, immediate. = TRUE)
  #   return(NULL)
  # }
  # if (!'variable' %in% names(x = hvf.info)) {
  #   stop(msg, call. = FALSE)
  # }
  # vf <- rownames(x = hvf.info)[which(x = hvf.info$variable)]
  # if ('rank' %in% names(x = hvf.info)) {
  #   vf <- vf[order(hvf.info$rank[which(x = hvf.info$variable)])]
  # } else {
  #   warning(
  #     "No variable feature rank found, returning features in assay order",
  #     call. = FALSE,
  #     immediate. = TRUE
  #   )
  # }
  # return(vf)
}

#' @param simplify When pulling for multiple layers, combine into a single
#' vector and select a common set of variable features for all layers
#' @param nfeatures Maximum number of features to select when simplifying
#'
#' @rdname VariableFeatures
#' @method VariableFeatures Assay5
#' @export
#'
VariableFeatures.Assay5 <- VariableFeatures.StdAssay

#' @rdname VariableFeatures-StdAssay
#' @method VariableFeatures<- StdAssay
#' @export
#'
"VariableFeatures<-.StdAssay" <- function(
  object,
  method = 'custom',
  layer = NULL,
  ...,
  value
) {
  if (!length(x = value)) {
    return(object)
  }
  value <- intersect(x = value, y = rownames(x = object))
  if (!length(x = value)) {
    stop("None of the features specified are present in this assay", call. = FALSE)
  }
  object[['var.features']] <- value
  # add rank
  object[['var.features.rank']] <- NA
  object[[]][row.names(object[[]]) %in% value,]$var.features.rank <- match(row.names(object[[]])[row.names(object[[]]) %in% value], value)

  # layer <- Layers(object = object, search = layer)
  # df <- data.frame(TRUE, seq_along(along.with = value), row.names = value)
  # for (lyr in layer) {
  #   names(x = df) <- paste('vf', method, lyr, c('variable', 'rank'), sep = '_')
  #   object[] <- df
  # }
  return(object)
}

#' @rdname VariableFeatures
#' @method VariableFeatures<- Assay5
#' @export
#'
"VariableFeatures<-.Assay5" <- `VariableFeatures<-.StdAssay`

#' @rdname VariableFeatures
#' @export
#' @method SVFInfo StdAssay
#'
SVFInfo.StdAssay <- function(
  object,
  method = c("markvariogram", "moransi"),
  status = FALSE,
  selection.method = deprecated(),
  ...
) {
  CheckDots(...)
  if (is_present(arg = selection.method)) {
    .Deprecate(
      when = '5.0.0',
      what = 'SVFInfo(selection.method = )',
      with = 'SVFInfo(method = )'
    )
    method <- selection.method
  }
  method <- match.arg(arg = method)
  vars <- switch(
    EXPR = method,
    markvariogram = grep(
      pattern = "r.metric",
      x = colnames(x = object[[]]),
      value = TRUE
    ),
    moransi = grep(
      pattern = 'MoransI',
      x = colnames(x = object[[]]),
      value = TRUE
    ),
    abort(message = paste("Unknown method:", sQuote(x = method)))
  )

  tryCatch(
    expr = svf.info <- object[[vars]],
    error = function(e) {
      stop(
        "Unable to find spatially variable feature information for method '",
        method,
        "'",
        call. = FALSE
      )
    }
  )

  colnames(x = svf.info) <- vars
  if (status) {
    svf.info$variable <- object[[paste0(method, '.spatially.variable')]]
    svf.info$rank <- object[[paste0(method, '.spatially.variable.rank')]]
  }
  return(svf.info)
}

#' @rdname SpatiallyVariableFeatures-StdAssay
#' @method SpatiallyVariableFeatures StdAssay
#' @export
#'
SpatiallyVariableFeatures.StdAssay <- function(
  object,
  method = "moransi",
  decreasing = TRUE,
  selection.method = deprecated(),
  ...
) {
  CheckDots(...)
  if (is_present(arg = selection.method)) {
    .Deprecate(
      when = '5.0.0',
      what = 'SpatiallyVariableFeatures(selection.method = )',
      with = 'SpatiallyVariableFeatures(method = )'
    )
    method <- selection.method
  }
  vf <- SVFInfo(object = object, method = method, status = TRUE)
  vf <- vf[rownames(x = vf)[which(!is.na(x = vf[, "variable"]))], ]
  if (!is.null(x = decreasing)) {
    vf <- vf[order(x = vf[, "rank"][, 1], decreasing = !decreasing), ]
  }
  return(rownames(x = vf))
}

#' @rdname SpatiallyVariableFeatures
#' @method SpatiallyVariableFeatures Assay5
#' @export
#'
SpatiallyVariableFeatures.Assay5 <- SpatiallyVariableFeatures.StdAssay

#' @rdname VariableFeatures
#' @method SVFInfo Assay5
#' @export
#'
SVFInfo.Assay5 <- SVFInfo.StdAssay

#' @method WhichCells StdAssay
#' @export
#'
WhichCells.StdAssay <- WhichCells.Assay

# WhichCells.StdAssay <- function(
#   object,
#   cells = NULL,
#   expression = missing_arg(),
#   invert = FALSE,
#   ...
# ) {
#   cells <- cells %||% colnames(x = object)
#   if (!is_missing(x = expression) && !is.null(x = substitute(expr = expression))) {
#     key.pattern <- paste0('^', Key(object = object))
#     expr <- if (tryCatch(expr = is_quosure(x = expression), error = \(...) FALSE)) {
#       expression
#     } else if (is.call(x = enquo(arg = expression))) {
#       enquo(arg = expression)
#     } else {
#       parse(text = expression)
#     }
#     expr.char <- suppressWarnings(expr = as.character(x = expr))
#     expr.char <- unlist(x = lapply(X = expr.char, FUN = strsplit, split = ' '))
#     expr.char <- gsub(
#       pattern = key.pattern,
#       replacement = '',
#       x = expr.char,
#       perl = TRUE
#     )
#     expr.char <- gsub(
#       pattern = '(',
#       replacement = '',
#       x = expr.char,
#       fixed = TRUE
#     )
#     expr.char <- gsub(
#       pattern = '`',
#       replacement = '',
#       x = expr.char
#     )
#   }
#   if (isTRUE(x = invert)) {
#     cells <- setdiff(x = colnames(x = object), y = cells)
#   }
#   cells <- ''
#   return(as.character(x = cells))
# }

#' @method WhichCells Assay5
#' @export
#'
WhichCells.Assay5 <- WhichCells.StdAssay

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @inherit .DollarNames.Assay5 params return title description details sections
#'
#' @importFrom utils .DollarNames
#'
#' @keywords internal
#' @method .DollarNames StdAssay
#' @export
#'
#' @family stdassay
#'
.DollarNames.StdAssay <- function(x, pattern = '') {
  layers <- as.list(x = Layers(object = x))
  names(x = layers) <- unlist(x = layers)
  return(.DollarNames(x = layers, pattern = pattern))
}

#' Dollar-sign Autocompletion
#'
#' Autocompletion for \code{$} access on an \code{\link{Assay5}} object
#'
#' @inheritParams [.Assay5
#' @inheritParams utils::.DollarNames
#'
#' @return The layer name matches for \code{pattern}
#'
#' @importFrom utils .DollarNames
#'
#' @keywords internal
#'
#' @method .DollarNames Assay5
#' @export
#'
#' @concept assay5
#'
#' @seealso \code{\link[utils:.DollarNames]{utils::.DollarNames}}
#'
.DollarNames.Assay5 <- .DollarNames.StdAssay

#' @inherit $.Assay5 params return title description details sections
#'
#' @keywords internal
#' @method $ StdAssay
#' @export
#'
#' @family stdassay
#'
"$.StdAssay" <- function(x, i) {
  return(LayerData(object = x, layer = i))
}

#' Layer Data
#'
#' Get and set layer data
#'
#' @inheritParams [.Assay5
#'
#' @return {$}: Layer data for layer \code{i}
#'
#' @method $ Assay5
#' @export
#'
#' @family assay5
#'
"$.Assay5" <- `$.StdAssay`


#' @rdname cash-.StdAssay
#'
#' @method $<- StdAssay
#' @export
#'
"$<-.StdAssay" <- `$<-.Assay`

#' @return \code{$<-}: \code{x} with layer data \code{value} saved as \code{i}
#'
#' @rdname cash-.Assay5
#'
#' @method $<- Assay5
#' @export
#'
"$<-.Assay5" <- `$<-.StdAssay`

#' @inherit [.Assay5 params return title description details sections
#'
#' @keywords internal
#' @method [ StdAssay
#' @export
#'
#' @family stdassay
#'
"[.StdAssay" <- `[.Assay`

#' Layer Data
#'
#' Get and set layer data
#'
#' @inheritParams [[.Assay5
#' @param i Name of layer data to get or set
#' @param ... Arguments passed to \code{\link{LayerData}}
#'
#' @return \code{[}: The layer data for layer \code{i}
#'
#' @method [ Assay5
#' @export
#'
#' @family assay5
#'
#' @seealso \code{\link{LayerData}}
#'
#' @order 1
#'
"[.Assay5" <- `[.StdAssay`

#' @inherit [[.Assay5 params return title description details sections
#'
#' @keywords internal
#' @method [[ StdAssay
#' @export
#'
#' @family stdassay
#'
"[[.StdAssay" <- function(x, i, j, ..., drop = FALSE) {
  if (missing(x = i)) {
    i <- colnames(x = slot(object = x, name = 'meta.data'))
  }
  data.return <- slot(object = x, name = 'meta.data')[, i, drop = FALSE, ...]
  if (nrow(x = data.return) == 0) {
    return(data.return)
  }
  row.names(x = data.return) <- rownames(x = x)
  if (isTRUE(x = drop)) {
    data.return <- unlist(x = data.return, use.names = FALSE)
    names(x = data.return) <- rep.int(
      x = rownames(x = x),
      times = length(x = i)
    )
  }
  return(data.return)
}

#' Feature-Level Meta Data
#'
#' Get and set feature-level meta data
#'
#' @param x An \code{\link{Assay5}} object
#' @param i Name of feature-level meta data to fetch or add
#' @param j Ignored
#' @param drop See \code{\link{drop}}
#' @template param-dots-ignored
#'
#' @return \code{[[}: The feature-level meta data for \code{i}
#'
#' @method [[ Assay5
#' @export
#'
#' @family assay5
#'
#' @order 1
#'
"[[.Assay5" <- `[[.StdAssay`

#' @inherit dim.Assay5 params return title description details sections
#'
#' @keywords internal
#' @method dim StdAssay
#' @export
#'
#' @family stdassay
#'
dim.StdAssay <- function(x) {
  return(vapply(
    X = c('features', 'cells'),
    FUN = function(s) {
      return(nrow(x = slot(object = x, name = s)))
    },
    FUN.VALUE = numeric(length = 1L),
    USE.NAMES = FALSE
  ))
}

#' Feature and Cell Numbers
#'
#' @inheritParams [[.Assay5
#'
#' @return A two-length numeric vector with the total number of
#' features and cells in \code{x}
#'
#' @method dim Assay5
#' @export
#'
#' @family assay5
#'
dim.Assay5 <- dim.StdAssay

#' @inherit dimnames.Assay5 params return title description details sections
#'
#' @keywords internal
#' @method dimnames StdAssay
#' @export
#'
#' @seealso \code{\link{Cells}} \code{\link{Features}}
#' @family stdassay
#'
dimnames.StdAssay <- function(x) {
  return(list(Features(x = x, layer = NA), Cells(x = x, layer = NA)))
}

#' Assay-Level Feature and Cell Names
#'
#' Get and set feature and cell names in v5 Assays
#'
#' @inheritParams [[.Assay5
#'
#' @return \code{dimnames}: A two-length list with the following values:
#' \itemize{
#'  \item A character vector with all features in \code{x}
#'  \item A character vector with all cells in \code{x}
#' }
#'
#' @method dimnames Assay5
#' @export
#'
#' @family assay5
#' @family dimnames
#'
dimnames.Assay5 <- dimnames.StdAssay

#' @rdname dimnames.StdAssay
#'
#' @method dimnames<- StdAssay
#' @export
#'
"dimnames<-.StdAssay" <- function(x, value) {
  msg <- "Invalid 'dimnames' given for an assay"
  if (!is_bare_list(x = value, n = 2L)) {
    stop(msg, call. = FALSE)
  } else if (!all(sapply(X = value, FUN = length) == dim(x = x))) {
    stop(msg, call. = FALSE)
  }
  value <- lapply(X = value, FUN = as.character)
  rownames(x = slot(object = x, name = 'features')) <- value[[1L]]
  rownames(x = slot(object = x, name = 'cells')) <- value[[2L]]
  validObject(object = x)
  return(x)
}

#' @param value A two-length list with updated feature and/or cells names
#'
#' @return \code{dimnames<-}: \code{x} with the feature and/or cell
#' names updated to \code{value}
#'
#' @rdname dimnames.Assay5
#'
#' @method dimnames<- Assay5
#' @export
#'
"dimnames<-.Assay5" <- `dimnames<-.StdAssay`

#' @rdname sub-sub-.StdAssay
#'
#' @method head StdAssay
#' @export
#'
head.StdAssay <- head.Assay

#' @param n Number of meta data rows to show
#'
#' @return \code{head}: The first \code{n} rows of feature-level meta data
#'
#' @rdname sub-sub-.Assay5
#'
#' @method head Assay5
#' @export
#'
head.Assay5 <- head.StdAssay

#' @inherit merge.Assay5 params return title description details sections
#'
#' @note All assays must be of the same type; merging different v5 assays (eg.
#' \code{\link{Assay5}} and \code{\link{Assay5T}}) is currently unsupported
#'
#' @keywords internal
#' @method merge StdAssay
#' @export
#'
merge.StdAssay <- function(
  x,
  y,
  labels = NULL,
  add.cell.ids = NULL,
  collapse = FALSE,
  ...
) {
  assays <- c(x, y)
  for (i in seq_along(assays)) {
    if (inherits(x = assays[[i]], what = 'Assay')) {
      assays[[i]] <- as(object = assays[[i]], Class = "Assay5") # TODO: support Assay5T
    }
      }
  labels <- labels %||% as.character(x = seq_along(along.with = assays))
  # add.cell.ids <- add.cell.ids %||% labels
  # TODO: Support collapsing layers
  if (isTRUE(x = collapse)) {
    abort(message = "Collapsing layers is not yet supported")
  }
  for (i in seq_along(along.with = assays)) {
    if (is_na(x = labels[i])) {
      labels[i] <- as.character(x = i)
    }
    if (is_na(x = add.cell.ids[i])) {
      add.cell.ids[i] <- as.character(x = i)
    }
    if (!is.null(x = add.cell.ids[i])) {
      colnames(x = assays[[i]]) <- paste(
        colnames(x = assays[[i]]),
        add.cell.ids[i], sep = '.'
      )
    }
  }
  features.all <- LogMap(y = Reduce(
    f = union,
    x = lapply(X = assays, FUN = rownames)
  ))
  combined <- new(
    Class = class(x = x),
    layers = list(),
    cells = LogMap(y = Reduce(
      f = union,
      x = lapply(X = assays, FUN = colnames)
    )),
    features = features.all,
    meta.data = EmptyDF(n = nrow(x = features.all)),
    misc = list(),
    key = Key(object = x) %||% character(length = 0L)
  )

  # Add layers
  # TODO: Support collapsing layers
  if (isTRUE(x = collapse)) {
    abort(message = "Collapsing layers is not yet supported")
  } else {
    # Get default layer as default of first assay
    default <- DefaultLayer(assays[[1]])
    for (i in seq_along(along.with = assays)) {
      for (lyr in Layers(object = assays[[i]])) {
        LayerData(
          object = combined,
          layer = paste(lyr, labels[i], sep = '.'),
          features = Features(x = assays[[i]], layer = lyr),
          cells = Cells(x = assays[[i]], layer = lyr)
        ) <- LayerData(object = assays[[i]], layer = lyr, fast = TRUE)
      }
    }
  }

  # Add feature-level metadata
  for (i in seq_along(along.with = assays)) {
    # Rename HVF columns
    mf <- assays[[i]][[]]
    if (!ncol(x = mf)) {
      next
    }
    for (type in c('vf')) {
      vf.idx <- grep(pattern = paste0('^', type, '_'), x = names(x = mf))
      if (length(x = vf.idx)) {
        names(x = mf)[vf.idx] <- vapply(
          X = names(x = mf)[vf.idx],
          FUN = function(vf) {
            vf <- unlist(x = strsplit(x = vf, split = '_'))
            vf <- paste(
              paste(vf[1:2], collapse = '_'),
              paste(
                paste(vf[3:(length(x = vf) - 1L)], collapse = '_'),
                labels[i],
                sep = '.'
              ),
              vf[length(x = vf)],
              sep = '_'
            )
          },
          FUN.VALUE = character(length = 1L)
        )
      }
    }
    combined[[]] <- mf
  }
  # TODO: Add misc
  DefaultLayer(combined) <- Layers(object = combined, search = default)
  validObject(object = combined)
  return(combined)
}

#' Merge Assays
#'
#' Merge one or more v5 assays together
#'
#' \strong{Note}: collapsing layers is currently not supported
#'
#' @inheritParams [.Assay5
#' @template param-dots-ignored
#' @param y One or more \code{\link{Assay5}} objects
#' @param labels A character vector equal to the number of objects; defaults to
#' \code{as.character(seq_along(c(x, y)))}
#' @param add.cell.ids A character vector equal to the number of objects
#' provided to append to all cell names; if \code{TRUE}, uses \code{labels} as
#' \code{add.cell.ids}
#' @param collapse If \code{TRUE}, merge layers of the same name together; if
#' \code{FALSE}, appends \code{labels} to the layer name
#'
#' @return A new v5 assay with data merged from \code{c(x, y)}
#'
#' @method merge Assay5
#' @export
#'
#' @family assay5
#'
merge.Assay5 <- merge.StdAssay

#' @inherit split.Assay5 params return title description details sections
#'
#' @keywords internal
#' @method split StdAssay
#' @export
#'
#' @family stdassay
#'
split.StdAssay <- function(
  x,
  f,
  drop = FALSE,
  layers = c("counts", "data"),
  ret = c('assay', 'multiassays', 'layers'),
  ...
) {
  op <- options(Seurat.object.assay.brackets = 'v5')
  on.exit(expr = options(op))
  ret <- ret[1L]
  ret <- match.arg(arg = ret)
  layers.to.split <- Layers(object = x, search = layers)
  if (!identical(Layers(object = x), layers.to.split)) {
     message(
       'Splitting ',
       paste(sQuote(x = layers.to.split), collapse = ', '),
       ' layers. Not splitting ',
       paste(
         sQuote(x = setdiff(Layers(object = x), layers.to.split)),
         collapse = ', '
       ),
       '. If you would like to split other layers, set in `layers` argument.'
     )
  }
  layers <- Layers(object = x, search = layers)
  layers.split <- list()
  for (i in seq_along(along.with = layers)) {
    if (length(x = colnames(x = x[layers[i]])) != length(x = colnames(x = x))) {
      layers.split[[i]] <- layers[i]
    }
  }
  layers.split <- unlist(x = layers.split)
  if (length(x = layers.split)) {
    abort(message = paste(
      strwrap(x = paste(
        "The following layers are already split:",
        paste(sQuote(x = layers.split), collapse = ', '),
        "\nPlease join before splitting"
      ))
    ))
  }
  default <- ifelse(
    test = DefaultLayer(object = x) %in% layers,
    yes = DefaultLayer(object = x),
    no = layers[1L]
  )
  cells <- Cells(x = x, layer = layers)
  if (is_named(x = f)) {
    f <- f[cells]
  }
  if (length(x = f) != length(x = cells)) {
    abort(message = "Not enough splits for this assay")
  }
  if (any(is.na(x = f))) {
    f <- factor(x = f, levels = c(unique(as.character(f)), 'na'))
    f[is.na(x = f)] <- 'na'
  } else {
    f <- factor(x = f, levels = unique(x = as.character(x = f)))
  }
  splits <- split(x = cells, f = f, drop = drop)
  names(x = splits) <- .MakeNames(x = names(x = splits))
  return(switch(
    EXPR = ret,
    assay = {
      for (lyr in layers) {
        p <- progressor(steps = length(x = splits))
        p(
          message = paste(
            'Splitting layer',
            sQuote(x = lyr),
            'into',
            length(x = splits),
            'splits'
          ),
          class = 'sticky',
          amount = 0
        )
        lcells <- Cells(x = x, layer = lyr)
        for (i in seq_along(along.with = splits)) {
          p(
            message = paste(
              'Creating split for',
              sQuote(x = names(x = splits)[i])
            ),
            class = 'sticky',
            amount = 0
          )
          group <- paste(lyr, names(x = splits)[i], sep = '.')
          xcells <- intersect(x = splits[[i]], y = lcells)
          LayerData(object = x, layer = group, cells = xcells) <- LayerData(
            object = x,
            layer = lyr,
            cells = xcells
          )
          p()
        }
        p(type = 'finish')
        suppressWarnings(expr = LayerData(object = x, layer = lyr) <- NULL)
        DefaultLayer(object = x) <- default
      }
      x
    },
    multiassays = {
      value <- vector(mode = 'list', length = length(x = splits))
      names(x = value) <- names(x = splits)
      for (group in names(x = splits)) {
        value[[group]] <- subset(
          x = x,
          cells = splits[[group]],
          layers = layers
        )
        Key(object = value[[group]]) <- Key(object = group, quiet = TRUE)
      }
      value
    },
    layers = {
      groups <- apply(
        X = expand.grid(layers, names(x = splits)),
        MARGIN = 1L,
        FUN = paste,
        collapse = '.'
      )
      value <- vector(mode = 'list', length = length(x = groups))
      names(x = value) <- groups
      for (lyr in layers) {
        lcells <- Cells(x = x, layer = lyr)
        for (i in seq_along(along.with = splits)) {
          group <- paste(lyr, names(x = splits)[i], sep = '.')
          xcells <- intersect(x = splits[[i]], y = lcells)
          value[[group]] <- LayerData(object = x, layer = lyr, cells = xcells)
        }
      }
      value
    },
    abort(message = paste("Unknown split return type", sQuote(x = ret)))
  ))
}

#' Split an Assay
#'
#' @inheritParams [.Assay5
#' @inheritParams base::split
#' @param layers Names of layers to include in the split; pass \code{NA} for
#' all layers; pass \code{NULL} for the \link[=DefaultLayer]{default layer}
#' @param ret Type of return value; choose from:
#' \itemize{
#'  \item \dQuote{\code{assay}}: a single \code{\link{Assay5}} object
#'  \item \dQuote{\code{multiassay}}: a list of \code{\link{Assay5}} objects
#'  \item \dQuote{\code{layers}}: a list of layer matrices
#' }
#' @template param-dots-ignored
#'
#' @return Depends on the value of \code{ret}:
#' \itemize{
#'  \item \dQuote{\code{assay}}: \code{x} with the layers requested in
#'  \code{layers} split based on \code{f}; all other layers are left as-is
#'  \item \dQuote{\code{multiassay}}: a list of \code{\link{Assay5}} objects;
#'  the list contains one value per split and each assay contains only the
#'  layers requested in \code{layers} with the \link[=Key]{key} set to the split
#'  \item \dQuote{\code{layers}}: a list of matrices of length
#'  \code{length(assays) * length(unique(f))}; the list is named as
#'  \dQuote{\code{layer.split}}
#' }
#'
#' @method split Assay5
#' @export
#'
#' @family assay5
#'
#' @template section-progressr
#'
split.Assay5 <- split.StdAssay

#' @inherit subset.Assay5 params return title description details sections
#'
#' @keywords internal
#' @method subset StdAssay
#' @export
#'
#' @family stdassay
#'
subset.StdAssay <- function(
  x,
  cells = NULL,
  features = NULL,
  layers = NULL,
  ...
) {
  # define an inner function to validate the `cells` and `features` params
  .validate_param <- function(name, values, allowed) {
    # if `values` is null or contains only null values, keep all allowed values
    if (all(is.na(values))) {
      values <- allowed
    } else if (any(is.na(x = values))) {
      # if any values are NA, issue a warning and remove NAs
      warning(
        paste0("NAs passed in ", name, " vector, removing NAs"),
        call. = FALSE,
        immediate. = TRUE
      )
      # and drop null values from `values`
      values <- values[!is.na(x = values)]
    }
    # if `values` is numeric, treat them as indices
    if (is.numeric(values)) {
      values <- allowed[values]
    }
    # ensure `values` are in the allowed set
    values <- intersect(values, allowed)
    # if no valid values remain, stop execution with an error
    if (!length(values)) {
      stop(paste0("None of the ", name, " provided found in this assay"), call. = FALSE)
    }
    return(values)
  }

  # if no subsetting is specified, return the original object
  if (is.null(cells) && is.null(features) && is.null(layers)) {
    return(x)
  }

  # validate and filter cells
  all_cells <- Cells(x)
  cells <- .validate_param("cells", cells, all_cells)
  # validate and filter features
  all_features <- Features(x = x, layer = NA)
  features <- .validate_param("features", features, all_features)
  # validate and filter layers
  all_layers <- Layers(object = x)
  layers <- layers %||% all_layers
  layers <- match.arg(
    arg = layers,
    choices = all_layers,
    several.ok = TRUE
  )
  
  # subset cells and features layer by layer
  for (layer_name in all_layers) {
    # maybe drop the layer
    if (!layer_name %in% layers) {
      LayerData(x, layer = layer_name) <- NULL
      next
    }
    # otherwise, filter the the layer's cells and features
    # `MatchCells` is a bit of a misnomer - assuming that `new` is a 
    # subset of `old`, the function returns a list of indices mapping
    # the values of `new` to their order in `orig`
    layer_cells <- MatchCells(
      new = Cells(x = x, layer = layer_name),
      orig = cells,
      ordered = TRUE
    )
    layer_features <- MatchCells(
      new = Features(x = x, layer = layer_name),
      orig = features,
      ordered = TRUE
    )
    # if no valid cells or features, drop the layer data
    if (is.null(layer_cells) || is.null(layer_features)) {
      LayerData(object = x, layer = layer_name) <- NULL
      next
    } 
    # otherwise, apply the subset
    LayerData(object = x, layer = layer_name) <- LayerData(
      object = x,
      layer = layer_name,
      cells = layer_cells,
      features = layer_features
    )
  }

  # clean up the cells and features slots
  slot(x, name = "cells") <- droplevels(slot(x, name = "cells"))
  slot(x, name = "features") <- droplevels(slot(x, name = "features"))

  # in case any features were found in a only one layer and it was dropped
  # in the previous loop, we need to make sure our feature list is updated
  features <- intersect(features, Features(x = x, layer = NA))
  # update the features to match the valid list - see note above on `MatchCells`
  mfeatures <- MatchCells(
    new = all_features,
    orig = features,
    ordered = TRUE
  )
  # subset the meta.data slot accordingly
  slot(x, name = "meta.data") <- slot(x, name = "meta.data")[mfeatures, , drop = FALSE]

  # ensure the object is valid
  validObject(x)

  return(x)
}

#' Subset an Assay
#'
#' @inheritParams [[.Assay5
#' @param cells Cell names
#' @param features Feature names
#' @param layers Layer to keep; defaults to all layers
#'
#' @return \code{x} with just the cells and features specified by
#' \code{cells} and \code{features} for the layers specified by \code{layers}
#'
#' @method subset Assay5
#' @export
#'
#' @family assay5
#'
subset.Assay5 <- subset.StdAssay

#' @rdname sub-sub-.StdAssay
#'
#' @method tail StdAssay
#' @export
#'
tail.StdAssay <- tail.Assay

#' @return \code{tail}: the last \code{n} rows of feature-level meta data
#'
#' @rdname sub-sub-.Assay5
#'
#' @method tail Assay5
#' @export
#'
tail.Assay5 <- tail.StdAssay

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.VFLayers <- function(
  object,
  type = c('hvf', 'svf'),
  layers = NA,
  missing = FALSE
) {
  type <- type[1L]
  type <- match.arg(arg = type)
  pattern <- switch(
    EXPR = type,
    'hvf' = '^vf_',
    stop("Unknown type: '", type, "'", call. = FALSE)
  )
  vf.cols <- grep(
    pattern = paste0(pattern, '[[:alnum:]]+_'),
    x = colnames(x = object[[]]),
    value = TRUE
  )
  vf.layers <- unique(x = unlist(x = lapply(
    X = strsplit(x = vf.cols, split = '_'),
    FUN = function(x) {
      return(paste(x[3L:(length(x = x) - 1L)], collapse = '_'))
    }
  )))

  if (!isTRUE(x = missing)) {
    vf.layers <- intersect(
      x = vf.layers,
      y = Layers(object = object, search = layers)
    )
  }
  if (!length(x = vf.layers)) {
    vf.layers <- NULL
  }
  return(vf.layers)
}

#' @param object A \code{\link{StdAssay}} object
#' @param type Type of variable feature method to pull; choose from:
#' \itemize{
#'  \item \dQuote{\code{hvf}}: highly variable features
#'  \item \dQuote{\code{svf}}: spatially variable features
#' }
#' @param layers Vector of layers to restrict methods for, or a search pattern
#' for multiple layers
#'
#' @return A vector of variable feature methods found in \code{object}
#'
#' @noRd
#'
.VFMethods <- function(
  object,
  type = c('hvf', 'svf'),
  layers = NA,
  missing = FALSE
) {
  type <- type[1L]
  type <- match.arg(arg = type)
  pattern <- switch(
    EXPR = type,
    'hvf' = '^vf_',
    abort(message = paste("Unknown type:", sQuote(x = type)))
  )
  vf.cols <- grep(
    pattern = paste0(pattern, '[[:alnum:]]+_'),
    x = colnames(x = object[[]]),
    value = TRUE
  )
  # layers <- Layers(object = object, search = layers)
  layers <- .VFLayers(
    object = object,
    type = type,
    layers = layers,
    missing = missing
  )
  vf.cols <- Filter(
    f = function(x) {
      x <- unlist(x = strsplit(x = x, split = '_'))
      x <- paste(x[3:(length(x = x) - 1L)], collapse = '_')
      return(x %in% layers)
    },
    x = vf.cols
  )
  vf.methods <- unique(x = unlist(x = lapply(
    X = strsplit(x = vf.cols, split = '_'),
    FUN = '[[',
    2L
  )))
  if (!length(x = vf.methods)) {
    vf.methods <- NULL
  }
  return(vf.methods)
}

#' @param object A \code{\link{StdAssay}} object
#' @param type Type of variable feature method to pull; choose from:
#' \itemize{
#'  \item \dQuote{\code{hvf}}: highly variable features
#'  \item \dQuote{\code{svf}}: spatially variable features
#' }
#' @param layers Vector of layers to restrict methods for, or a search pattern
#' for multiple layers
#'
#' @return A vector of variable feature methods and corresponding layers found in \code{object}
#'
#' @importFrom stats setNames
#' @importFrom utils modifyList
#'
#' @noRd
#'
.VFMethodsLayers <- function(
  object,
  type = c('hvf', 'svf'),
  layers = NA,
  missing = FALSE
) {
  type <- type[1L]
  type <- match.arg(arg = type)
  pattern <- switch(
    EXPR = type,
    'hvf' = '^vf_',
    abort(message = paste("Unknown type:", sQuote(x = type)))
  )
  vf.cols <- grep(
    pattern = paste0(pattern, '[[:alnum:]]+_'),
    x = colnames(x = object[[]]),
    value = TRUE
  )
  # layers <- Layers(object = object, search = layers)
  layers <- .VFLayers(
    object = object,
    type = type,
    layers = layers,
    missing = missing
  )
  vf.cols <- Filter(
    f = function(x) {
      x <- unlist(x = strsplit(x = x, split = '_'))
      x <- paste(x[3:(length(x = x) - 1L)], collapse = '_')
      return(x %in% layers)
    },
    x = vf.cols
  )
  # Extract methods and layers
  vf.methods.layers <- lapply(vf.cols, function(col) {
    components <- strsplit(col, split = "_")[[1]]
    method <- components[2]
    layer <- paste(components[3:(length(components) - 1)], collapse = "_")
    return(c(method = method, layer = layer))
  })

  # Combine into a list
  vf.list <- lapply(unique(unlist(lapply(vf.methods.layers, `[[`, "method"))), function(method) {
    layers <- unique(unlist(lapply(vf.methods.layers, function(x) {
      if (x['method'] == method)
        return(x['layer'])
    })))
    return(setNames(list(layers), method))
  })
  vf.list <- Reduce(modifyList, vf.list)
  if (!length(x = vf.list)) {
    vf.list <- NULL
  }
  return(vf.list)
}

CalcN5 <- function(object) {
  if (IsMatrixEmpty(x = LayerData(object = object))) {
    return(NULL)
  }
  return(list(
    nCount = colSums(x = object),
    nFeature = colSums(x = LayerData(object = object) > 0)
  ))
}

# Join single layers
#
JoinSingleLayers <- function(
  object,
  layers = NULL,
  new = NULL,
  default = TRUE,
  nfeatures = Inf,
  ...
) {
  if (is.null(x = layers)) {
    stop('Layers cannot be NULL')
  }
  if (length(x = layers) > 1L) {
    stop('The length of input layers should be 1')
  }
  layers <- Layers(object = object, search = layers)
  new <- new %||% 'newlayer'
  if (length(x = layers) == 1L) {
    LayerData(object = object, layer = new) <- LayerData(object = object, layer = layers)
    return(object)
  }
  if (length(x = layers) == 0L) {
    return(object)
  }
  # Stitch the layers together
  ldat <- StitchMatrix(
    x = LayerData(object = object, layer = layers[1L]),
    y = lapply(X = layers[2:length(x = layers)], FUN = LayerData, object = object),
    rowmap = slot(object = object, name = 'features')[, layers],
    colmap = slot(object = object, name = 'cells')[, layers]
  )
  LayerData(object = object, layer = new) <- ldat
  # Set the new layer as default
  if (isTRUE(x = default)) {
    DefaultLayer(object = object) <- new
  }
  # Remove the old layers
  for (lyr in layers) {
    object[lyr] <- NULL
  }
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setAs(
  from = 'Assay',
  to = 'Assay5',
  def = function(from) {
    # Initialize the new object
    to <- new(
      Class = 'Assay5',
      cells = LogMap(y = colnames(x = from)),
      features = LogMap(y = rownames(x = from)),
      assay.orig = DefaultAssay(object = from) %||% character(length = 0L),
      meta.data = EmptyDF(n = nrow(x = from)),
      key = Key(object = from)
    )
    # browser()
    # Add the expression matrices
    for (i in c('counts', 'data', 'scale.data')) {
      adata <- GetAssayData(object = from, layer = i)
      if (IsMatrixEmpty(x = adata)) {
        next
      }
      LayerData(object = to, layer = i) <- adata
    }
    # Set the default layer
    DefaultLayer(object = to) <- ifelse(
      test = 'counts' %in% Layers(object = to) && !'scale.data' %in% Layers(object = to),
      yes = 'counts',
      no = 'data'
    )
    # Add feature-level meta data
    to[[]] <- from[[]]
    # Set Variable features
    VariableFeatures(object = to) <- VariableFeatures(object = from)
    # Add miscellaneous data
    mdata <- Misc(object = from)
    for (i in names(x = mdata)) {
      Misc(object = to, slot = i) <- mdata[[i]]
    }

    return(to)
  }
)

setAs(
  from = 'Assay5',
  to = 'Assay',
  def = function(from) {
    data.list <- c()
    original.layers <- Layers(object = from)
    layers.saved <- c()
    for (i in c('counts', 'data', 'scale.data')) {
      layers.saved <- c(layers.saved, Layers(object = from, search = i))
      if (length(Layers(object = from, search = i)) > 1) {
        warning("Joining '", i, "' layers. If you have the same cells in multiple layers, ",
                "the expression value for the cell in the '",
                i, "' slot will be the value from the '",
                Layers(object = from, search = i)[1], "' layer.",
                call. = FALSE,
                immediate. = TRUE)
        from <- JoinLayers(object = from,
                           layers = i,
                           new = i)
      }
      if (i == "data") {
        if (isTRUE(Layers(object = from, search = i) == "scale.data")) {
          warning("No counts or data slot in object. Setting 'data' slot using",
                  " data from 'scale.data' slot. To recreate 'data' slot, you",
                  " must set and normalize data from a 'counts' slot.",
                  call. = FALSE)
        }
      }
      adata <- LayerData(object = from, layer = i)
      if(inherits(x = adata, what = "IterableMatrix")) {
        warning("Converting IterableMatrix to sparse dgCMatrix",
                call. = FALSE)
        adata <- as(object = adata, Class = "dgCMatrix")
      }
      data.list[[i]] <- adata
    }
    if (IsMatrixEmpty(x = data.list[["data"]])){
      data.list[["data"]] <- data.list[["counts"]]
    }
    if (any(!(original.layers %in% layers.saved))){
      layers.remove <- original.layers[!(original.layers %in% layers.saved)]
      warning("Layers ", paste0(layers.remove, collapse = ', '),
              " will be removed from the object as v3 assays only support",
              " 'counts', 'data', or 'scale.data' slots.",
              call. = FALSE,
              immediate. = TRUE)
    }
    to <- new(
      Class = 'Assay',
      counts = data.list[["counts"]],
      data = data.list[["data"]],
      scale.data = data.list[["scale.data"]],
      assay.orig = DefaultAssay(object = from) %||% character(length = 0L),
      meta.features = data.frame(row.names = rownames(x = data.list[["data"]])),
      key = Key(object = from)
    )
    # Add feature-level meta data
    suppressWarnings(to[[]] <- from[[]])
    # set variable features
    VariableFeatures(object = to) <- VariableFeatures(object = from)
    mdata <- Misc(object = from)
    for (i in names(x = mdata)) {
      Misc(object = to, slot = i) <- mdata[[i]]
    }
    return(to)
  }
)

#' @rdname sub-.StdAssay
#'
setMethod(
  f = '[<-',
  signature = c(x = 'StdAssay', i = 'character'),
  definition = function(x, i, ..., value) {
    LayerData(object = x, layer = i, ...) <- value
    return(x)
  }
)

#' @param value A matrix-like object to add as a new layer
#'
#' @return \code{[<-}: \code{x} with layer data \code{value} saved as \code{i}
#'
#' @rdname sub-.Assay5
#'
setMethod(
  f = '[<-',
  signature = c(x = 'Assay5', i = 'character'),
  definition = function(x, i, ..., value) {
    return(callNextMethod(x = x, i = i, ..., value = value))
  }
)

#' @rdname sub-sub-.StdAssay
#'
setMethod(
  f = '[[<-',
  signature = c(
    x = 'StdAssay',
    i = 'character',
    j = 'missing',
    value = 'data.frame'
  ),
  definition = function(x, i, ..., value) {
    if (!length(x = i) && !ncol(x = value)) {
      return(x)
    }
    i <- match.arg(arg = i, choices = colnames(x = value), several.ok = TRUE)
    names.intersect <- intersect(
      x = row.names(x = value),
      y = Features(x = x, layer = NA)
    )
    if (length(x = names.intersect)) {
      value <- value[names.intersect, , drop = FALSE]
    } else if (nrow(x = value) == nrow(x = x)) {
      row.names(x = value) <- Features(x = x, layer = NA)
    } else {
      abort(message = "Cannot add more or less meta data without feature names")
    }
    for (n in i) {
      v <- value[[n]]
      names(x = v) <- row.names(value)
      x[[n]] <- v
    }
    return(x)
  }
)

#' @rdname sub-sub-.StdAssay
#'
setMethod(
  f = '[[<-',
  signature = c(
    x = 'StdAssay',
    i = 'missing',
    j = 'missing',
    value = 'data.frame'
  ),
  definition = function(x, ..., value) {
    # Allow removing all meta data
    if (IsMatrixEmpty(x = value)) {
      x[[names(x = x[[]])]] <- NULL
      return(x)
    }
    if (is.null(names(x = value))) {
      warn(message = 'colnames of input cannot be NULL')
    } else {
      # If no `i` provided, use the column names from value
      x[[names(x = value)]] <- value
    }
    return(x)
  }
)

#' @importFrom methods selectMethod
#'
#' @rdname sub-sub-.StdAssay
#'
setMethod(
  f = '[[<-',
  signature = c(x = 'StdAssay', i = 'character', j = 'missing', value = 'factor'),
  definition = function(x, i, ..., value) {
    f <- slot(
      object = selectMethod(
        f = '[[<-',
        signature = c(
          x = 'StdAssay',
          i = 'character',
          j = 'missing',
          value = 'vector'
        )
      ),
      name = '.Data'
    )
    return(f(x = x, i = i, value = value))
  }
)

#' @rdname sub-sub-.StdAssay
#'
setMethod(
  f = '[[<-',
  signature = c(x = 'StdAssay', i = 'character', j = 'missing', value = 'NULL'),
  definition = function(x, i, ..., value) {
    for (name in i) {
      slot(object = x, name = 'meta.data')[[name]] <- NULL
    }
    return(x)
  }
)

#' @rdname sub-sub-.StdAssay
#'
setMethod(
  f = '[[<-',
  signature = c(x = 'StdAssay', i = 'character', j = 'missing', value = 'vector'),
  definition = function(x, i, ..., value) {
    # Add multiple bits of metadata
    if (length(x = i) > 1L) {
      value <- rep_len(x = value, length.out = length(x = i))
      for (idx in seq_along(along.with = i)) {
        x[[i[idx]]] <- value[[idx]]
      }
    } else {
      # Add a single column of metadata
      if (is.null(x = names(x = value))) {
        if (length(x = unique(x = value)) == 1) {
          value <- rep_len(x = value, length.out = nrow(x = x))
          names(x = value) <- Features(x = x, layer = NA)
        } else {
          names(x = value) <- value
        }
      }
      names.intersect <- intersect(
        x = names(x = value),
        y = Features(x = x, layer = NA)
      )
      if (!length(x = names.intersect)) {
        abort(message = "No feature overlap between new meta data and assay")
      }
      value <- value[names.intersect]
      df <- EmptyDF(n = nrow(x = x))
      rownames(x = df) <- Features(x = x, layer = NA)
      # df[[i]] <- if (i %in% names(x = x[[]])) {
      #   x[[i]]
      # } else {
      #   NA
      # }
      df[names(x = value), i] <- value
      if (nrow(x = slot(object = x, name = 'meta.data')) == 0) {
        slot(object = x, name = 'meta.data') <- EmptyDF(n = nrow(x = x))
      }
      slot(object = x, name = 'meta.data')[, i] <- df[[i]]
    }
    validObject(object = x)
    return(x)
  }
)

#' @rdname sub-sub-.StdAssay
#'
setMethod(
  f = '[[<-',
  signature = c(x = 'StdAssay', i = 'numeric', j = 'missing', value = 'ANY'),
  definition = function(x, i, ..., value) {
    if (ncol(x = x[[]])) {
      i <- colnames(x = x[[]])[as.integer(x = i)]
      i <- i[!is.na(x = i)]
      if (length(x = i)) {
        x[[i]] <- value
      }
    }
    return(x)
  }
)

#' @rdname sub-sub-.StdAssay
#'
setMethod(
  f = '[[<-',
  signature = c(x = 'StdAssay', i = 'missing', j = 'missing', value = 'NULL'),
  definition = function(x, ..., value) {
    slot(object = x, name = 'meta.data') <- EmptyDF(n = nrow(x = x))
    return(x)
  }
)

#' @param value Feature-level meta data to add
#'
#' @return \code{[[<-}: \code{x} with \code{value} added as \code{i}
#' in feature-level meta data
#'
#' @rdname sub-sub-.Assay5
#'
#' @order 2
#'
setMethod(
  f = '[[<-',
  signature = c(x = 'Assay5'),
  definition = function(x, i, ..., value) {
    return(callNextMethod(x = x, i = i, ..., value = value))
  }
)

#' V5 Assay Summaries
#'
#' Summary maths for \code{\link{StdAssay}} Objects
#'
#' @inheritParams base::colSums
#' @param layer Name of layer to run function on
#' @template param-dots-ignored
#'
#' @return The results of the summary math function for the layer specified
#'
#' @name v5-assay-summaries
#' @rdname v5-assay-summaries
#'
#' @keywords internal
#'
NULL

#' @rdname v5-assay-summaries
#'
setMethod(
  f = 'colMeans',
  signature = c(x = 'StdAssay'),
  definition = function(x, na.rm = FALSE, dims = 1, layer = NULL, ...) {
    return(Matrix::colMeans(
      x = LayerData(object = x, layer = layer),
      na.rm = na.rm,
      dims = dims
    ))
  }
)

#' @rdname v5-assay-summaries
#'
setMethod(
  f = 'colSums',
  signature = c(x = 'StdAssay'),
  definition = function(x, na.rm = FALSE, dims = 1, layer = NULL, ...) {
    return(Matrix::colSums(
      x = LayerData(object = x, layer = layer),
      na.rm = na.rm,
      dims = dims
    ))
  }
)

#' @rdname v5-assay-summaries
#'
setMethod(
  f = 'rowMeans',
  signature = c(x = 'StdAssay'),
  definition = function(x, na.rm = FALSE, dims = 1, layer = NULL, ...) {
    return(Matrix::rowMeans(
      x = LayerData(object = x, layer = layer),
      na.rm = na.rm,
      dims = dims
    ))
  }
)

#' @rdname v5-assay-summaries
#'
setMethod(
  f = 'rowSums',
  signature = c(x = 'StdAssay'),
  definition = function(x, na.rm = FALSE, dims = 1, layer = NULL, ...) {
    return(Matrix::rowSums(
      x = LayerData(object = x, layer = layer),
      na.rm = na.rm,
      dims = dims
    ))
  }
)

#' V5 Assay Overview
#'
#' Overview of a \code{\link{StdAssay}} object
#'
#' @param object A v5 Assay
#'
#' @template return-show
#'
#' @keywords internal
#'
#' @family stdassay
#'
setMethod(
  f = 'show',
  signature = 'StdAssay',
  definition = function(object) {
    # Basic assay info
    cat(
      .AssayClass(object = object),
      'data with',
      nrow(x = object),
      'features for',
      ncol(x = object),
      'cells\n'
    )
    # Feature information
    if (length(x = VariableFeatures(object = object))) {
      top.ten <- head(x = VariableFeatures(object = object), n = 10L)
      top <- 'Top'
      variable <- 'variable'
    } else {
      top.ten <- head(x = Features(x = object), n = 10L)
      top <- 'First'
      variable <- ''
    }
    features <- paste(
      variable,
      paste0(
        ifelse(test = length(x = top.ten) == 1L, yes = 'feature', no = 'features'),
        ":\n"
      )
    )
    features <- gsub(pattern = '^\\s+', replacement = '', x = features)
    cat(
      top,
      length(x = top.ten),
      features,
      paste(strwrap(x = paste(top.ten, collapse = ', ')), collapse = '\n'),
      '\n'
    )
    cat(
      "Layers:\n",
      paste(strwrap(x = paste(Layers(object = object), collapse = ', ')), collapse = '\n'),
      "\n"
    )
    return(invisible(x = NULL))
  }
)

#' @rdname split.StdAssay
#'
setMethod( # Because R is stupid
  f = 'split',
  signature = c(x = 'StdAssay'),
  definition = split.StdAssay
)

#' V5 Assay Validity
#'
#' @templateVar cls StdAssay
#' @template desc-validity
#'
#' @section Layer Validation:
#' blah
#'
#' @inheritSection Key-validity Key Validation
#'
#' @keywords internal
#'
#' @name StdAssay-validity
#'
#' @family stdassay
#'
#' @seealso \code{\link[methods]{validObject}}
#'
setValidity(
  Class = 'StdAssay',
  method = function(object) {
    if (isFALSE(x = getOption(x = "Seurat.object.validate", default = TRUE))) {
      warn(
        message = paste("Not validating", class(x = object)[1L], "objects"),
        class = 'validationWarning'
      )
      return(TRUE)
    }
    valid <- NULL
    # Check layers
    dorder <- c(
      features = .MARGIN(x = object, type = 'features'),
      cells = .MARGIN(x = object, type = 'cells')
    )
    adims <- dim(x = object) # c(features, cells)
    if (!IsNamedList(x = slot(object = object, name = 'layers'), pass.zero = TRUE)) {
      valid <- c(valid, "'layers' must be a named list")
    }
    for (layer in Layers(object = object)) {
      # Reorder dimensions of layer to c(features, cells)
      ldims <- dim(x = slot(object = object, name = 'layers')[[layer]])[dorder]
      if (length(x = ldims) != 2L) {
        valid <- c(valid, "Layers must be two-dimensional objects")
        break
      }
      # Check that we have the correct features and cells
      for (i in seq.int(from = 1L, to = 2L)) {
        if (ldims[i] > adims[i]) {
          valid <- c(
            valid,
            paste0(
              "Layers may not have more ",
              names(x = dorder)[i],
              " than present in the assay (offending layer",
              layer,
              ")"
            )
          )
        }
      }
      # Check that we've recorded the cells and features in the maps
      for (i in c('cells', 'features')) {
        didx <- c(features = 1L, cells = 2L)[i]
        if (!layer %in% colnames(x = slot(object = object, name = i))) {
          valid <- c(
            valid,
            paste0(
              "All layers must have a record in the ",
              i,
              " map (offending layer: ",
              layer,
              ")"
            )
          )
        } else {
          nmap <- length(x = slot(object = object, name = i)[[layer]])
          if (nmap != ldims[didx]) {
            valid <- c(
              valid,
              paste0(
                "Layers must have the same ",
                i,
                " as present in the map (offending layer: ",
                layer,
                ")"
              )
            )
          }
        }
      }
    }
    didx <- slot(object = object, name = 'default')
    if (length(x = didx)) {
      if (didx < 0 || didx > length(x = Layers(object = object))) {
        valid <- c(
          valid,
          "'default' must be between 0 and the number of layers present"
        )
      }
    }
    # TODO: Check variable features
    # TODO: Check meta features
    # TODO: Check key
    # TODO: Check misc
    return(valid %||% TRUE)
  }
)

#' @inherit StdAssay-validity title details sections
#'
#' @templateVar cls Assay5
#' @template desc-validity
#'
#' @name Assay5-validity
#'
#' @family assay5
#'
#' @seealso \code{\link[methods]{validObject}}
#'
NULL
