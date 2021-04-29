#' @include zzz.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Core Assay Infrastructure
#'
#' The \code{StdAssay} class is a virtual class that provides core
#' infrastructure for assay data in \pkg{Seurat}. Assays contain expression
#' data along with additional representations of the expression data (layers)
#' and associated metadata. Derived classes
#' (eg. \link[Assay-class]{the v5 Assay}) must set the storage mechanism for
#' expression data (eg. \code{\link[base]{matrix}} or
#' \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}) and may optionally define
#' additional functionality
#'
#' @slot data A two-dimensional object with expression data with cells as
#' columns and features as rows
#' @slot layers A named list containing alternate representations of
#' \code{data}; layers must have the same cells as \code{data} and either the
#' same or a subset of the features present in \code{data}
#' @slot key ...
#' @slot cells A vector of cell names present in the assay
#' @slot assay.orig ...
#' @slot features A matrix containing the presence/absence status of features
#' in the layers with the following properties;
#' \itemize{
#'  \item The dimensions of the matrix should be m\emph{f} rows and n\emph{l}
#'  columns where m\emph{f} is the number of features (rows) in \code{data}
#'  and n\emph{l} is the number of layers in \code{layers}
#'  \item The row names should be the features (in order) present in \code{data}
#'  \item The column names should be the names of each layer present in
#'  \code{layers}
#'  \item Each value should be either \code{TRUE} if a given feature is present
#'  in a given layer, otherwise \code{FALSE}
#' }
#' @slot meta.data ...
#' @slot misc ...
#'
#' @exportClass StdAssay
#'
#' @aliases StdAssay
#'
#' @seealso \code{\link{Assay5-class}}
#'
setClass(
  Class = 'StdAssay',
  contains = 'VIRTUAL',
  slots = c(
    data = 'ANY',
    layers = 'list',
    key = 'character',
    cells = 'character',
    assay.orig = 'character',
    features = 'matrix',
    # var.features = 'character',
    meta.data = 'data.frame',
    misc = 'list'
  )
)

#' The v5 \code{Assay} Object
#'
#' The v5 \code{Assay} is the typical \code{Assay} class used in \pkg{Seurat}
#' v5; ...
#'
#' @slot data A \code{\link[base]{matrix}} or
#' \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}} with single-cell expression
#' data
#'
#' @exportClass Assay5
#'
#' @aliases Assay5
#'
#' @seealso \code{\link{StdAssay-class}}
#'
setClass(
  Class = 'Assay5',
  contains = 'StdAssay',
  slots = c(
    data = 'AnyMatrix'
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CreateAssay5Object <- function(
  counts,
  min.cells = 0,
  min.features = 0,
  ...
) {
  # Filter based on min.features
  if (min.features > 0) {
    nfeatures <- Matrix::colSums(x = counts > 0)
    counts <- counts[, which(x = nfeatures >= min.features)]
  }
  # Filter features based on the number of cells expressing
  if (min.cells > 0) {
    num.cells <- Matrix::rowSums(x = counts > 0)
    counts <- counts[which(x = num.cells >= min.cells),]
  }
  features <- rownames(x = counts)
  cells <- colnames(x = counts)
  dimnames(x = counts) <- list(NULL, NULL)
  return(new(
    Class = 'Assay5',
    data = counts,
    layers = list(),
    features = matrix(
      nrow = nrow(x = counts),
      ncol = 0,
      dimnames = list(features, NULL)
    ),
    cells = cells,
    meta.data = as.data.frame(x = matrix(nrow = nrow(x = counts), ncol = 0))
  ))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname AddMetaData
#' @export
#' @method AddMetaData StdAssay
#'
AddMetaData.StdAssay <- .AddMetaData

#' @rdname Cells
#' @method Cells StdAssay
#' @export
#'
Cells.StdAssay <- function(x, ...) {
  return(slot(object = x, name = 'cells'))
}

#' @rdname DefaultAssay
#' @export
#' @method DefaultAssay StdAssay
#'
DefaultAssay.StdAssay <- function(object, ...) {
  # object <- UpdateSlots(object = object)
  return(slot(object = object, name = 'assay.orig'))
}

#' @rdname DefaultAssay
#' @export
#' @method DefaultAssay<- StdAssay
#'
"DefaultAssay<-.StdAssay" <- function(object, ..., value) {
  # object <- UpdateSlots(object = object)
  slot(object = object, name = 'assay.orig') <- value
  return(object)
}

#' @param layer
#'
#' @rdname Cells
#' @export
#' @method Features StdAssay
#'
Features.StdAssay <- function(x, layer = 'data', ...) {
  layer <- layer[1]
  layer <- match.arg(arg = layer, choices = Layers(object = x))
  fmat <- slot(object = x, name = 'features')
  return(switch(
    EXPR = layer,
    'data' = rownames(x = fmat),
    rownames(x = fmat)[fmat[, layer]]
  ))
}

#' @rdname AssayData
#' @export
#' @method GetAssayData StdAssay
#'
#' @examples
#' \dontrun{
#' # Get the data directly from an Assay object
#' GetAssayData(pbmc_small[["RNA"]], slot = "data")[1:5,1:5]
#' }
#'
GetAssayData.StdAssay <- function(
  object,
  slot = 'data',
  ...
) {
  CheckDots(...)
  slot <- slot[1]
  slot <- match.arg(arg = slot, choices = Layers(object = object))
  return(slot(object = object, name = slot))
}

#' @rdname VariableFeatures
#' @export
#' @method HVFInfo StdAssay
#'
#' @examples
#' \dontrun{
#' # Get the HVF info directly from an Assay object
#' HVFInfo(pbmc_small[["RNA"]], selection.method = 'vst')[1:5, ]
#' }
#'
HVFInfo.StdAssay <- function(object, selection.method, status = FALSE, ...) {
  .NotYetImplemented()
  CheckDots(...)
  disp.methods <- c('mean.var.plot', 'dispersion', 'disp')
  if (tolower(x = selection.method) %in% disp.methods) {
    selection.method <- 'mvp'
  }
  selection.method <- switch(
    EXPR = tolower(x = selection.method),
    'sctransform' = 'sct',
    selection.method
  )
  vars <- switch(
    EXPR = selection.method,
    'vst' = c('mean', 'variance', 'variance.standardized'),
    'mvp' = c('mean', 'dispersion', 'dispersion.scaled'),
    'sct' = c('gmean', 'variance', 'residual_variance'),
    stop("Unknown method: '", selection.method, "'", call. = FALSE)
  )
  tryCatch(
    expr = hvf.info <- object[[paste(selection.method, vars, sep = '.')]],
    error = function(e) {
      stop(
        "Unable to find highly variable feature information for method '",
        selection.method,
        "'",
        call. = FALSE
      )
    }
  )
  colnames(x = hvf.info) <- vars
  if (status) {
    hvf.info$variable <- object[[paste0(selection.method, '.variable')]]
  }
  return(hvf.info)
}

#' @rdname Key
#' @export
#' @method Key StdAssay
#'
#' @examples
#' \dontrun{
#' # Get an Assay key
#' Key(pbmc_small[["RNA"]])
#' }
#'
Key.StdAssay <- function(object, ...) {
  CheckDots(...)
  return(slot(object = object, name = 'key'))
}

#' @rdname Key
#' @export
#' @method Key<- StdAssay
#'
#' @examples
#' \dontrun{
#' # Set the key for an Assay
#' Key(pbmc_small[["RNA"]]) <- "newkey_"
#' Key(pbmc_small[["RNA"]])
#' }
#'
"Key<-.StdAssay" <- function(object, ..., value) {
  CheckDots(...)
  slot(object = object, name = 'key') <- value
  return(object)
}

#' @method LayerData StdAssay
#' @export
#'
LayerData.StdAssay <- function(object, layer = 'data', ...) {
  layer <- layer[1]
  layer <- match.arg(arg = layer, choices = Layers(object = object))
  return(switch(
    EXPR = layer,
    'data' = slot(object = object, name = layer),
    slot(object = object, name = 'layers')[[layer]]
  ))
}

#' @method LayerData Assay5
#' @export
#'
LayerData.Assay5 <- function(object, layer = 'data', dnames = TRUE, ...) {
  ldat <- NextMethod()
  if (isTRUE(x = dnames)) {
    dimnames(x = ldat) <- list(
      Features(x = object, layer = layer),
      Cells(x = object)
    )
  }
  return(ldat)
}

#' @method LayerData<- StdAssay
#' @export
#'
"LayerData<-.StdAssay" <- function(object, layer, features = NULL, ..., value) {
  features <- features %||% rownames(x = value)
  if (layer == 'data') {
    # Modifying data
    if (!identical(x = dim(x = value), y = dim(x = object))) {
      stop("'data' cannot change dimensions")
    }
    if (!is.null(x = features)) {
      fmatch <- if (is.numeric(x = features)) {
        features
      } else {
        match(x = features, table = Features(x = object))
      }
      if (any(is.na(x = fmatch))) {
        stop("Mismatched features between value and current")
      }
      value <- value[fmatch, ]
    }
    # TODO: check cell order
    try(
      expr = dimnames(x = value) <- list(NULL, NULL),
      silent = TRUE
    )
    slot(object = object, name = 'data') <- value
  } else if (is.null(x = value)) {
    # Removing a layer
    slot(object = object, name = 'layers')[[layer]] <- NULL
    fmat <- slot(object = object, name = 'features')
    fidx <- which(x = colnames(x = fmat) == layer)
    if (length(x = fidx)) {
      fmat <- fmat[, -fidx]
    }
    slot(object = object, name = 'features') <- fmat
  } else {
    # Adding a layer
    if (!identical(x = ncol(x = value), y = ncol(x = object))) {
      stop("Layers must have the same cells as currently present")
    }
    if (is.null(x = features)) {
      if (nrow(x = value) != nrow(x = object)) {
        stop("If features are not provided, then layers must have the same features as 'data'")
      }
      fmatch <- seq_len(length.out = nrow(x = value))
    } else {
      fmatch <- if (is.numeric(x = features)) {
        features
      } else {
        as.vector(x = na.omit(object = match(
          x = features,
          table = Features(x = object)
        )))
      }
      value <- value[order(fmatch), ]
    }
    # TODO: check cell order
    try(
      expr = dimnames(x = value) <- list(NULL, NULL),
      silent = TRUE
    )
    slot(object = object, name = 'layers')[[layer]] <- value
    fmat <- cbind(
      slot(object = object, name = 'features'),
      matrix(data = FALSE, nrow = nrow(x = object), dimnames = list(NULL, layer))
    )
    fmat[fmatch, layer] <- TRUE
    slot(object = object, name = 'features') <- fmat
  }
  return(object)
}

#' @method LayerData<- Assay5
#' @export
#'
"LayerData<-.Assay5" <- function(object, layer, ..., value) {
  if (!inherits(x = value, what = c('AnyMatrix', 'NULL'))) {
    value <- as.sparse(x = value)
  }
  object <- NextMethod()
  return(object)
}

#' @param data Include \dQuote{data} as a layer
#'
#' @rdname Layers
#' @method Layers StdAssay
#' @export
#'
Layers.StdAssay <- function(object, data = TRUE, ...) {
  layers <- names(x = slot(object = object, name = 'layers'))
  if (isTRUE(x = data)) {
    layers <- c('data', layers)
  }
  return(layers)
}

#' @param slot Name of specific bit of meta data to pull
#'
#' @rdname Misc
#' @export
#' @method Misc StdAssay
#'
Misc.StdAssay <- .Misc

#' @rdname Misc
#' @export
#' @method Misc<- StdAssay
#'
"Misc<-.StdAssay" <- `.Misc<-`

#' @importFrom stats na.omit
#'
#' @rdname AssayData
#' @export
#' @method SetAssayData StdAssay
#'
#' @examples
#' \dontrun{
#' # Set an Assay slot directly
#' count.data <- GetAssayData(pbmc_small[["RNA"]], slot = "counts")
#' count.data <- as.matrix(x = count.data + 1)
#' new.assay <- SetAssayData(pbmc_small[["RNA"]], slot = "counts", new.data = count.data)
#' }
#'
SetAssayData.StdAssay <- function(
  object,
  slot,
  new.data,
  ...
) {
  .NotYetImplemented()
  LayerData(object = object, layer = slot) <- new.data
  return(object)
  CheckDots(...)
  slot <- slot[1]
  # slot <- match.arg(arg = slot, choices = '')
  if (!IsMatrixEmpty(x = new.data)) {
    if (any(grepl(pattern = '_', x = rownames(x = new.data)))) {
      warning(
        "Feature names cannot have underscores ('_'), replacing with dashes ('-')",
        call. = FALSE,
        immediate. = TRUE
      )
      rownames(x = new.data) <- gsub(
        pattern = '_',
        replacement = '-',
        x = rownames(x = new.data)
      )
    }
    if (ncol(x = new.data) != ncol(x = object)) {
      stop(
        "The new data doesn't have the same number of cells as the current data",
        call. = FALSE
      )
    }
    num.counts <- nrow(x = object)
    counts.names <- rownames(x = object)
    if (slot == 'scale.data' && nrow(x = new.data) > num.counts) {
      warning(
        "Adding more features than present in current data",
        call. = FALSE,
        immediate. = TRUE
      )
    } else if (slot %in% c('counts', 'data') && nrow(x = new.data) != num.counts) {
      warning(
        "The new data doesn't have the same number of features as the current data",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    if (!all(rownames(x = new.data) %in% counts.names)) {
      warning(
        "Adding features not currently present in the object",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    new.features <- na.omit(object = match(
      x = counts.names,
      table = rownames(x = new.data)
    ))
    new.cells <- colnames(x = new.data)
    if (!all(new.cells %in% colnames(x = object))) {
      stop(
        "All cell names must match current cell names",
        call. = FALSE
      )
    }
    new.data <- new.data[new.features, colnames(x = object), drop = FALSE]
    if (slot %in% c('counts', 'data') && !all(dim(x = new.data) == dim(x = object))) {
      stop(
        "Attempting to add a different number of cells and/or features",
        call. = FALSE
      )
    }
  }
  if (!is.vector(x = rownames(x = new.data))) {
    rownames(x = new.data) <- as.vector(x = rownames(x = new.data))
  }
  if (!is.vector(x = colnames(x = new.data))) {
    colnames(x = new.data) <- as.vector(x = colnames(x = new.data))
  }
  slot(object = object, name = slot) <- new.data
  return(object)
}

#' @rdname VariableFeatures
#' @export
#' @method VariableFeatures StdAssay
#'
VariableFeatures.StdAssay <- function(object, selection.method = NULL, ...) {
  .NotYetImplemented()
  CheckDots(...)
  if (!is.null(x = selection.method)) {
    vf <- HVFInfo(
      object = object,
      selection.method = selection.method,
      status = TRUE
    )
    return(rownames(x = vf)[which(x = vf[, "variable"][, 1])])
  }
  return(slot(object = object, name = 'var.features'))
}

#' @rdname VariableFeatures
#' @export
#' @method VariableFeatures<- StdAssay
#'
"VariableFeatures<-.StdAssay" <- function(object, ..., value) {
  .NotYetImplemented()
  CheckDots(...)
  if (length(x = value) == 0) {
    slot(object = object, name = 'var.features') <- character(length = 0)
    return(object)
  }
  if (any(grepl(pattern = '_', x = value))) {
    warning(
      "Feature names cannot have underscores '_', replacing with dashes '-'",
      call. = FALSE,
      immediate = TRUE
    )
    value <- gsub(pattern = '_', replacement = '-', x = value)
  }
  value <- split(x = value, f = value %in% rownames(x = object))
  if (length(x = value[['FALSE']]) > 0) {
    if (length(x = value[['TRUE']]) == 0) {
      stop(
        "None of the features provided are in this Assay object",
        call. = FALSE
      )
    } else {
      warning(
        "Not all features provided are in this Assay object, removing the following feature(s): ",
        paste(value[['FALSE']], collapse = ', '),
        call. = FALSE,
        immediate. = TRUE
      )
    }
  }
  slot(object = object, name = 'var.features') <- value[['TRUE']]
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' \code{StdAssay} Methods
#'
#' Methods for \code{\link{StdAssay}} objects for generics defined in
#' other packages
#'
#' @param x,object An \code{\link{StdAssay}} object
#' @param i,features For \code{[[}: metadata names; for all other methods,
#' feature names or indices
#' @param j,cells Cell names or indices
#' @param ... Arguments passed to other methods
#'
#' @name StdAssay-methods
#' @rdname StdAssay-methods
#'
#' @concept assay
#'
NULL

#' @describeIn StdAssay-methods Get expression data from an \code{StdAssay}
#'
#' @return \code{[}: The \code{data} slot for features \code{i} and cells
#' \code{j}
#'
#' @export
#' @method [ StdAssay
#'
"[.StdAssay" <- function(x, i, j, ...) {
  if (missing(x = i)) {
    i <- seq_len(length.out = nrow(x = x))
  }
  if (missing(x = j)) {
    j <- seq_len(length.out = ncol(x = x))
  }
  return(GetAssayData(object = x)[i, j, ..., drop = FALSE])
}

#' @describeIn StdAssay-methods Get feature-level metadata
#'
#' @param drop See \code{\link[base]{drop}}
#'
#' @return \code{[[}: The feature-level metadata for \code{i}
#'
#' @export
#' @method [[ StdAssay
#'
"[[.StdAssay" <- function(x, i, ..., drop = FALSE) {
  if (missing(x = i)) {
    i <- colnames(x = slot(object = x, name = 'meta.features'))
  }
  data.return <- slot(object = x, name = 'meta.features')[, i, drop = FALSE, ...]
  if (drop) {
    data.return <- unlist(x = data.return, use.names = FALSE)
    names(x = data.return) <- rep.int(x = rownames(x = x), times = length(x = i))
  }
  return(data.return)
}

#' @describeIn StdAssay-methods Number of cells and features for an \code{StdAssay}
#'
#' @return \code{dim}: The number of features (\code{nrow}) and cells
#' (\code{ncol})
#'
#' @export
#' @method dim StdAssay
#'
dim.StdAssay <- function(x) {
  return(dim(x = LayerData(object = x)))
}

#' @describeIn StdAssay-methods Cell- and feature-names for an \code{StdAssay}
#'
#' @return \code{dimnames}: Feature (row) and cell (column) names
#'
#' @export
#' @method dimnames StdAssay
#'
dimnames.StdAssay <- function(x) {
  return(list(Features(x = x), Cells(x = x)))
}

#' @describeIn StdAssay-methods Get the first rows of feature-level metadata
#'
#' @return \code{head}: The first \code{n} rows of feature-level metadata
#'
#' @export
#' @method head StdAssay
#'
head.StdAssay <- .head

#' @describeIn StdAssay-methods Get the last rows of feature-level metadata
#'
#' @return \code{tail}: The last \code{n} rows of feature-level metadata
#'
#' @importFrom utils tail
#'
#' @export
#' @method tail StdAssay
#'
tail.StdAssay <- .tail

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @describeIn StdAssay-methods Add feature-level metadata
#'
#' @param value Additional metadata to add
#'
#' @return \code{[[<-}: \code{x} with metadata \code{value} added as \code{i}
#'
#' @export
#'
setMethod(
  f = '[[<-',
  signature = c('x' = 'StdAssay'),
  definition = function(x, i, ..., value) {
    meta.data <- x[[]]
    feature.names <- rownames(x = meta.data)
    if (is.data.frame(x = value)) {
      value <- lapply(
        X = 1:ncol(x = value),
        FUN = function(index) {
          v <- value[[index]]
          names(x = v) <- rownames(x = value)
          return(v)
        }
      )
    }
    err.msg <- "Cannot add more or fewer meta.features information without values being named with feature names"
    if (length(x = i) > 1) {
      # Add multiple bits of feature-level metadata
      value <- rep_len(x = value, length.out = length(x = i))
      for (index in 1:length(x = i)) {
        names.intersect <- intersect(x = names(x = value[[index]]), feature.names)
        if (length(x = names.intersect) > 0) {
          meta.data[names.intersect, i[index]] <- value[[index]][names.intersect]
        } else if (length(x = value) %in% c(nrow(x = meta.data), 1) %||% is.null(x = value)) {
          meta.data[i[index]] <- value[index]
        } else {
          stop(err.msg, call. = FALSE)
        }
      }
    } else {
      # Add a single column to feature-level metadata
      value <- unlist(x = value)
      if (length(x = intersect(x = names(x = value), y = feature.names)) > 0) {
        meta.data[, i] <- value[feature.names]
      } else if (length(x = value) %in% c(nrow(x = meta.data), 1) || is.null(x = value)) {
        meta.data[, i] <- value
      } else {
        stop(err.msg, call. = FALSE)
      }
    }
    slot(object = x, name = 'meta.features') <- meta.data
    return(x)
  }
)

#' @describeIn StdAssay-methods Overview of an \code{StdAssay} object
#'
#' @return \code{show}: Prints summary to \code{\link[base]{stdout}} and
#' invisibly returns \code{NULL}
#'
#' @importFrom utils head
#' @importFrom methods show
#'
#' @export
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
      ncol(x = object), 'cells\n'
    )
    # Feature information
    # if (length(x = VariableFeatures(object = object)) > 0) {
    if (FALSE) {
      top.ten <- head(x = VariableFeatures(object = object), n = 10L)
      top <- 'Top'
      variable <- 'variable'
    } else {
      top.ten <- head(x = Features(x = object), n = 10L)
      top <- 'First'
      variable <- ''
    }
    features <- paste0(
      variable,
      ' feature',
      if (length(x = top.ten) != 1) {
        's'
      },
      ":\n"
    )
    features <- gsub(pattern = '^\\s+', replacement = '', x = features)
    cat(
      top,
      length(x = top.ten),
      features,
      paste(strwrap(x = paste(top.ten, collapse = ', ')), collapse = '\n'),
      '\n'
    )
    # Layer information
    if (length(x = Layers(object = object, data = FALSE))) {
      cat(
        "Additional layers:\n",
        paste(
          strwrap(x = paste(
            Layers(object = object, data = FALSE),
            collapse = ', '
          )),
          collapse = '\n'
        )
      )
    }
    return(invisible(x = NULL))
  }
)

setValidity(
  Class = 'StdAssay',
  method = function(object) {
    valid <- NULL
    # Check data
    data.dims <- na.dims <- c(NA_integer_, NA_integer_)
    if (IsMatrixEmpty(x = GetAssayData(object = object, slot = 'data'))) {
      valid <- c(valid, "'data' cannot be empty")
    } else {
      data.dims <- dim(x = GetAssayData(object = object, slot = 'data'))
      if (length(x = data.dims) != 2) {
        valid <- c(valid, "'data' must be a two-dimensional object")
        data.dims <- na.dims
      }
    }
    # Check layers
    if (length(x = slot(object = object, name = 'layers'))) {
      if (!IsNamedList(x = slot(object = object, name = 'layers'))) {
        valid <- c(valid, "'layers' must be a named list")
      }
      for (layer in Layers(object = object, data = FALSE)) {
        ldat <- LayerData(object = object, which = layer)
        if (ncol(x = ldat) != data.dims[2]) {
          valid <- c(valid, "All layers must have the same cells as 'data'")
          break
        } else if (nrow(x = ldat) > data.dims[1]) {
          valid <- c(
            valid,
            "All layers must have the same or a subset of features as 'data'"
          )
          break
        }
      }
    }
    # Check features
    fmat <- slot(object = object, name = 'features')
    if (nrow(x = fmat) != data.dims[1]) {
      valid <- c(valid, "'features' must have the same features as 'data'")
    } else if (is.null(x = rownames(x = fmat))) {
      valid <- c(valid, "'features' must have rownames assigned")
    } else if (!is.logical(x = fmat)) {
      valid <- c(valid, "'features' must be a logical matrix")
    }
    for (layer in Layers(object = object, data = FALSE)) {
      if (!layer %in% colnames(x = fmat)) {
        valid <- c(valid, "All layers must have a column in the feature matrix")
        break
      } else if (sum(fmat[, layer], na.rm = TRUE) < 1L) {
        valid <- c(valid, "All layers must have at least one feature present")
        break
      }
    }
    # TODO: Check variable features
    # TODO: Check meta features
    # TODO: Check key
    # TODO: Check misc
    return(valid %||% TRUE)
  }
)

setValidity(
  Class = 'Assay5',
  method = function(object) {
    valid <- NULL
    # Check class of layers
    for (layer in Layers(object = object, data = FALSE)) {
      cls <- inherits(
        x = LayerData(object = object, layer = layer),
        what = 'AnyMatrix'
      )
      if (!isTRUE(x = cls)) {
        valid <- c(valid, "layers must be either a 'matrix' or 'dgCMatrix'")
        break
      }
    }
    return(valid %||% TRUE)
  }
)
