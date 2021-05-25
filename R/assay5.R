#' @include zzz.R
#' @include logmap.R
#' @importFrom methods setAs
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
    layers = 'list',
    key = 'character',
    cells = 'LogMap',
    assay.orig = 'character',
    features = 'LogMap',
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
  contains = 'StdAssay'
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Generic Assay Creation
#'
#' Create an assay object; runs a standardized filtering scheme that
#' works regardless of the direction of the data (eg. cells as columns
#' and features as rows or vice versa) and creates an assay object based
#' on the  initialization scheme defined for \code{\link{StdAssay}}-derived
#' class \code{type}
#'
#' @param counts A two-dimensional expression matrix
#' @param min.cells Include features detected in at least this many cells;
#' will subset the counts matrix as well. To reintroduce excluded features,
#' create a new object with a lower cutoff
#' @param min.features Include cells where at least this many features
#' are detected
#' @param cells Vector of cell names
#' @param features Vector of feature names
#' @param layer Name of layer to store \code{counts} as
#' @param type Type of assay object to create; must be the name of a class
#' that's derived from \code{\link{StdAssay}}
#' @param csum Function for calculating cell sums
#' @param fsum Function for calculating feature sums
#' @param ... Extra parameters passed to \code{\link[methods]{new}} for
#' assay creation; used to set slots not defined by \code{\link{StdAssay}}
#'
#' @return An object of class \code{type} with a layer named \code{layer}
#' containing the data found in \code{counts}
#'
#' @importFrom methods getClass
#' @importFrom utils getS3method methods
#'
#' @export
#'
#' @keywords internal
#'
.CreateStdAssay <- function(
  counts,
  min.cells = 0,
  min.features = 0,
  cells = NULL,
  features = NULL,
  layer = 'counts',
  type = 'Assay5',
  csum = Matrix::colSums,
  fsum = Matrix::rowSums,
  ...
) {
  cdef <- getClass(Class = type)
  if (!'StdAssay' %in% names(x = slot(object = cdef, name = 'contains'))) {
    stop("Class '", type, "' does not inherit from StdAssay")
  }
  fmargin <- getS3method(
    f = '.MARGIN',
    class = unlist(x = strsplit(
      x = methods(generic.function = '.MARGIN', class = type)[1],
      split = '\\.'
    ))[3]
  )
  cdim <- fmargin(object = type, type = 'cells')
  fdim <- fmargin(object = type, type = 'features')
  cells <- cells %||% paste0(
    'Cell_',
    seq_len(length.out = dim(x = counts)[cdim])
  )
  features <- features %||% paste0(
    'Feature',
    seq_len(length.out = dim(x = counts)[fdim])
  )
  # Filter based on min.features
  if (min.features > 0) {
    cells.use <- which(x = csum(counts > 0) >= min.features)
    counts <- if (cdim == 1) {
      counts[cells.use, ]
    } else {
      counts[, cells.use]
    }
    cells <- cells[cells.use]
  }
  # Filter based on min.cells
  if (min.cells > 0) {
    features.use <- which(x = fsum(counts > 0) >= min.cells)
    counts <- if (fdim == 1) {
      counts[features.use, ]
    } else {
      counts[, features.use]
    }
    features <- features[features.use]
  }
  # Create the object
  object <- new(
    Class = type,
    layers = list(),
    features = LogMap(y = features),
    cells = LogMap(y = cells),
    meta.data = EmptyDF(n = length(x = features)),
    ...
  )
  LayerData(object = object, layer = layer) <- counts
  return(object)
}

#' Create a v5 Assay object
#'
#' Create an \code{\link{Assay5}} object from a feature expression matrix;
#' the expected format of the matrix is features x cells
#'
#' @inheritParams .CreateStdAssay
#' @param ... Extra parameters passed to \code{\link{as.sparse}}
#'
#' @return An \code{\link{Assay5}} object
#'
#' @export
#'
CreateAssay5Object <- function(
  counts,
  min.cells = 0,
  min.features = 0,
  layer = 'counts',
  ...
) {
  if (!inherits(x = counts, what = 'AnyMatrix')) {
    counts <- as.sparse(x = counts, ...)
  }
  cells <- colnames(x = counts)
  features <- rownames(x = counts)
  dimnames(x = counts) <- list(NULL, NULL)
  return(.CreateStdAssay(
    counts = counts,
    min.cells =  min.cells,
    min.features = min.features,
    cells = cells,
    features = features,
    layer = layer
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
Cells.StdAssay <- function(x, layer = NULL, ...) {
  layer <- layer %||% DefaultLayer(object = x)
  return(slot(object = x, name = 'cells')[[layer]])
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

#' @rdname DefaultLayer
#' @export
#' @method DefaultLayer StdAssay
#'
DefaultLayer.StdAssay <- function(object, ...) {
  return(Layers(object = object)[1])
}

#' @rdname DefaultLayer
#' @export
#' @method DefaultLayer<- StdAssay
#'
"DefaultLayer<-.StdAssay" <- function(object, ..., value) {
  value <- value[1]
  layers <- Layers(object = object)
  value <- match.arg(arg = value, choices = layers)
  idx <- which(x = layers == value)
  slot(object = object, name = 'layers') <- c(
    slot(object = object, name = 'layers')[idx],
    slot(object = object, name = 'layers')[-idx]
  )
  return(object)
}

#' @param layer
#'
#' @rdname Cells
#' @export
#' @method Features StdAssay
#'
Features.StdAssay <- function(x, layer = NULL, ...) {
  layer <- layer %||% DefaultLayer(object = x)
  return(slot(object = x, name = 'features')[[layer]])
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
GetAssayData.StdAssay <- function(object, slot = 'data', ...) {
  CheckDots(..., fxns = LayerData)
  return(LayerData(object = object, layer = slot, ...))
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
LayerData.StdAssay <- function(object, layer = NULL, ...) {
  layer <- layer[1] %||% DefaultLayer(object = object)
  layer <- match.arg(arg = layer, choices = Layers(object = object))
  return(slot(object = object, name = 'layers')[[layer]])
}

#' @method LayerData Assay5
#' @export
#'
LayerData.Assay5 <- function(object, layer = NULL, dnames = TRUE, ...) {
  ldat <- NextMethod()
  if (isTRUE(x = dnames)) {
    dimnames(x = ldat) <- list(
      Features(x = object, layer = layer),
      Cells(x = object, layer = layer)
    )
  }
  return(ldat)
}

#' @param features,cells Vectors of features/cells to include ...
#'
#' @rdname Layers
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
  # Remove a layer
  if (is.null(x = value)) {
    if (length(x = Layers(object = object)) == 1L) {
      stop("Cannot remove only layer")
    } else if (layer == DefaultLayer(object = object)) {
      DefaultLayer(object = object) <- Layers(object = object)[2]
      warning(
        "Removing default layer, setting default to ",
        DefaultLayer(object = object),
        call. = FALSE,
        immediate. = TRUE
      )
    }
    slot(object = object, name = 'layers')[[layer]] <- NULL
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
  fdim <- .MARGIN(object = object)
  cdim <- .MARGIN(object = object, type = 'cells')
  features <- features %||% rownames(x = value) %||% seq_len(
    length.out = dim(x = value)[fdim]
  )
  cells <- cells %||% colnames(x = value) %||% seq_len(
    length.out = dim(x = value)[cdim]
  )
  # Check features and cells
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
    if (!identical(x = features, y = Features(x = object, layer = layer))) {
      warning(
        "Different features in new layer data than already exists for ",
        layer,
        call. = FALSE,
        immediate. = TRUE
      )
    }
    if (!identical(x = cells, y = Cells(x = object, layer = layer))) {
      warning(
        "Different cells in new layer data than already exists for ",
        layer,
        call. = FALSE,
        immediate. = TRUE
      )
    }
  }
  # Reorder the layer data
  vdims <- list(fmatch, cmatch)
  names(x = vdims) <- c('i', 'j')[c(fdim, cdim)]
  value <- do.call(what = '[', args = c(list(x = value), vdims))
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
"LayerData<-.Assay5" <- function(object, layer, ..., value) {
  if (!inherits(x = value, what = c('AnyMatrix', 'NULL'))) {
    value <- as.sparse(x = value)
  }
  object <- NextMethod()
  validObject(object = object)
  return(object)
}

#' @rdname Layers
#' @method Layers StdAssay
#' @export
#'
Layers.StdAssay <- function(object, ...) {
  layers <- names(x = slot(object = object, name = 'layers'))
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
SetAssayData.StdAssay <- function(object, slot, new.data, ...) {
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
  return(vapply(
    X = c('features', 'cells'),
    FUN = function(s) {
      return(nrow(x = slot(object = x, name = s)))
    },
    FUN.VALUE = numeric(length = 1L),
    USE.NAMES = FALSE
  ))
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

CalcN5 <- function(object) {
  if (IsMatrixEmpty(x = LayerData(object = object))) {
    return(NULL)
  }
  return(list(
    nCount = colSums(x = object),
    nFeature = colSums(x = LayerData(object = object) > 0)
  ))
}

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
    # Check layers
    dorder <- c(
      features = .MARGIN(object = object, type = 'features'),
      cells = .MARGIN(object = object, type = 'cells')
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
      if (ldims[1] > adims[1]) {
        valid <- c(
          valid,
          paste0(
            "Layers may not have more features than present in the assay ",
            "(offending layer: ",
            layer,
            ")"
          )
        )
      }
      if (ldims[2] != adims[2]) {
        valid <- c(
          valid,
          paste0(
            "Layers must have the same cells as present in the assay ",
            "(offending layer: ",
            layer,
            ")"
          )
        )
      }
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
    return(TRUE)
    # Check class of layers
    for (layer in Layers(object = object, data = FALSE)) {
      cls <- inherits(
        x = LayerData(object = object, layer = layer, dnames = FALSE),
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
