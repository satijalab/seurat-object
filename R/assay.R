#' @include zzz.R
#' @include generics.R
#' @include default.R
#' @include graph.R
#' @include keymixin.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setClassUnion(name = 'AnyMatrix', members = c("matrix", "dgCMatrix"))

#' The Assay Class
#'
#' The Assay object is the basic unit of Seurat; each Assay stores raw,
#' normalized, and scaled data as well as cluster information, variable
#' features, and any other assay-specific metadata. Assays should contain single
#' cell expression data such as RNA-seq, protein, or imputed expression data.
#'
#' @slot counts Unnormalized data such as raw counts or TPMs
#' @slot data Normalized expression data
#' @slot scale.data Scaled expression data
# @slot key Key for the Assay
#' @slot assay.orig Original assay that this assay is based off of. Used to
#' track assay provenance
#' @slot var.features Vector of features exhibiting high variance across
#' single cells
#' @slot meta.features Feature-level metadata
# @slot misc Utility slot for storing additional data associated with the assay
#' @template slot-misc
#' @template slot-key
#'
#' @name Assay-class
#' @rdname Assay-class
#' @exportClass Assay
#'
#' @family assay
#'
#' @aliases Assay
#'
setClass(
  Class = 'Assay',
  contains = 'KeyMixin',
  slots = c(
    counts = 'AnyMatrix',
    data = 'AnyMatrix',
    scale.data = 'matrix',
    # key = 'character',
    assay.orig = 'OptionalCharacter',
    var.features = 'vector',
    meta.features = 'data.frame',
    misc = 'OptionalList'
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Create an Assay object
#'
#' Create an Assay object from a feature (e.g. gene) expression matrix. The
#' expected format of the input matrix is features x cells.
#'
#' Non-unique cell or feature names are not allowed. Please make unique before
#' calling this function.
#'
#' @param counts Unnormalized data such as raw counts or TPMs
#' @param data Prenormalized data; if provided, do not pass \code{counts}
#' @param min.cells Include features detected in at least this many cells. Will
#' subset the counts matrix as well. To reintroduce excluded features, create a
#' new object with a lower cutoff
#' @param min.features Include cells where at least this many features are
#' detected
#' @param key Optional key to initialize assay with
#' @param check.matrix Check counts matrix for NA, NaN, Inf, and
#' non-integer values
#' @param ... Arguments passed to \code{\link{as.sparse}}
#'
#' @return A \code{\link{Assay}} object
#'
#' @importFrom methods as
#' @importFrom Matrix colSums rowSums
#'
#' @export
#'
#' @family assay
#'
#' @examples
#' \dontrun{
#' pbmc_raw <- read.table(
#'   file = system.file('extdata', 'pbmc_raw.txt', package = 'Seurat'),
#'   as.is = TRUE
#' )
#' pbmc_rna <- CreateAssayObject(counts = pbmc_raw)
#' pbmc_rna
#' }
#'
CreateAssayObject <- function(
  counts,
  data,
  min.cells = 0,
  min.features = 0,
  key = NULL,
  check.matrix = FALSE,
  ...
) {
  if (missing(x = counts) && missing(x = data)) {
    abort(message = "Must provide either 'counts' or 'data'")
  } else if (!missing(x = counts) && !missing(x = data)) {
    abort(message = "Either 'counts' or 'data' must be missing; both cannot be provided")
  } else if (!missing(x = counts)) {
    # check that dimnames of input counts are unique
    if (anyDuplicated(x = rownames(x = counts))) {
      warn(
        message = "Non-unique features (rownames) present in the input matrix, making unique"
      )
      rownames(x = counts) <- make.unique(names = rownames(x = counts))
    }
    if (anyDuplicated(x = colnames(x = counts))) {
      warn(
        message = "Non-unique cell names (colnames) present in the input matrix, making unique"
      )
      colnames(x = counts) <- make.unique(names = colnames(x = counts))
    }
    if (is.null(x = colnames(x = counts))) {
      abort(message = "No cell names (colnames) names present in the input matrix")
    }
    if (any(rownames(x = counts) == '')) {
      abort(message = "Feature names of counts matrix cannot be empty")
    }
    if (nrow(x = counts) > 0 && is.null(x = rownames(x = counts))) {
      abort(message = "No feature names (rownames) names present in the input matrix")
    }
    if (!inherits(x = counts, what = 'dgCMatrix')) {
      if (inherits(x = counts, what = "data.frame")) {
        counts <- as.sparse(x = counts, ...)
      } else {
        counts <- as.sparse(x = counts)
      }
    }
    if (isTRUE(x = check.matrix)) {
      CheckMatrix(object = counts)
    }
    # Filter based on min.features
    if (min.features > 0) {
      nfeatures <- Matrix::colSums(x = counts > 0)
      counts <- counts[, which(x = nfeatures >= min.features)]
    }
    # filter genes on the number of cells expressing
    if (min.cells > 0) {
      num.cells <- Matrix::rowSums(x = counts > 0)
      counts <- counts[which(x = num.cells >= min.cells), ]
    }
    data <- counts
  } else if (!missing(x = data)) {
    # check that dimnames of input data are unique
    if (anyDuplicated(x = rownames(x = data))) {
      warn(
        message = "Non-unique features (rownames) present in the input matrix, making unique"
      )
      rownames(x = data) <- make.unique(names = rownames(x = data))
    }
    if (anyDuplicated(x = colnames(x = data))) {
      warn(
        message = "Non-unique cell names (colnames) present in the input matrix, making unique"
      )
      colnames(x = data) <- make.unique(names = colnames(x = data))
    }
    if (is.null(x = colnames(x = data))) {
      abort(message = "No cell names (colnames) names present in the input matrix")
    }
    if (any(rownames(x = data) == '')) {
      abort(message = "Feature names of data matrix cannot be empty", call. = FALSE)
    }
    if (nrow(x = data) > 0 && is.null(x = rownames(x = data))) {
      abort(message = "No feature names (rownames) names present in the input matrix")
    }
    if (min.cells != 0 | min.features != 0) {
      warn(
        message = "No filtering performed if passing to data rather than counts"
      )
    }
    counts <- new(Class = 'matrix')
  }
  # Ensure row- and column-names are vectors, not arrays
  if (!is.vector(x = rownames(x = counts))) {
    rownames(x = counts) <- as.vector(x = rownames(x = counts))
  }
  if (!is.vector(x = colnames(x = counts))) {
    colnames(x = counts) <- as.vector(x = colnames(x = counts))
  }
  if (!is.vector(x = rownames(x = data))) {
    rownames(x = data) <- as.vector(x = rownames(x = data))
  }
  if (!is.vector(x = colnames(x = data))) {
    colnames(x = data) <- as.vector(x = colnames(x = data))
  }
  counts <- CheckFeaturesNames(data = counts)
  data <- CheckFeaturesNames(data = data)
  # Initialize meta.features
  init.meta.features <- data.frame(row.names = rownames(x = data))
  misc <- if (.GetSeuratCompat() < '5.0.0') {
    list()
  } else {
    calcN_option <- getOption(
      x = 'Seurat.object.assay.calcn',
      default =  Seurat.options$Seurat.object.assay.calcn
    )
    list(calcN = calcN_option %||% TRUE)
  }
  assay <- new(
    Class = 'Assay',
    counts = counts,
    data = data,
    scale.data = new(Class = 'matrix'),
    key = Key(object = key)[1L] %||% '',
    meta.features = init.meta.features,
    misc = misc
  )
  return(assay)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom Matrix colSums
#'
#' @method .CalcN Assay
#' @export
#'
.CalcN.Assay <- function(object, layer = 'counts', ...) {
  layer <- tryCatch(
    expr = Layers(object = object, search = layer),
    error = function(...) NULL
  )
  if (is.null(x = layer)) {
    return(NULL)
  }
  ldat <- LayerData(object = object, layer = layer)
  if (IsMatrixEmpty(x = ldat) || !inherits(x = ldat, what = 'Matrix')) {
    return(NULL)
  }
  cells_stat <- .CalcN.default(object = ldat)
  return(cells_stat)
}

#' @rdname AddMetaData
#'
# @templateVar fname AddMetaData
# @templateVar version 4
# @template name-oldv
#'
#' @export
#' @method AddMetaData Assay
#'
AddMetaData.Assay <- function(object, metadata, col.name = NULL) {
  if (is.null(x = col.name) && (is.atomic(x = metadata) && !is.matrix(x = metadata))) {
    abort(message = "'col.name' must be provided for atomic meta data")
  }
  if (inherits(x = metadata, what = c('matrix', 'Matrix'))) {
    metadata <- as.data.frame(x = metadata)
  }
  col.name <- col.name %||% names(x = metadata) %||% colnames(x = metadata)
  if (is.null(x = col.name)) {
    abort(message = "No metadata name provided and could not infer it from metadata object")
  }
  object[[col.name]] <- metadata
  return(object)
}


#' @rdname DefaultAssay
#' @export
#' @method DefaultAssay Assay
#'
DefaultAssay.Assay <- function(object, ...) {
  object <- UpdateSlots(object = object)
  return(slot(object = object, name = 'assay.orig'))
}

#' @rdname DefaultAssay
#' @export
#' @method DefaultAssay<- Assay
#'
"DefaultAssay<-.Assay" <- function(object, ..., value) {
  object <- UpdateSlots(object = object)
  slot(object = object, name = 'assay.orig') <- value
  return(object)
}

#' @rdname DefaultLayer
#' @method DefaultLayer Assay
#' @export
#'
DefaultLayer.Assay <- function(object, ...) {
  return('data')
}

#' @method Features Assay
#' @export
#'
Features.Assay <- function(
  x,
  layer = c('data', 'scale.data', 'counts'),
  slot = deprecated(),
  ...
) {
  if (is_present(arg = slot)) {
    .Deprecate(
      when = '5.0.0',
      what = 'Features(slot = )',
      with = 'Features(layer = )'
    )
    layer <- slot
  }
  layer <- layer[1L] %||% 'data'
  layer <- match.arg(arg = layer)
  features <- rownames(x = GetAssayData(object = x, layer = layer))
  if (!length(x = features)) {
    features <- NULL
  }
  return(features)
}

#' @method FetchData Assay
#' @export
#'
FetchData.Assay <- function(
  object,
  vars,
  cells = NULL,
  layer = NULL,
  slot = deprecated(),
  ...
) {
  if (is_present(arg = slot)) {
    .Deprecate(
      when = '5.0.0',
      what = 'FetchData(slot = )',
      with = 'FetchData(layer = )'
    )
    layer <- layer %||% slot
  }
  # Identify slot to use
  layer <- layer %||% 'data'
  layer <- match.arg(arg = layer, choices = c('counts', 'data', 'scale.data'))
  # Identify cells to use
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  cells.orig <- cells
  cells <- intersect(x = cells, y = colnames(x = object))
  if (length(x = cells) != length(x = cells.orig)) {
    warn(message = paste(
      "Removing",
      length(x = cells.orig) - length(x = cells),
      "cells not present in the assay"
    ))
  }
  # Check vars
  orig <- vars
  vars <- gsub(
    pattern = paste0('^', Key(object = object)),
    replacement = '',
    x = vars
  )
  # Pull expression information
  mat <- GetAssayData(object = object, layer = layer)
  if (IsMatrixEmpty(x = mat)) {
    abort(message = paste("Layer", sQuote(x = layer), "is empty in this assay"))
  }
  vars <- intersect(x = vars, y = rownames(x = mat))
  tf <- .GetMethod(fxn = 't', cls = class(x = mat))
  data.fetched <- as.data.frame(x = as.matrix(
    x = tf(x = mat[vars, cells, drop = FALSE])
  ))
  # Add keys to keyed vars
  keyed.features <- paste0(Key(object = object), names(x = data.fetched))
  keyed.idx <- which(x = keyed.features %in% orig)
  if (length(x = keyed.idx)) {
    names(x = data.fetched)[keyed.idx] <- keyed.features[keyed.idx]
  }
  # Check the final list of features
  missing <- setdiff(x = orig, y = names(x = data.fetched))
  if (length(x = missing) == length(x = orig)) {
    abort(message = "None of the requested features found")
  } else if (length(x = missing)) {
    warn(message = paste(
      "The following features could not be found",
      paste(missing, collapse = ', ')
    ))
  }
  return(data.fetched)
}

#' @rdname AssayData
#' @export
#' @method GetAssayData Assay
#'
#' @examples
#' # Get the data directly from an Assay object
#' GetAssayData(pbmc_small[["RNA"]], layer = "data")[1:5,1:5]
#'
GetAssayData.Assay <- function(
  object,
  layer = c('data', 'scale.data', 'counts'),
  slot = deprecated(),
  ...
) {
  CheckDots(...)
  if (is_present(arg = slot)) {
    .Deprecate(
      when = '5.0.0',
      what = 'GetAssayData(slot = )',
      with = 'GetAssayData(layer = )'
    )
    layer <- slot
  }
  layer <- layer[1L] %||% 'data'
  layer <- match.arg(arg = layer)
  return(methods::slot(object = object, name = layer))
}

#' @rdname VariableFeatures
#' @export
#' @method HVFInfo Assay
#'
#' @examples
#' # Get the HVF info directly from an Assay object
#' HVFInfo(pbmc_small[["RNA"]], method = 'vst')[1:5, ]
#'
HVFInfo.Assay <- function(
  object,
  method,
  status = FALSE,
  selection.method = deprecated(),
  ...
) {
  CheckDots(...)
  if (is_present(arg = selection.method)) {
    .Deprecate(
      when = '5.0.0',
      what = 'HVFInfo(selection.method = )',
      with = 'HVFInfo(method = )'
    )
    method <- selection.method
  }
  disp.methods <- c('mean.var.plot', 'dispersion', 'disp')
  if (tolower(x = method) %in% disp.methods) {
    method <- 'mvp'
  }
  method <- switch(
    EXPR = tolower(x = method),
    sctransform = 'sct',
    method
  )
  vars <- switch(
    EXPR = method,
    vst = c('mean', 'variance', 'variance.standardized'),
    mvp = c('mean', 'dispersion', 'dispersion.scaled'),
    sct = c('gmean', 'variance', 'residual_variance'),
    abort(message = paste("Unknown method:", sQuote(x = method)))
  )
  tryCatch(
    expr = hvf.info <- object[[paste(method, vars, sep = '.')]],
    error = function(e) {
      stop(
        "Unable to find highly variable feature information for method '",
        method,
        "'",
        call. = FALSE
      )
    }
  )
  colnames(x = hvf.info) <- vars
  if (status) {
    hvf.info$variable <- object[[paste0(method, '.variable')]]
  }
  return(hvf.info)
}

#' @rdname Key
#' @export
#' @method Key Assay
#'
#' @examples
#' # Get an Assay key
#' Key(pbmc_small[["RNA"]])
#'
Key.Assay <- function(object, ...) {
  CheckDots(...)
  return(slot(object = object, name = 'key'))
}

#' @rdname Key
#' @export
#' @method Key<- Assay
#'
#' @examples
#' # Set the key for an Assay
#' Key(pbmc_small[["RNA"]]) <- "newkey_"
#' Key(pbmc_small[["RNA"]])
#'
"Key<-.Assay" <- function(object, ..., value) {
  CheckDots(...)
  slot(object = object, name = 'key') <- value
  return(object)
}

#' @rdname Layers
#' @method LayerData Assay
#' @export
#'
LayerData.Assay <- function(
  object,
  layer = NULL,
  cells = NULL,
  features = NULL,
  slot = deprecated(),
  ...
) {
  if (is_present(arg = slot)) {
    deprecate_stop(
      when = "5.0.0",
      what = "LayerData(slot = )",
      with = "LayerData(layer = )"
    )
  }
  # Figure out which matrix we're pulling
  layer <- layer[1L] %||% "data"

  # layer <- match.arg(
  #   arg = layer,
  #   choices = Layers(object = object, search = FALSE)
  # )
  # Handle empty layers
  if (IsMatrixEmpty(x = methods::slot(object = object, name = layer))) {
    msg <- paste("Layer", sQuote(x = layer), "is empty")
    opt <- getOption(
      x = 'Seurat.object.assay.v3.missing_layer',
      default = Seurat.options$Seurat.object.assay.v3.missing_layer
    )
    opt <- tryCatch(
      expr = arg_match0(arg = opt, values = c('matrix', 'null', 'error')),
      error = function(...) {
        return(Seurat.options$Seurat.object.assay.v3.missing_layer)
      }
    )
    if (opt == 'error') {
      abort(message = msg)
    }
    warn(message = msg)
    return(switch(
      EXPR = opt,
      matrix = switch(
        EXPR = layer,
        scale.data = new(Class = 'matrix'),
        new(Class = 'dgCMatrix')
      ),
      NULL
    ))
  }
  # Allow cell/feature subsets
  cells <- cells %||% colnames(x = object)
  features <- features %||% Features(x = object, layer = layer)
  if (is_bare_integerish(x = cells, finite = TRUE)) {
    cells <- colnames(x = object)[cells]
  }
  cells <- arg_match(
    arg = cells,
    values = colnames(x = object),
    multiple = TRUE
  )
  if (is_bare_integerish(x = features, finite = TRUE)) {
    features <- Features(x = object, layer = layer)[features]
  }
  features <- arg_match(
    arg = features,
    values = Features(x = object, layer = layer),
    multiple = TRUE
  )
  if (length(x = features) == 0) {
    stop('features are not found')
  }
  # Pull the matrix for the cells/features requested
  return(methods::slot(object = object, name = layer)[features, cells, drop = FALSE])
}

#' @rdname Layers
#' @method LayerData<- Assay
#' @export
#'
"LayerData<-.Assay" <- function(object, layer, ..., value) {
  # Check the layer name
  layer <- layer[1L]
  layer <- match.arg(
    arg = layer,
    choices = Layers(object = object, search = FALSE)
  )
  # Allow short-hand switch
  if (rlang::is_scalar_character(x = value)) {
    value <- arg_match0(arg = value, values = Layers(object = object))
    value <- LayerData(object = object, layer = value)
  }
  # Prepare an empty matrix if value is NULL
  value <- value %||% switch(
    EXPR = layer,
    scale.data = new(Class = 'matrix'),
    counts = new(Class = 'dgCMatrix'),
    data = {
      if (IsMatrixEmpty(x = suppressWarnings(expr = LayerData(object = object, layer = 'counts')))) {
        abort(message = "Cannot remove the data layer")
      }
      warn(message = "Resetting the data matrix to the raw counts")
      LayerData(object = object, layer = 'counts')
    }
  )
  # Check the class of the matrix
  if (!inherits(x = value, what = c('matrix', 'dgCMatrix'))) {
   abort(message = paste(
     "'value' must be a 'matrix' or 'dgCMatrix' in v3 Assays, not a",
     sQuote(x = class(x = value)[1L])
   ))
  }
  if (!IsMatrixEmpty(x = value)) {
    vnames <- dimnames(x = value)
    # Check presence of cell- and feature-names
    if (is.null(x = vnames)) {
      if (!all(dim(x = value) == dim(x = object))) {
        abort(message = "New data must have feature and cell names")
      }
      dimnames(x = value) <- dimnames(x = object)
    } else if (any(.IsNull(x = vnames)) || !all(unlist(x = lapply(X = vnames, FUN = nzchar)))) {
      abort(message = "New data must have feature and cell names")
    }
    # Remove underscores from feature names
    if (any(grepl(pattern = '_', x = rownames(x = value)))) {
      warn(
        message = "Feature names cannot have underscores ('_'), replacing with dashes ('-')"
      )
      rownames(x = value) <- gsub(
        pattern = '_',
        replacement = '-',
        x = rownames(x = value)
      )
    }
    # Check the the cells
    if (ncol(x = value) != ncol(x = object)) {
      abort(message = "The new data must have the same number of cells as the current data")
    } else if (!all(colnames(x = value) %in% colnames(x = object))) {
      abort(message = "The new data must have the same cells as the current data")
    }
    value <- value[, colnames(x = object), drop = FALSE]
    # Check the features
    if (!any(rownames(x = value) %in% rownames(x = object))) {
      abort(message = "None of the features provided are present in the existing data")
    } else if (!all(rownames(x = value) %in% rownames(x = object))) {
      warn(message = "Extra features present in the the new data compared to the existing data")
    }
    features <- intersect(x = rownames(x = object), y = rownames(x = value))
    value <- value[features, , drop = FALSE]
    if (layer %in% c('counts', 'data') && nrow(x = value) != nrow(x = object)) {
      abort(message = "The new data must have the same number of features as the current data")
    }
  }
  slot(object = object, name = layer) <- value
  validObject(object = object)
  return(object)
}

#' @rdname Layers
#' @method Layers Assay
#' @export
#'
Layers.Assay <- function(object, search = NA, ...) {
  layers <- c('counts', 'data', 'scale.data')
  if (isFALSE(x = search)) {
    return(layers)
  }
  layers <- Filter(
    f = function(x) {
      return(!IsMatrixEmpty(x = slot(object = object, name = x)))
    },
    x = layers
  )
  if (!length(x = layers)) {
    abort(message = "All matrices are empty in this Assay")
  }
  if (is.null(x = search)) {
    return(DefaultLayer(object = object))
  }
  if (!is_na(x = search)) {
    layers <- intersect(x = search, y = layers)
    if (length(x = layers) == 0) {
      warning(
        "Layer ",
        search,
        " isn't present in the assay ",
        deparse(expr = substitute(expr = object)),
        "; returning NULL",
        call. = FALSE,
        immediate. = TRUE
      )
      return(NULL)
    }
  }
  return(layers)
}

#' @param slot Name of specific bit of meta data to pull
#'
#' @rdname Misc
#' @export
#' @method Misc Assay
#'
Misc.Assay <- .Misc

#' @rdname Misc
#' @export
#' @method Misc<- Assay
#'
"Misc<-.Assay" <- `.Misc<-`

#' @param new.names vector of new cell names
#'
#' @rdname RenameCells
#' @export
#' @method RenameCells Assay
#'
#' @examples
#' # Rename cells in an Assay
#' head(x = colnames(x = pbmc_small[["RNA"]]))
#' renamed.assay <- RenameCells(
#'     pbmc_small[["RNA"]],
#'     new.names = paste0("A_", colnames(x = pbmc_small[["RNA"]]))
#' )
#' head(x = colnames(x = renamed.assay))
#'
RenameCells.Assay <- function(object, new.names = NULL, ...) {
  CheckDots(...)
  names(new.names) <- NULL
  for (data.slot in c("counts", "data", "scale.data")) {
    old.data <- GetAssayData(object = object, layer = data.slot)
    if (ncol(x = old.data) <= 1) {
      next
    }
    colnames(x = slot(object = object, name = data.slot)) <- new.names
  }
  return(object)
}

#' @importFrom stats na.omit
#'
#' @rdname AssayData
#' @export
#' @method SetAssayData Assay
#'
#' @examples
#' # Set an Assay layer directly
#' count.data <- GetAssayData(pbmc_small[["RNA"]], layer = "counts")
#' count.data <- as.matrix(x = count.data + 1)
#' new.assay <- SetAssayData(pbmc_small[["RNA"]], layer = "counts", new.data = count.data)
#'
SetAssayData.Assay <- function(
  object,
  layer = c('data', 'scale.data', 'counts'),
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
  CheckDots(...)
  layer <- layer[1]
  layer <- match.arg(arg = layer)
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
    if (layer == 'scale.data' && nrow(x = new.data) > num.counts) {
      warning(
        "Adding more features than present in current data",
        call. = FALSE,
        immediate. = TRUE
      )
    } else if (layer %in% c('counts', 'data') && nrow(x = new.data) != num.counts) {
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
    if (layer %in% c('counts', 'data') && !all(dim(x = new.data) == dim(x = object))) {
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
  slot(object = object, name = layer) <- new.data
  return(object)
}

#' @param decreasing Return features in decreasing order (most spatially
#' variable first).
#'
#' @rdname VariableFeatures
#' @export
#' @method SpatiallyVariableFeatures Assay
#'
SpatiallyVariableFeatures.Assay <- function(
  object,
  method = "moransi",
  decreasing = TRUE,
  selection.method = deprecated(),
  ...
) {
  CheckDots(...)
  if (is_present(arg = selection.method)) {
    .Deprecate(
      when = "5.0.0",
      what = "SpatiallyVariableFeatures(selection.method = )",
      with = "SpatiallyVariableFeatures(method = )"
    )
    method <- selection.method
  }
  vf <- SVFInfo(object = object, method = method, status = TRUE)
  vf <- vf[rownames(vf)[which(vf[, "variable"][, 1])], ]
  if (!is.null(x = decreasing)) {
    vf <- vf[order(vf[, "rank"][, 1], decreasing = !decreasing), ]
  }
  return(rownames(vf))
}

#' @rdname VariableFeatures
#' @export
#' @method SVFInfo Assay
#'
SVFInfo.Assay <- function(
  object,
  method = c("markvariogram", "moransi"),
  status = FALSE,
  selection.method = deprecated(),
  ...
) {
  CheckDots(...)
  if (is_present(arg = selection.method)) {
    .Deprecate(
      when = "5.0.0",
      what = "SVFInfo(selection.method = )",
      with = "SVFInfo(method = )"
    )
    method <- selection.method
  }
  method <- match.arg(method)
  vars <- switch(
    EXPR = method,
    markvariogram = grep(
      pattern = "r.metric",
      x = colnames(object[[]]),
      value = TRUE
    ),
    moransi = grep(
      pattern = "MoransI",
      x = colnames(object[[]]),
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
    svf.info$variable <- object[[paste0(method, ".spatially.variable")]]
    svf.info$rank <- object[[paste0(method, ".spatially.variable.rank")]]
  }
  return(svf.info)
}

#' @rdname VariableFeatures
#' @export
#' @method VariableFeatures Assay
#'
VariableFeatures.Assay <- function(
  object,
  method = NULL,
  selection.method = deprecated(),
  ...
) {
  suppressWarnings(CheckDots(...))
  if (is_present(arg = selection.method)) {
    .Deprecate(
      when = '5.0.0',
      what = 'VariableFeatures(selection.method = )',
      with = 'VariableFeatures(method = )'
    )
    method <- selection.method
  }
  if (!is.null(x = method)) {
    vf <- HVFInfo(
      object = object,
      method = method,
      status = TRUE
    )
    return(rownames(x = vf)[which(x = vf[, "variable"][, 1])])
  }
  return(slot(object = object, name = 'var.features'))
}

#' @rdname VariableFeatures
#' @export
#' @method VariableFeatures<- Assay
#'
"VariableFeatures<-.Assay" <- function(object, ..., value) {
  CheckDots(...)
  if (!length(x = value)) {
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
      abort(message = "None of the features provided are in this Assay object")
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

#' @param cells Subset of cell names
#' @param expression A predicate expression for feature/variable expression,
#' can evaluate anything that can be pulled by \code{FetchData}; please note,
#' you may need to wrap feature names in backticks (\code{``}) if dashes
#' between numbers are present in the feature name
#' @param invert Invert the selection of cells
#'
#' @importFrom stats na.omit
#'
#' @rdname WhichCells
#' @export
#' @method WhichCells Assay
#'
WhichCells.Assay <- function(
  object,
  cells = NULL,
  expression,
  invert = FALSE,
  ...
) {
  CheckDots(...)
  cells <- cells %||% colnames(x = object)
  if (!missing(x = expression) && !is.null(x = substitute(expr = expression))) {
    key.pattern <- paste0('^', Key(object = object))
    expr <- if (tryCatch(expr = is_quosure(x = expression), error = function(...) FALSE)) {
      expression
    } else if (is.call(x = enquo(arg = expression))) {
      enquo(arg = expression)
    } else {
      parse(text = expression)
    }
    expr.char <- suppressWarnings(expr = as.character(x = expr))
    expr.char <- unlist(x = lapply(X = expr.char, FUN = strsplit, split = ' '))
    expr.char <- gsub(
      pattern = key.pattern,
      replacement = '',
      x = expr.char,
      perl = TRUE
    )
    expr.char <- gsub(
      pattern = '(',
      replacement = '',
      x = expr.char,
      fixed = TRUE
    )
    expr.char <- gsub(
      pattern = '`',
      replacement = '',
      x = expr.char
    )
    vars.use <- which(x = expr.char %in% rownames(x = object))
    expr.char <- expr.char[vars.use]
    data.subset <- FetchData(object = object, vars = expr.char)
    cells <- rownames(x = data.subset)[eval_tidy(expr = expr, data = data.subset)]
  }
  if (invert) {
    cells <- colnames(x = object)[!colnames(x = object) %in% cells]
  }
  cells <- na.omit(object = unlist(x = cells, use.names = FALSE))
  return(as.character(x = cells))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @inherit .DollarNames.Assay5 return title details sections seealso
#'
#' @description  Autocompletion for \code{$} access on an
#' \code{\link{Assay}} object
#'
#' @inheritParams utils::.DollarNames
#' @param x An \code{\link{Assay}} object
#'
#' @importFrom utils .DollarNames
#'
#' @keywords internal
#'
#' @method .DollarNames Assay
#' @export
#'
#' @concept assay
#'
.DollarNames.Assay <- function(x, pattern = '') {
  slots.avial <- Layers(x)
  slots.avial <- as.list(slots.avial)
  names(slots.avial) <- unlist(slots.avial)
  return(.DollarNames(x = slots.avial, pattern = pattern))
}

#' @inherit $.Assay5 return title description details sections params
#'
#' @param x An \code{\link{Assay}} object
#'
#' @method $ Assay
#' @export
#'
#' @family assay
#'
#' @examples
#' rna <- pbmc_small[["RNA"]]
#'
#' # Fetch a layer with `$`
#' rna$data[1:10, 1:4]
#'
"$.Assay" <- function(x, i) {
  return(LayerData(object = x, layer = i))
}

#' @rdname cash-.Assay
#'
#' @method $<- Assay
#' @export
#'
#' @examples
#' # Add a layer with `$`
#' rna$data <- rna$counts
#' rna$data[1:10, 1:4]
#'
"$<-.Assay" <- function(x, i, value) {
  LayerData(object = x, layer = i) <- value
  return(x)
}

#' @inherit [.Assay5 return title description details sections
#'
#' @inheritParams [.Assay5
#' @param x An \code{\link{Assay}} object
#' @param j Ignored
#'
#' @method [ Assay
#' @export
#'
#' @order 1
#'
#' @seealso \code{\link{LayerData}}
#'
#' @family assay
#'
#' @examples
#' rna <- pbmc_small[["RNA"]]
#'
#' # Get a vector of layer names in this assay
#' rna[]
#'
#' # Fetch layer data
#' rna["data"][1:10, 1:4]
#'
"[.Assay" <- function(x, i = missing_arg(), j = missing_arg(), ...) {
  if (getOption(x = 'Seurat.object.assay.brackets', default = 'v5') == 'v3') {
    if (is_missing(x = i)) {
      i <- seq_len(length.out = nrow(x = x))
    }
    if (is_missing(x = j)) {
      j <- seq_len(length.out = ncol(x = x))
    }
    return(LayerData(
      object = x,
      layer = DefaultLayer(object = x)[1L],
      cells = j,
      features = i
    ))
  }
  if (is_missing(x = i)) {
    return(Layers(object = x))
  }
  return(LayerData(object = x, layer = i, ...))
}

#' @inherit [[.Assay5 return title description details sections
#'
#' @inheritParams [[.Assay5
#' @param x An \code{\link{Assay}} object
#'
#' @method [[ Assay
#' @export
#'
#' @family assay
#'
#' @order 1
#'
#' @examples
#' rna <- pbmc_small[["RNA"]]
#'
#' # Pull the entire feature-level meta data data frame
#' head(rna[[]])
#'
#' # Pull a specific column of feature-level meta data
#' head(rna[["vst.mean"]])
#' head(rna[["vst.mean", drop = TRUE]])
#'
"[[.Assay" <- function(x, i, ..., drop = FALSE) {
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

#' @inherit dim.Assay5 return title description details sections
#'
#' @inheritParams [.Assay
#'
#' @method dim Assay
#' @export
#'
#' @family assay
#'
#' @examples
#' rna <- pbmc_small[["RNA"]]
#' dim(rna)
#'
dim.Assay <- function(x) {
  return(dim(x = GetAssayData(object = x)))
}

#' @inherit dimnames.Assay5 title description details sections
#'
#' @inheritParams [.Assay
#'
#' @return \code{dimnames}: A two-length list with the following values:
#' \itemize{
#'  \item A character vector will all features in \code{x}
#'  \item A character vector will all cells in \code{x}
#' }
#'
#' @method dimnames Assay
#' @export
#'
#' @family assay
#' @family dimnames
#'
#' @examples
#' rna <- pbmc_small[["RNA"]]
#'
#' # Feature and cell names can be acquired with `rownames` and `colnames`
#' head(rownames(rna))
#' head(colnames(rna))
#'
dimnames.Assay <- function(x) {
  return(dimnames(x = GetAssayData(object = x)))
}

#' @param value A two-length list where the first entry is the existing feature
#' names for \code{x} and the second entry is the \emph{updated} cell names
#' for \code{x}
#'
#' @return \code{dimnames<-}: \code{x} with the cell names updated to those
#' in \code{value[[2L]]}
#'
#' @rdname dimnames.Assay
#'
#' @method dimnames<- Assay
#' @export
#'
#' @examples
#' # Cell names can be updated with `colnames<-`
#' colnames(rna)[1] <- "newcell"
#' head(colnames(rna))
#'
"dimnames<-.Assay" <- function(x, value) {
  op <- options(Seurat.object.validate = FALSE)
  on.exit(expr = options(op), add = TRUE)
  # Check the provided dimnames
  msg <- "Invalid 'dimnames' given for a Seurat object"
  if (!is_bare_list(x = value, n = 2L)) {
    abort(message = msg)
  } else if (!all(sapply(X = value, FUN = length) == dim(x = x))) {
    abort(message = msg)
  }
  value <- lapply(X = value, FUN = as.character)
  # Warn about changing features
  if (!all(value[[1L]] == rownames(x = slot(object = x, name = 'data')))) {
    warn(message = "Changing feature names in v3 Assays is not supported")
  }
  # Set cell names
  for (lyr in c('counts', 'data', 'scale.data')) {
    if (!IsMatrixEmpty(x = slot(object = x, name = lyr))) {
      suppressWarnings(expr = colnames(x = slot(object = x, name = lyr)) <- value[[2L]])
    }
  }
  # Validate and return the Seurat object
  options(op)
  validObject(object = x)
  return(x)
}

#' @rdname sub-sub-.Assay
#'
#' @method head Assay
#' @export
#'
#' @examples
#' # `head` and `tail` can be used to quickly view feature-level meta data
#' head(rna)
#'
head.Assay <- function(x, n = 10L, ...) {
  return(head(x[[]], n = 10L, ...))
}

#' Merge Assays
#'
#' Merge one or more v3 assays together
#'
#' @inheritParams [[.Assay
#' @param y One or more \code{\link{Assay}} objects
#' @param add.cell.ids A character vector of \code{length(x = c(x, y))};
#' appends the corresponding values to the start of each objects' cell names
#' @param merge.data Merge the data slots instead of just merging the counts
#' (which requires renormalization); this is recommended if the same
#' normalization approach was applied to all objects
#' @param labels,collapse Currently unused
#'
#' @return A new assay with data merged from \code{c(x, y)}
#'
#' @method merge Assay
#' @export
#'
#' @family assay
#'
merge.Assay <- function(
  x = NULL,
  y = NULL,
  add.cell.ids = NULL,
  merge.data = TRUE,
  labels = NULL,
  collapse = TRUE,
  ...
) {
  CheckDots(...)
  assays <- c(x, y)
  if (any(sapply(
    X = assays,
    FUN = function(assay.i) inherits(x = assay.i, what = "Assay5")
  ))) {
    return(merge(x = as(x, "Assay5"), y, ...))
  }
  if (!is.null(x = add.cell.ids)) {
    for (i in seq_along(along.with = assays)) {
      assays[[i]] <- RenameCells(
        object = assays[[i]],
        new.names = add.cell.ids[i]
      )
    }
  }
  # Merge the counts (if present)
  counts.mats <- lapply(X = assays, FUN = ValidateDataForMerge, slot = "counts")
  keys <- unlist(sapply(X = assays, FUN = Key))
  merged.counts <- RowMergeSparseMatrices(
    mat1 = counts.mats[[1]],
    mat2 = counts.mats[2:length(x = counts.mats)]
  )
  combined.assay <- CreateAssayObject(
    counts = merged.counts,
    min.cells = -1,
    min.features = -1
  )
  Key(object = combined.assay) <- keys[1]
  if (merge.data) {
    data.mats <- lapply(X = assays, FUN = ValidateDataForMerge, slot = "data")
    merged.data <- RowMergeSparseMatrices(
      mat1 = data.mats[[1]],
      mat2 = data.mats[2:length(x = data.mats)]
    )
    # only keep cells that made it through counts filtering params
    if (!all.equal(target = colnames(x = combined.assay), current = colnames(x = merged.data))) {
      merged.data <- merged.data[, colnames(x = combined.assay)]
    }
    combined.assay <- SetAssayData(
      object = combined.assay,
      layer = "data",
      new.data = merged.data
    )
  }
  return(combined.assay)
}

#' @inherit split.Assay5 title description details
#'
#' @inheritParams split.Assay5
#' @param x An \code{\link{Assay}} object
#'
#' @return Returns a v5 assay with splitted layers
#'
#' @method split Assay
#' @export
#'
#' @family assay
#'
split.Assay <- function(
  x,
  f,
  drop = FALSE,
  layers = NA,
  ...
) {
  warn(message = paste(
    strwrap(x = paste(
      "Input is a v3 assay and `split()` only works for v5 assays;",
      "converting to a v5 assay"
    ))
  ))
  x <- as(object = x, Class = 'Assay5')
  split.x <- split(
    x = x,
    f = f,
    drop = drop,
    layers = layers,
    ...
  )
  return(split.x)
}

#' @inherit subset.Assay5 title description details sections
#'
#' @inheritParams subset.Assay5
#' @param x An \code{\link{Assay}} object
#'
#' @return \code{x} with just the cells and features specified by
#' \code{cells} and \code{features}
#'
#' @importFrom stats na.omit
#'
#' @method subset Assay
#' @export
#'
#' @family assay
#'
#' @examples
#' rna <- pbmc_small[["RNA"]]
#' rna2 <- subset(rna, features = VariableFeatures(rna))
#' rna2
#'
subset.Assay <- function(x, cells = NULL, features = NULL, ...) {
  CheckDots(...)
  cells <- cells %||% colnames(x = x)
  if (all(is.na(x = cells))) {
    cells <- colnames(x = x)
  } else if (any(is.na(x = cells))) {
    warn(message = "NAs passed in cells vector, removing NAs")
    cells <- na.omit(object = cells)
  }
  cells <- intersect(x = colnames(x), y = cells)
  features <- features %||% rownames(x = x)
  if (all(is.na(x = features))) {
    features <- rownames(x = x)
  } else if (any(is.na(x = features))) {
    warn(message = "NAs passed in the features vector, removing NAs")
    features <- na.omit(object = features)
  }
  if (all(sapply(X = list(features, cells), FUN = length) == dim(x = x))) {
    return(x)
  }
  if (is.numeric(x = features)) {
    features <- rownames(x = x)[features]
  }
  features <- gsub(
    pattern = paste0('^', Key(object = x)),
    replacement = '',
    x = features
  )
  features <- intersect(x = features, y = rownames(x = x))
  if (length(x = features) == 0) {
    abort(message = "Cannot find features provided")
  }
  if (ncol(x = GetAssayData(object = x, layer = 'counts')) == ncol(x = x)) {
    slot(object = x, name = "counts") <- GetAssayData(object = x, layer = "counts")[features, cells, drop = FALSE]
  }
  slot(object = x, name = "data") <- GetAssayData(object = x, layer = "data")[features, cells, drop = FALSE]
  cells.scaled <- colnames(x = GetAssayData(object = x, layer = "scale.data"))
  cells.scaled <- cells.scaled[cells.scaled %in% cells]
  cells.scaled <- cells.scaled[na.omit(object = match(x = colnames(x = x), table = cells.scaled))]
  features.scaled <- rownames(x = GetAssayData(object = x, layer = 'scale.data'))
  features.scaled <- intersect(x = features, y = features.scaled)
  slot(object = x, name = "scale.data") <- if (length(x = cells.scaled) > 0 && length(x = features.scaled) > 0) {
    GetAssayData(object = x, layer = "scale.data")[features.scaled, cells.scaled, drop = FALSE]
  } else {
    new(Class = 'matrix')
  }
  VariableFeatures(object = x) <- VariableFeatures(object = x)[VariableFeatures(object = x) %in% features]
  slot(object = x, name = 'meta.features') <- x[[]][features, , drop = FALSE]
  validObject(object = x)
  return(x)
}

#' @rdname sub-sub-.Assay
#'
#' @method tail Assay
#' @export
#'
#' @examples
#' tail(rna)
#'
tail.Assay <- function(x, n = 10L, ...) {
  return(tail(x[[]], n = n, ...))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname sub-.Assay
#'
#' @examples
#' # Set layer data
#' rna["data"] <- rna["counts"]
#' rna["data"][1:10, 1:4]
#'
setMethod(
  f = '[<-',
  signature = c(x = 'Assay', i = 'character'),
  definition = function(x, i, ..., value) {
    LayerData(object = x, layer = i, ...) <- value
    return(x)
  }
)

#' @rdname sub-sub-.Assay
#'
#' @order 2
#'
setMethod(
  f = '[[<-',
  signature = c(x = 'Assay'),
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

#' @rdname sub-sub-.Assay
#'
setMethod(
  f = '[[<-',
  signature = c(
    x = 'Assay',
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

#' Row and Column Sums and Means
#'
#' Calculate \code{\link{rowSums}}, \code{\link{colSums}},
#' \code{\link{rowMeans}}, and \code{\link{colMeans}} on \code{Assay} objects
#'
#' @inheritParams [[.Assay
#' @inheritParams Matrix::colMeans
#' @param slot Name of assay expression matrix to calculate column/row
#' means/sums on
#'
#' @return \code{colMeans}: The column (cell-wise) means of \code{slot}
#'
#' @importFrom Matrix colMeans
#'
#' @keywords internal
#'
#' @export
#'
#' @concept assay
#'
#' @seealso \code{\link{Assay}}
#'
#' @examples
#' rna <- pbmc_small[["RNA"]]
#'
#' colMeans(rna)
#'
setMethod(
  f = 'colMeans',
  signature = c(x = 'Assay'),
  definition = function(x, na.rm = FALSE, dims = 1, ..., slot = 'data') {
    return(Matrix::colMeans(
      x = GetAssayData(object = x, layer = slot),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

#' @return \code{colSums}: The column (cell-wise) sums of \code{slot}
#'
#' @rdname colMeans-Assay-method
#'
#' @importFrom Matrix colSums
#'
#' @export
#'
#' @examples
#' colSums(rna)
#'
setMethod(
  f = 'colSums',
  signature = c(x = 'Assay'),
  definition = function(x, na.rm = FALSE, dims = 1, ..., slot = 'data') {
    return(Matrix::colSums(
      x = GetAssayData(object = x, layer = slot),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

#' @return \code{rowMeans}: The row (feature-wise) means of \code{slot}
#'
#' @rdname colMeans-Assay-method
#'
#' @importFrom Matrix rowMeans
#'
#' @export
#'
#' @examples
#' rowMeans(rna)
#'
setMethod(
  f = 'rowMeans',
  signature = c(x = 'Assay'),
  definition = function(x, na.rm = FALSE, dims = 1, ..., slot = 'data') {
    return(Matrix::rowMeans(
      x = GetAssayData(object = x, layer = slot),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

#' @return \code{rowSums}: The row (feature-wise) sums of \code{slot}
#'
#' @rdname colMeans-Assay-method
#'
#' @importFrom Matrix rowSums
#'
#' @export
#'
#' @examples
#' rowSums(rna)
#'
setMethod(
  f = 'rowSums',
  signature = c(x = 'Assay'),
  definition = function(x, na.rm = FALSE, dims = 1, ..., slot = 'data') {
    return(Matrix::rowSums(
      x = GetAssayData(object = x, layer = slot),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

#' V3 Assay Overview
#'
#' Overview of an \code{\link{Assay}} object
#'
#' @template return-show
#'
#' @keywords internal
#'
#' @concept assay
#'
#' @seealso \code{\link{Assay}}
#'
#' @examples
#' rna <- pbmc_small[["RNA"]]
#' rna
#'
setMethod(
  f = 'show',
  signature = 'Assay',
  definition = function(object) {
    cat(
      class(x = object)[1],
      'data with',
      nrow(x = object),
      'features for',
      ncol(x = object), 'cells\n'
    )
    if (length(x = VariableFeatures(object = object)) > 0) {
      top.ten <- head(x = VariableFeatures(object = object), n = 10L)
      top <- 'Top'
      variable <- 'variable'
    } else {
      top.ten <- head(x = rownames(x = object), n = 10L)
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
    return(invisible(x = NULL))
  }
)

#' V3 Assay Validity
#'
#' @templateVar cls Assay
#' @template desc-validity
#'
#' @section \code{data} Validation:
#' blah
#'
#' @section \code{counts} Validation:
#' blah
#'
#' @section \code{scale.data} Validation:
#' blah
#'
#' @section Feature-Level Meta Data Validation:
#' blah
#'
#' @section Variable Feature Validation:
#' blah
#'
#' @inheritSection Key-validity Key Validation
#'
#' @name Assay-validity
#'
#' @family assay
#' @seealso \code{\link[methods]{validObject}}
#'
#' @examples
#' rna <- pbmc_small[["RNA"]]
#' validObject(rna)
#'
setValidity(
  Class = 'Assay',
  method = function(object) {
    if (.GetSeuratCompat() < '5.0.0') {
      return(TRUE)
    }
    if (isFALSE(x = getOption(x = "Seurat.object.validate", default = TRUE))) {
      warn(
        message = paste("Not validating", class(x = object)[1L], "objects"),
        class = 'validationWarning'
      )
      return(TRUE)
    }
    valid <- NULL
    # Check matrices
    features <- rownames(x = slot(object = object, name = 'data'))
    if (anyDuplicated(x = features)) {
      valid <- c(valid, "duplicate feature names are not allowed")
    }
    cells <- colnames(x = slot(object = object, name = 'data'))
    if (anyDuplicated(x = cells)) {
      valid <- c(valid, "duplicate cell names are not allowed")
    }
    for (lyr in c('counts', 'scale.data')) {
      ldat <- slot(object = object, name = lyr)
      if (IsMatrixEmpty(x = ldat)) {
        next
      }
      if (!all(colnames(x = ldat) == cells)) {
        valid <- c(
          valid,
          paste0("'", lyr, "' must have the same cells as 'data'")
        )
      }
      if (lyr == 'counts' && !all(rownames(x = ldat) == features)) {
        valid <- c(
          valid,
          paste0("'", lyr, "' must have the same features as 'data'")
        )
      } else if (lyr == 'scale.data') {
        scaled <- rownames(x = ldat)
        if (!all(scaled %in% features)) {
          valid <- c(
            valid,
            "all features in 'scale.data' must be present in 'data'"
          )
        } else if (is.unsorted(x = MatchCells(new = scaled, orig = features, ordered = TRUE))) {
          valid <- c(
            valid,
            "features in 'scale.data' must be in the same order as in 'data'"
          )
        }
      }
    }
    # Check meta.features
    mf <- slot(object = object, name = 'meta.features')
    if (nrow(x = mf) != length(x = features)) {
      valid <- c(
        valid,
        "'meta.features' must have the same number of rows as 'data'"
      )
    } else if (!all(row.names(x = mf) == features)) {
      valid <- c(valid, "meta.features' must have the same features as 'data'")
    }
    # Check variable features
    vf <- slot(object = object, name = 'var.features')
    if (length(x = vf) && !all(vf %in% features)) {
      valid <- c(valid, "all 'var.features' must be present in")
    }
    # TODO: Check assay.orig
    return(valid %||% TRUE)
  }
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Calculate nCount and nFeature
#'
#' @param object An \code{\link{Assay}} object
#'
#' @return A named list with nCount and nFeature
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' \donttest{
#' calcn <- SeuratObject:::CalcN(pbmc_small[["RNA"]])
#' head(as.data.frame(calcn))
#' }
#'
CalcN <- .CalcN

#' Subset cells in vst data
#'
#' @param sct.info A vst.out list
#' @param cells vector of cells to retain
#' @param features vector of features to retain
#'
#' @keywords internal
#'
#' @noRd
#'
SubsetVST <- function(sct.info, cells, features) {
  cells.keep <- intersect(x = cells, y = rownames(x = sct.info$cell_attr))
  sct.info$cell_attr <- sct.info$cell_attr[cells.keep, ]
  # find which subset of features are in the SCT assay
  feat.keep <- intersect(x = features, y = rownames(x = sct.info$gene_attr))
  sct.info$gene_attr <- sct.info$gene_attr[feat.keep, ]
  return(sct.info)
}

#' Validate Assay Data for Merge
#'
#' Pulls the proper data matrix for merging assay data. If the slot is empty,
#' will return an empty matrix with the proper dimensions from one of the
#' remaining data slots.
#'
#' @param assay Assay to pull data from
#' @param slot Slot to pull from
#'
#' @return Returns the data matrix if present (i.e.) not 0x0. Otherwise,
#' returns an appropriately sized empty sparse matrix
#'
#' @importFrom methods as
#' @importFrom Matrix Matrix
#'
#' @keywords internal
#'
#' @noRd
#'
ValidateDataForMerge <- function(assay, slot) {
  mat <- GetAssayData(object = assay, layer = slot)
  if (any(dim(x = mat) == c(0, 0))) {
    slots.to.check <- setdiff(x = c("counts", "data", "scale.data"), y = slot)
    for (ss in slots.to.check) {
      data.dims <- dim(x = GetAssayData(object = assay, layer = ss))
      data.slot <- ss
      if (!any(data.dims == c(0, 0))) {
        break
      }
    }
    if (any(data.dims == c(0, 0))) {
      stop("The counts, data, and scale.data slots are all empty for the provided assay.")
    }
    mat <- Matrix(
      data = 0,
      nrow = data.dims[1],
      ncol = data.dims[2],
      dimnames = dimnames(x = GetAssayData(object = assay, layer = data.slot))
    )
    mat <- as.sparse(x = mat)
  }
  return(mat)
}
