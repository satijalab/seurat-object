#' @include zzz.R
#' @include generics.R
#' @include default.R
#' @include graph.R
#' @importFrom methods new setClass setValidity slot slot<-
#' @importFrom methods new slot slot<-
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
#' @slot key Key for the Assay
#' @slot assay.orig Original assay that this assay is based off of. Used to
#' track assay provenance
#' @slot var.features Vector of features exhibiting high variance across
#' single cells
#' @slot meta.features Feature-level metadata
#' @slot misc Utility slot for storing additional data associated with the assay
#'
#' @name Assay-class
#' @rdname Assay-class
#' @exportClass Assay
#'
#' @concept assay
#'
#' @seealso \code{\link{Assay-methods}}
#'
Assay <- setClass(
  Class = 'Assay',
  # contains = 'StdAssay',
  slots = c(
    counts = 'AnyMatrix',
    data = 'AnyMatrix',
    scale.data = 'matrix',
    key = 'character',
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
#' new object with a lower cutoff.
#' @param min.features Include cells where at least this many features are
#' detected.
#' @param check.matrix Check counts matrix for NA, NaN, Inf, and non-integer values
#' @param ... Arguments passed to \code{\link{as.sparse}}
#'
#' @return A \code{\link{Assay}} object
#'
#' @importFrom methods as
#' @importFrom Matrix colSums rowSums
#'
#' @export
#'
#' @concept assay
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
  check.matrix = FALSE,
  ...
) {
  if (missing(x = counts) && missing(x = data)) {
    stop("Must provide either 'counts' or 'data'")
  } else if (!missing(x = counts) && !missing(x = data)) {
    stop("Either 'counts' or 'data' must be missing; both cannot be provided")
  } else if (!missing(x = counts)) {
    # check that dimnames of input counts are unique
    if (anyDuplicated(x = rownames(x = counts))) {
      warning(
        "Non-unique features (rownames) present in the input matrix, making unique",
        call. = FALSE,
        immediate. = TRUE
      )
      rownames(x = counts) <- make.unique(names = rownames(x = counts))
    }
    if (anyDuplicated(x = colnames(x = counts))) {
      warning(
        "Non-unique cell names (colnames) present in the input matrix, making unique",
        call. = FALSE,
        immediate. = TRUE
      )
      colnames(x = counts) <- make.unique(names = colnames(x = counts))
    }
    if (is.null(x = colnames(x = counts))) {
      stop("No cell names (colnames) names present in the input matrix")
    }
    if (any(rownames(x = counts) == '')) {
      stop("Feature names of counts matrix cannot be empty", call. = FALSE)
    }
    if (nrow(x = counts) > 0 && is.null(x = rownames(x = counts))) {
      stop("No feature names (rownames) names present in the input matrix")
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
      warning(
        "Non-unique features (rownames) present in the input matrix, making unique",
        call. = FALSE,
        immediate. = TRUE
      )
      rownames(x = data) <- make.unique(names = rownames(x = data))
    }
    if (anyDuplicated(x = colnames(x = data))) {
      warning(
        "Non-unique cell names (colnames) present in the input matrix, making unique",
        call. = FALSE,
        immediate. = TRUE
      )
      colnames(x = data) <- make.unique(names = colnames(x = data))
    }
    if (is.null(x = colnames(x = data))) {
      stop("No cell names (colnames) names present in the input matrix")
    }
    if (any(rownames(x = data) == '')) {
      stop("Feature names of data matrix cannot be empty", call. = FALSE)
    }
    if (nrow(x = data) > 0 && is.null(x = rownames(x = data))) {
      stop("No feature names (rownames) names present in the input matrix")
    }
    if (min.cells != 0 | min.features != 0) {
      warning(
        "No filtering performed if passing to data rather than counts",
        call. = FALSE,
        immediate. = TRUE
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
  if (any(grepl(pattern = '_', x = rownames(x = counts))) || any(grepl(pattern = '_', x = rownames(x = data)))) {
    warning(
      "Feature names cannot have underscores ('_'), replacing with dashes ('-')",
      call. = FALSE,
      immediate. = TRUE
    )
    rownames(x = counts) <- gsub(
      pattern = '_',
      replacement = '-',
      x = rownames(x = counts)
    )
    rownames(x = data) <- gsub(
      pattern = '_',
      replacement = '-',
      x = rownames(x = data)
    )
  }
  if (any(grepl(pattern = '|', x = rownames(x = counts), fixed = TRUE)) || any(grepl(pattern = '|', x = rownames(x = data), fixed = TRUE))) {
    warning(
      "Feature names cannot have pipe characters ('|'), replacing with dashes ('-')",
      call. = FALSE,
      immediate. = TRUE
    )
    rownames(x = counts) <- gsub(
      pattern = '|',
      replacement = '-',
      x = rownames(x = counts),
      fixed = TRUE
    )
    rownames(x = data) <- gsub(
      pattern = '|',
      replacement = '-',
      x = rownames(x = data),
      fixed = TRUE
    )
  }
  # Initialize meta.features
  init.meta.features <- data.frame(row.names = rownames(x = data))
  assay <- new(
    Class = 'Assay',
    counts = counts,
    data = data,
    scale.data = new(Class = 'matrix'),
    meta.features = init.meta.features,
    misc = list()
  )
  return(assay)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# @rdname AddMetaData
#' @inherit AddMetaData
#'
#' @templateVar fname AddMetaData
#' @templateVar version 4
#' @template name-oldv
#'
#' @export
#' @method AddMetaData Assay
#'
AddMetaData.Assay <- .AddMetaData

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
    f <- if (.IsFutureSeurat(version = '5.1.0')) {
      deprecate_stop
    } else if (.IsFutureSeurat(version = '5.0.0')) {
      deprecate_warn
    } else {
      deprecate_soft
    }
    f(
      when = '5.0.0',
      what = 'Features(slot = )',
      with = 'Features(layer = )'
    )
    layer <- slot
  }
  layer <- layer[1L] %||% 'data'
  layer <- match.arg(arg = layer)
  features <- rownames(x = GetAssayData(object = x, slot = layer))
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
    deprecate_soft(
      when = '4.9.0',
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
  # Check vars
  orig <- vars
  vars <- gsub(
    pattern = paste0('^', Key(object = object)),
    replacement = '',
    x = vars
  )
  # Pull expression information
  mat <- GetAssayData(object = object, slot = layer)
  if (IsMatrixEmpty(x = mat)) {
    stop("Slot ", layer, " is empty in this assay", call. = FALSE)
  }
  vars <- intersect(x = vars, y = rownames(x = mat))
  data.fetched <- as.data.frame(x = as.matrix(
    x = Matrix::t(x = mat[vars, cells, drop = FALSE])
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
    stop("None of the requested features found", call. = FALSE)
  } else if (length(x = missing)) {
    warning(
      "The following features could not be found ",
      paste(missing, collapse = ', '),
      call. = FALSE,
      immediate. = TRUE
    )
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
    f <- if (.IsFutureSeurat(version = '5.1.0')) {
      deprecate_stop
    } else if (.IsFutureSeurat(version = '5.0.0')) {
      deprecate_warn
    } else {
      deprecate_soft
    }
    f(
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
    f <- if (.IsFutureSeurat(version = '5.1.0')) {
      deprecate_stop
    } else if (.IsFutureSeurat(version = '5.0.0')) {
      deprecate_warn
    } else {
      deprecate_soft
    }
    f(
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
  selection.method <- switch(
    EXPR = tolower(x = method),
    'sctransform' = 'sct',
    method
  )
  vars <- switch(
    EXPR = method,
    'vst' = c('mean', 'variance', 'variance.standardized'),
    'mvp' = c('mean', 'dispersion', 'dispersion.scaled'),
    'sct' = c('gmean', 'variance', 'residual_variance'),
    stop("Unknown method: '", method, "'", call. = FALSE)
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
  for (data.slot in c("counts", "data", "scale.data")) {
    old.data <- GetAssayData(object = object, slot = data.slot)
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
#' # Set an Assay slot directly
#' count.data <- GetAssayData(pbmc_small[["RNA"]], slot = "counts")
#' count.data <- as.matrix(x = count.data + 1)
#' new.assay <- SetAssayData(pbmc_small[["RNA"]], slot = "counts", new.data = count.data)
#'
SetAssayData.Assay <- function(
  object,
  slot = c('data', 'scale.data', 'counts'),
  new.data,
  ...
) {
  CheckDots(...)
  slot <- slot[1]
  slot <- match.arg(arg = slot)
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

#' @param decreasing Return features in decreasing order (most spatially
#' variable first).
#'
#' @rdname VariableFeatures
#' @export
#' @method SpatiallyVariableFeatures Assay
#'
SpatiallyVariableFeatures.Assay <- function(
  object,
  selection.method = "markvariogram",
  decreasing = TRUE,
  ...
) {
  CheckDots(...)
  vf <- SVFInfo(object = object, selection.method = selection.method, status = TRUE)
  vf <- vf[rownames(x = vf)[which(x = vf[, "variable"][, 1])], ]
  if (!is.null(x = decreasing)) {
    vf <- vf[order(x = vf[, "rank"], decreasing = !decreasing), ]
  }
  return(rownames(x = vf)[which(x = vf[, "variable"][, 1])])
}

#' @rdname VariableFeatures
#' @export
#' @method SVFInfo Assay
#'
SVFInfo.Assay <- function(
  object,
  selection.method = c("markvariogram", "moransi"),
  status = FALSE,
  ...
) {
  CheckDots(...)
  selection.method <- selection.method[1]
  selection.method <- match.arg(arg = selection.method)
  vars <- switch(
    EXPR = selection.method,
    'markvariogram' = grep(
      pattern = "r.metric",
      x = colnames(x = object[[]]),
      value = TRUE
    ),
    'moransi' = grep(
      pattern = 'moransi',
      x = colnames(x = object[[]]),
      value = TRUE
    ),
    stop("Unknown method: '", selection.method, "'", call. = FALSE)
  )
  tryCatch(
    expr = svf.info <- object[[vars]],
    error = function(e) {
      stop(
        "Unable to find highly variable feature information for method '",
        selection.method,
        "'",
        call. = FALSE
      )
    }
  )
  colnames(x = svf.info) <- vars
  if (status) {
    svf.info$variable <- object[[paste0(selection.method, '.spatially.variable')]]
    svf.info$rank <- object[[paste0(selection.method, '.spatially.variable.rank')]]
  }
  return(svf.info)
}

#' @rdname VariableFeatures
#' @export
#' @method VariableFeatures Assay
#'
VariableFeatures.Assay <- function(object, selection.method = NULL, ...) {
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
#' @method VariableFeatures<- Assay
#'
"VariableFeatures<-.Assay" <- function(object, ..., value) {
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

#' @param cells Subset of cell names
#' @param expression A predicate expression for feature/variable expression,
#' can evaluate anything that can be pulled by \code{FetchData}; please note,
#' you may need to wrap feature names in backticks (\code{``}) if dashes
#' between numbers are present in the feature name
#' @param invert Invert the selection of cells
#'
#' @importFrom stats na.omit
#' @importFrom rlang is_quosure enquo eval_tidy
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
    data.subset <- as.data.frame(x = t(x = as.matrix(x = object[expr.char, ])))
    colnames(x = data.subset) <- expr.char
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

#' \code{Assay} Methods
#'
#' Methods for \code{\link{Assay}} objects for generics defined in
#' other packages
#'
#' @param x,object An \code{\link{Assay}} object
#' @param i,features For \code{[[}: metadata names; for all other methods,
#' feature names or indices
#' @param j,cells Cell names or indices
#' @param ... Arguments passed to other methods
#'
#' @name Assay-methods
#' @rdname Assay-methods
#'
#' @concept assay
#'
NULL

#' @describeIn Assay-methods Get expression data from an \code{Assay}
#'
#' @return \code{[}: The \code{data} slot for features \code{i} and cells
#' \code{j}
#'
#' @export
#' @method [ Assay
#'
"[.Assay" <- function(x, i, j, ...) {
  if (missing(x = i)) {
    i <- seq_len(length.out = nrow(x = x))
  }
  if (missing(x = j)) {
    j <- seq_len(length.out = ncol(x = x))
  }
  return(GetAssayData(object = x)[i, j, ..., drop = FALSE])
}

#' @describeIn Assay-methods Get feature-level metadata
#'
#' @param drop See \code{\link[base]{drop}}
#'
#' @return \code{[[}: The feature-level metadata for \code{i}
#'
#' @export
#' @method [[ Assay
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

#' @describeIn Assay-methods Number of cells and features for an \code{Assay}
#'
#' @return \code{dim}: The number of features (\code{nrow}) and cells
#' (\code{ncol})
#'
#' @export
#' @method dim Assay
#'
dim.Assay <- function(x) {
  return(dim(x = GetAssayData(object = x)))
}

#' @describeIn Assay-methods Cell- and feature-names for an \code{Assay}
#'
#' @return \code{dimnames}: Feature (row) and cell (column) names
#'
#' @export
#' @method dimnames Assay
#'
dimnames.Assay <- function(x) {
  return(dimnames(x = GetAssayData(object = x)))
}

#' @method dimnames<- Assay
#' @export
#'
"dimnames<-.Assay" <- function(x, value) {
  op <- options(Seurat.object.validate = FALSE)
  on.exit(expr = options(op), add = TRUE)
  # Check the provided dimnames
  msg <- "Invalid 'dimnames' given for a Seurat object"
  if (!is_bare_list(x = value, n = 2L)) {
    stop(msg, call. = FALSE)
  } else if (!all(sapply(X = value, FUN = length) == dim(x = x))) {
    stop(msg, call. = FALSE)
  }
  value <- lapply(X = value, FUN = as.character)
  # Warn about changing features
  if (!all(value[[1L]] == rownames(x = slot(object = x, name = 'data')))) {
    warning(
      "Changing feature names in v3 Assays is not supported",
      call. = FALSE,
      immediate. = TRUE
    )
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

#' @describeIn Assay-methods Get the first rows of feature-level metadata
#'
#' @inheritParams utils::head
#'
#' @return \code{head}: The first \code{n} rows of feature-level metadata
#'
#' @export
#' @method head Assay
#'
head.Assay <- .head

#' @describeIn Assay-methods Merge \code{Assay} objects
#'
#' @param y A vector or list of one or more objects to merge
#' @param add.cell.ids A character vector of \code{length(x = c(x, y))};
#' appends the corresponding values to the start of each objects' cell names
#' @param merge.data Merge the data slots instead of just merging the counts
#' (which requires renormalization); this is recommended if the same
#' normalization approach was applied to all objects
#'
#' @return \code{merge}: Merged object
#'
#' @export
#' @method merge Assay
#'
merge.Assay <- function(
  x = NULL,
  y = NULL,
  add.cell.ids = NULL,
  merge.data = TRUE,
  ...
) {
  CheckDots(...)
  assays <- c(x, y)
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
  keys <- sapply(X = assays, FUN = Key)
  merged.counts <- RowMergeSparseMatrices(
    mat1 = counts.mats[[1]],
    mat2 = counts.mats[2:length(x = counts.mats)]
  )
  combined.assay <- CreateAssayObject(
    counts = merged.counts,
    min.cells = -1,
    min.features = -1
  )
  if (length(x = unique(x = keys)) == 1) {
    Key(object = combined.assay) <- keys[1]
  }
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
      slot = "data",
      new.data = merged.data
    )
  }
  return(combined.assay)
}

#' @describeIn Assay-methods Subset an \code{Assay}
#'
#' @return \code{subset}: A subsetted \code{Assay}
#'
#' @importFrom stats na.omit
#'
#' @export
#' @method subset Assay
#'
subset.Assay <- function(x, cells = NULL, features = NULL, ...) {
  CheckDots(...)
  cells <- cells %||% colnames(x = x)
  if (all(is.na(x = cells))) {
    cells <- colnames(x = x)
  } else if (any(is.na(x = cells))) {
    warning("NAs passed in cells vector, removing NAs")
    cells <- na.omit(object = cells)
  }
  features <- features %||% rownames(x = x)
  if (all(is.na(x = features))) {
    features <- rownames(x = x)
  } else if (any(is.na(x = features))) {
    warning("NAs passed in the features vector, removing NAs")
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
    stop("Cannot find features provided")
  }
  if (ncol(x = GetAssayData(object = x, slot = 'counts')) == ncol(x = x)) {
    slot(object = x, name = "counts") <- GetAssayData(object = x, slot = "counts")[features, cells, drop = FALSE]
  }
  slot(object = x, name = "data") <- GetAssayData(object = x, slot = "data")[features, cells, drop = FALSE]
  cells.scaled <- colnames(x = GetAssayData(object = x, slot = "scale.data"))
  cells.scaled <- cells.scaled[cells.scaled %in% cells]
  cells.scaled <- cells.scaled[na.omit(object = match(x = colnames(x = x), table = cells.scaled))]
  features.scaled <- rownames(x = GetAssayData(object = x, slot = 'scale.data'))
  features.scaled <- features.scaled[features.scaled %in% features]
  slot(object = x, name = "scale.data") <- if (length(x = cells.scaled) > 0 && length(x = features.scaled) > 0) {
    GetAssayData(object = x, slot = "scale.data")[features.scaled, cells.scaled, drop = FALSE]
  } else {
    new(Class = 'matrix')
  }
  VariableFeatures(object = x) <- VariableFeatures(object = x)[VariableFeatures(object = x) %in% features]
  slot(object = x, name = 'meta.features') <- x[[]][features, , drop = FALSE]
  return(x)
}

#' @describeIn Assay-methods Get the last rows of feature-level metadata
#'
#' @return \code{tail}: The last \code{n} rows of feature-level metadata
#'
#' @importFrom utils tail
#'
#' @export
#' @method tail Assay
#'
tail.Assay <- .tail

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @describeIn Assay-methods Add feature-level metadata
#'
#' @param value Additional metadata to add
#'
#' @return \code{[[<-}: \code{x} with metadata \code{value} added as \code{i}
#'
#' @export
#'
setMethod(
  f = '[[<-',
  signature = c('x' = 'Assay'),
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

#' @describeIn Assay-methods Calculate \code{\link[base]{colMeans}} on an
#' \code{Assay}
#'
#' @param slot Name of assay expression matrix to calculate column/row
#' means/sums on
#' @inheritParams Matrix::colMeans
#'
#' @return \code{colMeans}: The column (cell-wise) means of \code{slot}
#'
#' @importFrom Matrix colMeans
#'
#' @export
#'
setMethod(
  f = 'colMeans',
  signature = c('x' = 'Assay'),
  definition = function(x, na.rm = FALSE, dims = 1, ..., slot = 'data') {
    return(Matrix::colMeans(
      x = GetAssayData(object = x, slot = slot),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

#' @describeIn Assay-methods Calculate \code{\link[base]{colSums}} on an
#' \code{Assay}
#'
#' @return \code{colSums}: The column (cell-wise) sums of \code{slot}
#'
#' @importFrom Matrix colSums
#'
#' @export
#'
setMethod(
  f = 'colSums',
  signature = c('x' = 'Assay'),
  definition = function(x, na.rm = FALSE, dims = 1, ..., slot = 'data') {
    return(Matrix::colSums(
      x = GetAssayData(object = x, slot = slot),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

#' @describeIn Assay-methods Calculate \code{\link[base]{rowMeans}} on an
#' \code{Assay}
#'
#' @return \code{rowMeans}: The row (feature-wise) means of \code{slot}
#'
#' @importFrom Matrix rowMeans
#'
#' @export
#'
setMethod(
  f = 'rowMeans',
  signature = c('x' = 'Assay'),
  definition = function(x, na.rm = FALSE, dims = 1, ..., slot = 'data') {
    return(Matrix::rowMeans(
      x = GetAssayData(object = x, slot = slot),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

#' @describeIn Assay-methods Calculate \code{\link[base]{rowSums}} on an
#' \code{Assay}
#'
#' @return \code{rowSums}: The row (feature-wise) sums of \code{slot}
#'
#' @importFrom Matrix rowSums
#'
#' @export
#'
setMethod(
  f = 'rowSums',
  signature = c('x' = 'Assay'),
  definition = function(x, na.rm = FALSE, dims = 1, ..., slot = 'data') {
    return(Matrix::rowSums(
      x = GetAssayData(object = x, slot = slot),
      na.rm = na.rm,
      dims = dims,
      ...
    ))
  }
)

#' @describeIn Assay-methods Overview of an \code{Assay} object
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

setValidity(
  Class = 'Assay',
  method = function(object) {
    if (isFALSE(x = getOption(x = "Seurat.object.validate", default = TRUE))) {
      warning(
        "Not validating",
        class(x = object)[1L],
        "objects",
        call. = FALSE,
        immediate. = TRUE
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
    # TODO: Check key
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
#' @importFrom Matrix colSums
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
CalcN <- function(object) {
  if (IsMatrixEmpty(x = GetAssayData(object = object, slot = "counts"))) {
    return(NULL)
  }
  return(list(
    nCount = Matrix::colSums(x = object, slot = 'counts'),
    nFeature = Matrix::colSums(x = GetAssayData(object = object, slot = 'counts') > 0)
  ))
}

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
  mat <- GetAssayData(object = assay, slot = slot)
  if (any(dim(x = mat) == c(0, 0))) {
    slots.to.check <- setdiff(x = c("counts", "data", "scale.data"), y = slot)
    for (ss in slots.to.check) {
      data.dims <- dim(x = GetAssayData(object = assay, slot = ss))
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
      dimnames = dimnames(x = GetAssayData(object = assay, slot = data.slot))
    )
    mat <- as(object = mat, Class = "dgCMatrix")
  }
  return(mat)
}
