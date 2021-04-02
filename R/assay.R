#' @include zzz.R
#' @include generics.R
#' @importFrom methods new setClass setValidity slot slot<-
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @exportClass StdAssay
#'
setClass(
  Class = 'StdAssay',
  contains = 'VIRTUAL',
  slots = c(
    counts = 'ANY',
    data = 'ANY',
    scale.data = 'ANY',
    key = 'character',
    assay.orig = 'OptionalCharacter',
    var.features = 'character',
    meta.features = 'data.frame',
    misc = 'OptionalList'
  )
)

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
  contains = 'StdAssay',
  slots = c(
    counts = 'AnyMatrix',
    data = 'AnyMatrix',
    scale.data = 'matrix'
    # key = 'character',
    # assay.orig = 'OptionalCharacter',
    # var.features = 'vector',
    # meta.features = 'data.frame',
    # misc = 'OptionalList'
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
      counts <- as.sparse(x = counts, ...)
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

#' @rdname AddMetaData
#' @export
#' @method AddMetaData StdAssay
#'
AddMetaData.StdAssay <- .AddMetaData

#' @rdname DefaultAssay
#' @export
#' @method DefaultAssay StdAssay
#'
DefaultAssay.StdAssay <- function(object, ...) {
  object <- UpdateSlots(object = object)
  return(slot(object = object, name = 'assay.orig'))
}

#' @rdname DefaultAssay
#' @export
#' @method DefaultAssay<- StdAssay
#'
"DefaultAssay<-.StdAssay" <- function(object, ..., value) {
  object <- UpdateSlots(object = object)
  slot(object = object, name = 'assay.orig') <- value
  return(object)
}

#' @rdname AssayData
#' @export
#' @method GetAssayData StdAssay
#'
#' @examples
#' # Get the data directly from an Assay object
#' GetAssayData(pbmc_small[["RNA"]], slot = "data")[1:5,1:5]
#'
GetAssayData.StdAssay <- function(
  object,
  slot = c('data', 'scale.data', 'counts'),
  ...
) {
  CheckDots(...)
  slot <- slot[1]
  slot <- match.arg(arg = slot)
  return(slot(object = object, name = slot))
}

#' @rdname VariableFeatures
#' @export
#' @method HVFInfo StdAssay
#'
#' @examples
#' # Get the HVF info directly from an Assay object
#' HVFInfo(pbmc_small[["RNA"]], selection.method = 'vst')[1:5, ]
#'
HVFInfo.StdAssay <- function(object, selection.method, status = FALSE, ...) {
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
#' # Get an Assay key
#' Key(pbmc_small[["RNA"]])
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
#' # Set the key for an Assay
#' Key(pbmc_small[["RNA"]]) <- "newkey_"
#' Key(pbmc_small[["RNA"]])
#'
"Key<-.StdAssay" <- function(object, ..., value) {
  CheckDots(...)
  slot(object = object, name = 'key') <- value
  return(object)
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
#' @method SetAssayData StdAssay
#'
#' @examples
#' # Set an Assay slot directly
#' count.data <- GetAssayData(pbmc_small[["RNA"]], slot = "counts")
#' count.data <- as.matrix(x = count.data + 1)
#' new.assay <- SetAssayData(pbmc_small[["RNA"]], slot = "counts", new.data = count.data)
#'
SetAssayData.StdAssay <- function(
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
#' @method VariableFeatures StdAssay
#'
VariableFeatures.StdAssay <- function(object, selection.method = NULL, ...) {
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

#' @describeIn Assay-methods Get feature-level metadata
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

#' @describeIn Assay-methods Number of cells and features for an \code{Assay}
#'
#' @return \code{dim}: The number of features (\code{nrow}) and cells
#' (\code{ncol})
#'
#' @export
#' @method dim StdAssay
#'
dim.StdAssay <- function(x) {
  return(dim(x = GetAssayData(object = x)))
}

#' @describeIn Assay-methods Cell- and feature-names for an \code{Assay}
#'
#' @return \code{dimnames}: Feature (row) and cell (column) names
#'
#' @export
#' @method dimnames StdAssay
#'
dimnames.StdAssay <- function(x) {
  return(dimnames(x = GetAssayData(object = x)))
}

#' @describeIn Assay-methods Get the first rows of feature-level metadata
#'
#' @return \code{head}: The first \code{n} rows of feature-level metadata
#'
#' @export
#' @method head StdAssay
#'
head.StdAssay <- .head

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
#' @method tail StdAssay
#'
tail.StdAssay <- .tail

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
  signature = 'StdAssay',
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
    # Check counts
    if (!IsMatrixEmpty(x = GetAssayData(object = object, slot = 'counts'))) {
      counts.dims <- dim(x = GetAssayData(object = object, slot = 'counts'))
      if (length(x = counts.dims) != 2) {
        valid <- c(valid, "'counts' must be a two-dimensional object")
      } else if (!isTRUE(x = all.equal(target = counts.dims, data.dims))) {
        valid <- c(valid, "'counts' must be the same size as 'data'")
      }
    }
    # TODO: Check scale data
    # TODO: Check variable features
    # TODO: Check meta features
    # TODO: Check key
    # TODO: Check misc
    # if (!is.null(x = ))
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
