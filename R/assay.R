#' @include zzz.R
#' @include generics.R
#' @importFrom methods new
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The Assay Class
#'
#' The Assay object is the basic unit of Seurat; each Assay stores raw, normalized, and scaled data
#' as well as cluster information, variable features, and any other assay-specific metadata.
#' Assays should contain single cell expression data such as RNA-seq, protein, or imputed expression
#' data.
#'
#' @slot counts Unnormalized data such as raw counts or TPMs
#' @slot data Normalized expression data
#' @slot scale.data Scaled expression data
#' @slot key Key for the Assay
#' @slot assay.orig Original assay that this assay is based off of. Used to track
#' assay provenance
#' @slot var.features Vector of features exhibiting high variance across single cells
#' @slot meta.features Feature-level metadata
#' @slot misc Utility slot for storing additional data associated with the assay
#'
#' @name Assay-class
#' @rdname Assay-class
#' @exportClass Assay
#'
#' @seealso \code{\link{Assay-methods}}
#'
Assay <- setClass(
  Class = 'Assay',
  slots = c(
    counts = 'AnyMatrix',
    data = 'AnyMatrix',
    scale.data = 'matrix',
    key = 'character',
    assay.orig = 'OptionalCharacter',
    var.features = 'vector',
    meta.features = 'data.frame',
    misc = 'ANY'
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
#'
#' @importFrom methods as
#' @importFrom Matrix colSums rowSums
#'
#' @export
#'
#' @examples
#' pbmc_raw <- read.table(
#'   file = system.file('extdata', 'pbmc_raw.txt', package = 'Seurat'),
#'   as.is = TRUE
#' )
#' pbmc_rna <- CreateAssayObject(counts = pbmc_raw)
#' pbmc_rna
#'
CreateAssayObject <- function(
  counts,
  data,
  min.cells = 0,
  min.features = 0
) {
  if (missing(x = counts) && missing(x = data)) {
    stop("Must provide either 'counts' or 'data'")
  } else if (!missing(x = counts) && !missing(x = data)) {
    stop("Either 'counts' or 'data' must be missing; both cannot be provided")
  } else if (!missing(x = counts)) {
    # check that dimnames of input counts are unique
    if (anyDuplicated(rownames(x = counts))) {
      warning(
        "Non-unique features (rownames) present in the input matrix, making unique",
        call. = FALSE,
        immediate. = TRUE
      )
      rownames(x = counts) <- make.unique(names = rownames(x = counts))
    }
    if (anyDuplicated(colnames(x = counts))) {
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
      counts <- as(object = as.matrix(x = counts), Class = 'dgCMatrix')
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
    if (anyDuplicated(rownames(x = data))) {
      warning(
        "Non-unique features (rownames) present in the input matrix, making unique",
        call. = FALSE,
        immediate. = TRUE
      )
      rownames(x = data) <- make.unique(names = rownames(x = data))
    }
    if (anyDuplicated(colnames(x = data))) {
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
    meta.features = init.meta.features
  )
  return(assay)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname AddMetaData
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

#' @param slot Specific information to pull (i.e. counts, data, scale.data, ...)
#'
#' @rdname GetAssayData
#' @export
#' @method GetAssayData Assay
#'
#' @examples
#' # Get the data directly from an Assay object
#' GetAssayData(object = pbmc_small[["RNA"]], slot = "data")[1:5,1:5]
#'
GetAssayData.Assay <- function(
  object,
  slot = c('data', 'scale.data', 'counts'),
  ...
) {
  CheckDots(...)
  slot <- slot[1]
  slot <- match.arg(arg = slot)
  return(slot(object = object, name = slot))
}

#' @param selection.method Which method to pull; choose one from
#' \code{c('sctransform', 'sct')}
#' or \code{c('mean.var.plot', 'dispersion', 'mvp', 'disp')}
#' @param status Add variable status to the resulting data.frame
#'
#' @rdname HVFInfo
#' @export
#' @method HVFInfo Assay
#'
#' @examples
#' # Get the HVF info directly from an Assay object
#' HVFInfo(object = pbmc_small[["RNA"]], selection.method = 'vst')[1:5, ]
#'
HVFInfo.Assay <- function(object, selection.method, status = FALSE, ...) {
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
#' @method Key Assay
#'
#' @examples
#' # Get an Assay key
#' Key(object = pbmc_small[["RNA"]])
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
#' Key(object = pbmc_small[["RNA"]]) <- "newkey_"
#' Key(object = pbmc_small[["RNA"]])
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
#'     object = pbmc_small[["RNA"]],
#'     new.names = paste0("A_", colnames(x = pbmc_small[["RNA"]]))
#' )
#' head(x = colnames(x = renamed.assay))
#'
RenameCells.Assay <- function(object, new.names = NULL, ...) {
  CheckDots(...)
  if (IsSCT(assay = object)) {
    if (is.null(x = Misc(object = object, slot = 'vst.set'))) {
      suppressWarnings(Misc(object = object, slot = "vst.out")$cells_step1 <- new.names)
      suppressWarnings(rownames(x = Misc(object = object, slot = "vst.out")$cell_attr) <- new.names)
    } else{
      suppressWarnings(
        Misc(object, slot = "vst.set") <- lapply(
          X = Misc(object = object, slot = "vst.set"),
          FUN = function(x) {
            new.names.vst <- new.names[which(x = x$cells_step1 %in% Cells(x = object))]
            x$cells_step1 <- new.names.vst
            rownames(x = x$cell_attr) <- new.names.vst
            return(x)
          }
        )
      )
    }
  }
  for (data.slot in c("counts", "data", "scale.data")) {
    old.data <- GetAssayData(object = object, slot = data.slot)
    if (ncol(x = old.data) <= 1) {
      next
    }
    colnames(x = slot(object = object, name = data.slot)) <- new.names
  }
  return(object)
}

#' @param slot Where to store the new data
#' @param new.data New data to insert
#'
#'
#' @importFrom stats na.omit
#'
#' @rdname SetAssayData
#' @export
#' @method SetAssayData Assay
#'
#' @examples
#' # Set an Assay slot directly
#' count.data <- GetAssayData(object = pbmc_small[["RNA"]], slot = "counts")
#' count.data <- as.matrix(x = count.data + 1)
#' new.assay <- SetAssayData(object = pbmc_small[["RNA"]], slot = "counts", new.data = count.data)
#'
SetAssayData.Assay <- function(
  object,
  slot = c('counts', 'data', 'scale.data'),
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

#' @param selection.method Which method to pull. Options: markvariogram, moransi
#' @param status Add variable status to the resulting data.frame
#'
#' @rdname SVFInfo
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
#' @param expression A predicate expression for feature/variable expression, can
#' evalue anything that can be pulled by \code{FetchData}; please note, you may
#' need to wrap feature names in backticks (\code{``}) if dashes between numbers
#' are present in the feature name
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
#' Methods for \code{\link{Assay}} for generics defined in other packages
#'
#' @param object,x An \code{\link{Assay}} object
#' @param i,features For \code{[[}: metadata names; for all other methods,
#' feature names or indices
#' @param j,cells Cell names or indices
#' @param ... Arguments passed to other methods
#'
#' @name Assay-methods
#' @rdname Assay-methods
#'
NULL

#' @describeIn Assay-methods Data accessor
#'
#' @return \code{[}: The \code{data} slot for features \code{i} and cells
#' \code{j}
#'
#' @export
#' @method [ Assay
#'
"[.Assay" <- function(x, i, j, ...) {
  if (missing(x = i)) {
    i <- 1:nrow(x = x)
  }
  if (missing(x = j)) {
    j <- 1:ncol(x = x)
  }
  return(GetAssayData(object = x)[i, j, ..., drop = FALSE])
}

#' @describeIn Assay-methods Metadata accessor
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

#' @describeIn Assay-methods Number of cells and features
#'
#' @return \code{dim}: The number of features (\code{nrow}) and cells
#' \code{ncol}
#'
#' @export
#' @method dim Assay
#'
dim.Assay <- function(x) {
  return(dim(x = GetAssayData(object = x)))
}

#' @describeIn Assay-methods Cell- and feature-names
#'
#' @return \code{dimnames}: Feature (row) and cell (column) names
#'
#' @export
#' @method dimnames Assay
#'
dimnames.Assay <- function(x) {
  return(dimnames(x = GetAssayData(object = x)))
}

#' @describeIn Assay-methods Merge objects
#'
#' @param y A vector or list of one or more objects to merge
#' @param add.cell.ids A character vector of \code{length(x = c(x, y))}; appends
#' the corresponding values to the start of each objects' cell names
#' @param merge.data Merge the data slots instead of just merging the counts
#' (which requires renormalization); this is recommended if the same normalization
#' approach was applied to all objects
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
    for (i in 1:length(assays)) {
      assays[[i]] <- RenameCells(object = assays[[i]], new.names = add.cell.ids[i])
    }
  }
  # Merge the counts (if present)
  merged.counts <- ValidateDataForMerge(assay = assays[[1]], slot = "counts")
  keys <- Key(object = assays[[1]])
  for (i in 2:length(x = assays)) {
    merged.counts <- RowMergeSparseMatrices(
      mat1 = merged.counts,
      mat2 = ValidateDataForMerge(assay = assays[[i]], slot = "counts")
    )
    if (length(Key(object = assays[[i]]) > 0)) {
      keys[i] <- Key(object = assays[[i]])
    }
  }
  combined.assay <- CreateAssayObject(
    counts = merged.counts,
    min.cells = -1,
    min.features = -1
  )
  if (length(x = unique(x = keys)) == 1) {
    Key(object = combined.assay) <- keys[1]
  }
  if (merge.data) {
    merged.data <- ValidateDataForMerge(assay = assays[[1]], slot = "data")
    for (i in 2:length(x = assays)) {
      merged.data <- RowMergeSparseMatrices(
        mat1 = merged.data,
        mat2 = ValidateDataForMerge(assay = assays[[i]], slot = "data")
      )
    }
    # only keep cells that made it through counts filtering params
    merged.data <- merged.data[, colnames(x = combined.assay)]
    combined.assay <- SetAssayData(
      object = combined.assay,
      slot = "data",
      new.data = merged.data
    )
  }
  # merge SCT assay misc vst info and scale.data
  if (all(IsSCT(assay = assays))) {
    vst.set.new <- list()
    idx <- 1
    umi.assay.new <- list()
    for (i in 1:length(x = assays)) {
      vst.set.old <- Misc(object = assays[[i]], slot = "vst.set")
      umi.assay.old <- Misc(object = assays[[i]], slot = "umi.assay")
      if (!is.null(x = vst.set.old)) {
        for (j in 1:length(x = vst.set.old)) {
          vst.set.new[[idx]] <- vst.set.old[[j]]
          umi.assay.new[[idx]] <- umi.assay.old[[j]]
          idx <- idx + 1
        }
      } else if (!is.null(x = Misc(object = assays[[i]], slot = "vst.out"))) {
        vst.set.new[[idx]] <- Misc(object = assays[[i]], slot = "vst.out")
        umi.assay.new[[idx]] <- Misc(object = assays[[i]], slot = "umi.assay")
        idx <- idx + 1
      }
    }
    Misc(object = combined.assay, slot = "vst.set") <- vst.set.new
    Misc(object = combined.assay, slot = "umi.assay") <- umi.assay.new
    scale.data <- do.call(
      what = cbind,
      args = lapply(X = assays, FUN = function(x) GetAssayData(object = x, slot = "scale.data"))
    )
    combined.assay <- SetAssayData(
      object = combined.assay,
      slot = "scale.data",
      new.data = scale.data
    )
  }
  return(combined.assay)
}

#' @describeIn Assay-methods Subset object
#'
#' @return \code{subset}: A subsetted object
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
    slot(object = x, name = "counts") <- GetAssayData(
      object = x,
      slot = "counts"
    )[features, cells, drop = FALSE]
  }
  slot(object = x, name = "data") <- GetAssayData(
    object = x,
    slot = "data"
  )[features, cells, drop = FALSE]
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
  if (IsSCT(assay = x)) {
    # subset cells and genes in the SCT assay
    obj.misc <- Misc(object = x)
    if ("vst.set" %in% names(x = obj.misc)) {
      # set of vst.out objects
      vst.info <- obj.misc[["vst.set"]]
      for (i in seq_along(along.with = vst.info)) {
        vst.info[[i]] <- SubsetVST(
          sct.info = vst.info[[i]],
          cells = cells,
          features = features
        )
      }
      obj.misc[["vst.set"]] <- vst.info
    } else {
      # just one vst.out
      obj.misc[["vst.out"]] <- SubsetVST(
        sct.info = obj.misc[["vst.out"]],
        cells = cells,
        features = features
      )
    }
    slot(object = x, name = "misc") <- obj.misc
  }
  return(x)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname AddMetaData
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

#' @importFrom Matrix colMeans
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

#' @importFrom Matrix colSums
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

#' @importFrom Matrix rowMeans
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

#' @importFrom Matrix rowSums
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

#' @importFrom utils head
#'
setMethod(
  f = 'show',
  signature = 'Assay',
  definition = function(object) {
    cat('Assay data with', nrow(x = object), 'features for', ncol(x = object), 'cells\n')
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
      if (length(x = top.ten) != 1) {'s'}, ":\n"
    )
    features <- gsub(pattern = '^\\s+', replacement = '', x = features)
    cat(
      top,
      length(x = top.ten),
      features,
      paste(strwrap(x = paste(top.ten, collapse = ', ')), collapse = '\n'),
      '\n'
    )
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
CalcN <- function(object) {
  if (IsMatrixEmpty(x = GetAssayData(object = object, slot = "counts"))) {
    return(NULL)
  }
  return(list(
    nCount = colSums(x = object, slot = 'counts'),
    nFeature = colSums(x = GetAssayData(object = object, slot = 'counts') > 0)
  ))
}

#' Subset cells in vst data
#'
#' @param sct.info A vst.out list
#' @param cells vector of cells to retain
#' @param features vector of features to retain
#'
#' @keywords internal
SubsetVST <- function(sct.info, cells, features) {
  cells.keep <- intersect(x = cells, y = rownames(x = sct.info$cell_attr))
  sct.info$cell_attr <- sct.info$cell_attr[cells.keep, ]
  # find which subset of features are in the SCT assay
  feat.keep <- intersect(x = features, y = rownames(x = sct.info$gene_attr))
  sct.info$gene_attr <- sct.info$gene_attr[feat.keep, ]
  return(sct.info)
}
