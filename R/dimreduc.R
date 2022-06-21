#' @include zzz.R
#' @include generics.R
#' @include jackstraw.R
#' @include keymixin.R
#' @importFrom methods new slot slot<- slotNames
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The Dimensional Reduction Class
#'
#' The DimReduc object stores a dimensionality reduction taken out in Seurat;
#' each DimReduc consists of a cell embeddings matrix, a feature loadings
#' matrix, and a projected feature loadings matrix.
#'
#' @slot cell.embeddings Cell embeddings matrix (required)
#' @slot feature.loadings Feature loadings matrix (optional)
#' @slot feature.loadings.projected Projected feature loadings matrix (optional)
#' @slot assay.used Name of assay used to generate \code{DimReduc} object
#' @slot global Is this \code{DimReduc} global/persistent? If so, it will not be
#' removed when removing its associated assay
#' @slot stdev A vector of standard deviations
#' @slot key Key for the \code{DimReduc}, must be alphanumeric characters
#' followed by an underscore
#' @slot jackstraw A \code{\link{JackStrawData-class}} object associated with
#' this \code{DimReduc}
#' @slot misc Utility slot for storing additional data associated with the
#' \code{DimReduc} (e.g. the total variance of the PCA)
#'
#' @name DimReduc-class
#' @rdname DimReduc-class
#' @exportClass DimReduc
#'
DimReduc <- setClass(
  Class = 'DimReduc',
  contains = 'KeyMixin',
  slots = c(
    cell.embeddings = 'matrix',
    feature.loadings = 'matrix',
    feature.loadings.projected = 'matrix',
    assay.used = 'character',
    global = 'logical',
    stdev = 'numeric',
    # key = 'character',
    jackstraw = 'JackStrawData',
    misc = 'list'
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Create a DimReduc object
#'
#' @param embeddings A matrix with the cell embeddings
#' @param loadings A matrix with the feature loadings
#' @param projected A matrix with the projected feature loadings
#' @param assay Assay used to calculate this dimensional reduction
#' @param stdev Standard deviation (if applicable) for the dimensional reduction
#' @param key A character string to facilitate looking up features from a
#' specific DimReduc
#' @param global Specify this as a global reduction (useful for visualizations)
#' @param jackstraw Results from the JackStraw function
#' @param misc list for the user to store any additional information associated
#' with the dimensional reduction
#'
#' @return A \code{\link{DimReduc}} object
#'
#' @aliases SetDimReduction
#'
#' @export
#'
#' @concept dimreduc
#'
#' @examples
#' data <- GetAssayData(pbmc_small[["RNA"]], slot = "scale.data")
#' pcs <- prcomp(x = data)
#' pca.dr <- CreateDimReducObject(
#'   embeddings = pcs$rotation,
#'   loadings = pcs$x,
#'   stdev = pcs$sdev,
#'   key = "PC",
#'   assay = "RNA"
#' )
#'
CreateDimReducObject <- function(
  embeddings = new(Class = 'matrix'),
  loadings = new(Class = 'matrix'),
  projected = new(Class = 'matrix'),
  assay = NULL,
  stdev = numeric(),
  key = NULL,
  global = FALSE,
  jackstraw = NULL,
  misc = list()
) {
  if (is.null(x = assay)) {
    warning(
      "No assay specified, setting assay as RNA by default.",
      call. = FALSE,
      immediate. = TRUE
    )
    assay <- "RNA"
  }
  # Try to infer key from column names
  if (is.null(x = key) && is.null(x = colnames(x = embeddings))) {
    stop("Please specify a key for the DimReduc object")
  } else if (is.null(x = key)) {
    key <- regmatches(
      x = colnames(x = embeddings),
      m = regexec(pattern = '^[[:alnum:]]+_', text = colnames(x = embeddings))
    )
    key <- unique(x = unlist(x = key, use.names = FALSE))
  }
  if (length(x = key) != 1) {
    stop("Please specify a key for the DimReduc object")
  } else if (!grepl(pattern = '^[[:alnum:]]+_$', x = key)) {
    old.key  <- key
    key <- UpdateKey(key = old.key)
    colnames(x = embeddings) <- gsub(
      x = colnames(x = embeddings),
      pattern = old.key,
      replacement = key
    )
    warning(
      "All keys should be one or more alphanumeric characters followed by an underscore '_', setting key to ",
      key,
      call. = FALSE,
      immediate. = TRUE
    )
  }
  # ensure colnames of the embeddings are the key followed by a numeric
  if (is.null(x = colnames(x = embeddings))) {
    warning(
      "No columnames present in cell embeddings, setting to '",
      key,
      "1:",
      ncol(x = embeddings),
      "'",
      call. = FALSE,
      immediate. = TRUE
    )
    colnames(x = embeddings) <- paste0(key, 1:ncol(x = embeddings))
  } else if (!all(grepl(pattern = paste0('^', key, "[[:digit:]]+$"), x = colnames(x = embeddings)))) {
    digits <- unlist(x = regmatches(
      x = colnames(x = embeddings),
      m = regexec(pattern = '[[:digit:]]+$', text = colnames(x = embeddings))
    ))
    if (length(x = digits) != ncol(x = embeddings)) {
      stop("Please ensure all column names in the embeddings matrix are the key plus a digit representing a dimension number")
    }
    colnames(x = embeddings) <- paste0(key, digits)
  }
  if (!IsMatrixEmpty(x = loadings)) {
    if (any(rownames(x = loadings) == '')) {
      stop("Feature names of loadings matrix cannot be empty", call. = FALSE)
    }
    colnames(x = loadings) <- colnames(x = embeddings)
  }
  if (!IsMatrixEmpty(x = projected)) {
    if (any(rownames(x = loadings) == '')) {
      stop("Feature names of projected loadings matrix cannot be empty", call. = FALSE)
    }
    colnames(x = projected) <- colnames(x = embeddings)
  }
  jackstraw <- jackstraw %||% new(Class = 'JackStrawData')
  dim.reduc <- new(
    Class = 'DimReduc',
    cell.embeddings = embeddings,
    feature.loadings = loadings,
    feature.loadings.projected = projected,
    assay.used = assay,
    global = global,
    stdev = stdev,
    key = key,
    jackstraw = jackstraw,
    misc = misc
  )
  return(dim.reduc)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname Cells
#' @export
#' @method Cells DimReduc
#'
Cells.DimReduc <- function(x, ...) {
  return(rownames(x = x))
}

#' @rdname DefaultAssay
#' @export
#' @method DefaultAssay DimReduc
#'
DefaultAssay.DimReduc <- function(object, ...) {
  CheckDots(...)
  return(slot(object = object, name = 'assay.used'))
}

#' @rdname DefaultAssay
#' @export
#' @method DefaultAssay<- DimReduc
#'
"DefaultAssay<-.DimReduc" <- function(object, ..., value) {
  CheckDots(...)
  slot(object = object, name = 'assay.used') <- value
  return(object)
}

#' @method Features DimReduc
#' @export
#'
Features.DimReduc <- function(x, projected = NULL, ...) {
  projected <- isTRUE(x = projected %||% Projected(object = x))
  features <- rownames(x = Loadings(object = x, projected = projected))
  if (!length(x = features)) {
    features <- NULL
  }
  return(features)
}

#' @rdname FetchData
#' @export
#' @method FetchData DimReduc
#'
FetchData.DimReduc <- function(
  object,
  vars,
  cells = NULL,
  slot = c('embeddings', 'loadings', 'projected'),
  ...
) {
  slot <- slot[1L]
  slot <- match.arg(arg = slot)
  cells <- cells %||% Cells(x = object)
  if (is.numeric(x = cells)) {
    cells <- Cells(x = object)[cells]
  }
  pattern <- paste0('^(', Key(object = object), ')?[[:digit:]]+$')
  vars <- grep(pattern = pattern, x = vars, value = TRUE)
  if (!length(x = 'vars')) {
    stop("invalid vars")
  }
  vars <- gsub(pattern = Key(object = object), replacement = '', x = vars)
  vars <- as.integer(x = vars)
  vars <- paste0(Key(object = object), vars)
  data <- switch(
    EXPR = slot,
    'embeddings' = Embeddings(object = object),
    Loadings(object = object, projected = slot == 'projected')
  )
  missing <- setdiff(x = vars, y = colnames(x = data))
  if (length(x = missing) == length(x = vars)) {
    stop("Cannot find any of the requested dimensions")
  } else if (length(x = missing)) {
    warning(
      "Cannot find the following dimensions: ", paste0(missing, collapse = ', '),
      call. = FALSE,
      immediate. = TRUE
    )
    vars <- setdiff(x = vars, y = missing)
  }
  return(as.data.frame(x = data)[cells, vars, drop = FALSE])
}

#' @rdname Embeddings
#' @export
#' @method Embeddings DimReduc
#'
#' @examples
#' # Get the embeddings directly from a DimReduc object
#' Embeddings(object = pbmc_small[["pca"]])[1:5, 1:5]
#'
Embeddings.DimReduc <- function(object, ...) {
  CheckDots(...)
  return(slot(object = object, name = 'cell.embeddings'))
}

#' @method FetchData DimReduc
#' @export
#'
FetchData.DimReduc <- function(
  object,
  vars,
  cells = NULL,
  layer = c('embeddings', 'loadings', 'projected'),
  ...
) {
  layer <- layer[1L]
  layer <- match.arg(arg = layer)
  cells <- cells %||% Cells(x = object)
  if (is.numeric(x = cells)) {
    cells <- Cells(x = object)[cells]
  }
  key <- Key(object = object)
  ovars <- vars
  vars <- grep(
    pattern = paste0('^(', key, ')?[[:digit:]]+$'),
    x = vars,
    value = TRUE
  )
  if (!length(x = vars)) {
    stop(
      "None of the vars provided are valid for reduced dimensions",
      call. = FALSE
    )
  } else if (length(x = vars) != length(x = ovars)) {
    warning(
      "The following requested vars are not valid: ",
      paste(setdiff(x = ovars, y = vars), collapse = ', '),
      call. = FALSE,
      immediate. = TRUE
    )
  }
  vars <- paste0(
    key,
    as.integer(x = gsub(pattern = key, replacement = '', x = vars))
  )
  data <- switch(
    EXPR = layer,
    'embeddings' = Embeddings(object = object),
    Loadings(object = object, projected = layer == 'projected')
  )
  missing <- setdiff(x = vars, y = colnames(x = data))
  if (length(x = missing) == length(x = vars)) {
    stop("Cannot find any of the requested dimensions", call. = FALSE)
  } else if (length(x = missing)) {
    warning(
      "Cannot find the following dimensions: ", paste0(missing, collapse = ', '),
      call. = FALSE,
      immediate. = TRUE
    )
    vars <- setdiff(x = vars, y = missing)
  }
  return(as.data.frame(x = data)[cells, vars, drop = FALSE])
}

#' @rdname IsGlobal
#' @export
#' @method IsGlobal DimReduc
#'
IsGlobal.DimReduc <- function(object, ...) {
  object <- UpdateSlots(object = object)
  return(slot(object = object, name = 'global'))
}

#' @param slot Name of slot to store JackStraw scores to
#' Can shorten to 'empirical', 'fake', 'full', or 'overall'
#'
#' @rdname JS
#' @export
#' @method JS DimReduc
#'
JS.DimReduc <- function(object, slot = NULL, ...) {
  CheckDots(...)
  jackstraw <- slot(object = object, name = 'jackstraw')
  if (!is.null(x = slot)) {
    jackstraw <- JS(object = jackstraw, slot = slot)
  }
  return(jackstraw)
}

#' @rdname JS
#' @export
#' @method JS<- DimReduc
#'
"JS<-.DimReduc" <- function(object, slot = NULL, ..., value) {
  CheckDots(...)
  if (inherits(x = value, what = 'JackStrawData')) {
    slot(object = object, name = 'jackstraw') <- value
  } else if (is.null(x = NULL)) {
    stop("A slot must be specified")
  } else {
    JS(object = JS(object = object), slot = slot) <- value
  }
  return(object)
}

#' @rdname Key
#' @export
#' @method Key DimReduc
#'
#' @examples
#' # Get a DimReduc key
#' Key(object = pbmc_small[["pca"]])
#'
Key.DimReduc <- function(object, ...) {
  CheckDots(...)
  return(NextMethod())
}

#' @rdname Key
#' @export
#' @method Key<- DimReduc
#'
#' @examples
#' # Set the key for DimReduc
#' Key(object = pbmc_small[["pca"]]) <- "newkey2_"
#' Key(object = pbmc_small[["pca"]])
#'
"Key<-.DimReduc" <- function(object, ..., value) {
  CheckDots(...)
  object <- UpdateSlots(object = object)
  object <- NextMethod()
  return(object)
}

#' @param projected Pull the projected feature loadings?
#'
#' @rdname Loadings
#' @export
#' @method Loadings DimReduc
#'
#' @examples
#' # Get the feature loadings for a given DimReduc
#' Loadings(object = pbmc_small[["pca"]])[1:5,1:5]
#'
Loadings.DimReduc <- function(object, projected = FALSE, ...) {
  CheckDots(...)
  projected <- projected %||% Projected(object = object)
  slot <- ifelse(
    test = projected,
    yes = 'feature.loadings.projected',
    no = 'feature.loadings'
  )
  return(slot(object = object, name = slot))
}

#' @rdname Loadings
#' @export
#' @method Loadings<- DimReduc
#'
#' @examples
#' # Set the feature loadings for a given DimReduc
#' new.loadings <- Loadings(object = pbmc_small[["pca"]])
#' new.loadings <- new.loadings + 0.01
#' Loadings(object = pbmc_small[["pca"]]) <- new.loadings
#'
"Loadings<-.DimReduc" <- function(object, projected = TRUE, ..., value) {
  CheckDots(...)
  slot.use <- ifelse(
    test = projected,
    yes = 'feature.loadings.projected',
    no = 'feature.loadings'
  )
  if (ncol(x = value) != length(x = object)) {
    stop("New feature loadings must have the dimensions as currently calculated")
  }
  slot(object = object, name = slot.use) <- value
  return(object)
}

#' @rdname Misc
#' @export
#' @method Misc DimReduc
#'
Misc.DimReduc <- .Misc

#' @rdname Misc
#' @export
#' @method Misc<- DimReduc
#'
"Misc<-.DimReduc" <- `.Misc<-`

#' @rdname RenameCells
#' @export
#' @method RenameCells DimReduc
#'
#' @examples
#' # Rename cells in a DimReduc
#' head(x = Cells(x = pbmc_small[["pca"]]))
#' renamed.dimreduc <- RenameCells(
#'     object = pbmc_small[["pca"]],
#'     new.names = paste0("A_", Cells(x = pbmc_small[["pca"]]))
#' )
#' head(x = Cells(x = renamed.dimreduc))
#'
RenameCells.DimReduc <- function(object, new.names = NULL, ...) {
  CheckDots(...)
  old.data <- Embeddings(object = object)
  rownames(x = old.data) <- new.names
  slot(object = object, name = "cell.embeddings") <- old.data
  validObject(object = object)
  return(object)
}

#' @rdname Stdev
#' @export
#' @method Stdev DimReduc
#'
#' @examples
#' # Get the standard deviations for each PC from the DimReduc object
#' Stdev(object = pbmc_small[["pca"]])
#'
Stdev.DimReduc <- function(object, ...) {
  CheckDots(...)
  return(slot(object = object, name = 'stdev'))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' \code{DimReduc} Methods
#'
#' Methods for \code{\link{DimReduc}} objects for generics defined in
#' other packages
#'
#' @inheritParams Assay-methods
#' @param x,object A \code{\link{DimReduc}} object
#' @param i For \code{[}: feature names or indices; for \code{[[}: cell names
#' or indices
#' @param j Dimensions to pull for
#' @param ... Arguments passed to other methods
#'
#' @name DimReduc-methods
#' @rdname DimReduc-methods
#'
#' @concept dimreduc
#'
NULL

#' @describeIn DimReduc-methods Pull feature loadings
#'
#' @return \code{[}: Feature loadings for features \code{i} and dimensions
#' \code{j}
#'
#' @export
#' @method [ DimReduc
#'
"[.DimReduc" <- function(x, i, j, drop = FALSE, ...) {
  loadings <- Loadings(object = x)
  if (missing(x = i)) {
    i <- 1:nrow(x = loadings)
  }
  if (missing(x = j)) {
    j <- names(x = x)
  } else if (is.numeric(x = j)) {
    j <- names(x = x)[j]
  }
  bad.j <- j[!j %in% colnames(x = loadings)]
  j <- j[!j %in% bad.j]
  if (length(x = j) == 0) {
    stop("None of the requested loadings are present.")
  }
  if (length(x = bad.j) > 0) {
    warning(
      "The following loadings are not present: ",
      paste(bad.j, collapse = ", ")
    )
  }
  return(Loadings(object = x)[i, j, drop = drop, ...])
}

#' @describeIn DimReduc-methods Pull cell embeddings
#'
#' @return \code{[[}: Cell embeddings for cells \code{i} and dimensions \code{j}
#'
#' @export
#' @method [[ DimReduc
#'
"[[.DimReduc" <- function(x, i, j, drop = FALSE, ...) {
  if (missing(x = i)) {
    i <- 1:nrow(x = x)
  }
  if (missing(x = j)) {
    j <- names(x = x)
  } else if (is.numeric(x = j)) {
    j <- names(x = x)[j]
  }
  embeddings <- Embeddings(object = x)
  bad.j <- j[!j %in% colnames(x = embeddings)]
  j <- j[!j %in% bad.j]
  if (length(x = j) == 0) {
    stop("None of the requested embeddings are present.")
  }
  if (length(x = bad.j) > 0) {
    warning(
      "The following embeddings are not present: ",
      paste(bad.j, collapse = ", ")
    )
  }
  return(embeddings[i, j, drop = drop, ...])
}

#' @describeIn DimReduc-methods The number of cells and dimensions for a
#' \code{DimReduc}
#'
#' @return \code{dim}: The number of cells (\code{nrow}) and dimensions
#' (\code{ncol})
#'
#' @export
#' @method dim DimReduc
#'
dim.DimReduc <- function(x) {
  return(dim(x = Embeddings(object = x)))
}

#' @describeIn DimReduc-methods The cell and dimension names for a
#' \code{DimReduc} object
#'
#' @return \code{dimnames}: The cell (row) and dimension (column) names
#'
#' @export
#' @method dimnames DimReduc
#'
dimnames.DimReduc <- function(x) {
  return(dimnames(x = Embeddings(object = x)))
}

#' @describeIn DimReduc-methods The number of dimensions for a \code{DimReduc}
#' object
#'
#' @return \code{length}: The number of dimensions
#'
#' @export
#' @method length DimReduc
#'
length.DimReduc <- function(x) {
  return(ncol(x = Embeddings(object = x)))
}

#' @describeIn DimReduc-methods Merge two or more \code{DimReduc} objects
#' together
#'
#' @export
#' @method merge DimReduc
#'
merge.DimReduc <- function(
  x = NULL,
  y = NULL,
  add.cell.ids = NULL,
  ...
) {
  CheckDots(...)
  drs <- c(x, y)
  if (!is.null(x = add.cell.ids)) {
    for (i in 1:length(x = drs)) {
      drs[[i]] <- RenameCells(object = drs[[i]], new.names = add.cell.ids[i])
    }
  }
  embeddings.mat <- list()
  min.dim <- c()
  for (i in 1:length(x = drs)) {
    embeddings.mat[[i]] <- Embeddings(object = drs[[i]])
    min.dim <- c(min.dim, ncol(x = embeddings.mat[[i]]))
  }
  if (length(x = unique(x = min.dim)) > 1) {
    min.dim <- min(min.dim)
    warning(
      "Reductions contain differing numbers of dimensions, merging first ",
      min.dim,
      call. = FALSE,
      immediate. = TRUE
    )
    embeddings.mat <- lapply(
      X = embeddings.mat,
      FUN = function(x) {
        return(x[, 1:min.dim])
      }
    )
  }
  embeddings.mat <- do.call(what = rbind, args = embeddings.mat)
  merged.dr <- CreateDimReducObject(
    embeddings = embeddings.mat,
    loadings = Loadings(object = drs[[1]], projected = FALSE),
    projected = Loadings(object = drs[[1]], projected = TRUE),
    assay = DefaultAssay(object = drs[[1]]),
    key = Key(object = drs[[1]]),
    global = IsGlobal(object = drs[[1]])
  )
  return(merged.dr)
}

#' @describeIn DimReduc-methods The dimension names for a \code{DimReduc} object
#'
#' @return \code{names}: The names for the dimensions (eg. \dQuote{PC_1})
#'
#' @export
#' @method names DimReduc
#'
names.DimReduc <- function(x) {
  return(colnames(x = Embeddings(object = x)))
}

#' @describeIn DimReduc-methods Prints a set of features that most strongly
#' define a set of components; \strong{note}: requires feature loadings to be
#' present in order to work
#'
#' @param dims Number of dimensions to display
#' @param nfeatures Number of genes to display
#' @param projected Use projected slot
#' @param ... Arguments passed to other methods
#'
#' @return \code{print}: Displays set of features defining the components and
#' invisibly returns \code{x}
#'
#' @aliases print
#' @seealso \code{\link[base]{cat}}
#'
#' @export
#' @method print DimReduc
#'
print.DimReduc <- function(
  x,
  dims = 1:5,
  nfeatures = 20,
  projected = FALSE,
  ...
) {
  CheckDots(...)
  loadings <- Loadings(object = x, projected = projected)
  if (!IsMatrixEmpty(x = loadings)) {
    nfeatures <- min(nfeatures, nrow(x = loadings))
    if (ncol(x = loadings) == 0) {
      warning("Dimensions have not been projected. Setting projected = FALSE")
      projected <- FALSE
      loadings <- Loadings(object = x, projected = projected)
    }
    if (min(dims) > ncol(x = loadings)) {
      stop("Cannot print dimensions greater than computed")
    }
    if (max(dims) > ncol(x = loadings)) {
      warning("Only ", ncol(x = loadings), " dimensions have been computed.")
      dims <- intersect(x = dims, y = seq_len(length.out = ncol(x = loadings)))
    }
    for (dim in dims) {
      # features <- TopFeatures(
      #   object = x,
      #   dim = dim,
      #   nfeatures = nfeatures * 2,
      #   projected = projected,
      #   balanced = TRUE
      # )
      features <- Top(
        data = loadings[, dim, drop = FALSE],
        num = nfeatures * 2,
        balanced = TRUE
      )
      cat(Key(object = x), dim, '\n')
      pos.features <- split(
        x = features$positive,
        f = ceiling(x = seq_along(along.with = features$positive) / 10)
      )
      cat("Positive: ", paste(pos.features[[1]], collapse = ", "), '\n')
      pos.features[[1]] <- NULL
      if (length(x = pos.features) > 0) {
        for (i in pos.features) {
          cat("\t  ", paste(i, collapse = ", "), '\n')
        }
      }
      neg.features <- split(
        x = features$negative,
        f = ceiling(x = seq_along(along.with = features$negative) / 10)
      )
      cat("Negative: ", paste(neg.features[[1]], collapse = ", "), '\n')
      neg.features[[1]] <- NULL
      if (length(x = neg.features) > 0) {
        for (i in neg.features) {
          cat("\t  ", paste(i, collapse = ", "), '\n')
        }
      }
    }
  }
  return(invisible(x = x))
}

#' @describeIn DimReduc-methods Subset a \code{DimReduc} object
#'
#' @param cells,features Cells and features to keep during the subset
#'
#' @return \code{subset}: \code{x} for cells \code{cells} and features
#' \code{features}
#'
#' @export
#' @method subset DimReduc
#'
subset.DimReduc <- function(x, cells = NULL, features = NULL, ...) {
  CheckDots(...)
  cells <- Cells(x = x) %iff% cells %||% Cells(x = x)
  if (all(is.na(x = cells))) {
    cells <- Cells(x = x)
  } else if (any(is.na(x = cells))) {
    warning("NAs passed in cells vector, removing NAs")
    cells <- na.omit(object = cells)
  }
  # features <- rownames(x = x) %iff% features %||% rownames(x = x)
  features <- rownames(x = Loadings(object = x)) %iff% features %||% rownames(x = Loadings(object = x))
  if (all(sapply(X = list(features, cells), FUN = length) == dim(x = x))) {
    return(x)
  }
  slot(object = x, name = 'cell.embeddings') <- if (is.null(x = cells)) {
    new(Class = 'matrix')
  } else {
    if (is.numeric(x = cells)) {
      cells <- Cells(x = x)[cells]
    }
    cells <- intersect(x = cells, y = Cells(x = x))
    if (length(x = cells) == 0) {
      stop("Cannot find cell provided", call. = FALSE)
    }
    x[[cells, , drop = FALSE]]
  }
  slot(object = x, name = 'feature.loadings') <- if (is.null(x = features)) {
    new(Class = 'matrix')
  } else {
    if (is.numeric(x = features)) {
      features <- rownames(x = x)[features]
    }
    features.loadings <- intersect(
      x = rownames(x = Loadings(object = x, projected = FALSE)),
      y = features
    )
    if (length(x = features.loadings) == 0) {
      stop("Cannot find features provided", call. = FALSE)
    }
    Loadings(object = x, projected = FALSE)[features.loadings, , drop = FALSE]
  }
  slot(object = x, name = 'feature.loadings.projected') <- if (is.null(x = features) || !Projected(object = x)) {
    new(Class = 'matrix')
  } else {
    features.projected <- intersect(
      x = rownames(x = Loadings(object = x, projected = TRUE)),
      y = features
    )
    if (length(x = features.projected) == 0) {
      stop("Cannot find features provided", call. = FALSE)
    }
    Loadings(object = x, projected = TRUE)[features.projected, , drop = FALSE]
  }
  slot(object = x, name = 'jackstraw') <- new(Class = 'JackStrawData')
  return(x)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @describeIn DimReduc-methods Show basic summary of a \code{DimReduc} object
#'
#' @return \code{show}: Prints summary to \code{\link[base]{stdout}} and
#' invisibly returns \code{NULL}
#'
#' @importFrom methods show
#'
#' @export
#'
setMethod(
  f = 'show',
  signature = 'DimReduc',
  definition = function(object) {
    cat(
      "A dimensional reduction object with key", Key(object = object), '\n',
      'Number of dimensions:', length(x = object), '\n',
      'Projected dimensional reduction calculated: ', Projected(object = object), '\n',
      'Jackstraw run:', as.logical(x = JS(object = object)), '\n',
      'Computed using assay:', DefaultAssay(object = object), '\n'
    )
    return(invisible(x = NULL))
  }
)

setValidity(
  Class = 'DimReduc',
  method = function(object) {
    valid <- NULL
    return(valid %||% TRUE)
  }
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.RenameFeatures <- function(object, old.names = NULL, new.names = NULL) {
  if (is.null(x = old.names) && is.null(x = new.names)) {
    return(object)
  }
  # Checks
  op <- options(Seurat.object.validate = FALSE)
  on.exit(expr = options(op))
  stopifnot(length(x = old.names) == length(x = new.names))
  stopifnot(all(nzchar(x = old.names)))
  stopifnot(all(nzchar(x = new.names)))
  if (is.null(x = Features(x = object)) && length(x = new.names)) {
    warning("No features present in this DimReduc", call. = FALSE, immediate. = TRUE)
  }
  # Rename features
  names(x = new.names) <- old.names
  features <- Features(x = object, projected = FALSE)
  ldat <- Loadings(object = object, projected = FALSE)
  rownames(x = ldat) <- unname(obj = new.names[features])
  Loadings(object = object, projected = FALSE) <- ldat
  if (isTRUE(x = Projected(object = object))) {
    pdat <- Loadings(object = object, projected = TRUE)
    pfeatures <- Features(x = object, projected = TRUE)
    rownames(x = pdat) <- unname(obj = new.names[pfeatures])
    Loadings(object = object, projected = TRUE) <- pdat
  }
  # Validate and return
  options(op)
  validObject(object = object)
  return(object)
}

#' Check to see if projected loadings have been set
#'
#' @param object a DimReduc object
#'
#' @return TRUE if projected loadings have been set, else FALSE
#'
#' @keywords internal
#'
#' @noRd
#'
Projected <- function(object) {
  return(!IsMatrixEmpty(x = Loadings(object = object, projected = TRUE)))
}
