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
#' @slot jackstraw A \code{\link{JackStrawData-class}} object associated with
#' this \code{DimReduc}
#' @template slot-misc
#' @template slot-key
#'
#' @exportClass DimReduc
#'
#' @aliases DimReduc
#'
#' @family dimreduc
#'
setClass(
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
#' @family dimreduc
#'
#' @examples
#' data <- GetAssayData(pbmc_small[["RNA"]], layer = "scale.data")
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
    warn(message = "No assay specified, setting assay as RNA by default.")
    assay <- 'RNA'
  }
  # Try to infer key from column names
  if (is.null(x = key) && is.null(x = colnames(x = embeddings))) {
    abort(message = "Please specify a key for the DimReduc object")
  } else if (is.null(x = key)) {
    key <- regmatches(
      x = colnames(x = embeddings),
      m = regexec(pattern = '^[[:alnum:]]+_', text = colnames(x = embeddings))
    )
    key <- unique(x = unlist(x = key, use.names = FALSE))
  }
  if (length(x = key) != 1) {
    abort(message = "Please specify a key for the DimReduc object")
  } else if (!grepl(pattern = .KeyPattern(), x = key)) {
    old.key  <- key
    key <- Key(object = key)
    colnames(x = embeddings) <- gsub(
      x = colnames(x = embeddings),
      pattern = old.key,
      replacement = key
    )
  }
  # ensure colnames of the embeddings are the key followed by a numeric
  if (is.null(x = colnames(x = embeddings))) {
    warn(message = paste0(
      "No columnames present in cell embeddings, setting to '",
      key,
      "1:",
      ncol(x = embeddings),
      "'"
    ))
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
      abort(message = "Feature names of loadings matrix cannot be empty")
    }
    colnames(x = loadings) <- colnames(x = embeddings)
  }
  if (!IsMatrixEmpty(x = projected)) {
    if (any(rownames(x = loadings) == '')) {
      abort(message = "Feature names of projected loadings matrix cannot be empty")
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
  # layer = c('embeddings', 'loadings', 'projected'),
  # layer = 'embeddings',
  ...
) {
  layer <- 'embeddings'
  layer <- arg_match0(arg = layer, values = 'embeddings')
  cells <- cells %||% Cells(x = object)
  if (is.numeric(x = cells)) {
    cells <- Cells(x = object)[cells]
  }
  cells <- intersect(x = cells, y = Cells(x = object))
  if (!length(x = cells)) {
    abort(message = "None of the cells requested found in this dimensional reduction")
  }
  key <- Key(object = object)
  ovars <- vars
  vars <- grep(
    pattern = paste0('^(', key, ')?[[:digit:]]+$'),
    x = vars,
    value = TRUE
  )
  if (!length(x = vars)) {
    abort(message = "None of the vars provided are valid for reduced dimensions")
  } else if (length(x = vars) != length(x = ovars)) {
    warn(message = paste(
      "The following requested vars are not valid:",
      paste(setdiff(x = ovars, y = vars), collapse = ', '),
    ))
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
    abort(message = "Cannot find any of the requested dimensions")
  } else if (length(x = missing)) {
    warn(message = paste(
      "Cannot find the following dimensions:",
      paste0(missing, collapse = ', ')
    ))
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
Key.DimReduc <- .Key

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
  op <- options(Seurat.object.validate = FALSE)
  on.exit(expr = options(op), add = TRUE)
  old <- Key(object = object)
  suppressWarnings(expr = object <- NextMethod(), classes = 'validationWarning')
  for (i in c("cell.embeddings", "feature.loadings", "feature.loadings.projected")) {
    mat <- slot(object = object, name = i)
    if (IsMatrixEmpty(x = mat)) {
      next
    }
    colnames(x = mat) <- gsub(
      pattern = paste0('^', old),
      replacement = Key(object = object),
      x = colnames(x = mat)
    )
    slot(object = object, name = i) <- mat
  }
  options(op)
  validObject(object = object)
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

#' Get Feature Loadings
#'
#' Pull feature loadings from a \link[=DimReduc]{dimensional reduction}
#'
#' \code{[} does not distinguish between projected and unprojected feature
#' loadings; to select whether projected or unprojected loadings should be
#' pulled, please use \code{\link{Loadings}}
#'
#' @param x A \code{\link{DimReduc}} object
#' @param i Feature identifiers or indices
#' @param j Dimension identifiers or indices
#' @param drop Coerce the result to the lowest possible dimension; see
#' \code{\link{drop}} for further details
#' @template param-dots-method
#'
#' @return Feature loadings for features \code{i} and dimensions \code{j}
#'
#' @method [ DimReduc
#' @export
#'
#' @family dimreduc
#'
#' @seealso \code{\link{Loadings}}
#'
#' @examples
#' pca <- pbmc_small[["pca"]]
#' pca[1:10, 1:5]
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

#' Get Cell Embeddings
#'
#' Pull cell embeddings from a \link[=DimReduc]{dimensional reduction}
#'
#' @inheritParams [.DimReduc
#' @param i Cell names or indices
#'
#' @return Cell embeddings for cells \code{i} and dimensions \code{j}
#'
#' @method [[ DimReduc
#' @export
#'
#' @family dimreduc
#'
#' @seealso \code{\link{Embeddings}}
#'
#' @examples
#' pca <- pbmc_small[["pca"]]
#' pca[[1:10, 1:5]]
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

#' Dimensional Reduction Meta-Information
#'
#' Pull meta-information about cells and dimensions for a given
#' \link[=DimReduc]{dimensional reduction}; cell meta-information is stored
#' as row meta-information (eg. \code{nrow}, \code{rownames}) and dimension
#' meta-information is stored as column meta-information (eg. \code{ncol},
#' \code{colnames})
#'
#' @inheritParams [.DimReduc
#'
#' @return \code{dim}: The number of cells (\code{nrow}) and dimensions
#' (\code{ncol})
#'
#' @method dim DimReduc
#' @export
#'
#' @family dimreduc
#'
#' @examples
#' pca <- pbmc_small[["pca"]]
#' pca
#' dim(pca)
#'
#' # nrow is number of cells
#' nrow(pca)
#'
#' # rownames pulls cell names
#' head(rownames(pca))
#'
#' # ncol and length are number of dimensions
#' ncol(pca)
#' length(pca)
#'
#' # colnames and names pull dimension identifiers
#' head(colnames(pca))
#' head(names(pca))
#'
dim.DimReduc <- function(x) {
  return(dim(x = Embeddings(object = x)))
}

#' @return \code{dimnames}: The cell (row) and dimension (column) names
#'
#' @rdname dim.DimReduc
#'
#' @method dimnames DimReduc
#' @export
#'
#' @seealso \code{Cells}
#'
dimnames.DimReduc <- function(x) {
  return(dimnames(x = Embeddings(object = x)))
}

#' @return \code{length}: The number of dimensions
#'
#' @rdname dim.DimReduc
#'
#' @export
#' @method length DimReduc
#'
length.DimReduc <- function(x) {
  return(ncol(x = x))
}

#' Merge Dimensional Reductions
#'
#' Merge two or more \link[=DimReduc]{dimensional reduction} together
#'
#' @inheritParams [.DimReduc
#' @inheritParams merge.Assay5
#' @param y One or more \code{\link{DimReduc}} objects
#' @template param-dots-ignored
#'
#' @return A new \code{DimReduc} object with data merged from \code{c(x, y)}
#'
#' @method merge DimReduc
#' @export
#'
#' @family dimreduc
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
    add.cell.ids <- unique(x = add.cell.ids)
    if (!is_bare_character(x = add.cell.ids, n = length(x = drs))) {
      abort(
        message = "'add.cell.ids' must be unique for every dimensional reduction"
      )
    }
    for (i in seq_along(along.with = drs)) {
      drs[[i]] <- RenameCells(object = drs[[i]], new.names = add.cell.ids[i])
    }
  }
  all.cells <- unlist(x = lapply(X = drs, FUN = Cells))
  if (anyDuplicated(x = all.cells)) {
    abort(message = "Duplicate cells in provided dimensional reductions")
  }
  embeddings.mat <- lapply(X = drs, FUN = Embeddings)
  min.dim <- vapply(
    X = embeddings.mat,
    FUN = ncol,
    FUN.VALUE = integer(length = 1L),
    USE.NAMES = FALSE
  )
  # embeddings.mat <- list()
  # min.dim <- c()
  # for (i in 1:length(x = drs)) {
  #   embeddings.mat[[i]] <- Embeddings(object = drs[[i]])
  #   min.dim <- c(min.dim, ncol(x = embeddings.mat[[i]]))
  # }
  if (length(x = unique(x = min.dim)) > 1) {
    min.dim <- min(min.dim)
    warn(message = paste(
      "Reductions contain differing numbers of dimensions, merging first",
      min.dim
    ))
    # warning(
    #   "Reductions contain differing numbers of dimensions, merging first ",
    #   min.dim,
    #   call. = FALSE,
    #   immediate. = TRUE
    # )
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

#' @return \code{names}: The dimension identifiers
#'
#' @rdname dim.DimReduc
#'
#' @method names DimReduc
#' @export
#'
names.DimReduc <- function(x) {
  # return(colnames(x = Embeddings(object = x)))
  return(colnames(x = x))
}

#' Print Top Feature Loadings
#'
#' Prints a set of features that most strongly define a set of components;
#' \strong{note}: requires feature loadings to be present in order to work
#'
#' @inheritParams [.DimReduc
#' @param dims Number of dimensions to display
#' @param nfeatures Number of genes to display
#' @param projected Use projected slot
#' @template param-dots-ignored
#'
#' @return Displays set of features defining the components and
#' invisibly returns \code{x}
#'
#' @method print DimReduc
#' @export
#'
#' @aliases print
#'
#' @family dimreduc
#'
#' @seealso \code{\link[base]{cat}}
#'
#' @examples
#' pca <- pbmc_small[["pca"]]
#' print(pca)
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

#' Subset a Dimensional Reduction
#'
#' Subset a \code{\link{DimReduc}} object
#'
#' @inheritParams [.DimReduc
#' @param cells,features Cells and features to keep during the subset
#' @template param-dots-ignored
#'
#' @return \code{x} for cells \code{cells} and features \code{features}
#'
#' @method subset DimReduc
#' @export
#'
#' @family dimreduc
#'
subset.DimReduc <- function(x, cells = NULL, features = NULL, ...) {
  CheckDots(...)
  cells <- Cells(x = x) %iff% cells %||% Cells(x = x)
  if (all(is.na(x = cells))) {
    cells <- Cells(x = x)
  } else if (any(is.na(x = cells))) {
    warn(message = "NAs passed in cells vector, removing NAs")
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
    cells <- intersect(x = Cells(x = x), y = cells)
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Dimensional Reduction Overview
#'
#' Overview of a \code{\link{DimReduc}} object
#'
#' @param object A dimensional reduction
#'
#' @template return-show
#'
#' @keywords internal
#'
#' @seealso \code{\link{DimReduc}}
#'
#' @examples
#' pca <- pbmc_small[["pca"]]
#' pca
#'
setMethod(
  f = 'show',
  signature = 'DimReduc',
  definition = function(object) {
    cat(
      "A dimensional reduction object with key", Key(object = object), '\n',
      'Number of dimensions:', length(x = object), '\n',
      'Number of cells:', length(x = Cells(x = object)), '\n',
      'Projected dimensional reduction calculated: ', Projected(object = object), '\n',
      'Jackstraw run:', as.logical(x = JS(object = object)), '\n',
      'Computed using assay:', DefaultAssay(object = object), '\n'
    )
    return(invisible(x = NULL))
  }
)

#' Dimensional Reduction Validity
#'
#' @templateVar cls DimReduc
#' @template desc-validity
#'
#' @section Cell Embeddings Validation:
#' The cell embeddings matrix must be a numeric matrix of dimensions
#' \eqn{n_{cells}} by \eqn{d_{dimensions}}; row names must be the cell names
#' and column names must be the dimension identifier. The dimension identifier
#' must be \dQuote{\code{key_dimension}} (eg. \dQuote{\code{PC_1}}). Dimension
#' identifiers must be in order and cannot be skipped
#'
#' @section Feature and Projected Feature Loadings Validation:
#' blah
#'
#' @inheritSection Key-validity Key Validation
#'
#' @section Standard Deviations Validation:
#' blah
#'
#' @name DimReduc-validity
#'
#' @family dimreduc
#'
setValidity(
  Class = 'DimReduc',
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
    ValidColnames <- function(
      mat,
      ref = NULL,
      type = c('embeddings', 'loadings', 'projected')
    ) {
      ret <- NULL
      if (IsMatrixEmpty(x = mat)) {
        return(ret)
      }
      type <- match.arg(arg = type)
      type <- switch(
        EXPR = type,
        embeddings = 'cell.embeddings',
        loadings = 'feature.loadings',
        projected = 'feature.loadings.projected'
      )
      mat.names <- colnames(x = mat)
      key <- paste0('^', Key(object = object))
      if (is.null(x = mat.names)) {
        ret <- c(ret, paste("colnames must be present in", sQuote(x = type)))
      } else if (!all(grepl(pattern = key, x = mat.names))) {
        ret <- c(
          ret,
          paste(
            "colnames for",
            sQuote(x = type),
            "must start with reduction key",
            paste0("(", Key(object = object), ")")
          )
        )
      } else {
        dims <- as.numeric(x = gsub(pattern = key, replacement = '', x = mat.names))
        if (!is_bare_integerish(x = dims, n = ncol(x = mat), finite = TRUE) || any(dims < 1L)) {
          ret <- c(
            ret,
            paste(
              "dimension names for",
              sQuote(x = type),
              "must be positive integers"
            )
          )
        } else if (is.unsorted(x = dims)) {
          ret <- c(
            ret,
            paste("dimensions for", sQuote(x = type), "must be in order")
          )
        }
      }
      if (!is.null(x = ref)) {
        if (length(x = mat.names) != length(x = ref)) {
          ret <- c(
            ret,
            paste(sQuote(x = type), "must have", length(x = ref), "dimensions")
          )
        } else if (!all(mat.names == ref)) {
          ret <- c(
            ret,
            paste(
              "dimensions in",
              sQuote(x = type),
              "do not match dimensions in reduction"
            )
          )
        }
      }
      return(ret)
    }
    # Validate cell embeddings
    emb <- Embeddings(object = object)
    if (!is.numeric(x = emb)) {
      valid <- c(valid, "'cell.embeddings' must be a numeric matrix")
    }
    if (is.null(x = rownames(x = emb)) || !all(nzchar(x = rownames(x = emb)))) {
      valid <- c(valid, "rownames must be present in 'cell.embeddings'")
    }
    valid <- c(valid, ValidColnames(mat = emb, type = 'embeddings'))
    if (!is.null(x = valid)) {
      return(valid)
    }
    dims <- colnames(x = emb)
    # if (is.null(x = colnames(x = emb))) {
    #   valid <- c(valid, "colnames must be present in 'cell.embeddings'")
    # } else {
    #   emb.names <- colnames(x = emb)
    #   if (!all(grepl(pattern = paste0('^', Key(object = object)), x = emb.names))) {
    #     valid <- c(
    #       valid,
    #       paste0(
    #         "colnames for 'cell.embeddings' must start with reduction key (",
    #         Key(object = object),
    #         ")"
    #       )
    #     )
    #   }
    # }
    # if (!is.null(x = valid)) {
    #   return(valid)
    # }
    # TODO: Validate feature loadings
    lds <- Loadings(object = object, projected = FALSE)
    valid <- c(valid, ValidColnames(mat = lds, type = 'loadings'))
    # TODO: Validate projected loadings
    prj <- Loadings(object = object, projected = TRUE)
    valid <- c(valid, ValidColnames(mat = prj, type = 'projected'))
    # TODO: Validate assay used
    if (!rlang::is_scalar_character(x = DefaultAssay(object = object))) {
      valid <- c(valid, "'assay.orig' must be a 1-length character")
    }
    # Validate globalness
    if (!rlang::is_scalar_logical(x = IsGlobal(object = object))) {
      valid <- c(valid, "'global' must be a 1-length logical")
    } else if (is_na(x = IsGlobal(object = object))) {
      valid <- c(valid, "'global' must be TRUE or FALSE")
    }
    # TODO: Validate standard deviations
    # TODO: Validate JackStraw data
    # TODO: Validate misc
    return(valid %||% TRUE)
  }
)
