#' @include zzz.R
#' @include assay.R
#' @include layers.R
#' @include logmap.R
#' @include keymixin.R
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
#' data (layers) and associated feature-level metadata. Derived classes
#' (eg. \link[Assay-class]{the v5 Assay}) may optionally define additional
#' functionality
#'
#' @slot layers A named list containing alternate representations of
#' \code{data}; layers must have the same cells as \code{data} and either the
#' same or a subset of the features present in \code{data}
#' @slot cells A \code{\link{LogMap}} describing the cell membership for
#' each layer
#' @slot features A \code{\link{LogMap}} describing the feature membership
#' for each layer
#' @slot assay.orig ...
#' @slot meta.data ...
#' @slot misc ...
#' @slot key A one-length character vector with the object's key; keys must
#' be one or more alphanumeric characters followed by an underscore
#' \dQuote{\code{_}} (regex pattern \dQuote{\code{^[[:alnum:]]+_$}})
#'
#' @exportClass StdAssay
#'
#' @aliases StdAssay
#'
#' @seealso \code{\link{Assay5-class}}
#'
setClass(
  Class = 'StdAssay',
  contains = c('VIRTUAL', 'KeyMixin'),
  slots = c(
    layers = 'list',
    cells = 'LogMap',
    features = 'LogMap',
    default = 'integer',
    # var.features = 'character',
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

setClass(
  Class = 'Assay5T',
  contains = 'StdAssay'
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @method .AssayClass Assay5T
#' @export
#'
.AssayClass.Assay5T <- function(object) {
  return('Transposed Assay (v5)')
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
  .NotYetImplemented()
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
  # Check layer names
  if (is.null(x = names(x = counts))) {
    names(x = counts) <- paste0('counts', seq_along(along.with = counts))
  } else if (any(!nzchar(x = names(x = counts)))) {
    names(x = counts)[!nzchar(x = names(x = counts))] <- 'counts'
  }
  names(x = counts) <- make.unique(names = names(x = counts), sep = '')
  # Check cell/feature names for all layers
  if (is.atomic(x = cells)) {
    cells <- rep_len(x = list(cells), length.out = length(x = counts))
  }
  if (!is_bare_list(x = cells) || length(x = cells) != length(x = counts)) {
    stop("Not enough cells for the counts matrices provided", call. = FALSE)
  }
  cells <- .CheckNames(x = cells, n = names(x = counts))
  if (is.atomic(x = features)) {
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
      cells.use <- which(x = csum(counts[[layer]] > 0) >= min.features)
      counts[[layer]] <- if (cdim == 1L) {
        counts[[layer]][cells.use, ]
      } else {
        counts[[layer]][, cells.use]
      }
      cells[[layer]] <- cells[[layer]][cells.use]
    }
  }
  # Filter based on min.cells
  if (min.cells > 0) {
    for (layer in names(x = counts)) {
      features.use <- which(x = fsum(counts[[layer]] > 0) >= min.cells)
      counts[[layer]] <- if (fdim == 1L) {
        counts[[layer]][features.use, ]
      } else {
        counts[[layer]][, features.use]
      }
      features[[layer]] <- features[[layer]][features.use]
    }
  }
  features.all <- Reduce(f = union, x = features)
  # Create the object
  object <- new(
    Class = type,
    layers = list(),
    default = 0L,
    features = LogMap(y = features.all),
    cells = LogMap(y = Reduce(f = union, x = cells)),
    meta.data = EmptyDF(n = length(x = features.all)),
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
.MARGIN.Assay5T <- function(object, type = c('features', 'cells'), ...) {
  type <- type[1]
  type <- match.arg(arg = type)
  return(unname(obj = c(features = 2L, cells = 1L)[type]))
}

#' @rdname AddMetaData
#' @method AddMetaData StdAssay
#' @export
#'
AddMetaData.StdAssay <- .AddMetaData

#' @importFrom methods as
#'
#' @method CastAssay StdAssay
#' @export
#'
CastAssay.StdAssay <- function(object, to, layers = NULL, verbose = TRUE, ...) {
  layers <- Layers(object = object, search = layers)
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
      warning(
        "Unable to cast layer ",
        lyr,
        ": ",
        e$message,
        call. = FALSE,
        immediate. = TRUE
      )
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

#' @rdname Cells
#' @method Cells StdAssay
#' @export
#'
Cells.StdAssay <- function(x, layer = NULL, ...) {
  layer <- layer[1L] %||% DefaultLayer(object = x)
  if (is_na(x = layer)) {
    return(rownames(x = slot(object = x, name = 'cells')))
  }
  layer <- match.arg(arg = layer, choices = Layers(object = x))
  return(slot(object = x, name = 'cells')[[layer]])
}

#' @rdname CreateAssay5Object
#' @method CreateAssay5Object default
#' @export
#'
CreateAssay5Object.default <- function(
  counts,
  min.cells = 0,
  min.features = 0,
  layer = 'counts',
  transpose = FALSE,
  csum = NULL,
  fsum = NULL,
  ...
) {
  if (isTRUE(x = transpose)) {
    type <- 'Assay5T'
    csum <- csum %||% Matrix::rowSums
    fsum <- fsum %||% Matrix::colSums
  } else {
    type <- 'Assay5'
    csum <- csum %||% Matrix::colSums
    fsum <- fsum %||% Matrix::rowSums
  }
  return(.CreateStdAssay(
    counts = counts,
    min.cells = min.cells,
    min.features = min.features,
    transpose = transpose,
    type = type,
    layer = layer,
    csum = csum,
    fsum = fsum,
    ...
  ))
}

#' @rdname CreateAssay5Object
#' @method CreateAssay5Object list
#' @export
#'
CreateAssay5Object.list <- function(
  counts,
  min.cells = 0,
  min.features = 0,
  transpose = FALSE,
  csum = NULL,
  fsum = NULL,
  ...
) {
  if (any(sapply(X = counts, FUN = inherits, what = 'spam'))) {
    colsums <- spam::colSums
    rowsums <- spam::rowSums
  } else {
    colsums <- Matrix::colSums
    rowsums <- Matrix::rowSums
  }
  if (isTRUE(x = transpose)) {
    type <- 'Assay5T'
    csum <- csum %||% rowsums
    fsum <- fsum %||% colsums
  } else {
    type <- 'Assay5'
    csum <- csum %||% colsums
    fsum <- fsum %||% rowsums
  }
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

#' @rdname CreateAssay5Object
#' @method CreateAssay5Object Matrix
#' @export
#'
CreateAssay5Object.Matrix <- function(
  counts,
  min.cells = 0,
  min.features = 0,
  transpose = FALSE,
  layer = 'counts',
  ...
) {
  return(.CreateStdAssay(
    counts = counts,
    min.cells = min.cells,
    min.features = min.features,
    transpose = transpose,
    type = ifelse(test = isTRUE(x = transpose), yes = 'Assay5T', no = 'Assay5'),
    layer = layer,
    ...
  ))
}

#' @rdname CreateAssay5Object
#' @method CreateAssay5Object matrix
#' @export
#'
CreateAssay5Object.matrix <- CreateAssay5Object.Matrix

#' @rdname CreateAssay5Object
#' @method CreateAssay5Object spam
#' @export
#'
CreateAssay5Object.spam <- function(
  counts,
  min.cells = 0,
  min.features = 0,
  transpose = FALSE,
  layer = 'counts',
  ...
) {
  if (isTRUE(x = transpose)) {
    type <- 'Assay5T'
    csum <- spam::rowSums
    fsum <- spam::colSums
  } else {
    type <- 'Assay5'
    csum <- spam::colSums
    fsum <- spam::rowSums
  }
  return(.CreateStdAssay(
    counts = counts,
    min.cells = min.cells,
    min.features = min.features,
    transpose = transpose,
    type = type,
    layer = layer,
    csum = csum,
    fsum = fsum,
    ...
  ))
}

#' @rdname DefaultAssay
#' @export
#' @method DefaultAssay StdAssay
#'
DefaultAssay.StdAssay <- function(object, ...) {
  return(slot(object = object, name = 'assay.orig'))
}

#' @rdname DefaultAssay
#' @export
#' @method DefaultAssay<- StdAssay
#'
"DefaultAssay<-.StdAssay" <- function(object, ..., value) {
  slot(object = object, name = 'assay.orig') <- value
  return(object)
}

#' @rdname DefaultLayer
#' @export
#' @method DefaultLayer StdAssay
#'
DefaultLayer.StdAssay <- function(object, ...) {
  idx <- slot(object = object, name = 'default')
  if (!length(x = idx) || idx == 0L) {
    idx <- 1L
  }
  return(Layers(object = object)[seq_len(length.out = idx)])
}

#' @rdname DefaultLayer
#' @export
#' @method DefaultLayer<- StdAssay
#'
"DefaultLayer<-.StdAssay" <- function(object, ..., value) {
  # value <- value[1]
  layers <- Layers(object = object)
  # value <- match.arg(arg = value, choices = layers, several.ok = TRUE)
  value <- Layers(object = object, search = value)
  # idx <- which(x = layers == value)
  idx <- MatchCells(new = layers, orig = value, ordered = TRUE)
  slot(object = object, name = 'layers') <- c(
    slot(object = object, name = 'layers')[idx],
    slot(object = object, name = 'layers')[-idx]
  )
  slot(object = object, name = 'default') <- length(x = value)
  return(object)
}

#' @param layer
#'
#' @rdname Cells
#' @export
#' @method Features StdAssay
#'
Features.StdAssay <- function(x, layer = NULL, ...) {
  layer <- layer[1L] %||% DefaultLayer(object = x)[1L]
  if (is_na(x = layer)) {
    return(rownames(x = slot(object = x, name = 'features')))
  }
  layer <- match.arg(arg = layer, choices = Layers(object = x))
  return(slot(object = x, name = 'features')[[layer]])
}

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
  layer <- rev(x = Layers(
    object = object,
    search = layer %||% DefaultLayer(object = object)
  ))
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
    data.fetched[lcells, lvars] <- t(x = LayerData(
      object = object,
      layer = lyr,
      cells = lcells,
      features = lvars
    ))
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
      warning(
        "Removing ",
        length(x = no.data),
        " cells missing data for features requested",
        call. = FALSE,
        immediate. = TRUE
      )
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
    stop("None of the requested variables found", call. = FALSE)
  } else if (length(x = missing)) {
    warning(
      "The following variables could not be found: ",
      paste(missing, collapse = ', '),
      call. = FALSE,
      immediate. = TRUE
    )
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

#' @importFrom utils adist
#'
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
HVFInfo.StdAssay <- function(
  object,
  method = NULL,
  status = FALSE,
  layer = NULL,
  strip = TRUE,
  ...
) {
  # Find available HVF methods and layers
  vf.methods <- .VFMethods(object = object, type = 'hvf')
  vf.layers <- .VFLayers(object = object, type = 'hvf')
  # Determine which method and layer to use
  method <- (method %||% vf.methods)[1L]
  method <- tryCatch(
    expr = match.arg(arg = method, choices = vf.methods),
    error = function(...) {
      return(NULL)
    }
  )
  # If no methods found, return NULL
  if (is.null(x = method)) {
    return(method)
  }
  layer <- (layer %||% vf.layers)[1L]
  layer <- vf.layers[which.min(x = adist(x = layer, y = vf.layers))]
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
  return(NextMethod())
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
  object <- NextMethod()
  return(object)
}

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
#' @method LayerData StdAssay
#' @export
#'
LayerData.StdAssay <- function(
  object,
  layer = NULL,
  cells = NULL,
  features = NULL,
  fast = FALSE,
  ...
) {
  # Figure out the layer we're pulling
  layer <- layer[1L] %||% DefaultLayer(object = object)[1L]
  layer <- match.arg(arg = layer, choices = Layers(object = object))
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
  # Pull the layer data
  ldat <- if (.MARGIN(object = object) == 1L) {
    slot(object = object, name = 'layers')[[layer]][features, cells, drop = FALSE]
  } else {
    slot(object = object, name = 'layers')[[layer]][cells, features, drop = FALSE]
  }
  # Add dimnames and transpose if requested
  ldat <- if (isTRUE(x = fast)) {
    ldat
  } else if (is.na(x = fast)) {
    .GetLayerData(
      x = ldat,
      dnames = dnames,
      fmargin = 1L,
      ...
    )
  } else {
    .GetLayerData(
      x = ldat,
      dnames = dnames,
      fmargin = .MARGIN(object = object, type = 'features'),
      ...
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
    if (slot(object = object, name = 'default') > length(x = Layers(object = object))) {
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
  fdim <- .MARGIN(object = object, type = 'features')
  cdim <- .MARGIN(object = object, type = 'cells')
  # Assume input matrix is features x cells
  dnames <- list(
    # features %||% dimnames(x = value)[[fdim]],
    features %||% dimnames(x = value)[[1L]],
    # cells %||% dimnames(x = value)[[cdim]]
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
  value <- .PrepLayerData(
    x = value,
    target = dim(x = object),
    dnames = dnames,
    fmargin = fdim,
    ...
  )
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

#' @param search A pattern to search layer names for
#'
#' @rdname Layers
#' @method Layers StdAssay
#' @export
#'
Layers.StdAssay <- function(object, search = NULL, ...) {
  layers <- names(x = slot(object = object, name = 'layers'))
  if (!is.null(x = search)) {
    layers <- unique(x = unlist(x = lapply(
      X = search,
      FUN = function(lyr) {
        if (lyr %in% layers) {
          return(lyr)
        }
        return(grep(pattern = lyr, x = layers, value = TRUE, ...))
      }
    )))
    if (!length(x = layers)) {
      stop("No layers found matching search pattern provided", call. = FALSE)
    }
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
VariableFeatures.StdAssay <- function(object, method = NULL, layer = NULL, ...) {
  hvf.info <- HVFInfo(
    object = object,
    method = method,
    layer = layer,
    status = TRUE,
    strip = TRUE
  )
  msg <- 'No variable features found'
  if (is.null(x = hvf.info)) {
    warning(msg, call. = FALSE, immediate. = TRUE)
    return(NULL)
  }
  if (!'variable' %in% names(x = hvf.info)) {
    stop(msg, call. = FALSE)
  }
  vf <- rownames(x = hvf.info)[which(x = hvf.info$variable)]
  if ('rank' %in% names(x = hvf.info)) {
    vf <- vf[order(hvf.info$rank[which(x = hvf.info$variable)])]
  } else {
    warning(
      "No variable feature rank found, returning features in assay order",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  return(vf)
}

#' @rdname VariableFeatures
#' @export
#' @method VariableFeatures<- StdAssay
#'
"VariableFeatures<-.StdAssay" <- function(object, method = 'custom', layer = NULL, ..., value) {
  value <- intersect(x = value, y = rownames(x = object))
  if (!length(x = value)) {
    stop("None of the features specified are present in this assay", call. = FALSE)
  }
  layer <- Layers(object = object, search = layer)
  df <- data.frame(TRUE, seq_along(along.with = value), row.names = value)
  for (lyr in layer) {
    names(x = df) <- paste('vf', method, lyr, c('variable', 'rank'), sep = '_')
    object[[]] <- df
  }
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
#' @param value For \code{[[} Feature-level metadata to add to the assay; for
#' \code{dimnames<-}, a list of two character vectors with the first entry being
#' new feature names and the second being new cell names for the assay
#' @param ... Arguments passed to other methods
#'
#' @details The following methods are provided for interacting with a
#' \code{StdAssay} object
#'
#' @name StdAssay-methods
#' @rdname StdAssay-methods
#'
#' @concept assay
#'
NULL

#' @details \code{[}: Get expression data from an \code{StdAssay}
#'
#' @return \code{[}: The \code{data} slot for features \code{i} and cells
#' \code{j}
#'
#' @rdname StdAssay-methods
#'
#' @method [ StdAssay
#' @export
#'
"[.StdAssay" <- function(x, i, j, ...) {
  if (missing(x = i)) {
    i <- seq_len(length.out = nrow(x = x))
  }
  if (missing(x = j)) {
    j <- seq_len(length.out = ncol(x = x))
  }
  return(LayerData(object = x, cells = j, features = i, ...))
}

#' @details \code{[[}: Get feature-level metadata
#'
#' @inheritParams base::`[[.data.frame`
#'
#' @return \code{[[}: The feature-level metadata for \code{i}
#'
#' @rdname StdAssay-methods
#'
#' @method [[ StdAssay
#' @export
#'
"[[.StdAssay" <- function(x, i, ..., drop = FALSE) {
  if (missing(x = i)) {
    i <- colnames(x = slot(object = x, name = 'meta.data'))
  }
  data.return <- slot(object = x, name = 'meta.data')[, i, drop = FALSE, ...]
  # row.names(x = data.return) <- Features(x = x, layer = NA)
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

#' @details \code{dim}: Number of cells and features for an \code{StdAssay}
#'
#' @return \code{dim}: The number of features (\code{nrow}) and cells
#' (\code{ncol})
#'
#' @rdname StdAssay-methods
#'
#' @method dim StdAssay
#' @export
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

#' @details \code{dimnames}: Get the feature and cell names
#'
#' @return \code{dimnames}: A list with feature (row) and cell (column) names
#'
#' @rdname StdAssay-methods
#'
#' @method dimnames StdAssay
#' @export
#'
#' @seealso \code{\link{Cells}} \code{\link{Features}}
#'
dimnames.StdAssay <- function(x) {
  return(list(Features(x = x, layer = NA), Cells(x = x, layer = NA)))
}

#' @details \code{dimnames<-}: Set the feature and cell names
#'
#' @return \code{dimnames<-}: \code{x} with the cell and/or feature
#' names updated to \code{value}
#'
#' @rdname StdAssay-methods
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

#' @details \code{head} and \code{tail}: Get the first or last rows of feature
#' level meta data
#'
#' @return \code{head}: The first \code{n} rows of feature-level metadata
#'
#' @importFrom utils head
#'
#' @rdname StdAssay-methods
#'
#' @method head StdAssay
#' @export
#'
head.StdAssay <- .head

#' @details \code{merge}: Merge multiple assays together
#'
#' @return \code{merge}: A new assay ...
#'
#' @rdname StdAssay-methods
#'
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
  # TODO: Support multiple types of assays
  if (length(x = unique(x = sapply(X = assays, FUN = class))) != 1L) {
    stop("Multiple types of assays provided", call. = FALSE)
  }
  labels <- labels %||% as.character(x = seq_along(along.with = assays))
  # add.cell.ids <- add.cell.ids %||% labels
  # TODO: Support collapsing layers
  if (isTRUE(x = collapse)) {
    stop("Collapsing layers is not yet supported", call. = FALSE)
  }
  for (i in seq_along(along.with = assays)) {
    if (is.na(x = labels[i])) {
      labels[i] <- as.character(x = i)
    }
    if (isTRUE(x = is.na(x = add.cell.ids[i]))) {
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
    stop("Collapsing layers is not yet supported", call. = FALSE)
  } else {
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
  validObject(object = combined)
  return(combined)
}

#' @details \code{subset}: Subset an assay to a given set of cells
#' and/or features. \strong{Note}: reordering of cells/features is
#' \emph{not} permitted
#'
#' @return \code{subset}: The assay subsetted to the cells and/or features given
#'
#' @rdname StdAssay-methods
#'
#' @method subset StdAssay
#' @export
#'
subset.StdAssay <- function(
  x,
  cells = NULL,
  features = NULL,
  layers = NULL,
  ...
) {
  if (is.null(x = cells) && is.null(x = features)) {
    return(x)
  }
  # Check the cells vector
  if (all(is.na(x = cells))) {
    cells <- Cells(x = x, layer = NA)
  } else if (any(is.na(x = cells))) {
    warning(
      "NAs passed in cells vector, removing NAs",
      call. = FALSE,
      immediate. = TRUE
    )
    cells <- cells[!is.na(x = cells)]
  }
  if (is.numeric(x = cells)) {
    cells <- Cells(x = x, layer = NA)[cells]
  }
  cells <- intersect(x = cells, y = Cells(x = x, layer = NA))
  if (!length(x = cells)) {
    stop("None of the cells provided found in this assay", call. = FALSE)
  }
  # Check the features vector
  if (all(is.na(x = features))) {
    features <- Features(x = x, layer = NA)
  } else if (any(is.na(x = features))) {
    warning(
      "NAs passed in features vector, removing NAs",
      call. = FALSE,
      immediate. = TRUE
    )
    features <- features[!is.na(x = features)]
  }
  if (is.numeric(x = features)) {
    features <- Features(x = x, layer = NA)[features]
  }
  features <- intersect(x = features, y = Features(x = x, layer = NA))
  if (!length(x = features)) {
    stop("None of the features provided found in this assay", call. = FALSE)
  }
  # Check the layers
  layers.all <- Layers(object = x)
  layers <- layers %||% layers.all
  layers <- match.arg(
    arg = layers,
    choices = layers.all,
    several.ok = TRUE
  )
  # Remove unused layers
  for (lyr in setdiff(x = layers.all, y = layers)) {
    LayerData(object = x, layer = lyr) <- NULL
  }
  # Perform the subsets
  for (l in layers) {
    lcells <- MatchCells(
      new = Cells(x = x, layer = l),
      orig = cells,
      ordered = TRUE
    )
    lfeatures <- MatchCells(
      new = Features(x = x, layer = l),
      orig = features,
      ordered = TRUE
    )
    if (is.null(x = lcells) || is.null(x = features)) {
      LayerData(object = x, layer = l) <- NULL
    } else {
      LayerData(object = x, layer = l) <- LayerData(
        object = x,
        layer = l,
        cells = lcells,
        features = lfeatures
      )
    }
  }
  slot(object = x, name = 'cells') <- droplevels(x = slot(
    object = x,
    name = 'cells'
  ))
  # TODO: Subset feature-level metadata
  mfeatures <- MatchCells(
    new = Features(x = x, layer = NA),
    orig = features,
    ordered = TRUE
  )
  slot(object = x, name = 'meta.data') <- slot(
    object = x,
    name = 'meta.data'
  )[mfeatures, , drop = FALSE]
  # Update the cell/feature maps
  for (i in c('cells', 'features')) {
    slot(object = x, name = i) <- droplevels(x = slot(object = x, name = i))
  }
  validObject(object = x)
  return(x)
}

#' @return \code{tail}: The last \code{n} rows of feature-level metadata
#'
#' @importFrom utils tail
#'
#' @rdname StdAssay-methods
#'
#' @method tail StdAssay
#' @export
#'
tail.StdAssay <- .tail

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.VFLayers <- function(
  object,
  type = c('hvf', 'svf'),
  layers = NULL,
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
      return(x[3L:(length(x = x) - 1L)])
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
  layers = NULL,
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
    # Add the expression matrices
    for (i in c('counts', 'data', 'scale.data')) {
      adata <- GetAssayData(object = from, slot = i)
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
    # Add miscellaneous data
    mdata <- Misc(object = from)
    for (i in names(x = mdata)) {
      Misc(object = to, slot = i) <- mdata[[i]]
    }
    return(to)
  }
)

#' @details \code{[[<-}: Add or remove pieces of meta data
#'
#' @return \code{[[<-}: \code{x} with the metadata updated
#'
#' @rdname StdAssay-methods
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
      stop(
        "Cannot add more or less meta data without feature names",
        call. = FALSE
      )
    }
    for (n in i) {
      v <- value[[n]]
      names(x = v) <- row.names(value)
      x[[n]] <- v
    }
    return(x)
  }
)

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
    # If no `i` provided, use the column names from value
    x[[names(x = value)]] <- value
    return(x)
  }
)

#' @importFrom methods selectMethod
#'
#' @rdname StdAssay-methods
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

#' @rdname StdAssay-methods
#'
setMethod(
  f = '[[<-',
  signature = c(x = 'StdAssay', i = 'character', j = 'missing', value = 'NULL'),
  definition = function(x, i, ..., value) {
    browser()
    for (name in i) {
      slot(object = x, name = 'meta.data')[[name]] <- NULL
    }
    return(x)
  }
)

#' @rdname StdAssay-methods
#'
setMethod(
  f = '[[<-',
  signature = c(x = 'StdAssay', i = 'character', j = 'missing', value = 'vector'),
  definition = function(x, i, ..., value) {
    # Add multiple bits of metadata
    if (length(x = i) > 1L) {
      value <- rep_len(x = value, length.out = length(x = i))
      for (idx in seq_along(along.with = i)) {
        x[i[idx]] <- value[[idx]]
      }
    } else {
      # Add a single column of metadata
      if (is.null(x = names(x = value))) {
        value <- rep_len(x = value, length.out = nrow(x = x))
        names(x = value) <- Features(x = x, layer = NA)
      } else {
        names.intersect <- intersect(
          x = names(x = value),
          y = Features(x = x, layer = NA)
        )
        if (!length(x = names.intersect)) {
          stop(
            "No feature overlap between new meta data and assay",
            call. = FALSE
          )
        }
        value <- value[names.intersect]
      }
      df <- EmptyDF(n = nrow(x = x))
      rownames(x = df) <- Features(x = x, layer = NA)
      df[[i]] <- if (i %in% names(x = x[[]])) {
        x[[i]]
      } else {
        NA
      }
      df[names(x = value), i] <- value
      slot(object = x, name = 'meta.data')[, i] <- df[[i]]
    }
    validObject(object = x)
    return(x)
  }
)

#' @rdname StdAssay-methods
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

#' @rdname StdAssay-methods
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
    # If no `i` provided, use the column names from value
    x[[names(x = value)]] <- value
    return(x)
  }
)

#' @rdname StdAssay-methods
#'
setMethod(
  f = '[[<-',
  signature = c(x = 'StdAssay', i = 'missing', j = 'missing', value = 'NULL'),
  definition = function(x, ..., value) {
    slot(object = x, name = 'meta.data') <- EmptyDF(n = nrow(x = x))
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

#' @details \code{show}: Overview of an \code{StdAssay} object
#'
#' @return \code{show}: Prints summary to \code{\link[base]{stdout}} and
#' invisibly returns \code{NULL}
#'
#' @importFrom utils head
#' @importFrom methods show
#'
#' @rdname StdAssay-methods
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
    cat("Default layer:", DefaultLayer(object = object) %||% "NULL")
    # Layer information
    layers <- setdiff(
      x = Layers(object = object),
      y = DefaultLayer(object = object)
    )
    if (length(x = layers)) {
      cat(
        "\nAdditional layers:\n",
        paste(strwrap(x = paste(layers, collapse = ', ')), collapse = '\n'),
        "\n"
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
