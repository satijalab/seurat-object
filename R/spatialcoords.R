#' @include zzz.R
#' @include generics.R
#' @include centroids.R
#' @include molecules.R
#' @include segmentation.R
#' @include spatial.R
#' @include logmap.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The Spirula Coordinates Object
#'
#' A modern container for storing coordinates of spatially-resolved single
#' cells. Capable of storing multiple layers of centroid/voxel- and
#' polygon/segmentation-based coordinate systems. Supports coordinates for
#' spatially-resolved molecule (FISH) data. Compatible with Seurat's
#' \code{\link{SpatialImage}}
#'
#' @slot cells (\code{\link{LogMap}}) A map of the cells present
#' in each segmentation layer
#' @slot molecules (\code{\link[base]{list}}) A named list of
#' \code{\link{Molecules}} objects defining spatially-resolved
#' molecular coordinates
#' @slot segmentations (\code{[named]\link[base]{list}}
#' \{\code{\link{Segmentation}}, \code{\link{Centroids}}\}) A named list of
#' \code{\link{Segmentation}} and \code{\link{Centroids}} objects defining
#' spatially-resolved cellular segmentations or centroids
#' @slot assay (\code{\link[base:character]{character [1L]}}) A character naming
#' the associated assay of the spatial coordinates
#' @slot key (\code{\link[base:character]{character [1L]}}) The key for the
#' spatial coordinates
#'
#' @exportClass SpatialCoords
#'
#' @aliases SpatialCoords
#'
#' @seealso \code{\link{SpatialCoords-methods}}
#'
setClass(
  Class = 'SpatialCoords',
  contains = 'SpatialImage',
  slots = list(
    cells = 'LogMap',
    molecules = 'list',
    segmentations = 'list'
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' \code{SpatialCoords} Methods
#'
#' Methods for \code{\link{SpatialCoords}} objects
#'
#' @details The following methods are defined for interacting with a
#' \code{SpatialCoords} object:
#'
#' @param x,object A \code{\link{SpatialCoords}} object
#' @param layer Name of segmentation or molecular layer to extract cell or
#' feature names for; pass \code{NA} to return all cells or feature names
#' @param i,cells For \code{[[} and \code{[[<-}, the name of a segmentation or
#' \dQuote{molecules}; for \code{FetchData}, \code{subset}. and \code{[}, a
#' vector of cells to keep
#' @param j,features For \code{subset} and \code{[}, a vector of features to
#' keep; for \code{[[<-}, not used
#' @param value For \code{[[<-}, a replacement \code{\link{Molecules}},
#' \code{\link{Centroids}}, or \code{\link{Segmentations}} object or \code{NULL}
#' to remove the layer stored at \code{i}
#' @param ... Arguments passed to other methods
#'
#' @name SpatialCoords-methods
#' @rdname SpatialCoords-methods
#'
#' @seealso \code{\link{SpatialCoords-class}}
#'
NULL

#' @template method-cells
#'
#' @rdname SpatialCoords-methods
#' @method Cells SpatialCoords
#' @export
#'
Cells.SpatialCoords <- function(x, layer = NULL, ...) {
  layer <- layer[1L] %||% DefaultSegmentation(object = x)
  if (is.na(x = layer)) {
    return(Reduce(
      f = union,
      x = lapply(X = slot(object = x, name = 'segmentations'), FUN = Cells)
    ))
  }
  layer <- match.arg(arg = layer, choices = Segmentations(object = x))
  return(Cells(x = x[[layer]]))
}

#' @rdname CreateSpatialCoords
#' @method CreateSpatialCoords Centroids
#' @export
#'
CreateSpatialCoords.Centroids <- function(
  coords,
  molecules = NULL,
  assay = 'Spatial',
  key = NULL,
  name = NULL,
  ...
) {
  name <- name %||% as.character(x = tolower(x = class(x = coords)[1L]))
  coords <- list(coords)
  names(x = coords) <- name
  return(CreateSpatialCoords(
    coords = coords,
    molecules = molecules,
    assay = assay,
    key = key
  ))
}

#' @inheritParams CreateCentroids
#' @param type When providing a \code{\link[base]{data.frame}}, specify if the
#' coordinates represent a cell segmentation or voxel centroids
#' @param molecules A \code{\link[base]{data.frame}} with spatially-resolved
#' molecule information or a \code{\link{Molecules}} object
#' @param assay Name of associated assay
#' @param key Key for these spatial coordinates
#' @param name When \code{coords} is a \code{\link[base]{data.frame}},
#' \code{\link{Centroids}}, or \code{\link{Segmentation}}, name to store
#' coordinates as
#'
#' @rdname CreateSpatialCoords
#' @method CreateSpatialCoords data.frame
#' @export
#'
CreateSpatialCoords.data.frame <- function(
  coords,
  type = c('segmentation', 'centroids'),
  nsides = Inf,
  radius = NULL,
  theta = 0L,
  molecules = NULL,
  assay = 'Spatial',
  key = NULL,
  name = NULL,
  ...
) {
  type <- match.arg(arg = type)
  name <- name %||% type
  coords <- switch(
    EXPR = type,
    'segmentation' = CreateSegmentation(coords = coords),
    'centroids' = CreateCentroids(
      coords = coords,
      nsides = nsides,
      radius = radius,
      theta = theta
    )
  )
  return(CreateSpatialCoords(
    coords = coords,
    molecules = molecules,
    assay = assay,
    key = key
  ))
}

#'
#' @rdname CreateSpatialCoords
#' @method CreateSpatialCoords list
#' @export
#'
CreateSpatialCoords.list <- function(
  coords,
  molecules = NULL,
  assay = 'Spatial',
  key = NULL,
  ...
) {
  cells <- LogMap(y = unname(obj = unique(x = unlist(x = lapply(
    X = coords,
    FUN = Cells
  )))))
  for (i in names(x = coords)) {
    cells[[i]] <- Cells(x = coords[[i]])
    coords[[i]] <- subset(x = coords[[i]], cells = cells[[i]])
  }
  obj <- new(
    Class = 'SpatialCoords',
    cells = cells,
    segmentations = coords,
    molecules = molecules %iff% list(molecules = CreateMolecules(
      coords = molecules,
      key = 'mols_'
    )) %||% list(),
    assay = assay,
    key = key %||% Key(object = assay, quiet = TRUE)
  )
  validObject(object = obj)
  return(obj)
}

#' @rdname CreateSpatialCoords
#' @method CreateSpatialCoords Segmentation
#' @export
#'
CreateSpatialCoords.Segmentation <- CreateSpatialCoords.Centroids

#' @rdname Crop
#' @method Crop SpatialCoords
#' @export
#'
Crop.SpatialCoords <- function(
  object,
  x = NULL,
  y = NULL,
  coords = c("plot", "tissue"),
  ...
) {
  if (is.null(x = x) && is.null(x = y)) {
    return(object)
  }
  for (s in Segmentations(object = object, molecules = TRUE)) {
    object[[s]] <- Crop(object = object[[s]], x = x, y = y, coords = coords)
  }
  return(object)
}

#' @rdname Segmentations
#' @method DefaultSegmentation SpatialCoords
#' @export
#'
DefaultSegmentation.SpatialCoords <- function(object) {
  return(names(x = object)[1])
}

#' @rdname Segmentations
#' @method DefaultSegmentation<- SpatialCoords
#' @export
#'
"DefaultSegmentation<-.SpatialCoords" <- function(object, ..., value) {
  value <- match.arg(arg = value, choices = Segmentations(object = object))
  idx <- which(x = Segmentations(object = object) == value)
  norder <- c(
    idx,
    setdiff(x = seq_len(length.out = length(x = object)), y = idx)
  )
  slot(object = object, name = 'segmentations') <- slot(
    object = object,
    name = 'segmentations'
  )[norder]
  return(object)
}

#' @template method-features
#'
#' @rdname SpatialCoords-methods
#' @method Features SpatialCoords
#' @export
#'
Features.SpatialCoords <- function(x, layer = NULL, ...) {
  if (!length(x = Molecules(object = x))) {
    return(NULL)
  }
  layer <- layer[1L] %||% Molecules(object = x)[1L]
  if (is.na(x = layer)) {
    return(Reduce(
      f = union,
      x = lapply(X = slot(object = x, name = 'molecules'), FUN = Features)
    ))
  }
  layer <- match.arg(arg = layer, choices = Molecules(object = x))
  return(Features(x = x[[layer]]))
}

#' @param vars A vector of variables to fetch; can be the name of a
#' segmentation layer, to get tissue coordinates, or molecule names,
#' to get molecule coordinates
#' @param simplify If only returning either tissue or molecule coordinates,
#' return a single data frame instead of a list
#'
#' @details \code{FetchData}: Fetch tissue and/or molecule coordinates from
#' a \code{SpatialCoords} object
#'
#' @return \code{FetchData}: If both molecule and tissue coordinates are
#' requested, then a two-length list:
#' \itemize{
#'  \item \dQuote{\code{molecules}}: A data frame with the molecule coordinates
#'   requested. If molecules requested are keyed, the keys are preserved in the
#'   data frame
#'  \item \dQuote{\code{coordinates}}: A data frame with coordinates from the
#'   segmentation layers requested
#' }
#' If \code{simplify} is \code{TRUE} and only one data frame is generated, then
#' only the data frame is returned. Otherwise, a one-length list is returned
#' with the single data frame generated
#'
#' @rdname SpatialCoords-methods
#' @method FetchData SpatialCoords
#' @export
#'
FetchData.SpatialCoords <- function(
  object,
  vars,
  cells = NULL,
  simplify = TRUE,
  ...
) {
  vars.orig <- vars
  if (is.numeric(x = cells)) {
    cells <- Cells(x = object)[cells]
  }
  # Find keyed molecules
  object.keys <- Keys(object = object)
  keyed.mols <- sapply(
    X = object.keys,
    FUN = function(key) {
      return(grep(pattern = paste0('^', key), x = vars, value = TRUE))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  keyed.mols <- Filter(f = length, x = keyed.mols)
  mols.fetched <- sapply(
    X = names(x = keyed.mols),
    FUN = function(x) {
      df <- FetchData(object = object[[x]], vars = keyed.mols[[x]], ...)
      df$molecule <- paste0(Key(object = object[[x]]), df$molecule)
      return(df)
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  vars <- setdiff(
    x = vars,
    y = unique(x = lapply(
      X = mols.fetched,
      FUN = function(df) {
        return(unique(x = df$molecule))
      }
    ))
  )
  # Find all other molecules
  unkeyed.mols <- Filter(
    f = function(x) {
      return(x %in% Features(x = object, layer = NA))
    },
    x = vars
  )
  if (length(x = unkeyed.mols)) {
    mols.default <- Molecules(object = object)[1L]
    unkeyed.fetched <- FetchData(
      object = object[[mols.default]],
      vars = unkeyed.mols,
      ...
    )
    if (mols.default %in% names(x = mols.fetched)) {
      unkeyed.fetched$molecule <- paste0(
        Key(object = object[[mols.default]]),
        unkeyed.fetched$molecule
      )
      vars <- setdiff(x = vars, y = unique(x = unkeyed.mols))
    }
    mols.fetched <- append(x = mols.fetched, values = list(unkeyed.fetched))
  }
  # Assembled the molecules data frame
  mols.fetched <- do.call(what = 'rbind', args = mols.fetched)
  rownames(x = mols.fetched) <- NULL
  vars <- setdiff(x = vars, y = unique(x = mols.fetched$molecule))
  # Find all coordinates for the cells requested
  coords <- Filter(
    f = function(x) {
      return(x %in% Segmentations(object = object, molecules = FALSE))
    },
    x = vars
  )
  coords.fetched <- sapply(
    X = coords,
    FUN = function(x) {
      if (!is.null(x = cells) && !any(cells %in% Cells(x = object, layer = coords))) {
        return(NULL)
      }
      df <- GetTissueCoordinates(object = subset(x = object[[x]], cells = cells))
      df$layer <- x
      return(df)
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  coords.fetched <- do.call(what = 'rbind', args = coords.fetched)
  rownames(x = coords.fetched) <- NULL
  vars <- setdiff(x = vars, y = unique(x = coords.fetched$layer))
  # Warn/error about missing vars
  if (identical(x = vars, y = vars.orig)) {
    stop("Unable to find any of the provided vars", call. = FALSE)
  } else if (length(x = vars)) {
    warning(
      "The following vars were not found: ",
      paste(vars, collapse = ', '),
      call. = FALSE,
      immediate. = TRUE
    )
  }
  # Return fetched data
  data.fetched <- list(molecules = mols.fetched, coordinates = coords.fetched)
  data.fetched <- Filter(f = Negate(f = is.null), x = data.fetched)
  if (length(x = data.fetched) == 1L && isTRUE(x = simplify)) {
    return(data.fetched[[1L]])
  }
  return(data.fetched)
}

#' @details \code{GetTissueCoordinates}: Get cell or molecule  coordinates from
#' a \code{SpatialCoords} object
#'
#' @return \code{GetTissueCoordinates}: ...
#'
#' @rdname SpatialCoords-methods
#' @method GetTissueCoordinates SpatialCoords
#' @export
#'
GetTissueCoordinates.SpatialCoords <- function(object, which = NULL, ...) {
  which <- which %||% DefaultSegmentation(object = object)
  which <- match.arg(
    arg = which,
    choices = Segmentations(object = object, molecules = TRUE)
  )
  return(GetTissueCoordinates(object = object[[which]], ...))
}

#' @details \code{Keys}: Get the keys of molecule sets contained within a
#' \code{SpatialCoords} object
#'
#' @return \code{Keys}: A named vector of molecule set keys; names are the
#' names of the molecule sets and values are the keys for the respective
#' molecule set
#'
#' @rdname SpatialCoords-methods
#' @method Keys SpatialCoords
#' @export
#'
Keys.SpatialCoords <- function(object, ...) {
  return(sapply(X = slot(object = object, name = 'molecules'), FUN = Key))
}

#' @rdname Segmentations
#' @method Molecules SpatialCoords
#' @export
#'
Molecules.SpatialCoords <- function(object, ...) {
  return(names(x = slot(object = object, name = 'molecules')))
}

#' @param molecules Include \dQuote{molecules} as a segmentation layer, if
#' molecule information is present
#'
#' @rdname Segmentations
#' @method Segmentations SpatialCoords
#' @export
#'
Segmentations.SpatialCoords <- function(object, molecules = FALSE, ...) {
  segmentations <- names(x = object)
  if (isTRUE(x = molecules)) {
    segmentations <- c(segmentations, Molecules(object = object))
  }
  return(segmentations)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom utils .DollarNames
#'
#' @method .DollarNames SpatialCoords
#' @export
#'
.DollarNames.SpatialCoords <- function(x, pattern = '') {
  segmentations <- as.list(x = Segmentations(object = x, molecules = TRUE))
  names(x = segmentations) <- unlist(x = segmentations)
  return(.DollarNames(x = segmentations, pattern = pattern))
}

#' @rdname SpatialCoords-methods
#' @method $ SpatialCoords
#' @export
#'
"$.SpatialCoords" <- function(x, i, ...) {
  return(x[[i]])
}

#' @rdname SpatialCoords-methods
#' @method [ SpatialCoords
#' @export
#'
"[.SpatialCoords" <- function(x, i, j, ...) {
  if (missing(x = i)) {
    i <- NULL
  }
  if (missing(x = j)) {
    j <- NULL
  }
  return(subset(x = x, cells = i, features = j, ...))
}

#' @details \code{$}, \code{[[}: Extract a segmentation layer
#'
#' @return \code{$}, \code{[[}: The segmentation layer or spatially-resolved
#' molecule information stored at \code{i}
#'
#' @rdname SpatialCoords-methods
#' @method [[ SpatialCoords
#' @export
#'
"[[.SpatialCoords" <- function(x, i, ...) {
  i <- match.arg(arg = i, choices = Segmentations(object = x, molecules = TRUE))
  slot.use <- ifelse(
    test = i %in% Molecules(object = x),
    yes = 'molecules',
    no = 'segmentations'
  )
  return(slot(object = x, name = slot.use)[[i]])
}

#' @method dim SpatialCoords
#' @export
#'
dim.SpatialCoords <- function(x) {
  return(c(0, 0))
}

#' @details \code{length}: Get the number of segmentation layers in a
#' \code{SpatialCoords} object
#'
#' @return \code{length}: The number of segmentation layers
#' (\code{\link{Segmentation}} or \code{\link{Centroids}} objects)
#'
#' @rdname SpatialCoords-methods
#' @method length SpatialCoords
#' @export
#'
length.SpatialCoords <- function(x) {
  return(length(x = slot(object = x, name = 'segmentations')))
}

#' @details \code{names}: Get the names of segmentation layers, equivalent to
#' \code{\link[Spirula:Segmentations]{Segmentations(x, molecules = FALSE)}}
#'
#' @return \code{names}: A vector of segmentation layer names
#'
#' @rdname SpatialCoords-methods
#' @method names SpatialCoords
#' @export
#'
names.SpatialCoords <- function(x) {
  return(names(x = slot(object = x, name = 'segmentations')))
}

#' @details \code{subset}, \code{[}: Subset a \code{SpatialCoords} object
#'
#' @return \code{subset}: \code{x} with just the cells and features specified
#'
#' @rdname SpatialCoords-methods
#' @method subset SpatialCoords
#' @export
#'
#' @noRd
#'
subset.SpatialCoords <- function(x, cells = NULL, features = NULL, ...) {
  features <- Features(x = x) %iff% features
  if (is.null(x = cells) && is.null(x = features)) {
    return(x)
  }
  for (i in Molecules(object = x)) {
    x[[i]] <- subset(x = x[[i]], features = features)
  }
  for (i in Segmentations(object = x, molecules = FALSE)) {
    x[[i]] <- subset(x = x[[i]], cells = cells)
  }
  slot(object = x, name = 'cells') <- .PruneLogMap(x = slot(
    object = x,
    name = 'cells'
  ))
  validObject(object = x)
  return(x)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Add a Segmentation Layer
#'
#' @param x A \code{\link{SpatialCoords}} object
#' @param i Name to store segmentation layer as
#' @param ... Ignored
#' @param value A \code{\link{Segmentation}} or \code{\link{Centroids}} object
#' to add
#'
#' @return \code{x} with \code{value} saved as \code{i}
#'
#' @importFrom methods as
#'
#' @keywords internal
#'
#' @noRd
#'
.AddSegmentation <- function(x, i, ..., value) {
  if (i %in% Molecules(object = x)) {
    stop("'", i, "' already present as molecules")
  }
  # Check bounding box
  if (!.BboxIntersect(i = bbox(obj = value), j = bbox(obj = x), constraint = 'overlap')) {
    stop("New segmentation layer does not overlap with existing bounds")
  }
  # Reorder cells
  vcells <- MatchCells(
    new = Cells(x = value),
    orig = Cells(x = x),
    ordered = TRUE
  )
  vcells <- c(
    vcells,
    setdiff(
      x = seq.int(from = 1L, to = length(x = Cells(x = value))),
      y = vcells
    )
  )
  value <- value[vcells]
  # Check class
  if (i %in% names(x = x)) {
    same.class <- vapply(
      X = list(x[[i]], value),
      FUN = inherits,
      FUN.VALUE = logical(length = 1L),
      what = 'Segmentation'
    )
    if (length(x = unique(x = same.class)) != 1L) {
      warning(
        "Replacement value for ",
        i,
        " not of class ",
        class(x = x[[i]]),
        call. = FALSE,
        immediate. = TRUE
      )
    }
  }
  # Add segmentation layer
  slot(object = x, name = 'segmentations')[[i]] <- value
  # Update cells LogMap
  if (i %in% names(x = x)) {
    slot(object = x, name = 'cells')[[i]] <- Cells(x = value)
    slot(object = x, name = 'segmentations')[[i]] <- subset(
      x = x[[i]],
      cells = slot(object = x, name = 'cells')[[i]]
    )
  } else {
    vmap <- LogMap(y = Cells(x = value))
    vmap[[i]] <- Cells(x = value)
    cells <- merge(
      x = slot(object = x, name = 'cells'),
      y = vmap,
      by = 0,
      all = TRUE
    )
    rownames(x = cells) <- cells$Row.names
    cells <- cells[, setdiff(x = colnames(x = cells), y = 'Row.names')]
    cells <- as(object = as.matrix(x = cells), Class = 'LogMap')
    slot(object = x, name = 'cells') <- cells
  }
  # Validate and return
  validObject(object = x)
  return(x)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @details \code{[[<-}: Add or remove segmentation layers and molecule
#' information to/from a \code{SpatialCoords} object
#'
#' @return \code{[[<-}: Varies depending on the class of \code{value}:
#' \itemize{
#'  \item If \code{value} is \code{NULL}, returns \code{x} with the layer
#'  \code{i} removed; also allows removing \code{molecules}; does not allow
#'  removing the default segmentation
#'  \item If \code{value} is a \code{Molecules}, returns \code{x} with
#'  \code{value} stored in \code{molecules}; requires that \code{i} is
#'  \dQuote{molecules}
#'  \item Otherwise, stores \code{value} as a segmentation layer named \code{i}
#' }
#'
#' @rdname SpatialCoords-methods
#'
setMethod(
  f = '[[<-',
  signature = c(
    x = 'SpatialCoords',
    i = 'character',
    j = 'missing',
    value = 'Centroids'
  ),
  definition = .AddSegmentation
)

#' @rdname SpatialCoords-methods
#'
setMethod(
  f = '[[<-',
  signature = c(
    x = 'SpatialCoords',
    i = 'character',
    j = 'missing',
    value = 'Molecules'
  ),
  definition = function(x, i, ..., value) {
    if (i %in% Segmentations(object = x, molecules = FALSE)) {
      stop("'", i, "' already present as a segmentation layer")
    }
    check.key <- TRUE
    # Check bounding box for incoming molecules
    if (!.BboxIntersect(i = bbox(obj = value), j = bbox(obj = x), constraint = 'overlap')) {
      stop("New molecules do not overlap with existing bounds")
    }
    # TODO: Check replacement molecules
    if (i %in% Molecules(object = x)) {
      check.key <- Key(object = value) != Key(object = x[[i]])
    }
    if (isTRUE(x = check.key)) {
      if (Key(object = value) %in% Keys(object = x)) {
        key <- Key(object = i, quiet = TRUE)
        while (key %in% Keys(object = x)) {
          key <- Key(object = RandomName(), quiet = TRUE)
        }
        warning(
          "Duplicate moleculecular keys, changing to '",
          key, "'",
          call. = FALSE,
          immediate. = TRUE
        )
        Key(object = value) <- key
      }
    }
    # Add incoming molecules
    slot(object = x, name = 'molecules')[[i]] <- value
    # Validate and return
    validObject(object = x)
    return(x)
  }
)

#' @importFrom methods as
#'
#' @rdname SpatialCoords-methods
#'
setMethod(
  f = '[[<-',
  signature = c(
    x = 'SpatialCoords',
    i = 'character',
    j = 'missing',
    value = 'NULL'
  ),
  definition = function(x, i, ..., value) {
    i <- match.arg(arg = i, choices = Segmentations(object = x, molecules = TRUE))
    if (inherits(x = x[[i]], what = 'Molecules')) {
      slot(object = x, name = 'molecules')[[i]] <- NULL
    } else if (i == DefaultSegmentation(object = x)) {
      stop("Cannot remove default segmentation", call. = FALSE)
    } else {
      slot(object = x, name = 'segmentations')[[i]] <- NULL
      slot(object = x, name = 'cells')[[i]] <- NULL
      slot(object = x, name = 'cells') <- .PruneLogMap(x = slot(
        object = x,
        name = 'cells'
      ))
    }
    validObject(object = x)
    return(x)
  }
)

#' @rdname SpatialCoords-methods
#'
setMethod(
  f = '[[<-',
  signature = c(
    x = 'SpatialCoords',
    i = 'character',
    j = 'missing',
    value = 'Segmentation'
  ),
  definition = .AddSegmentation
)

setMethod(
  f = 'bbox',
  signature = 'SpatialCoords',
  definition = function(obj) {
    boxes <- lapply(X = slot(object = obj, name = 'segmentations'), FUN = bbox)
    boxes <- do.call(what = 'cbind', args = boxes)
    return(bbox(obj = t(x = boxes)))
  }
)

#' @template method-show
#'
#' @rdname SpatialCoords-methods
#'
setMethod(
  f = 'show',
  signature = c(object = 'SpatialCoords'),
  definition = function(object) {
    # Show cell information
    cat(
      "Spatial coordinates for",
      length(x = Cells(x = object, layer = NA)),
      "cells"
    )
    # Show molecule information
    if (length(x = Features(x = object, layer = NA))) {
      cat(" and", length(x = Features(x = object, layer = NA)), "molecules\n")
      cat(
        " First 10 molecules:",
        strwrap(x = paste(
          head(x = Features(x = object, layer = NA)),
          collapse = ', '
        ))
      )
    }
    cat("\n")
    # Show segmentation information
    cat("Default segmentation:", DefaultSegmentation(object = object), "\n")
    if (length(x = Segmentations(object = object)) > 1L) {
      segs <- setdiff(
        x = Segmentations(object = object),
        y = DefaultSegmentation(object = object)
      )
      cat(
        character(),
        length(x = segs),
        "other",
        ifelse(
          test = length(x = segs) == 1L,
          yes = "segmentation",
          no = "segmentations"
        ),
        "present:",
        strwrap(x = paste(segs, collapse = ', ')),
        "\n"
      )
    }
    # Show associated assay
    cat("Associated assay:", DefaultAssay(object = object), "\n")
    # Show key
    cat("Key:", Key(object = object), "\n")
    return(invisible(x = NULL))
  }
)

setValidity(
  Class = 'SpatialCoords',
  method = function(object) {
    valid <- NULL
    # Check cells
    cmap <- slot(object = object, name = 'cells')
    if (ncol(x = cmap) != length(x = slot(object = object, name = 'segmentations'))) {
      valid <- c(valid, "'cells' must have the same length as 'segmentations'")
    } else if (all(sort(x = colnames(x = cmap)) != sort(x = names(x = object)))) {
      valid <- c(valid, "'cells' must have the same names as 'segmentations'")
    }
    # Check segmentations
    nlist <- IsNamedList(
      x = slot(object = object, name = 'segmentations'),
      pass.zero = TRUE
    )
    if (!isTRUE(x = nlist)) {
      valid <- c(valid, "'segmentations' must be a named list")
    } else {
      for (s in Segmentations(object = object)) {
        if (!inherits(x = object[[s]], what = c('Segmentation', 'Centroids'))) {
          valid <- c(
            valid,
            "All segmentations must be either either a 'Segmentation' or 'Centroids' object"
          )
          break
        } else if (!identical(x = Cells(x = object[[s]]), y = cmap[[s]])) {
          valid <- c(
            valid,
            "All segmentations must have the cells specified in 'cells'"
          )
          break
        }
      }
    }
    # Check molecules
    nlist <- IsNamedList(
      x = slot(object = object, name = 'molecules'),
      pass.zero = TRUE
    )
    if (!isTRUE(x = nlist)) {
      valid <- c(valid, "'molecules' must be a named list")
    } else {
      for (m in Molecules(object = object)) {
        if (!inherits(x = object[[m]], what = 'Molecules')) {
          valid <- c(valid, "All molecules must inherit from 'Molecules'")
          break
        }
      }
    }
    return(valid %||% TRUE)
  }
)
