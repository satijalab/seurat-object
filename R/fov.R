#' @include zzz.R
#' @include generics.R
#' @include centroids.R
#' @include spatial.R
#' @include molecules.R
#' @include segmentation.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The Field of View Object
#'
#' A modern container for storing coordinates of spatially-resolved single
#' cells. Capable of storing multiple cell segmentation boundary masks.
#' Supports coordinates for spatially-resolved molecule (FISH) data.
#' Compatible with \code{\link{SpatialImage}}
#'
#' @slot molecules A named list of
#' \code{\link[SeuratObject:Molecules-class]{Molecules}} objects defining
#' spatially-resolved molecular coordinates
#' @slot boundaries  A named list of
#' \code{\link[SeuratObject:Segmentation-class]{Segmentation}} and
#' \code{\link[SeuratObject:Centroids-class]{Centroids}} objects defining
#' spatially-resolved boundaries
#' @slot coords_x_orientation A character indicating which axis 
#' \code{x} coordinates are associated with in spatial plots. 
#' Currently only applies to Visium objects. Ensures consistency in 
#' plotting spatial data across versions, as objects prior to the 
#' addition of this slot had \code{x} coordinates mapped to the vertical axis.
#' @slot assay A character naming the associated assay
#' of the spatial coordinates
#' @template slot-key
#'
#' @exportClass FOV
#'
#' @aliases FOV
#'
#' @concept fov
#'
#' @seealso \code{\link{FOV-methods}}
#'
setClass(
  Class = 'FOV',
  contains = 'SpatialImage',
  slots = list(
    molecules = 'list',
    boundaries = 'list',
    coords_x_orientation = 'character'
  ),
  prototype = list(
    coords_x_orientation = character(0)
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' \code{FOV} Methods
#'
#' Methods for \code{\link{FOV}} objects
#'
#' @details The following methods are defined for interacting with a
#' \code{FOV} object:
#'
#' @param x,object A \code{\link{FOV}} object
#' @param boundary,set Name of segmentation boundary or molecule set  to
#' extract cell or feature names for; pass \code{NA} to return all
#' cells or feature names
#' @param i,cells For \code{[[} and \code{[[<-}, the name of a segmentation or
#' \dQuote{molecules}; for \code{FetchData}, \code{subset}. and \code{[}, a
#' vector of cells to keep
#' @param j,features For \code{subset} and \code{[}, a vector of features to
#' keep; for \code{[[<-}, not used
#' @param value For \code{[[<-}, a replacement
#' \code{\link[SeuratObject:Molecules-class]{Molecules}},
#' \code{\link[SeuratObject:Centroids-class]{Centroids}}, or
#' \code{\link[SeuratObject:Segmentation-class]{Segmentation}} object;
#' otherwise \code{NULL} to remove the boundary stored at \code{i}
#' @param ... Arguments passed to other methods
#'
#' @name FOV-methods
#' @rdname FOV-methods
#'
#' @concept fov
#'
#' @seealso \code{\link{FOV-class}}
#'
NULL

#' @rdname Boundaries
#' @method Boundaries FOV
#' @export
#'
Boundaries.FOV <- function(object, ...) {
  return(names(x = slot(object = object, name = 'boundaries')))
}

#' @template method-cells
#'
#' @rdname FOV-methods
#' @method Cells FOV
#' @export
#'
Cells.FOV <- function(x, boundary = NULL, ...) {
  boundary <- boundary[1L] %||% DefaultBoundary(object = x)
  if (is.na(x = boundary)) {
    return(Reduce(
      f = union,
      x = lapply(X = slot(object = x, name = 'boundaries'), FUN = Cells)
    ))
  }
  boundary <- match.arg(arg = boundary, choices = Boundaries(object = x))
  return(Cells(x = x[[boundary]]))
}

#' @rdname CreateFOV
#' @method CreateFOV Centroids
#' @export
#'
CreateFOV.Centroids <- function(
  coords,
  molecules = NULL,
  assay = 'Spatial',
  key = NULL,
  misc = NULL,
  name = NULL,
  ...
) {
  name <- name %||% as.character(x = tolower(x = class(x = coords)[1L]))
  coords <- list(coords)
  names(x = coords) <- name
  return(CreateFOV(
    coords = coords,
    molecules = molecules,
    assay = assay,
    key = key,
    misc = misc %||% list()
  ))
}

#' @inheritParams CreateCentroids
#' @param type When providing a \code{\link[base]{data.frame}}, specify if
#' the coordinates represent a cell segmentation or voxel centroids
#' @param molecules A \code{\link[base]{data.frame}} with spatially-resolved
#' molecule information or a
#' \code{\link[SeuratObject:Molecules-class]{Molecules}} object
#' @param assay Name of associated assay
#' @param key Key for these spatial coordinates
#' @param name When \code{coords} is a \code{\link[base]{data.frame}},
#' \code{\link[SeuratObject:Centroids-class]{Centroids}}, or
#' \code{\link[SeuratObject:Segmentation-class]{Segmentation}}, name
#' to store coordinates as
#' @param misc A list of miscellaneous information to store with the object
#'
#' @rdname CreateFOV
#' @method CreateFOV data.frame
#' @export
#'
CreateFOV.data.frame <- function(
  coords,
  type = c('segmentation', 'centroids'),
  nsides = Inf,
  radius = NULL,
  theta = 0L,
  molecules = NULL,
  assay = 'Spatial',
  key = NULL,
  misc = NULL,
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
  return(CreateFOV(
    coords = coords,
    molecules = molecules,
    assay = assay,
    key = key,
    misc = misc %||% list()
  ))
}

#'
#' @rdname CreateFOV
#' @method CreateFOV list
#' @export
#'
CreateFOV.list <- function(
  coords,
  molecules = NULL,
  assay = 'Spatial',
  key = NULL,
  misc = NULL,
  ...
) {
  # Create a list of Molecules objects if provided; otherwise use an empty list
  molecules <- molecules %iff% list(molecules = CreateMolecules(
    coords = molecules,
    key = 'mols_'
  )) %||% list()
  # Create and validate the FOV object
  obj <- new(
    Class = 'FOV',
    boundaries = coords,
    molecules = molecules,
    assay = assay,
    key = key %||% Key(object = assay, quiet = TRUE),
    misc = misc %||% list()
  )
  return(obj)
}

#' @rdname CreateFOV
#' @method CreateFOV Segmentation
#' @export
#'
CreateFOV.Segmentation <- CreateFOV.Centroids

#' @rdname Crop
#' @method Crop FOV
#' @export
#'
Crop.FOV <- function(
  object,
  x = NULL,
  y = NULL,
  coords = c("plot", "tissue"),
  ...
) {
  if (is.null(x = x) && is.null(x = y)) {
    return(object)
  }
  for (s in names(x = object)) {
    object[[s]] <- Crop(object = object[[s]], x = x, y = y, coords = coords)
  }
  return(object)
}

#' @rdname Boundaries
#' @method DefaultBoundary FOV
#' @export
#'
DefaultBoundary.FOV <- function(object) {
  return(Boundaries(object = object)[1])
}

#' @rdname Boundaries
#' @method DefaultBoundary<- FOV
#' @export
#'
"DefaultBoundary<-.FOV" <- function(object, ..., value) {
  value <- match.arg(arg = value, choices = Boundaries(object = object))
  idx <- which(x = Boundaries(object = object) == value)
  norder <- c(
    idx,
    setdiff(x = seq_len(length.out = length(x = object)), y = idx)
  )
  slot(object = object, name = 'boundaries') <- slot(
    object = object,
    name = 'boundaries'
  )[norder]
  return(object)
}

#' @template method-features
#'
#' @rdname FOV-methods
#' @method Features FOV
#' @export
#'
Features.FOV <- function(x, set = NULL, ...) {
  if (!length(x = Molecules(object = x))) {
    return(NULL)
  }
  set <- set[1L] %||% Molecules(object = x)[1L]
  if (is.na(x = set)) {
    return(Reduce(
      f = union,
      x = lapply(X = slot(object = x, name = 'molecules'), FUN = Features)
    ))
  }
  set <- match.arg(arg = set, choices = Molecules(object = x))
  return(Features(x = x[[set]]))
}

#' @param vars A vector of variables to fetch; can be the name of a
#' segmentation boundary, to get tissue coordinates, or molecule names,
#' to get molecule coordinates
#' @param simplify If only returning either boundary or molecule coordinates,
#' return a single data frame instead of a list
#'
#' @details \code{FetchData}: Fetch boundary and/or molecule coordinates from
#' a \code{FOV} object
#'
#' @return \code{FetchData}: If both molecule and boundary coordinates are
#' requested, then a two-length list:
#' \itemize{
#'  \item \dQuote{\code{molecules}}: A data frame with the molecule coordinates
#'   requested. If molecules requested are keyed, the keys are preserved in the
#'   data frame
#'  \item \dQuote{\code{coordinates}}: A data frame with coordinates from the
#'   segmentation boundaries requested
#' }
#' If \code{simplify} is \code{TRUE} and only one data frame is generated, then
#' only the data frame is returned. Otherwise, a one-length list is returned
#' with the single data frame generated
#'
#' @rdname FOV-methods
#' @method FetchData FOV
#' @export
#'
FetchData.FOV <- function(
  object,
  vars,
  cells = NULL,
  simplify = TRUE,
  ...
) {
  vars.orig <- vars
  if (is.numeric(x = cells)) {
    cells <- Cells(x = object)[cells]
  } else if (is.null(cells)) {
    cells <- Cells(x = object)
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
      return(x %in% Features(x = object, set = NA))
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
      return(x %in% Boundaries(object = object))
    },
    x = vars
  )
  coords.fetched <- sapply(
    X = coords,
    FUN = function(x) {
      if (!is.null(x = cells) && !any(cells %in% Cells(x = object, boundary = coords))) {
        return(NULL)
      }
      df <- GetTissueCoordinates(object = subset(x = object[[x]], cells = cells))
      df$boundary <- x
      return(df)
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  coords.fetched <- do.call(what = 'rbind', args = coords.fetched)
  rownames(x = coords.fetched) <- NULL
  vars <- setdiff(x = vars, y = unique(x = coords.fetched$boundary))
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

#' @param which Name of segmentation boundary or molecule set to retrieve coordinates for;
#' if NULL, will retrieve coordinates for the default boundary
#'
#' @examples
#' \dontrun{
#' GetTissueCoordinates(object, which = "centroids")
#' }
#' @rdname GetTissueCoordinates
#' @method GetTissueCoordinates FOV
#' @export
#'
GetTissueCoordinates.FOV <- function(object, which = NULL, ...) {
  which <- which %||% DefaultBoundary(object = object)
  which <- match.arg(arg = which, choices = names(x = object))
  return(GetTissueCoordinates(object = object[[which]], ...))
}

#' @details \code{Keys}: Get the keys of molecule sets contained within a
#' \code{FOV} object
#'
#' @return \code{Keys}: A named vector of molecule set keys; names are the
#' names of the molecule sets and values are the keys for the respective
#' molecule set
#'
#' @rdname FOV-methods
#' @method Keys FOV
#' @export
#'
Keys.FOV <- function(object, ...) {
  return(sapply(X = slot(object = object, name = 'molecules'), FUN = Key))
}

#' @rdname Boundaries
#' @method Molecules FOV
#' @export
#'
Molecules.FOV <- function(object, ...) {
  return(names(x = slot(object = object, name = 'molecules')))
}

#' @details \code{RenameCells}: Update cell names
#'
#' @inheritParams RenameCells
#'
#' @return \code{RenameCells}: \code{object} with the cells renamed to
#' \code{new.names}
#'
#' @rdname FOV-methods
#' @method RenameCells FOV
#' @export
#'
RenameCells.FOV <- function(object, new.names = NULL, ...) {
  if (is.null(x = new.names)) {
    return(object)
  }
  new.names <- make.unique(names = new.names)
  all.cells <- Cells(x = object, boundary = NA)
  if (length(x = new.names) != length(x = all.cells)) {
    stop("Cannot partially rename cells", call. = FALSE)
  }
  for (boundary in Boundaries(object = object)) {
    idx <- MatchCells(
      new = all.cells,
      orig = Cells(x = object[[boundary]]),
      ordered = TRUE
    )
    if (!length(x = idx)) {
      next
    }
    object[[boundary]] <- RenameCells(
      object = object[[boundary]],
      new.names = new.names[idx]
    )
  }
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom utils .DollarNames
#'
#' @method .DollarNames FOV
#' @export
#'
.DollarNames.FOV <- function(x, pattern = '') {
  layers <- as.list(x = names(x = x))
  names(x = layers) <- unlist(x = layers)
  return(.DollarNames(x = layers, pattern = pattern))
}

#' @rdname FOV-methods
#' @method $ FOV
#' @export
#'
"$.FOV" <- function(x, i, ...) {
  return(x[[i]])
}

#' @rdname FOV-methods
#' @method [ FOV
#' @export
#'
"[.FOV" <- function(x, i, j, ...) {
  if (missing(x = i)) {
    i <- NULL
  }
  if (missing(x = j)) {
    j <- NULL
  }
  return(subset(x = x, cells = i, features = j, ...))
}

#' @details \code{$}, \code{[[}: Extract a segmentation boundary
#'
#' @return \code{$}, \code{[[}: The segmentation boundary or spatially-resolved
#' molecule information stored at \code{i}
#'
#' @rdname FOV-methods
#' @method [[ FOV
#' @export
#'
"[[.FOV" <- function(x, i, ...) {
  i <- match.arg(arg = i, choices = names(x = x))
  slot.use <- ifelse(
    test = i %in% Molecules(object = x),
    yes = 'molecules',
    no = 'boundaries'
  )
  return(slot(object = x, name = slot.use)[[i]])
}

#' Aggregate Molecules into an Expression Matrix
#'
#' @param x An object with spatially-resolved molecule information
#' @param by Name of a
#' \code{\link[SeuratObject:Segmentation-class]{Segmentation}} within
#' \code{object} or a
#' \code{\link[SeuratObject:Segmentation-class]{Segmentation}} object
#' @param set Name of molecule set to aggregate
#' @param drop Drop molecules not present in a segmentation; if \code{FALSE},
#' adds a column called \dQuote{\code{boundless}} consisting of molecule counts
#' not in a segmentation
#' @param ... Arguments passed to other methods
#'
#' @return An expression matrix
#'
#' @importFrom stats aggregate
#'
#' @name aggregate
#' @rdname aggregate
#'
#' @keywords internal
#'
#' @method aggregate FOV
#' @export
#'
#' @template section-progressr
#' @template section-future
#'
#' @order 1
#'
aggregate.FOV <- function(x, by = NULL, set = NULL, drop = TRUE, ...) {
  # Check molecules
  set <- set[1L] %||% Molecules(object = x)[1L]
  if (is.null(x = set)) {
    stop("No molecules present in this FOV", call. = FALSE)
  }
  set <- match.arg(arg = set, choices = Molecules(object = x))
  # Check segmentation boundaries
  by <- by[1L] %||% Filter(
    f = function(b) {
      return(inherits(x = x[[b]], what = 'Segmentation'))
    },
    x = Boundaries(object = x)
  )[1L]
  if (is.character(x = by)) {
    by <- x[[by]]
  }
  if (!inherits(x = by, what = 'SpatialPolygons')) {
    stop("'by' is not a segmentation boundary", call. = FALSE)
  }
  # TODO: Check bbox intersect
  # Aggregate
  return(aggregate(x = x[[set]], by = by, drop = drop, ...))
}

#' @method dim FOV
#' @export
#'
dim.FOV <- function(x) {
  return(c(0, 0))
}

#' @details \code{length}: Get the number of segmentation layers in a
#' \code{FOV} object
#'
#' @return \code{length}: The number of segmentation layers
#' (\code{\link[SeuratObject:Segmentation-class]{Segmentation}} or
#' \code{\link[SeuratObject:Centroids-class]{Centroids}} objects)
#'
#' @rdname FOV-methods
#' @method length FOV
#' @export
#'
length.FOV <- function(x) {
  return(length(x = slot(object = x, name = 'boundaries')))
}

#' @details \code{names}: Get the names of segmentation layers and molecule sets
#'
#' @return \code{names}: A vector of segmentation boundary and molecule set names
#'
#' @rdname FOV-methods
#' @method names FOV
#' @export
#'
names.FOV <- function(x) {
  return(c(Boundaries(object = x), Molecules(object = x)))
}

#' @details \code{subset}, \code{[}: Subset a \code{FOV} object
#'
#' @return \code{subset}: \code{x} with just the cells and features specified
#'
#' @rdname FOV-methods
#' @method subset FOV
#' @export
#'
subset.FOV <- function(x, cells = NULL, features = NULL, ...) {
  features <- Features(x = x) %iff% features
  if (is.null(x = cells) && is.null(x = features)) {
    return(x)
  }
  for (i in Molecules(object = x)) {
    x[[i]] <- subset(x = x[[i]], features = features)
  }
  if (is.numeric(x = cells)) {
    cells <- Cells(x = x, boundary = NA)[cells]
  }
  for (i in Boundaries(object = x)) {
    x[[i]] <- subset(x = x[[i]], cells = cells)
  }
  safeValidityCheck(object = x)
  return(x)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Add a Segmentation Boundary
#'
#' @param x A \code{\link{FOV}} object
#' @param i Name to store segmentation boundary as
#' @param ... Ignored
#' @param value A \code{\link[SeuratObject:Segmentation-class]{Segmentation}}
#' or [SeuratObject:Centroids-class]\code{\link{Centroids}} object
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
    stop("'", i, "' already present as molecules", call. = FALSE)
  }
  # Check bounding box
  if (!.BboxIntersect(i = bbox(obj = value), j = bbox(obj = x), constraint = 'overlap')) {
    stop(
      "New segmentation boundary does not overlap with existing bounds",
      call. = FALSE
    )
  }
  # # Reorder cells
  # vcells <- MatchCells(
  #   new = Cells(x = value),
  #   orig = Cells(x = x, boundary = NA),
  #   ordered = TRUE
  # )
  # vcells <- c(
  #   vcells,
  #   setdiff(
  #     x = seq.int(from = 1L, to = length(x = Cells(x = value))),
  #     y = vcells
  #   )
  # )
  # value <- value[vcells]
  # Check class
  if (i %in% Boundaries(object = x)) {
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
  # Add segmentation boundary
  slot(object = x, name = 'boundaries')[[i]] <- value
  # Reorder cells
  x <- .OrderCells(object = x)
  # Validate and return
  safeValidityCheck(object = x)
  return(x)
}

#' Order cells in an FOV
#'
#' @param object An \code{\link[SeuratObject:FOV-class]{FOV}} object
#'
#' @return \code{object} with the cells in each boundary ordered
#'
#' @keywords internal
#'
#' @noRd
#'
.OrderCells <- function(object) {
  all.cells <- Cells(x = object, boundary = NA)
  for (b in Boundaries(object = object)) {
    bcells <- MatchCells(
      new = Cells(x = object[[b]]),
      orig = all.cells,
      ordered = TRUE
    )
    bcells <- c(
      bcells,
      setdiff(x = seq_along(along.with = Cells(x = object[[b]])), y = bcells)
    )
    slot(object = object, name = 'boundaries')[[b]] <- object[[b]][bcells]
  }
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @details \code{[[<-}: Add or remove segmentation layers and molecule
#' information to/from a \code{FOV} object
#'
#' @return \code{[[<-}: Varies depending on the class of \code{value}:
#' \itemize{
#'  \item If \code{value} is \code{NULL}, returns \code{x} with the boundary
#'  \code{i} removed; also allows removing \code{molecules}; does not allow
#'  removing the default segmentation
#'  \item If \code{value} is a \code{Molecules}, returns \code{x} with
#'  \code{value} stored in \code{molecules}; requires that \code{i} is
#'  \dQuote{molecules}
#'  \item Otherwise, stores \code{value} as a segmentation boundary named \code{i}
#' }
#'
#' @rdname FOV-methods
#'
setMethod(
  f = '[[<-',
  signature = c(
    x = 'FOV',
    i = 'character',
    j = 'missing',
    value = 'Centroids'
  ),
  definition = .AddSegmentation
)

#' @rdname FOV-methods
#'
setMethod(
  f = '[[<-',
  signature = c(
    x = 'FOV',
    i = 'character',
    j = 'missing',
    value = 'Molecules'
  ),
  definition = function(x, i, ..., value) {
    if (i %in% Boundaries(object = x)) {
      stop("'", i, "' already present as a segmentation boundary")
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
    safeValidityCheck(object = x)
    return(x)
  }
)

#' @importFrom methods as
#'
#' @rdname FOV-methods
#'
setMethod(
  f = '[[<-',
  signature = c(
    x = 'FOV',
    i = 'character',
    j = 'missing',
    value = 'NULL'
  ),
  definition = function(x, i, ..., value) {
    i <- match.arg(arg = i, choices = names(x = x))
    if (inherits(x = x[[i]], what = 'Molecules')) {
      slot(object = x, name = 'molecules')[[i]] <- NULL
    } else if (i == DefaultBoundary(object = x)) {
      stop("Cannot remove default boundary", call. = FALSE)
    } else {
      slot(object = x, name = 'boundaries')[[i]] <- NULL
    }
    safeValidityCheck(object = x)
    return(x)
  }
)

#' @rdname FOV-methods
#'
setMethod(
  f = '[[<-',
  signature = c(
    x = 'FOV',
    i = 'character',
    j = 'missing',
    value = 'Segmentation'
  ),
  definition = .AddSegmentation
)

setMethod(
  f = 'bbox',
  signature = 'FOV',
  definition = function(obj) {
    boxes <- lapply(X = slot(object = obj, name = 'boundaries'), FUN = bbox)
    boxes <- do.call(what = 'cbind', args = boxes)
    return(bbox(obj = t(x = boxes)))
  }
)

#' @importFrom methods initialize
#'
setMethod(
  f = 'initialize',
  signature = 'FOV',
  definition = function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
    .Object <- .OrderCells(object = .Object)
    safeValidityCheck(object = .Object)
    return(.Object)
  }
)

#' @importClassesFrom sp Spatial
#' @rdname Overlay
#'
setMethod(
  f = 'Overlay',
  signature = c(x = 'FOV', y = 'Spatial'),
  definition = .OverBbox
)

#' @rdname Overlay
#'
setMethod(
  f = 'Overlay',
  signature = c(x = 'FOV', y = 'SpatialPolygons'),
  definition = function(x, y, invert = FALSE, ...) {
    for (i in names(x = x)) {
      x[[i]] <- Overlay(x = x[[i]], y = y, invert = invert, ...)
    }
    return(x)
  }
)

#' @rdname Overlay
#'
setMethod(
  f = 'Overlay',
  signature = c(x = 'FOV', y = 'FOV'),
  definition = .OverBbox
)

#' @template method-show
#'
#' @rdname FOV-methods
#'
setMethod(
  f = 'show',
  signature = c(object = 'FOV'),
  definition = function(object) {
    # Show cell information
    cat(
      "Spatial coordinates for",
      length(x = Cells(x = object, boundary = NA)),
      "cells"
    )
    # Show molecule information
    if (length(x = Features(x = object, boundary = NA))) {
      cat(" and", length(x = Features(x = object, boundary = NA)), "molecules\n")
      cat(
        " First 10 molecules:",
        strwrap(x = paste(
          head(x = Features(x = object, boundary = NA)),
          collapse = ', '
        ))
      )
    }
    cat("\n")
    # Show segmentation information
    cat(
      "Default segmentation boundary:",
      DefaultBoundary(object = object),
      "\n"
    )
    if (length(x = Boundaries(object = object)) > 1L) {
      segs <- setdiff(
        x = Boundaries(object = object),
        y = DefaultBoundary(object = object)
      )
      cat(
        character(),
        length(x = segs),
        "other segmentation boundaries present:",
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

#' FOV Validity
#'
#' @templateVar cls FOV
#' @template desc-validity
#'
#' @section Boundary Validation:
#' blah
#'
#' @section Molecule Validation:
#' blah
#'
#' @name FOV-validity
#'
#' @family fov
#'
#' @seealso \code{\link[methods]{validObject}}
#'
setValidity(
  Class = 'FOV',
  method = function(object) {
    if (isFALSE(x = getOption(x = "Seurat.object.validate", default = TRUE))) {
      warn(
        message = paste("Not validating", class(x = object)[1L], "objects"),
        class = 'validationWarning'
      )
      return(TRUE)
    }
    valid <- NULL
    # Check boundaries
    nlist <- IsNamedList(
      x = slot(object = object, name = 'boundaries'),
      pass.zero = TRUE
    )
    if (!isTRUE(x = nlist)) {
      valid <- c(valid, "'boundaries' must be a named list")
    } else {
      all.cells <- Cells(x = object, boundary = NA)
      for (s in Boundaries(object = object)) {
        if (!inherits(x = object[[s]], what = c('Segmentation', 'Centroids'))) {
          valid <- c(
            valid,
            "All segmentation boundaries must be either either a 'Segmentation' or 'Centroids' object"
          )
          break
        } else {
          cells <- Cells(x = object[[s]])
          if (!is.null(cells)) {
            matched.cells <- MatchCells(
              new = all.cells,
              orig = cells,
              ordered = TRUE
            )
            if (length(x = matched.cells) != length(x = Cells(x = object[[s]]))) {
              valid <- c(
                valid,
                "All segmentation boundaries must have cells"
              )
              break
            }
          } else {
            valid <- c(
              valid,
              paste(s, "contains 0 cells")
            )
            break
          }
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
