#' @include zzz.R
#' @include generics.R
#' @include keymixin.R
#' @importFrom future.apply future_lapply
#' @importClassesFrom sp CRS SpatialPoints
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The Spatial Molecules Class
#'
#' @slot .Data A list of \code{\link[sp]{SpatialPoints}} objects
#' @slot key The key for the \code{Molecules}
#'
#' @family segmentation
#' @templateVar cls Molecules
#' @template seealso-methods
#'
setClass(
  Class = 'Molecules',
  contains = c('KeyMixin', 'list')
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' \code{Molecules} Methods
#'
#' Methods for \code{\link{Molecules}} objects
#'
#' @inheritParams Centroids-methods
#' @param x,object A \code{\link{Molecules}} object
#' @param features A vector of molecule names to keep; if \code{NULL}, defaults
#' to all molecules
#'
#' @name Molecules-methods
#' @rdname Molecules-methods
#'
#' @seealso \code{\link{Molecules-class}}
#'
NULL

#' @param key A key to set for the molecules
#'
#' @importFrom sp SpatialPoints
#'
#' @rdname CreateMolecules
#' @method CreateMolecules data.frame
#' @export
#'
CreateMolecules.data.frame <- function(coords, key = '', ...) {
  idx <- NameIndex(x = coords, names = c('x', 'y', 'gene'), MARGIN = 2L)
  xy <- idx[c('x', 'y')]
  coords <- split(x = coords, f = coords[[idx[['gene']]]])
  p <- progressor(steps = length(x = coords))
  coords <- sapply(
    X = coords,
    FUN = function(x) {
      p()
      mat <- as.matrix(x = x[, xy])
      rownames(x = mat) <- NULL
      return(SpatialPoints(coords = mat))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  obj <- new(
    Class = 'Molecules',
    .Data = coords,
    key = ifelse(test = nchar(x = key), yes = Key(object = key), no = key)
  )
  validObject(object = obj)
  return(obj)
}

#' @rdname CreateMolecules
#' @method CreateMolecules Molecules
#' @export
#'
CreateMolecules.Molecules <- function(coords, ...) {
  return(coords)
}

#' @rdname CreateMolecules
#' @method CreateMolecules NULL
#' @export
#'
CreateMolecules.NULL <- function(coords, ...) {
  return(NULL)
}

#' @method Crop Molecules
#' @export
#'
Crop.Molecules <- function(
  object,
  x = NULL,
  y = NULL,
  coords = c('plot',' tissue'),
  ...
) {
  object <- .Crop(object = object, x = x, y = y, coords = coords, ...)
  for (i in names(x = object)) {
    if (!nrow(x = slot(object = object[[i]], name = 'coords'))) {
      object[[i]] <- NULL
    }
  }
  if (!length(x = object)) {
    return(NULL)
  }
  return(object)
}

#' @template method-features
#'
#' @rdname Molecules-methods
#' @method Features Molecules
#' @export
#'
Features.Molecules <- function(x, ...) {
  return(names(x = x))
}

#' @method FetchData Molecules
#' @export
#'
FetchData.Molecules <- function(
  object,
  vars,
  nmols = NULL,
  seed = NA_integer_,
  ...
) {
  vars <- gsub(
    pattern = paste0('^', Key(object = object)),
    replacement = '',
    x = vars
  )
  coords <- GetTissueCoordinates(object = object, features = vars)
  if (!is.null(x = nmols)) {
    if (!is.na(x = seed)) {
      set.seed(seed = seed)
    }
    coords <- lapply(
      X = unique(x = coords$molecule),
      FUN = function(m) {
        df <- coords[coords$molecule == m, , drop = FALSE]
        if (nrow(x = df) > nmols) {
          idx <- sample(x = seq_len(length.out = nrow(x = df)), size = nmols)
          df <- df[idx, , drop = FALSE]
        }
        return(df)
      }
    )
    coords <- do.call(what = 'rbind', args = coords)
  }
  return(coords)
  # return(fortify(model = object, data = vars, ...))
}

#' @details \code{GetTissueCoordinates}: Get spatially-resolved
#' molecule coordinates
#'
#' @return \code{GetTissueCoordinates}: A data frame with three columns:
#' \itemize{
#'  \item \dQuote{\code{x}}: the x-coordinate of a molecule
#'  \item \dQuote{\code{y}}: the y-coordinate of a molecule
#'  \item \dQuote{\code{molecule}}: the molecule name
#' }
#'
#' @importFrom sp coordinates
#'
#' @rdname Molecules-methods
#' @method GetTissueCoordinates Molecules
#' @export
#'
GetTissueCoordinates.Molecules <- function(object, features = NULL, ...) {
  features <- features %||% Features(x = object)
  coords <- lapply(
    X = features,
    FUN = function(f) {
      fcoords <- object[[f]]
      if (is.null(x = fcoords)) {
        return(NULL)
      }
      fcoords <- as.data.frame(x = coordinates(obj = fcoords))
      rownames(x = fcoords) <- NULL
      fcoords$molecule <- f
      return(fcoords)
    }
  )
  return(do.call(what = 'rbind', args = coords))
}

#' @method Simplify Molecules
#' @export
#'
Simplify.Molecules <- function(coords, tol, topologyPreserve = TRUE) {
  .NotYetImplemented()
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom stats aggregate
#' @rdname aggregate
#' @method aggregate Molecules
#' @export
#'
aggregate.Molecules <- function(x, by, drop = TRUE, ...) {
  if (!inherits(x = by, what = 'SpatialPolygons')) {
    stop(
      "Aggregation of molecules works only with spatial polygons",
      call. = FALSE
    )
  }
  oob <- 'boundless'
  cells <- unname(obj = names(x = by))
  if (isFALSE(x = drop)) {
    if (oob %in% cells) {
      oob <- RandomName(length = 7L)
      warning(
        "'boundless' already present in cell names, changing to '",
        oob,
        "'",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    cells <- c(cells, oob)
  }
  idx <- over(x = x, y = by)
  m <- Matrix(
    data = 0,
    nrow = length(x = idx),
    ncol = length(x = cells),
    dimnames = list(Features(x = x), cells),
    sparse = TRUE
  )
  p <- progressor(along = idx)
  p(message = "Creating expression matrix", class = 'sticky', amount = 0)
  for (i in seq_along(along.with = idx)) {
    x <- idx[[i]]
    f <- names(x = idx)[i]
    if (isFALSE(x = drop)) {
      m[f, oob] <- sum(is.na(x = x))
    }
    x <- sort(x = unname(obj = x[!is.na(x = x)]))
    x <- rle(x = x)
    m[f, x$values] <- x$lengths
    p()
  }
  return(m)
}


#' @details \code{subset}: Subset a \code{Molecules} object to certain molecules
#'
#' @return \code{subset}: \code{x} subsetted to the features specified
#' by \code{features}
#'
#' @rdname Molecules-methods
#' @method subset Molecules
#' @export
#'
subset.Molecules <- function(x, features = NULL, ...) {
  features <- Features(x = x) %iff% features
  if (is.null(x = features)) {
    return(x)
  } else if (length(x = features) == 1L && is.na(x = features)) {
    return(CreateMolecules(coords = NULL))
  }
  if (is.numeric(x = features)) {
    features <- names(x = x)[features]
    features <- features[!is.na(x = features)]
  }
  features <- MatchCells(new = Features(x = x), orig = features, ordered = TRUE)
  if (!length(x = features)) {
    warning("None of the requested features found", immediate. = TRUE)
    return(CreateMolecules(coords = NULL))
  }
  x <- x[features]
  return(as(object = x, Class = 'Molecules'))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setMethod(
  f = 'bbox',
  signature = 'Molecules',
  definition = function(obj) {
    bounds <- lapply(X = obj, FUN = bbox)
    bounds <- do.call(what = 'cbind', args = bounds)
    bounds <- apply(X = bounds, MARGIN = 1L, FUN = range)
    rownames(x = bounds) <- c('min', 'max')
    return(t(x = bounds))
  }
)

setMethod(
  f = 'over',
  signature = c(x = 'Molecules', y = 'SpatialPolygons'),
  definition = function(x, y, returnList = FALSE, fn = NULL, ...) {
    f <- Features(x = x)
    p <- progressor(steps = length(x = f))
    out <- future_lapply(
      X = seq_along(along.with = f),
      FUN = function(i) {
        fi <- f[[i]]
        p()
        return(over(x = x[[fi]], y = y, returnList = returnList, fn = fn, ...))
      }
    )
    names(x = out) <- f
    return(out)
  }
)

#' @rdname Overlay
#' @export
#'
setMethod(
  f = 'Overlay',
  signature = c(x = 'Molecules', y = 'SpatialPolygons'),
  definition = function(x, y, invert = FALSE, ...) {
    idx <- over(x = x, y = y)
    for (f in names(x = idx)) {
      select <- idx[[f]]
      select <- select[!is.na(x = select)]
      cells <- as.integer(x = names(x = select))
      if (isTRUE(x = invert)) {
        cells <- -cells
      }
      x[[f]] <- x[[f]][cells]
    }
    for (i in names(x = x)) {
      if (!nrow(x = slot(object = x[[i]], name = 'coords'))) {
        x[[i]] <- NULL
      }
    }
    if (!length(x = x)) {
      return(NULL)
    }
    validObject(object = x)
    return(x)
  }
)

#' @template method-show
#'
#' @rdname Molecules-methods
#'
setMethod(
  f = 'show',
  signature = c(object = 'Molecules'),
  definition = function(object) {
    nmols <- length(x = object)
    molword <- ifelse(test = nmols == 1L, yes = 'molecule', no = 'molecules')
    cat("Coordinates for", nmols, molword, "\n")
    hmols <- min(nmols, 10)
    cat(
      "First",
      hmols,
      paste0(molword, ':'),
      paste(
        strwrap(x = paste(
          head(x = names(x = object), n = hmols),
          collapse = ', '
        )),
        collapse = '\n '
      ),
      '\n'
    )
  }
)

setValidity(
  Class = 'Molecules',
  method = function(object) {
    if (isFALSE(x = getOption(x = "Seurat.object.validate", default = TRUE))) {
      warn(
        message = paste("Not validating", class(x = object)[1L], "objects"),
        class = 'validationWarning'
      )
      return(TRUE)
    } else if (!length(x = object)) {
      return(TRUE)
    }
    valid <- NULL
    # Check names
    unnamed <- is.null(x = names(x = object)) ||
      length(x = Filter(f = nchar, x = names(x = object))) != length(x = object)
    if (isTRUE(x = unnamed)) {
      valid <- c(valid, 'All entries must be named')
    }
    # Check SpatialPoints
    points <- vapply(
      X = object,
      FUN = inherits,
      FUN.VALUE = logical(length = 1L),
      what = 'SpatialPoints'
    )
    if (!all(points)) {
      valid <- c(valid, 'All entries must inherit from sp::SpatialPoints')
    }
    # Check that coordinates have no rownames
    rnames <- vapply(
      X = object,
      FUN = function(x) {
        return(is.null(x = rownames(x = slot(object = x, name = 'coords'))))
      },
      FUN.VALUE = logical(length = 1L)
    )
    if (!all(rnames)) {
      valid <- c(valid, 'All entries must have unnamed coordinates')
    }
    return(valid %||% TRUE)
  }
)
