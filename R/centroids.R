#' @include zzz.R
#' @include generics.R
#' @importFrom methods as callNextMethod
#' @importClassesFrom sp SpatialPoints
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The \code{Centroids} Class
#'
#' @slot cells (\code{\link[base:character]{character [n]}}) A vector of cell
#' names; there should be as many cell names as there are points and no
#' duplicate names
#' @slot nsides (\code{\link[base:integer]{integer [1L]}}) The number of sides
#' to draw when plotting centroids; must be either \code{0L} for circles or
#' greater than 3
#' @slot radius (\code{\link[base:numeric]{numeric [1L]}}) The radius of the
#' shape when plotting the centroids
#' @slot theta (\code{\link[base:numeric]{numeric [1L]}}) The angle in degrees
#' to adjust the shape when plotting the centroids
#'
#' @family segmentation
#' @templateVar cls Centroids
#' @template seealso-methods
#'
setClass(
  Class = 'Centroids',
  contains = 'SpatialPoints',
  slots = list(
    cells = 'character',
    nsides = 'integer',
    radius = 'numeric',
    theta = 'numeric'
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' \code{Centroids} Methods
#'
#' Methods for \code{\link[SeuratObject:Centroids-class]{Centroids}} objects
#'
#' @param x,object A \code{\link[SeuratObject:Centroids-class]{Centroids}}
#' object
#' @param i,cells A vector of cells to keep; if \code{NULL}, defaults
#' to all cells
#' @param j,drop Ignored
#' @param ... Arguments passed to other methods
#'
#' @name Centroids-methods
#' @rdname Centroids-methods
#'
#' @seealso \code{\link{Centroids-class}}
#'
#' @family segmentation
#'
NULL

#'
#' @rdname Centroids-methods
#' @method Cells Centroids
#' @export
#'
Cells.Centroids <- function(x, ...) {
  return(slot(object = x, name = 'cells'))
}

#' @importFrom sp SpatialPoints
#'
#' @method CreateCentroids default
#' @export
#'
CreateCentroids.default <- function(
  coords,
  nsides = Inf,
  radius = NULL,
  theta = 0L
) {
  if (inherits(x = coords, what = 'sf')) {
    # Set the attribute-geometry relationship to constant
    # See https://r-spatial.github.io/sf/reference/sf.html#details
    sf::st_agr(coords) <- "constant"

    # Extract centroids from sf object
    centroids <- sf::st_centroid(coords)

    # Convert to data frame and format
    centroid_coords <- sf::st_coordinates(centroids)
    centroids_df <- data.frame(
      x = centroid_coords[, "X"],
      y = centroid_coords[, "Y"],
      row.names = coords$barcodes,
      stringsAsFactors = FALSE
    )
    coords <- centroids_df
  }
  cnames <- c('x', 'y')
  if (ncol(x = coords) >= 3) {
    cnames <- append(x = cnames, values = 'cell')
  }
  idx <- NameIndex(x = coords, names = cnames, MARGIN = 2L)
  cells <- if ('cell' %in% names(x = idx)) {
    as.character(x = coords[, idx[['cell']], drop = TRUE])
  } else {
    rownames(x = coords)
  }
  coords <- as.matrix(x = coords[, idx[c('x', 'y'), drop = FALSE]])
  colnames(x = coords) <- c('x', 'y')
  rownames(x = coords) <- NULL
  if (is.infinite(x = nsides)) {
    nsides <- 0L
  }
  radius <- radius %||% .AutoRadius(coords = coords)
  obj <- as(object = SpatialPoints(coords = coords), Class = 'Centroids')
  slot(object = obj, name = 'cells') <- cells
  slot(object = obj, name = 'nsides') <- as.integer(x = nsides)
  slot(object = obj, name = 'radius') <- as.numeric(x = radius)
  slot(object = obj, name = 'theta') <- as.numeric(x = theta)
  validObject(object = obj)
  return(obj)
}

#' @method CreateCentroids Centroids
#' @export
#'
CreateCentroids.Centroids <- function(
  coords,
  nsides = NULL,
  radius = NULL,
  theta = NULL
) {
  return(CreateCentroids(
    coords = GetTissueCoordinates(object = coords),
    nsides = nsides %||% length(x = coords),
    radius = radius %||% Radius(object = coords),
    theta = theta %||% Theta(object = coords)
  ))
}

#' @method Crop Centroids
#' @export
#'
Crop.Centroids <- .Crop

#' @details \code{GetTissueCoordinates}: Get cell spatial coordinates
#'
#' @param full Expand the coordinates to the full polygon
#'
#' @return \code{GetTissueCoordinates}: A data frame with three columns:
#' \itemize{
#'  \item \dQuote{\code{x}}: the x-coordinate
#'  \item \dQuote{\code{y}}: the y-coordinate
#'  \item \dQuote{\code{cell}}: the cell name
#' }
#' If \code{full} is \code{TRUE}, then each coordinate will indicate a vertex
#' for the cell polygon (created based on nsides, radius, and theta); 
#' otherwise, each coordinate will indicate a centroid for the cell.
#'
#' @importFrom sp coordinates
#'
#' @rdname Centroids-methods
#' @method GetTissueCoordinates Centroids
#' @export
#'
GetTissueCoordinates.Centroids <- function(object, full = TRUE, ...) {
  coords <- as.data.frame(x = coordinates(obj = object))
  colnames(x = coords) <- c('x', 'y')
  coords$cell <- Cells(x = object)
  if (isTRUE(x = full) && is.finite(x = object)) {
    ct <- mapply(
      FUN = PolyVtx,
      xc = coords$x,
      yc = coords$y,
      MoreArgs = list(
        n = length(x = object),
        r = Radius(object = object),
        t1 = Theta(object = object)
      )
    )
    xt <- vector(mode = 'list', length = length(x = ct) / 2L)
    xi <- 1L
    for (i in seq.int(from = 1L, to = length(x = ct), by = 2L)) {
      xt[[xi]] <- data.frame(
        x = ct[[i]],
        y = ct[[i + 1L]],
        cell = coords$cell[xi],
        stringsAsFactors = FALSE
      )
      xi <- xi + 1L
    }
    coords <- do.call(what = 'rbind', args = xt)
  }
  return(coords)
}

#' @details \code{Radius}: Get the centroid radius
#'
#' @return \code{Radius} The radius of the centroids
#'
#' @rdname Centroids-methods
#' @method Radius Centroids
#' @export
#'
Radius.Centroids <- function(object, ...) {
  return(slot(object = object, name = 'radius'))
}

#' @details \code{RenameCells}: Update cell names
#'
#' @inheritParams RenameCells
#'
#' @return \code{RenameCells}: \code{object} with the cells renamed to
#' \code{new.names}
#'
#' @rdname Centroids-methods
#' @method RenameCells Centroids
#' @export
#'
RenameCells.Centroids <- function(object, new.names = NULL, ...) {
  if (is.null(x = new.names)) {
    return(object)
  }
  new.names <- make.unique(names = new.names)
  if (length(x = new.names) != length(x = Cells(x = object))) {
    stop("Cannot partially rename centroid cells", call. = FALSE)
  }
  slot(object = object, name = 'cells') <- new.names
  return(object)
}


#' @details \code{Theta}: Get the offset angle
#'
#' @return \code{Theta}: The offset angle in degrees
#'
#' @rdname Centroids-methods
#' @method Theta Centroids
#' @export
#'
Theta.Centroids <- function(object) {
  return(slot(object = object, name = 'theta'))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @details \code{is.finite}, \code{is.infinite}: Test to see if the centroids
#' are circular or polygonal
#'
#' @return \code{is.finite}: \code{TRUE} if the centroids are polygonal,
#' \code{FALSE} if circular
#'
#' @rdname Centroids-methods
#' @method is.finite Centroids
#' @export
#'
is.finite.Centroids <- .FiniteCentroids

#' @return \code{is.infinite}: The opposite of \code{is.finite}
#'
#' @rdname Centroids-methods
#' @method is.infinite Centroids
#' @export
#'
is.infinite.Centroids <- Negate(f = .FiniteCentroids)

#' @details \code{length}: Get the number of sides for the polygonal centroid
#'
#' @return \code{length}: \code{0} if the centroids are circular, otherwise the
#' number of sides of the polygonal centroid
#'
#' @rdname Centroids-methods
#' @method length Centroids
#' @export
#'
length.Centroids <- function(x) {
  return(slot(object = x, name = 'nsides'))
}

#' @template method-lengths
#'
#' @rdname Centroids-methods
#' @method lengths Centroids
#' @export
#'
lengths.Centroids <- function(x, use.names = TRUE) {
  return(rle(x = Cells(x = x)))
}

#' @details \code{subset}, \code{[}: Subset a \code{Centroids} object to
#' certain cells
#'
#' @return \code{subset}, \code{[}: \code{x} subsetted to the cells specified
#' by \code{cells}/\code{i}
#'
#' @rdname Centroids-methods
#' @method subset Centroids
#' @export
#'
subset.Centroids <- function(x, cells = NULL, ...) {
  args <- list(...)
  if (is.null(x = cells)) {
    return(x)
  }
  if (is.numeric(x = cells)) {
    cells <- Cells(x = x)[cells]
    cells <- cells[!is.na(x = cells)]
  }
  cells <- MatchCells(new = Cells(x = x), orig = cells, ordered = TRUE)
  if (!length(x = cells)) {
    stop("None of the requested cells found")
  }
  return(CreateCentroids(
    coords = GetTissueCoordinates(object = x)[cells, ],
    nsides = args$nsides %||% length(x = x),
    radius = args$radius %||% Radius(object = x),
    theta = args$theta %||% Theta(object = x)
  ))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname Centroids-methods
#'
setMethod(
  f = '[',
  signature = c(x = 'Centroids', i = 'character'),
  definition = function(x, i, j, ..., drop = TRUE) {
    i <- MatchCells(new = Cells(x = x), orig = i, ordered = TRUE)
    return(x[i, drop = drop, ...])
  }
)

#' @rdname Centroids-methods
#'
setMethod(
  f = '[',
  signature = c(x = 'Centroids', i = 'numeric'),
  definition = function(x, i, j, ..., drop = TRUE) {
    info <- list(
      nsides = length(x = x),
      radius = Radius(object = x),
      theta = Theta(object = x)
    )
    cells <- Cells(x = x)[i]
    cells <- cells[!is.na(x = cells)]
    x <- callNextMethod()
    for (n in names(x = info)) {
      slot(object = x, name = n) <- info[[n]]
    }
    slot(object = x, name = 'cells') <- cells
    validObject(object = x)
    return(x)
  }
)

setMethod(
  f = 'over',
  signature = c(x = 'Centroids', y = 'SpatialPolygons'),
  definition = function(x, y, returnList = FALSE, fn = NULL, ...) {
    deprecate_stop(
      when = '5.0.0',
      what = 'over()',
      details = "Future integration with `sf` is on the roadmap with no current ETA"
    )
    check_installed(pkg = 'sf')
    return(over(
      x = as(object = x, Class = 'sf'),
      y = as(object = y, Class = 'sf'),
      returnList = returnList,
      fn = fn,
      ...
    ))
  }
)

#' @rdname Overlay
#' @export
#'
setMethod(
  f = 'Overlay',
  signature = c(x = 'Centroids', y = 'SpatialPolygons'),
  definition = function(x, y, invert = FALSE, ...) {
    check_installed(pkg = 'sf', reason = 'to overlay spatial information')
    idx <- sf::st_intersects(
      x = as(object = x, Class = 'sf'),
      y = as(object = y, Class = 'sf'),
      sparse = FALSE
    )
    idx <- which(idx)
    names_in_sf_object1 <- if (!is.null(x = row.names(x = x))) {
      row.names(x = x)[idx]
    } else {
      x$id[idx]
    }
    idx <- setNames(
      object = rep(x = TRUE, length(x = idx)),
      nm = names_in_sf_object1
    )
    if (!length(x = idx)) {
      warn(message = "The selected region does not contain any cell centroids")
      return(NULL)
    }
    idx <- sort(x = as.integer(x = names(x = idx)))
    if (isTRUE(x = invert)) {
      idx <- -idx
    }
    return(x[idx])
  }
)

#' @template method-show
#'
#' @rdname Centroids-methods
#'
setMethod(
  f = 'show',
  signature = c(object = 'Centroids'),
  definition = function(object) {
    cat("Spatial centroids for", length(x = Cells(x = object)), "cells\n")
    cat(
      "",
      ifelse(
        test = length(x = object) == 0L,
        yes = "Circular",
        no = paste0(length(x = object), '-sided')
      ),
      "spots\n"
    )
    cat(" Radius:", Radius(object = object), '\n')
    cat(" Offset angle:", Theta(object = object), "degrees\n")
    return(invisible(x = NULL))
  }
)

setValidity(
  Class = 'Centroids',
  method = function(object) {
    if (isFALSE(x = getOption(x = "Seurat.object.validate", default = TRUE))) {
      warn(
        message = paste("Not validating", class(x = object)[1L], "objects"),
        class = 'validationWarning'
      )
      return(TRUE)
    }
    valid <- NULL
    # Check cell names
    cells <- Cells(x = object)
    ucells <- Filter(f = nchar, x = unique(x = cells))
    if (length(x = ucells) != length(x = cells)) {
      valid <- c(
        valid,
        "'cells' must be a vector of unique cell names with one name for every coordinate"
      )
    }
    if (length(x = cells) != nrow(x = slot(object = object, name = 'coords'))) {
      valid <- c(
        valid,
        "the length of 'cells' must equal the number of rows in 'coords'"
      )
    }
    if (!is.null(x = rownames(x = slot(object = object, name = 'coords')))) {
      valid <- c(valid, "'coords' must not have any rownames")
    }
    # Check nsides
    nsides <- length(x = object)
    if (nsides < 0L || nsides %in% seq.int(from = 1L, to = 2L)) {
      valid <- c(
        valid,
        "'nsides' must be either 0 or greater than or equal to 3"
      )
    }
    # Check radius
    if (Radius(object = object) <= 0) {
      valid <- c(valid, "'radius' must be greater than 0")
    }
    # Check theta
    theta <- Theta(object = object)
    if (theta < 0 || theta > 360) {
      valid <- c(valid, "'theta must be between 0 and 360")
    }
    return(valid %||% TRUE)
  }
)
