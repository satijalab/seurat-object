#' @include zzz.R
#' @include generics.R
#' @importFrom sp coordinates
#' @importFrom methods as callNextMethod
#' @importClassesFrom sp SpatialPolygons

NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The \code{Segmentation} Class
#'
#' A container for cell segmentation boundaries.
#' Inherits from \code{\link[sp:SpatialPolygons-class]{SpatialPolygons}}.
#' Supports storing boundaries in objects of class \code{\link[sf]{sf}}.
#'
#' @slot sf.data Segmentation boundaries in \code{\link[sf]{sf}} format
#' @slot is.lightweight Logical indicating whether ot not the object only stores 
#' segmentation information in the \code{sf.data} slot
#'
#' @family segmentation
#' @templateVar cls Segmentation
#' @template seealso-methods
#'
setClass(
  Class = 'Segmentation',
  contains = 'SpatialPolygons',
  slots = list(
    sf.data = 'ANY',
    is.lightweight = 'OptionalLogical'
  ),
  prototype = list(
    is.lightweight = FALSE
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' \code{Segmentation} Methods
#'
#' Methods for \code{\link[SeuratObject:Segmentation-class]{Segmentation}}
#' objects
#'
#' @inheritParams Centroids-methods
#' @param x,object,obj A
#' \code{\link[SeuratObject:Segmentation-class]{Segmentation}} object
#'
#' @name Segmentation-methods
#' @rdname Segmentation-methods
#'
#' @section Progress Updates with \pkg{progressr}:
#' The following methods use
#' \href{https://cran.r-project.org/package=progressr}{\pkg{progressr}} to
#' render status updates and progress bars:
#' \itemize{
#'  \item \code{RenameCells}
#' }
#' To enable progress updates, wrap
#' the function call in \code{\link[progressr]{with_progress}} or run
#' \code{\link[progressr:handlers]{handlers(global = TRUE)}} before running
#' this function. For more details about \pkg{progressr}, please read
#' \href{https://progressr.futureverse.org/articles/progressr-intro.html}{\code{vignette("progressr-intro")}}

#'
#' @section Parallelization with \pkg{future}:
#' The following methods use
#' \href{https://cran.r-project.org/package=future}{\pkg{future}} to enable
#' parallelization:
#' \itemize{
#'  \item \code{RenameCells}
#' }
#' Parallelization strategies can be set using
#' \code{\link[future]{plan}}. Common plans include \dQuote{\code{sequential}}
#' for non-parallelized processing or \dQuote{\code{multisession}} for parallel
#' evaluation using multiple \R sessions; for other plans, see the
#' \dQuote{Implemented evaluation strategies} section of
#' \code{\link[future:plan]{?future::plan}}. For a more thorough introduction
#' to \pkg{future}, see
#' \href{https://future.futureverse.org/articles/future-1-overview.html}{\code{vignette("future-1-overview")}}
#'
#' @concept future
#'
#' @seealso \code{\link{Segmentation-class}}
#'
#' @family segmentation
#'
NULL

#' @template method-cells
#'
#' @rdname Segmentation-methods
#' @method Cells Segmentation
#' @export
#'
Cells.Segmentation <- function(x, ...) {
  if (slot(object = x, name = 'is.lightweight')) {
    sf_data <- slot(object = x, name = 'sf.data')
    return(unname(obj = sf_data$cell))
  }
  return(unname(obj = names(x = x)))
}

#' @importFrom sp Polygon Polygons SpatialPolygons
#'
#' @rdname CreateSegmentation
#' @method CreateSegmentation data.frame
#' @export
#'
CreateSegmentation.data.frame <- function(coords, lightweight = FALSE) {
  if (lightweight) {
    # Create minimal valid SpatialPolygons structure to satisfy inheritance
    minimal_coords <- matrix(c(0, 1, 1, 1, 1, 0, 0, 0, 0, 1), ncol = 2, byrow = TRUE)
    minimal_polygon <- Polygons(
      srl = list(Polygon(coords = minimal_coords)),
      ID = "placeholder"
    )

    sp_base <- SpatialPolygons(list(minimal_polygon))

    # Now create Segmentation object with valid SpatialPolygons inheritance
    obj <- new(
      Class = 'Segmentation',
      sp_base,
      sf.data = coords
    )

    # Override with empty polygons for lightweight mode
    slot(obj, 'polygons') <- list()
    slot(obj, 'plotOrder') <- integer(0)
    slot(obj, 'proj4string') <- CRS(as.character(NA))

    # Get bbox from coordinate ranges in the dataframe
    x_range <- range(coords$x, na.rm = TRUE)
    y_range <- range(coords$y, na.rm = TRUE)
    slot(obj, 'bbox') <- matrix(c(x_range[1], y_range[1], x_range[2], y_range[2]), 
                                nrow = 2,
                                ncol = 2, 
                                dimnames = list(c("x", "y"), c("min", "max")))

    slot(obj, 'is.lightweight') <- TRUE
    return(obj)
  }
  idx <- NameIndex(x = coords, names = c('cell', 'x', 'y'), MARGIN = 2L)
  xy <- idx[c('x', 'y')]
  cell.idx <- idx[['cell']]
  coords <- split(x = coords, f = coords[[cell.idx]])
  coords <- sapply(
    X = coords,
    FUN = function(x) {
      cx <- as.matrix(x = x[, xy])
      colnames(x = cx) <- c('x', 'y')
      return(Polygons(
        srl = list(Polygon(coords = cx)),
        ID = unique(x = as.character(x = x[[cell.idx]]))
      ))
    }
  )
  coords <- SpatialPolygons(Srl = coords)
  CheckGC()
  return(as(object = coords, Class = 'Segmentation'))
}

#' @rdname CreateSegmentation
#' @method CreateSegmentation Segmentation
#' @export
#'
CreateSegmentation.Segmentation <- function(coords) {
  return(coords)
}

#' @rdname CreateSegmentation
#' @method CreateSegmentation sf
#' @export
#'
CreateSegmentation.sf <- function(coords, lightweight = FALSE) {
  # Method is called when creating Segmentation from an sf object
  if (lightweight) {
    # Create minimal valid SpatialPolygons structure to satisfy inheritance
    minimal_coords <- matrix(c(0, 1, 1, 1, 1, 0, 0, 0, 0, 1), ncol = 2, byrow = TRUE)
    minimal_polygon <- Polygons(
      srl = list(Polygon(coords = minimal_coords)),
      ID = "placeholder"
    )

    sp_base <- SpatialPolygons(list(minimal_polygon))

    # Now create Segmentation object with valid SpatialPolygons inheritance
    obj <- new(
      Class = 'Segmentation',
      sp_base,
      sf.data = coords
    )
    
    # Override with empty polygons for lightweight mode
    slot(obj, 'polygons') <- list()
    slot(obj, 'plotOrder') <- integer(0)
    slot(obj, 'proj4string') <- CRS(as.character(NA))

    # Get bbox from sf data (can be helpful to see coord range)
    slot(obj, 'bbox') <- matrix(sf::st_bbox(coords), 
                                nrow = 2,
                                ncol = 2, 
                                dimnames = list(c("x", "y"), c("min", "max")))

    slot(obj, 'is.lightweight') <- TRUE
    return(obj)
  } else {
    # Convert sf object to SpatialPolygons first
    sf_obj <- sf::st_as_sf(coords)
    sp_obj <- as(object = sf_obj, Class = 'Spatial')
    
    obj <- new(
      Class = 'Segmentation',
      sp_obj,
      sf.data = coords, # Store sf data in its original format
      is.lightweight = FALSE
    )
    
    # Set the cell IDs properly by updating the polygon IDs
    if ("cell" %in% names(coords)) {
      polygons <- slot(object = obj, name = 'polygons')
      for (i in seq_along(polygons)) {
        slot(object = polygons[[i]], name = 'ID') <- coords$cell[i]
      }
      # Update the names of the polygons list
      names(polygons) <- coords$cell
      slot(object = obj, name = 'polygons') <- polygons
    }
  }
  return(obj)
}

#' @method Crop Segmentation
#' @export
#'
Crop.Segmentation <- function(object, ...) {
  if (slot(object = object, name = 'is.lightweight')) {
    warn("Cropping is not yet supported for lightweight Segmentation objects")
    return(object)
  }
  return(.Crop(object, ...))
}

#' @details \code{GetTissueCoordinates}, \code{coordinates}: Get
#' tissue coordinates
#'
#' @inheritParams Centroids-methods
#'
#' @return \code{GetTissueCoordinates}, \code{coordinates}: A data frame with
#' three columns:
#' \itemize{
#'  \item \dQuote{\code{x}}: the x-coordinate
#'  \item \dQuote{\code{y}}: the y-coordinate
#'  \item \dQuote{\code{cell}} or \dQuote{\code{ID}}: the cell name
#' }
#' If \code{full} is \code{TRUE}, then each coordinate will indicate a vertex
#' for the cell polygon; otherwise, each coordinate will indicate a centroid
#' for the cell. Note: \code{GetTissueCoordinates} ....
#'
#' @rdname Segmentation-methods
#' @method GetTissueCoordinates Segmentation
#' @export
#'
GetTissueCoordinates.Segmentation <- function(object, full = TRUE, ...) {
  coords <- coordinates(obj = object, full = full, ...)
  colnames(x = coords) <- c('x', 'y', 'cell')
  rownames(x = coords) <- NULL
  return(coords)
}

#' @details \code{RenameCells}: Update cell names
#'
#' @inheritParams RenameCells
#'
#' @return \code{RenameCells}: \code{object} with the cells renamed to
#' \code{new.names}
#'
#' @importFrom future.apply future_mapply
#'
#' @rdname Segmentation-methods
#' @method RenameCells Segmentation
#' @export
#'
RenameCells.Segmentation <- function(object, new.names = NULL, ...) {
  if (is.null(x = new.names)) {
    return(object)
  }
  new.names <- make.unique(names = new.names)
  if (length(x = new.names) != length(x = Cells(x = object))) {
    stop("Cannot partially rename segmentation cells", call. = FALSE)
  }
  if (slot(object = object, name = 'is.lightweight')) {
    sf_data <- slot(object = object, name = 'sf.data')
    sf_data$cell <- new.names
    slot(object = object, name = 'sf.data') <- sf_data
  } else {
    names(x = slot(object = object, name = 'polygons')) <- new.names
    p <- progressor(along = slot(object = object, name = 'polygons'))
    slot(object = object, name = 'polygons') <- future_mapply(
      FUN = function(polygon, name) {
        slot(object = polygon, name = 'ID') <- name
        p()
        return(polygon)
      },
      polygon = slot(object = object, name = 'polygons'),
      name = new.names,
      SIMPLIFY = FALSE,
      USE.NAMES = TRUE
    )
    sf <- slot(object = object, name = 'sf.data')
    if (!is.null(x = sf)) {
      sf$cell <- new.names
      slot(object = object, name = 'sf.data') <- sf
    }
  }
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @template method-lengths
#'
#' @rdname Segmentation-methods
#' @method lengths Segmentation
#' @export
#'
lengths.Segmentation <- function(x, use.names = TRUE) {
  return(rle(x = GetTissueCoordinates(object = x, full = TRUE)$cell))
}

#' @details \code{subset}, \code{[}: Subset a \code{Segmentation} object to
#' certain cells
#'
#' @return \code{subset}, \code{[}: \code{x} subsetted to the cells specified
#' by \code{cells}/\code{i}
#'
#' @rdname Segmentation-methods
#' @method subset Segmentation
#' @export
#'
subset.Segmentation <- function(x, cells = NULL, ...) {
  if (is.null(x = cells)) {
    return(x)
  }
  if (slot(object = x, name = 'is.lightweight')) {
    sf_data <- slot(object = x, name = 'sf.data')
    sf_data <- sf_data[sf_data$cell %in% cells, ]
    sf_data <- sf_data[order(as.numeric(row.names(sf_data))), ]
    slot(object = x, name = 'sf.data') <- sf_data
    return(x)
  } else {
    if (is.numeric(x = cells)) {
      cells <- Cells(x = x)[cells]
      cells <- MatchCells(new = Cells(x = x), orig = cells, ordered = TRUE)
    } else {
      cells <- intersect(Cells(x), cells)
    }
    if (!length(x = cells)) {
      stop("None of the requested cells found")
    }
    x <- x[cells]
    result <- as(object = x, Class = 'Segmentation')
    # If sf.data is present, subset it as well
    if (!is.null(x = sf_data)) {
      sf_data <- sf_data[sf_data$cell %in% cells, ]
      sf_data <- sf_data[order(as.numeric(row.names(sf_data))), ]
      slot(object = result, name = 'sf.data') <- sf_data
    }
    return(result)
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @details \code{[[<-}: Attach or remove \code{sf}-derived data to/from a \code{Segmentation} object
#'
#' @return \code{[[<-}: 
#' \itemize{
#'  \item If \code{value} is an \code{data.frame} object,
#'  returns \code{x} with \code{value} stored in \code{sf.data};
#'  requires that \code{i} is \dQuote{sf.data}.
#'  \item If \code{value} is \code{NULL}, returns \code{x} with \code{sf.data} removed.
#' }
#' @param value The value to assign to the slot specified by \code{i} in the \code{Segmentation} object.
#' @rdname Segmentation-methods
#'
setMethod(
  f = '[[<-',
  signature = c(
    x = 'Segmentation',
    i = 'character',
    j = 'missing',
    value = 'ANY'
  ),
  definition = function(x, i, ..., value) {
    if (i == "sf.data") {
      if (!is.null(x = value) && !inherits(x = value, what = 'data.frame')) {
        stop("Value assigned to 'sf.data' must inherit from the data.frame class", call. = FALSE)
      }
      # Update sf.data slot
      slot(object = x, name = 'sf.data') <- value
      validObject(x)
      return(x)
    } else {
      stop("Cannot assign value to slot '", i, "' in Segmentation object", call. = FALSE)
    }
  }
)

#' @importFrom methods as
#'
#' @param value The value to assign to the slot specified by \code{i} in the \code{Segmentation} object.
#' @rdname Segmentation-methods
#'
setMethod(
  f = '[[<-',
  signature = c(
    x = 'Segmentation',
    i = 'character',
    j = 'missing',
    value = 'NULL'
  ),
  definition = function(x, i, ..., value) {
    if (i == "sf.data") {
      if (is.null(x = slot(object = x, name = 'sf.data'))) {
        warning("The 'sf.data' slot is already NULL", call. = FALSE)
        return(x)
      } else {
        slot(object = x, name = 'sf.data') <- NULL
      }
    }
    validObject(object = x)
    return(x)
  }
)

#' @rdname Segmentation-methods
#'
setMethod(
  f = '[',
  signature = c(x = 'Segmentation'),
  definition = function(x, i, j, ..., drop = TRUE) {
    if (slot(object = x, name = 'is.lightweight')) {
      sf_data <- slot(object = x, name = 'sf.data')
      sf_data <- sf_data[i, , drop = drop]
      sf_data <- sf_data[order(as.numeric(row.names(sf_data))), ]
      slot(object = x, name = 'sf.data') <- sf_data
      return(x)
    } else {
      # Ensure that subsetting preserves sf.data
      sf_data <- slot(object = x, name = 'sf.data')
      sf_data <- sf_data[i, , drop = drop]
      sf_data <- sf_data[order(as.numeric(row.names(sf_data))), ]
      slot(object = x, name = 'sf.data') <- sf_data
      x <- callNextMethod()
      result <- as(object = x, Class = 'Segmentation')
      # Update the sf.data slot with the subsetted sf data, if it exists
      if (!is.null(x = sf_data)) {
        sf_data <- sf::st_as_sf(sf_data)
        slot(object = result, name = 'sf.data') <- sf_data
        
      }
      return(result)
    }
  }
)

#' @rdname Segmentation-methods
#'
setMethod(
  f = 'coordinates',
  signature = c(obj = 'Segmentation'),
  definition = function(obj, full = TRUE, ...) {
    if (slot(object = obj, name = 'is.lightweight')) {
      coords <- slot(object = obj, name = 'sf.data')
      if (isTRUE(x = full)) {
        return(coords)
      }
      sf_obj <- sf::st_as_sf(coords)
      centroids <- sf::st_centroid(sf_obj)
      centroids_coords <- as.data.frame(sf::st_coordinates(centroids))
      centroids_coords$cell <- coords$cell[centroids_coords[, 'L2']]
      return(centroids_coords)
    } else {
      if (!isTRUE(x = full)) {
        coords <- as.data.frame(x = callNextMethod(obj = obj))
        coords$cell <- Cells(x = obj)
        return(coords)
      }
      coords <- lapply(
        X = slot(object = obj, name = 'polygons'),
        FUN = function(x) {
          polys <- lapply(
            X = slot(object = x, name = 'Polygons'),
            FUN = slot,
            name = 'coords'
          )
          polys <- as.data.frame(x = do.call(what = 'rbind', args = polys))
          colnames(x = polys) <- c('x', 'y')
          polys$ID <- slot(object = x, name = 'ID')
          return(polys)
        }
      )
      coords <- do.call(what = 'rbind', args = coords)
      rownames(x = coords) <- NULL
      return(coords)
    }
  }
)

setMethod(
  f = 'over',
  signature = c(x = 'Segmentation', y = 'SpatialPolygons'),
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
      sparse = FALSE,
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
  signature = c(x = 'Segmentation', y = 'SpatialPolygons'),
  definition = function(x, y, invert = FALSE, ...) {
    check_installed(pkg = 'sf', reason = 'to overlay spatial information')
    is_lightweight <- slot(object = x, name = 'is.lightweight')
    if (is_lightweight) {
      x_sf <- sf::st_as_sf(slot(object = x, name = 'sf.data'))
    } else {
      x_sf <- as(object = x, Class = 'sf')
    }
    idx <- sf::st_intersects(
      x = x_sf,
      y = as(object = y, Class = 'sf'),
      sparse = FALSE
    )
    idx <- which(x = idx)
    if (is_lightweight) {
      names_in_sf_object1 <- x_sf$cell[idx]
    } else {
      names_in_sf_object1 <- if (!is.null(x = row.names(x = x))) {
        row.names(x = x)[idx]
      } else {
        x$id[idx]
      }
    }
    idx <- setNames(
      object = rep.int(x = TRUE, times = length(x = idx)),
      nm = names_in_sf_object1
    )
    if (!length(x = idx)) {
      warn("The selected region does not contain any cell segmentations")
      return(NULL)
    }
    if (is_lightweight) {
      cells <- names(x = idx)
    } else {
      names(x = idx) <- vapply(
        X = strsplit(x = names(x = idx), split = '\\.'),
        FUN = '[[',
        FUN.VALUE = character(length = 1L),
        1L,
        USE.NAMES = FALSE
      )
      cells <- names(x = idx)
    }
    cells <- if (isTRUE(x = invert)) {
      setdiff(x = Cells(x = x), y = cells)
    } else {
      names(x = idx)
    }
    x <- x[cells]
    return(x)
  }
)

#' @template method-show
#'
#' @rdname Segmentation-methods
#'
setMethod(
  f = 'show',
  signature = c(object = 'Segmentation'),
  definition = function(object) {
    cat("A spatial segmentation for", length(x = object), "cells\n")
  }
)

#' Segmentation Validity
#'
#' @templateVar cls Segmentation
#' @template desc-validity
#'
#' @section sf.data Validation:
#' Validates that the sf.data slot contains an object of class \code{sf}.
#'
#' @name Segmentation-validity
#'
#' @family segmentation
#'
#' @seealso \code{\link[methods]{validObject}}
#'
setValidity(
  Class = 'Segmentation',
  method = function(object) {
    if (isFALSE(x = getOption(x = "Seurat.object.validate", default = TRUE))) {
      warn(
        message = paste("Not validating", class(x = object)[1L], "objects"),
        class = 'validationWarning'
      )
      return(TRUE)
    }
    valid <- NULL
    # Check sf.data slot
    sf_data <- slot(object = object, name = 'sf.data')
    if (!is.null(x = sf_data)) {
      # If sf.data is populated, it should inherit from 'sf'
      if (!inherits(x = sf_data, 'data.frame')) {
        valid <- c(
          valid,
          "'sf.data' slot must inherit from 'data.frame' class"
        )
      }
    }
    return(valid %||% TRUE)
  }
)