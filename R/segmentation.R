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
#'
#' @family segmentation
#' @templateVar cls Segmentation
#' @template seealso-methods
#'
setClass(
  Class = 'Segmentation',
  contains = 'SpatialPolygons',
  slots = list(
    sf.data = 'ANY'
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
  return(unname(obj = names(x = x)))
}

#' @importFrom sp Polygon Polygons SpatialPolygons
#'
#' @rdname CreateSegmentation
#' @method CreateSegmentation data.frame
#' @export
#'
CreateSegmentation.data.frame <- function(coords) {
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
CreateSegmentation.sf <- function(coords) {
  # Method is called when creating Segmentation from an sf object
  # Convert sf object to SpatialPolygons first
  sp_obj <- as(object = coords, Class = 'Spatial')
  
  obj <- new(
    Class = 'Segmentation',
    sp_obj,
    sf.data = coords # Store sf data in its original format
  )
  
  # Set the cell IDs properly by updating the polygon IDs
  if ("barcodes" %in% names(coords)) {
    polygons <- slot(object = obj, name = 'polygons')
    for (i in seq_along(polygons)) {
      slot(object = polygons[[i]], name = 'ID') <- coords$barcodes[i]
    }
    # Update the names of the polygons list
    names(polygons) <- coords$barcodes
    slot(object = obj, name = 'polygons') <- polygons
  }
  
  return(obj)
}

#' @method Crop Segmentation
#' @export
#'
Crop.Segmentation <- .Crop

#' Get coordinates for polygon plotting from an \code{sf} object
#'
#' Extract coordinates from an \code{sf} object in a format suitable
#' for use with \code{ggplot2::geom_polygon()}.
#'
#' @param sf_data A \code{sf} object
#' @return A data.frame with columns x, y, and cell (barcode), such that each row
#' represents a single vertex of the segmentation corresponding to that barcode.
#'
#' @export
#'
GetPolygonCoordinates <- function(sf_data) {
  if (!is.null(sf_data) && inherits(sf_data, 'sf')) {
    # Extract coordinates as a dataframe from sf object
    coords_mat <- sf::st_coordinates(sf_data)
    coords_df <- as.data.frame(coords_mat)
    
    n_coords <- nrow(coords_mat)
    coords_df <- data.frame(x = coords_mat[, 1],
                            y = coords_mat[, 2],
                            cell = character(n_coords))

    # L2 column corresponds to polygon (cell) index
    coords_df$cell <- sf_data$barcodes[coords_mat[, "L2"]]

    # Remove geometry, centroids, and cell id (numeric); keep all other columns
    plot_data <- sf_data
    sf::st_geometry(plot_data) <- NULL
    data_columns <- setdiff(names(plot_data), c("x", "y", "cell_id"))
    plot_data <- plot_data[, data_columns, drop = FALSE]
    
    # Check if there are any additional columns to merge
    if (ncol(plot_data) > 0 && "barcodes" %in% names(plot_data)) {
      # Reattach data (ex. metadata, expression data, if any) to coordinates dataframe
      coords_df <- merge(
        coords_df, plot_data, 
        by.x = "cell", by.y = "barcodes", 
        all.x = TRUE, sort = FALSE
      )
    }
    return(coords_df)
  } else {
    stop("Unable to extract polygon coordinates; the 'sf.data' slot does not contain an sf object.")
  }
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
    sf$barcodes <- new.names
    slot(object = object, name = 'sf.data') <- sf
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
  sf_data <- slot(object = x, name = 'sf.data')
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
    sf_data <- sf_data[sf_data$barcodes %in% cells, ]
    sf_data <- sf::st_as_sf(sf_data)
    slot(object = result, name = 'sf.data') <- sf_data
  }
  return(result)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @details \code{[[<-}: Attach or remove \code{sf} object to/from a \code{Segmentation} object
#'
#' @return \code{[[<-}: 
#' \itemize{
#'  \item If \code{value} is an \code{sf} object,
#'  returns \code{x} with \code{value} stored in \code{sf};
#'  requires that \code{i} is \dQuote{sf.data}.
#'  \item If \code{value} is \code{NULL}, returns \code{x} with \code{sf} removed.
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
      if (!is.null(x = value) && !inherits(x = value, what = 'sf')) {
        stop("Value assigned to 'sf.data' must inherit from the sf class", call. = FALSE)
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
    # Ensure that subsetting preserves sf.data
    sf_data <- slot(object = x, name = 'sf.data')
    if (!is.null(x = sf_data)) {
      sf_data <- sf_data[i, , drop = drop]
    }
    x <- callNextMethod()
    result <- as(object = x, Class = 'Segmentation')
    # Update the sf.data slot with the subsetted sf data, if it exists
    if (!is.null(x = sf_data)) {
      sf_data <- sf::st_as_sf(sf_data)
      slot(object = result, name = 'sf.data') <- sf_data
    }
    return(result)
  }
)

#' @rdname Segmentation-methods
#'
setMethod(
  f = 'coordinates',
  signature = c(obj = 'Segmentation'),
  definition = function(obj, full = TRUE, ...) {
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
    idx <- sf::st_intersects(
      x = as(object = x, Class = 'sf'),
      y = as(object = y, Class = 'sf'),
      sparse = FALSE
    )
    idx <- which(x = idx)
    names_in_sf_object1 <- if (!is.null(x = row.names(x = x))) {
      row.names(x = x)[idx]
    } else {
      x$id[idx]
    }
    idx <- setNames(
      object = rep.int(x = TRUE, times = length(x = idx)),
      nm = names_in_sf_object1
    )
    if (!length(x = idx)) {
      warn("The selected region does not contain any cell segmentations")
      return(NULL)
    }
    names(x = idx) <- vapply(
      X = strsplit(x = names(x = idx), split = '\\.'),
      FUN = '[[',
      FUN.VALUE = character(length = 1L),
      1L,
      USE.NAMES = FALSE
    )
    cells <- if (isTRUE(x = invert)) {
      setdiff(x = Cells(x = x), y = names(x = idx))
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
      if (!inherits(x = sf_data, 'sf')) {
        valid <- c(
          valid,
          "'sf.data' slot must inherit from 'sf' class"
        )
      }
    }
    return(valid %||% TRUE)
  }
)