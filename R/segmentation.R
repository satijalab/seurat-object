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
#' @slot compact Logical indicating whether or not the object only stores 
#' segmentation information in the \code{sf.data} slot, rather than also in the
#' standard \code{SpatialPolygons} slots, to save memory and processing time.
#' Currently only relevant for Visium data.
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
    compact = 'OptionalLogical'
  ),
  prototype = list(
    compact = FALSE
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
#' \href{https://progressr.futureverse.org/articles/progressr-01-intro.html}{\code{vignette("progressr-intro")}}

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
  compact <- .hasSlot(object = x, name = 'compact') && slot(object = x, name = 'compact')
  if (compact) {
    sf_data <- slot(object = x, name = 'sf.data')
    return(unique(sf_data$cell))
  }
  return(unname(obj = names(x = x)))
}

#' @importFrom sp Polygon Polygons SpatialPolygons CRS
#' @param compact Logical indicating whether or not the object should only store segmentation data
#' in the \code{sf.data} slot; see \link{Segmentation-class} for details.
#'
#' @rdname CreateSegmentation
#' @method CreateSegmentation data.frame
#' @export
#'
CreateSegmentation.data.frame <- function(coords, compact = FALSE) {
  if (compact) {
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

    # Override with empty polygons for compact mode
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

    slot(obj, 'compact') <- TRUE
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
CreateSegmentation.Segmentation <- function(coords, compact = FALSE) {
  return(coords)
}

#' @rdname CreateSegmentation
#' @method CreateSegmentation sf
#'
#' @param compact Logical indicating whether or not the object should only store segmentation data
#' in the \code{sf.data} slot; see \link{Segmentation-class} for details.
#'
#' @export
#'
CreateSegmentation.sf <- function(coords, compact = TRUE) {
  # Method is called when creating Segmentation from an sf object

  # Set the attribute-geometry relationship to constant
  # See https://r-spatial.github.io/sf/reference/sf.html#details
  sf::st_agr(coords) <- "constant"

  # Extract coordinates as a dataframe from sf object
  coords_mat <- sf::st_coordinates(coords)
  l2_indices <- coords_mat[, "L2"] # L2 column corresponds to polygon (cell) index
  
  coords_df <- data.frame(x = coords_mat[, 1],
                          y = coords_mat[, 2],
                          cell = coords$barcodes[l2_indices],
                          stringsAsFactors = FALSE)

  obj <- CreateSegmentation.data.frame(coords = coords_df, compact = compact)

  if (!compact) { # When compact = FALSE, make sure to store the sf dataframe
    slot(object = obj, name = 'sf.data') <- coords_df
  }

  return(obj)
}

#' @method Crop Segmentation
#' @export
#'
Crop.Segmentation <- function(object, ...) {
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
#' for the cell. Note: to compute centroids for segmentations stored in compact mode,
#' call \code{GetTissueCoordinates} on the associated \code{Centroids} object instead.
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
  compact <- .hasSlot(object = object, name = 'compact') && slot(object = object, name = 'compact')
  if (length(x = new.names) != length(x = Cells(x = object))) {
    stop("Cannot partially rename segmentation cells", call. = FALSE)
  }
  sf_data <- if (.hasSlot(object = object, name = 'sf.data')) slot(object = object, name = 'sf.data') else NULL
  id_map <- setNames(new.names, Cells(x = object))
  if (!compact) {
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
  }
  if (!is.null(x = sf_data)) {
    sf_data$cell <- id_map[ sf_data$cell ]
  }
  slot(object = object, name = 'sf.data') <- sf_data
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
  if (is.numeric(x = cells)) { # Account for the case when cells to subset are given as indices and not IDs
    cells <- Cells(x = x)[cells]
    cells <- MatchCells(new = Cells(x = x), orig = cells, ordered = TRUE)
  } else {
    cells <- intersect(Cells(x), cells)
  }
  if (!length(x = cells)) {
    stop("None of the requested cells found")
  }
  compact <- .hasSlot(object = x, name = 'compact') && slot(object = x, name = 'compact')
  sf_data <- if (.hasSlot(object = x, name = 'sf.data')) slot(object = x, name = 'sf.data') else NULL
  
  sf_data_subset <- NULL

  if (!is.null(sf_data)) {
    # Maintain original order of cells in subsetted object
    # Order is important when plotting polygons to avoid the sides (lines between vertices) being drawn incorrectly
    matching_indices <- which(sf_data$cell %in% cells)
    sf_data_subset <- sf_data[matching_indices, ]
  }
  
  if (!compact) {
    x <- x[cells]
    x <- as(object = x, Class = 'Segmentation')
  }
  slot(object = x, name = 'sf.data') <- sf_data_subset
  return(x)
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
    compact <- .hasSlot(object = x, name = 'compact') && slot(object = x, name = 'compact')
    sf_data <- if (.hasSlot(object = x, name = 'sf.data')) slot(object = x, name = 'sf.data') else NULL
    
    if (is.numeric(i)) { # Handle numeric indexing
        cells <- Cells(x = x)[i]
      } else {
        cells <- intersect(Cells(x), i)
      }
    
    sf_data_subset <- NULL
    if (!is.null(sf_data)) {
      matching_indices <- which(sf_data$cell %in% cells)
      sf_data_subset <- sf_data[matching_indices, ]
    }
    
    if (!compact) {
      x <- callNextMethod()
      x <- as(object = x, Class = 'Segmentation')
    }
    
    slot(object = x, name = 'sf.data') <- sf_data_subset
    return(x)
  }
)

#' @rdname Segmentation-methods
#'
setMethod(
  f = 'coordinates',
  signature = c(obj = 'Segmentation'),
  definition = function(obj, full = TRUE, ...) {
    compact <- .hasSlot(object = obj, name = 'compact') && slot(object = obj, name = 'compact')
    if (compact) {
      coords <- slot(object = obj, name = 'sf.data')
      if (isTRUE(x = full)) {
        return(coords)
      }
      message("Centroids cannot be computed for compact Segmentation objects; returning full coordinates.")
      return(coords)
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
    compact <- .hasSlot(object = x, name = 'compact') && slot(object = x, name = 'compact')
    sf_data <- if (.hasSlot(object = x, name = 'sf.data')) slot(object = x, name = 'sf.data') else NULL

    if (compact) {
      message("Overlaying compact Segmentation objects is currently not supported.")
      return(NULL)
    }
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
    sf_data <- if (!is.null(x = sf_data)) {
      matching_indices <- which(sf_data$cell %in% cells)
      sf_data[matching_indices, ]
    } else {
      NULL
    }
    slot(object = x, name = 'sf.data') <- sf_data
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
    compact <- .hasSlot(object = object, name = 'compact') && slot(object = object, name = 'compact')
    sf_data <- if (.hasSlot(object = object, name = 'sf.data')) slot(object = object, name = 'sf.data') else NULL
    if (compact && !is.null(x = sf_data)) {
      cat("A spatial segmentation for", length(unique(sf_data$cell)), "cells\n")
    } else {
      cat("A spatial segmentation for", length(x = object), "cells\n")
    }
    cat("Is compact:", compact, "\n")
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
      # If sf.data is populated, it should inherit from 'data.frame'
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