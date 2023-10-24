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
#' @family segmentation
#' @templateVar cls Segmentation
#' @template seealso-methods
#'
setClass(
  Class = 'Segmentation',
  contains = 'SpatialPolygons'
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

#' @method Crop Segmentation
#' @export
#'
Crop.Segmentation <- .Crop

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
  return(as(object = x, Class = 'Segmentation'))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname Segmentation-methods
#'
setMethod(
  f = '[',
  signature = c(x = 'Segmentation'),
  definition = function(x, i, j, ..., drop = TRUE) {
    x <- callNextMethod()
    return(as(object = x, Class = 'Segmentation'))
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
