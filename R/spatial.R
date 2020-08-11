#' @include zzz.R
#' @include generics.R
#' @importFrom methods setOldClass setClass slot slot<- new
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @note \code{scalefactors} objects can be created with \code{scalefactors()}
#'
#' @param spot Spot full resolution scale factor
#' @param fiducial Fiducial full resolution scale factor
#' @param hires High resolutoin scale factor
#' @param lowres Low resolution scale factor
#'
#' @rdname ScaleFactors
#' @export
#'
scalefactors <- function(spot, fiducial, hires, lowres) {
  object <- list(
    spot = spot,
    fiducial = fiducial,
    hires = hires,
    lowres = lowres
  )
  object <- sapply(X = object, FUN = as.numeric, simplify = FALSE, USE.NAMES = TRUE)
  return(structure(.Data = object, class = 'scalefactors'))
}

setOldClass(Classes = c('scalefactors'))

#' The SpatialImage class
#'
#' The \code{SpatialImage} class is a virtual class representing spatial
#' information for Seurat. All spatial image information must inherit from this
#' class for use with \code{Seurat} objects
#'
#' @slot assay Name of assay to associate image data with; will give this image
#' priority for visualization when the assay is set as the active/default assay
#' in a \code{Seurat} object
#' @slot key Key for the image
#'
#' @section Provided methods:
#' These methods are defined on the \code{SpatialImage} object and should not be
#' overwritten without careful thought
#' \itemize{
#'   \item \code{\link{DefaultAssay}} and \code{\link{DefaultAssay<-}}
#'   \item \code{\link{Key}} and \code{\link{Key<-}}
#'   \item \code{\link{IsGlobal}}
#'   \item \code{\link{Radius}}; this method \emph{can} be overridden to provide
#'   a spot radius for image objects
#' }
#'
#' @section Required methods:
#' All subclasses of the \code{SpatialImage} class must define the following
#' methods; simply relying on the \code{SpatialImage} method will result in
#' errors. For required parameters and their values, see the \code{Usage} and
#' \code{Arguments} sections
#' \describe{
#'   \item{\code{\link{Cells}}}{
#'    Return the cell/spot barcodes associated with each position
#'   }
#'   \item{\code{\link{dim}}}{
#'    Return the dimensions of the image for plotting in \code{(Y, X)} format
#'   }
#'   \item{\code{\link{GetImage}}}{
#'    Return image data; by default, must return a grob object
#'   }
#'   \item{\code{\link{GetTissueCoordinates}}}{
#'    Return tissue coordinates; by default, must return a two-column data.frame
#'    with x-coordinates in the first column and y-coordiantes in the second
#'   }
#'   \item{\code{\link{Radius}}}{
#'    Return the spot radius; returns \code{NULL} by default for use with
#'    non-spot image technologies
#'   }
#'   \item{\code{\link{RenameCells}}}{
#'    Rename the cell/spot barcodes for this image
#'   }
#'   \item{\code{\link{subset}} and \code{[}}{
#'    Subset the image data by cells/spots; \code{[} should only take \code{i}
#'    for subsetting by cells/spots
#'   }
#' }
#' These methods are used throughout Seurat, so defining them and setting the
#' proper defaults will allow subclasses of \code{SpatialImage} to work
#' seamlessly
#'
#' @name SpatialImage-class
#' @rdname SpatialImage-class
#'
#' @exportClass SpatialImage
#'
SpatialImage <- setClass(
  Class = 'SpatialImage',
  contains = 'VIRTUAL',
  slots = list(
    'assay' = 'character',
    'key' = 'character'
  )
)

#' The SlideSeq class
#'
#' The SlideSeq class represents spatial information from the Slide-seq platform
#'
#' @slot coordinates ...
#' @slot ...
#'
#' @name SlideSeq-class
#' @rdname SlideSeq-class
#'
#' @seealso For slots defined by the \code{SpatialImage} class, see
#' \code{\link{SpatialImage-class}}
#'
SlideSeq <- setClass(
  Class = 'SlideSeq',
  contains = 'SpatialImage',
  slots = list(
    'coordinates' = 'data.frame'
  )
)

#' The STARmap class
#'
#' The STARmap class represents spatial information from the STARmap platform
#'
#' @slot ... ...
#'
#' @seealso For slots defined by the \code{SpatialImage} class, see
#' \code{\link{SpatialImage-class}}
#'
STARmap <- setClass(
  Class = 'STARmap',
  contains = 'SpatialImage',
  slots = list(
    'coordinates' = 'data.frame',
    'qhulls' = 'data.frame'
  )
)

#' The VisiumV1 class
#'
#' The VisiumV1 class represents spatial information from the 10X Genomics
#' Visium platform
#'
#' @slot image A three-dimensional array with PNG image data, see
#' \code{\link[png]{readPNG}} for more details
#' @slot scale.factors An object of class \code{\link{scalefactors}}; see
#' \code{\link{scalefactors}} for more information
#' @slot coordinates A data frame with tissue coordinate information
#' @slot spot.radius Single numeric value giving the radius of the spots
#'
#' @name VisiumV1-class
#' @rdname VisiumV1-class
#' @exportClass VisiumV1
#'
#' @seealso For slots defined by the \code{SpatialImage} class, see
#' \code{\link{SpatialImage-class}}
#'
VisiumV1 <- setClass(
  Class = 'VisiumV1',
  contains = 'SpatialImage',
  slots = list(
    'image' = 'array',
    'scale.factors' = 'scalefactors',
    'coordinates' = 'data.frame',
    'spot.radius' = 'numeric'
  )
)

setClass(Class = 'SliceImage', contains = 'VisiumV1')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname Cells
#' @method Cells SlideSeq
#' @export
#'
Cells.SlideSeq <- function(x) {
  return(rownames(x = GetTissueCoordinates(object = x)))
}

#' @param x,object An object inheriting from \code{SpatialImage}
#'
#' @name SpatialImage-class
#' @rdname SpatialImage-class
#'
#' @method Cells SpatialImage
#' @export
#'
Cells.SpatialImage <- function(x) {
  stop(
    "'Cells' must be implemented for all subclasses of 'SpatialImage'",
    call. = FALSE
  )
}

#' @rdname Cells
#' @method Cells STARmap
#' @export
#'
Cells.STARmap <- function(x) {
  return(rownames(x = GetTissueCoordinates(object = x)))
}

#' @rdname Cells
#' @method Cells VisiumV1
#' @export
#'
Cells.VisiumV1 <- function(x) {
  return(rownames(x = GetTissueCoordinates(object = x, scale = NULL)))
}

#' @rdname DefaultAssay
#' @method DefaultAssay SpatialImage
#' @export
#'
DefaultAssay.SpatialImage <- function(object, ...) {
  CheckDots(...)
  return(slot(object = object, name = 'assay'))
}

#' @rdname DefaultAssay
#' @method DefaultAssay<- SpatialImage
#' @export
#'
"DefaultAssay<-.SpatialImage" <- function(object, ..., value) {
  CheckDots(...)
  slot(object = object, name = 'assay') <- value
  return(object)
}

#' @method GetImage SlideSeq
#' @export
#'
GetImage.SlideSeq <- function(
  object,
  mode = c('grob', 'raster', 'plotly', 'raw'),
  ...
) {
  mode <- match.arg(arg = mode)
  return(NullImage(mode = mode))
}

#' @inheritParams GetImage
#'
#' @rdname SpatialImage-class
#' @name SpatialImage-class
#'
#' @method GetImage SpatialImage
#' @export
#'
GetImage.SpatialImage <- function(
  object,
  mode = c('grob', 'raster', 'plotly', 'raw'),
  ...
) {
  mode <- match.arg(arg = mode)
  stop(
    "'GetImage' must be implemented for all subclasses of 'SpatialImage'",
    call. = FALSE
  )
}

#' @method GetImage STARmap
#' @export
#'
GetImage.STARmap <- function(
  object,
  mode = c('grob', 'raster', 'plotly', 'raw'),
  ...
) {
  mode <- match.arg(arg = mode)
  return(NullImage(mode = mode))
}

#' @importFrom plotly raster2uri
#' @importFrom grDevices as.raster
#' @importFrom grid rasterGrob unit
#'
#' @rdname GetImage
#' @method GetImage VisiumV1
#' @export
#'
GetImage.VisiumV1 <- function(
  object,
  mode = c('grob', 'raster', 'plotly', 'raw'),
  ...
) {
  mode <- match.arg(arg = mode)
  image <- slot(object = object, name = 'image')
  image <- switch(
    EXPR = mode,
    'grob' = rasterGrob(
      image = image,
      width = unit(x = 1, units = 'npc'),
      height = unit(x = 1, units = 'npc')
    ),
    'raster' = as.raster(x = image),
    'plotly' = list(
      source = raster2uri(r = GetImage(object = object, mode = 'raster')),
      xref = 'x',
      yref = 'y',
      # x = -7,
      # y = -7,
      sizex = ncol(x = object),
      sizey = nrow(x = object),
      sizing = 'stretch',
      opacity = 1,
      layer = 'below'
    ),
    'raw' = image,
    stop("Unknown image mode: ", mode, call. = FALSE)
  )
  return(image)
}

#' @rdname GetTissueCoordinates
#' @method GetTissueCoordinates SlideSeq
#' @export
#'
GetTissueCoordinates.SlideSeq <- function(object, ...) {
  coords <- slot(object = object, name = 'coordinates')
  colnames(x = coords) <- c('x', 'y')
  # coords$y <- -rev(x = coords$y) + 1
  # coords$y <- FlipCoords(x = coords$y)
  coords$cells <- rownames(x = coords)
  return(coords)
}

#' @inheritParams GetTissueCoordinates
#'
#' @rdname SpatialImage-class
#' @name SpatialImage-class
#'
#' @method GetTissueCoordinates SpatialImage
#' @export
#'
GetTissueCoordinates.SpatialImage <- function(object, ...) {
  stop(
    "'GetTissueCoordinates' must be implemented for all subclasses of 'SpatialImage'",
    call. = FALSE
  )
}

#' @param qhulls return qhulls instead of centroids
#'
#' @rdname GetTissueCoordinates
#' @method GetTissueCoordinates STARmap
#' @export
#'
GetTissueCoordinates.STARmap <- function(object, qhulls = FALSE, ...) {
  if (qhulls) {
    return(slot(object = object, name = 'qhulls'))
  }
  return(slot(object = object, name = 'coordinates'))
}

#' @param scale A factor to scale the coordinates by; choose from: 'tissue',
#' 'fiducial', 'hires', 'lowres', or \code{NULL} for no scaling
#' @param cols Columns of tissue coordinates data.frame to pull
#'
#' @rdname GetTissueCoordinates
#' @method GetTissueCoordinates VisiumV1
#' @export
#'
GetTissueCoordinates.VisiumV1 <- function(
  object,
  scale = 'lowres',
  cols = c('imagerow', 'imagecol'),
  ...
) {
  cols <- cols %||% colnames(x = slot(object = object, name = 'coordinates'))
  if (!is.null(x = scale)) {
    coordinates <- slot(object = object, name = 'coordinates')[, c('imagerow', 'imagecol')]
    scale <- match.arg(arg = scale, choices = c('spot', 'fiducial', 'hires', 'lowres'))
    scale.use <- ScaleFactors(object = object)[[scale]]
    coordinates <- coordinates * scale.use
  } else {
    coordinates <- slot(object = object, name = 'coordinates')[, cols]
  }
  return(coordinates)
}

#' @rdname IsGlobal
#' @method IsGlobal SpatialImage
#' @export
#'
IsGlobal.SpatialImage <- function(object) {
  return(TRUE)
}

#' @rdname Key
#' @method Key SpatialImage
#' @export
#'
Key.SpatialImage <- function(object, ...) {
  CheckDots(...)
  object <- UpdateSlots(object = object)
  return(slot(object = object, name = 'key'))
}

#' @rdname Key
#' @method Key<- SpatialImage
#' @export
#'
"Key<-.SpatialImage" <- function(object, ..., value) {
  CheckDots(...)
  object <- UpdateSlots(object = object)
  value <- UpdateKey(key = value)
  slot(object = object, name = 'key') <- value
  return(object)
}

#' @rdname Radius
#' @method Radius SlideSeq
#' @export
#'
Radius.SlideSeq <- function(object) {
  return(0.005)
}

#' @rdname SpatialImage-class
#' @name SpatialImage-class
#'
#' @method Radius SpatialImage
#' @export
#'
Radius.SpatialImage <- function(object) {
  return(NULL)
}

#' @rdname Radius
#' @method Radius VisiumV1
#' @export
#'
Radius.VisiumV1 <- function(object) {
  return(slot(object = object, name = 'spot.radius'))
}

#' @method RenameCells SlideSeq
#' @export
#'
RenameCells.SlideSeq <- function(object, new.names = NULL, ...) {
  return(RenameCells.VisiumV1(object = object, new.names = new.names))
}

#' @inheritParams  RenameCells
#'
#' @rdname SpatialImage-class
#' @name SpatialImage-class
#'
#' @method RenameCells SpatialImage
#' @export
#'
RenameCells.SpatialImage <- function(object, new.names = NULL, ...) {
  stop(
    "'RenameCells' must be implemented for all subclasses of 'SpatialImage'",
    call. = FALSE
  )
}

#' @method RenameCells STARmap
#' @export
#'
RenameCells.STARmap <- function(object, new.names = NULL, ...) {
  names(x = new.names) <- Cells(x = object)
  object <- RenameCells.VisiumV1(object = object, new.names = new.names)
  qhulls <- GetTissueCoordinates(object = object, qhull = TRUE)
  qhulls$cell <- new.names[qhulls$cell]
  slot(object = object, name = "qhulls") <- qhulls
  return(object)
}

#' @rdname RenameCells
#' @method RenameCells VisiumV1
#' @export
#'
RenameCells.VisiumV1 <- function(object, new.names = NULL, ...) {
  if (is.null(x = new.names)) {
    return(object)
  } else if (length(x = new.names) != length(x = Cells(x = object))) {
    stop("Wrong number of cell/spot names", call. = FALSE)
  }
  names(x = new.names) <- Cells(x = object)
  coordinates <- GetTissueCoordinates(object = object, scale = NULL, cols = NULL)
  rownames(x = coordinates) <- new.names[rownames(x = coordinates)]
  slot(object = object, name = 'coordinates') <- coordinates
  return(object)
}

#' @rdname ScaleFactors
#' @method ScaleFactors VisiumV1
#' @export
#'
ScaleFactors.VisiumV1 <- function(object, ...) {
  return(slot(object = object, name = 'scale.factors'))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @method [ SlideSeq
#' @export
#'
"[.SlideSeq" <- function(x, i, ...) {
  return(subset(x = x, cells = i, ...))
}

#' @param i,cells A vector of cells to keep
#'
#' @rdname SpatialImage-class
#' @name SpatialImage-class
#'
#' @method [ SpatialImage
#' @export
#'
"[.SpatialImage" <- function(x, i, ...) {
  stop(
    "'[' must be implemented for all subclasses of 'SpatialImage'",
    call. = FALSE
  )
}

#' @method [ VisiumV1
#' @export
#'
"[.VisiumV1" <- function(x, i, ...) {
  return(subset(x = x, cells = i))
}

#' @method dim SlideSeq
#' @export
#'
dim.SlideSeq <- function(x) {
  # return(dim(x = GetImage(object = x, mode = 'raw')))
  return(c(599, 600))
}

#' @rdname SpatialImage-class
#' @name SpatialImage-class
#'
#' @method dim SpatialImage
#' @export
#'
dim.SpatialImage <- function(x) {
  stop(
    "'dim' must be implemented for all subclasses of 'SpatialImage'",
    call. = FALSE
  )
}

#' @method dim STARmap
#' @export
#'
dim.STARmap <- function(x) {
  coords <- GetTissueCoordinates(object = x)
  return(c(
    max(coords[, 1]) - min(coords[, 1]),
    max(coords[, 2]) - min(coords[, 2])
  ))
}

#' @method dim VisiumV1
#' @export
#'
dim.VisiumV1 <- function(x) {
  return(dim(x = GetImage(object = x)$raster))
}

#' @method subset SlideSeq
#' @export
#'
subset.SlideSeq <- function(x, cells, ...) {
  x <- subset.VisiumV1(x = x, cells = cells, ...)
  return(x)
}

#' @rdname SpatialImage-class
#' @name SpatialImage-class
#'
#' @method subset SpatialImage
#' @export
#'
subset.SpatialImage <- function(x, cells, ...) {
  stop("'subset' must be implemented for all subclasses of 'SpatialImage'")
}

#' @method subset STARmap
#' @export
#'
subset.STARmap <- function(x, cells, ...) {
  x <- subset.VisiumV1(x = x, cells = cells, ...)
  qhulls <- GetTissueCoordinates(object = x, qhulls = TRUE)
  qhulls <- qhulls[qhulls$cell %in% cells, ]
  slot(object = x, name = 'qhulls') <- qhulls
  return(x)
}

#' @method subset VisiumV1
#' @export
#'
subset.VisiumV1 <- function(x, cells, ...) {
  coordinates <- GetTissueCoordinates(object = x, scale = NULL, cols = NULL)
  cells <- cells[cells %in% rownames(x = coordinates)]
  coordinates <- coordinates[cells, ]
  slot(object = x, name = 'coordinates') <- coordinates
  return(x)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setMethod(
  f = 'show',
  signature = 'SpatialImage',
  definition = function(object) {
    object <- UpdateSlots(object = object)
    cat(
      "Spatial data from the",
      class(x = object),
      "technology for",
      length(x = Cells(x = object)),
      "samples\n"
    )
    cat("Associated assay:", DefaultAssay(object = object), "\n")
    cat("Image key:", Key(object = object), "\n")
    return(invisible(x = NULL))
  }
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Return a null image
#'
#' @param mode Image representation to return
#' see \code{\link{GetImage}} for more details
#'
#' @importFrom grid nullGrob
#' @importFrom grDevices as.raster
#'
#' @keywords internal
#'
NullImage <- function(mode) {
  image <- switch(
    EXPR = mode,
    'grob' = nullGrob(),
    'raster' = as.raster(x = new(Class = 'matrix')),
    'plotly' = list('visible' = FALSE),
    'raw' = NULL,
    stop("Unknown image mode: ", mode, call. = FALSE)
  )
  return(image)
}
