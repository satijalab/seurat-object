#' @include zzz.R
#' @include generics.R
#' @importFrom methods new slot slot<-
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The Neighbor class
#'
#' The Neighbor class is used to store the results of neighbor finding
#' algorithms
#'
#' @slot nn.idx Matrix containing the nearest neighbor indices
#' @slot nn.dist Matrix containing the nearest neighbor distances
#' @slot alg.idx The neighbor finding index (if applicable). E.g. the annoy
#' index
#' @slot alg.info Any information associated with the algorithm that may be
#' needed downstream (e.g. distance metric used with annoy is needed when
#' reading in from stored file).
#' @slot cell.names Names of the cells for which the neighbors have been
#' computed.
#'
#' @name Neighbor-class
#' @rdname Neighbor-class
#' @exportClass Neighbor
#'
Neighbor <- setClass(
  Class = 'Neighbor',
  slots = c(
    nn.idx = 'matrix',
    nn.dist = 'matrix',
    alg.idx = 'ANY',
    alg.info = 'list',
    cell.names = 'character'
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname as.Neighbor
#' @export
#' @method as.Neighbor Graph
#'
as.Neighbor.Graph <- function(x, ...) {
  nn.mats <- GraphToNeighborHelper(mat = x)
  return(Neighbor(
    nn.idx = nn.mats[[1]],
    nn.dist = nn.mats[[2]],
    cell.names = rownames(x = x)
  ))
}

#' @rdname Cells
#' @method Cells Neighbor
#' @export
#'
Cells.Neighbor <- function(x) {
  return(slot(object = x, name = "cell.names"))
}

#' @rdname Distances
#' @export
#' @method Distances Neighbor
#'
Distances.Neighbor <- function(object, ...) {
  object <- UpdateSlots(object = object)
  distances <- slot(object = object, name = "nn.dist")
  rownames(x = distances) <- slot(object = object, name = "cell.names")
  return(distances)
}

#' @rdname Index
#' @export
#' @method Index Neighbor
#'
Index.Neighbor <- function(object, ...) {
  object <- UpdateSlots(object = object)
  index <- slot(object = object, name = "alg.idx")
  if (is.null(x = index)) {
    return(NULL)
  } else if (IsNullPtr(x = index$.pointer)) {
    return(NULL)
  }
  return(index)
}

#' @rdname Index
#' @export
#' @method Index<- Neighbor
#'
"Index<-.Neighbor" <- function(object, ..., value) {
  CheckDots(...)
  slot(object = object, name = "alg.idx") <- value
  return(object)
}

#' @rdname Indices
#' @export
#' @method Indices Neighbor
#'
Indices.Neighbor <- function(object, ...) {
  object <- UpdateSlots(object = object)
  indices <- slot(object = object, name = "nn.idx")
  rownames(x = indices) <- slot(object = object, name = "cell.names")
  return(indices)
}

#' @param old.names vector of old cell names
#' @rdname RenameCells
#' @export
#' @method RenameCells Neighbor
#'
RenameCells.Neighbor <- function(
  object,
  old.names = NULL,
  new.names = NULL,
  ...
) {
  CheckDots(...)
  neighbor.names <- Cells(x = object)
  names(x = new.names) <- old.names
  slot(object = object, name = "cell.names") <- unname(obj = new.names[neighbor.names])
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' \code{Neighbor} Methods
#'
#' Methods for \code{\link{Neighbor}} objects for generics defined in
#' other packages
#'
#' @param x,object A \code{\link{Neighbor}} object
#'
#' @name Neighbor-methods
#' @rdname Neighbor-methods
#'
#' @concept neighbor
#'
NULL

#' @describeIn Neighbor-methods Dimensions of the neighbor indices
#'
#' @return \code{dim} Dimensions of the indices matrix
#'
#' @export
#' @method dim Neighbor
#'
dim.Neighbor <- function(x) {
  return(dim(x = Indices(object = x)))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @describeIn Neighbor-methods Overview of a \code{Neighbor} object
#'
#' @return \code{show}: Prints summary to \code{\link[base]{stdout}} and
#' invisibly returns \code{NULL}
#'
#' @importFrom methods show
#'
#' @export
#'
setMethod(
  f = 'show',
  signature = 'Neighbor',
  definition = function(object) {
    cat(
      "A Neighbor object containing the",
      ncol(x = object),
      "nearest neighbors for",
      nrow(x = object),
      "cells"
    )
  }
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
