#' @include zzz.R
#' @include generics.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The Graph Class
#'
#' The Graph class inherits from \code{\link[Matrix:sparseMatrix]{dgCMatrix}}.
#' We do this to enable future expandability of graphs.
#'
#' @slot assay.used Optional name of assay used to generate \code{Graph} object
#'
#' @name Graph-class
#' @rdname Graph-class
#' @exportClass Graph
#'
#' @seealso \code{\link[Matrix]{dgCMatrix-class}}
#'
Graph <- setClass(
  Class = 'Graph',
  contains = "dgCMatrix",
  slots = list(
    assay.used = 'OptionalCharacter'
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom methods as
#'
#' @rdname as.Graph
#' @export
#' @method as.Graph Matrix
#'
#' @examples
#' # converting sparse matrix
#' mat <- Matrix::rsparsematrix(nrow = 10, ncol = 10, density = 0.1)
#' rownames(x = mat) <- paste0("feature_", 1:10)
#' colnames(x = mat) <- paste0("cell_", 1:10)
#' g <- as.Graph(x = mat)
#'
as.Graph.Matrix <- function(x, ...) {
  CheckDots(...)
  x <- as.sparse(x = x)
  if (is.null(x = rownames(x = x))) {
    stop("Please provide rownames to the matrix before converting to a Graph.")
  }
  if (is.null(x = colnames(x = x))) {
    stop("Please provide colnames to the matrix before converting to a Graph.")
  }
  return(as(object = x, Class = "Graph"))
}

#' @rdname as.Graph
#' @export
#' @method as.Graph matrix
#'
#' @examples
#' # converting dense matrix
#' mat <- matrix(data = 1:16, nrow = 4)
#' rownames(x = mat) <- paste0("feature_", 1:4)
#' colnames(x = mat) <- paste0("cell_", 1:4)
#' g <- as.Graph(x = mat)
#'
as.Graph.matrix <- as.Graph.Matrix

#' @param weighted If TRUE, fill entries in Graph matrix with value from the
#' nn.dist slot of the Neighbor object
#'
#' @rdname as.Graph
#' @export
#' @method as.Graph Neighbor
#'
as.Graph.Neighbor <- function(x, weighted = TRUE, ...) {
  CheckDots(...)
  j <- as.integer(x = Indices(object = x) - 1)
  i <- as.integer(x = rep(x = (1:nrow(x = x)) - 1, times = ncol(x = x)))
  vals <- if (weighted) {
    as.vector(x = Distances(object = x))
  } else {
    1
  }
  graph <- new(
    Class = "dgTMatrix",
    i = i,
    j = j,
    x = vals,
    Dim = as.integer(x = c(nrow(x = x), nrow(x = x)))
  )
  colnames(x = graph) <- rownames(x = graph) <- Cells(x = x)
  graph <- as.Graph.Matrix(x = graph)
  return(graph)
}

#' @rdname DefaultAssay
#' @export
#' @method DefaultAssay Graph
#'
DefaultAssay.Graph <- function(object, ...) {
  object <- UpdateSlots(object = object)
  return(slot(object = object, name = 'assay.used'))
}

#' @rdname DefaultAssay
#' @export
#' @method DefaultAssay<- Graph
#'
"DefaultAssay<-.Graph" <- function(object, ..., value) {
  object <- UpdateSlots(object = object)
  slot(object = object, name = 'assay.used') <- value
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
