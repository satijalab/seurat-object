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
