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
#'
#' @exportClass Graph
#'
#' @seealso \code{\link[Matrix]{dgCMatrix-class}}
#'
#' @family graph
#'
Graph <- setClass(
  Class = 'Graph',
  contains = "dgCMatrix",
  slots = list(
    assay.used = 'character'
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
    abort(message = "Please provide rownames to the matrix before converting to a Graph")
  }
  if (is.null(x = colnames(x = x))) {
    abort(message = "Please provide colnames to the matrix before converting to a Graph")
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

#' @method Cells Graph
#' @export
#'
Cells.Graph <- function(x, margin = 1L, ...) {
  margin <- as.integer(x = margin)
  if (is_na(x = margin)) {
    return(Reduce(f = union, x = dimnames(x = x)))
  }
  if (!isTRUE(margin %in% c(1L, 2L))) {
    stop("'margin' must be either 1 or 2", call. = FALSE)
  }
  return(dimnames(x = x)[[margin]])
}

#' @rdname DefaultAssay
#' @export
#' @method DefaultAssay Graph
#'
DefaultAssay.Graph <- function(object, ...) {
  # object <- UpdateSlots(object = object)
  assay <- slot(object = object, name = 'assay.used')
  if (!length(x = assay)) {
    assay <- NULL
  }
  return(assay)
}

#' @rdname DefaultAssay
#' @export
#' @method DefaultAssay<- Graph
#'
"DefaultAssay<-.Graph" <- function(object, ..., value) {
  op <- options(Seurat.object.validate = FALSE)
  on.exit(expr = options(op))
  object <- suppressWarnings(expr = UpdateSlots(object = object))
  if (!length(x = value) || !isTRUE(x = nzchar(x = value))) {
    value <- character(length = 0L)
  }
  slot(object = object, name = 'assay.used') <- value
  options(op)
  validObject(object = object)
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setValidity(
  Class = 'Graph',
  method = function(object) {
    if (.GetSeuratCompat() < '5.0.0') {
      return(TRUE)
    }
    if (isFALSE(x = getOption(x = "Seurat.object.validate", default = TRUE))) {
      warn(
        message = paste("Not validating", class(x = object)[1L], "objects"),
        class = 'validationWarning'
      )
      return(TRUE)
    }
    valid <- NULL
    # Check dimnames
    dnames <- dimnames(x = object)
    for (i in seq_along(along.with = dnames)) {
      type <- c('row', 'column')[i]
      if (is.null(x = dnames[[i]])) {
        valid <- c(valid, paste(type, "names must be provided"))
      } else if (any(!nzchar(x = dnames[[i]]))) {
        valid <- c(valid, paste(type, "names must not be empty strings"))
      } else if (anyDuplicated(x = dnames[[i]])) {
        valid <- c(valid, paste(type, "names may not contain duplicates"))
      }
    }
    # Check default assay
    assay <- DefaultAssay(object = object)
    if (length(x = assay) && !nzchar(x = assay)) {
      valid <- c(valid, "'assay.used' may not be an empty character ('')")
    }
    return(valid %||% TRUE)
  }
)

#' Graph Object Overview
#'
#' Overview of a \code{\link{Graph}} Object
#'
#' @template return-show
#'
#' @keywords internal
#'
#' @concept graph
#'
#' @examples
#' pbmc_small[["RNA_snn"]]
#'
setMethod(
  f = 'show',
  signature = 'Graph',
  definition = function(object) {
    cat(
      "A Graph object containing",
      nrow(x = object),
      "cells"
    )
  }
)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
