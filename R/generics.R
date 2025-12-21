#' @include zzz.R
#'
NULL

#' Assay Class Label
#'
#' @param object A \code{\link{StdAssay}} object
#'
#' @return The assay class label for \code{object}
#'
#' @keywords internal
#'
#' @export .AssayClass
#'
.AssayClass <- function(object) {
  UseMethod(generic = '.AssayClass', object = object)
}

#' Calculate nCount and nFeature
#'
#' @template param-dots-method
#' @param object An assay-like object
#'
#' @return A named list with ...
#'
#' @keywords internal
#'
#' @export .CalcN
#'
#' @examples
#' calcn <- .CalcN(pbmc_small[["RNA"]])
#' head(as.data.frame(calcn))
#'
.CalcN <- function(object, ...) {
  UseMethod(generic = '.CalcN', object = object)
}

#' Get the Package that Defines a Class
#'
#' @param object An object
#'
#' @return The package that defines the class of \code{object}
#'
#' @keywords internal
#'
#' @export .ClassPkg
#'
#' @examples
#' .ClassPkg(pbmc_small)
#'
.ClassPkg <- function(object) {
  UseMethod(generic = '.ClassPkg', object = object)
}

#' Generic Assay Creation
#'
#' Create an assay object; runs a standardized filtering scheme that
#' works regardless of the direction of the data (eg. cells as columns
#' and features as rows or vice versa) and creates an assay object based
#' on the  initialization scheme defined for \code{\link{StdAssay}}-derived
#' class \code{type}
#'
#' @param counts A two-dimensional expression matrix
#' @param min.cells Include features detected in at least this many cells;
#' will subset the counts matrix as well. To reintroduce excluded features,
#' create a new object with a lower cutoff
#' @param min.features Include cells where at least this many features
#' are detected
#' @param cells Vector of cell names
#' @param features Vector of feature names
#' @param type Type of assay object to create; must be the name of a class
#' that's derived from \code{\link{StdAssay}}
#' @param ... Extra parameters passed to \code{\link[methods]{new}} for
#' assay creation; used to set slots not defined by \code{\link{StdAssay}}
#'
#' @return An object of class \code{type} with a layer named \code{layer}
#' containing the data found in \code{counts}
#'
#' @keywords internal
#'
#' @export .CreateStdAssay
#'
#' @concept assay
#'
.CreateStdAssay <- function(
  counts,
  min.cells = 0,
  min.features = 0,
  cells = NULL,
  features = NULL,
  transpose = FALSE,
  type = 'Assay5',
  ...
) {
  UseMethod(generic = '.CreateStdAssay', object = counts)
}

#' Disk Loading Function
#'
#' Generate a function to load a matrix from an on-disk file
#'
#' @inheritParams .FilePath
#'
#' @return A one-length character that defines a function to load a matrix from
#' a file
#'
#' @keywords internal
#'
#' @export .DiskLoad
#'
.DiskLoad <- function(x) {
  UseMethod(generic = '.DiskLoad', object = x)
}

#' Find a File Path
#'
#' @param x A file-backed object
#'
#' @return The path to the file that backs \code{x}; if \code{x} is not a
#' file-backed object, returns \code{NULL}
#'
#' @keywords internal
#'
#' @export .FilePath
#'
.FilePath <- function(x) {
  UseMethod(generic = '.FilePath', object = x)
}

#' Get the Margin of an Object
#'
#' @param x An object
#'
#' @return The margin, eg. \code{1} for rows or \code{2} for columns
#'
#' @keywords internal
#'
#' @export .MARGIN
#'
.MARGIN <- function(x, ...) {
  UseMethod(generic = '.MARGIN', object = x)
}

#' Add in metadata associated with either cells or features.
#'
#' Adds additional data to the object. Can be any piece of information
#' associated with a cell (examples include read depth, alignment rate,
#' experimental batch, or subpopulation identity) or feature (ENSG name,
#' variance). To add cell level information, add to the Seurat object. If adding
#' feature-level metadata, add to the Assay object (e.g. \code{object[["RNA"]]})
#'
#' @param object An object
#' @param metadata A vector, list, or data.frame with metadata to add
#' @param col.name A name for meta data if not a named list or data.frame
#'
#' @return \code{object} with metadata added
#'
#' @rdname AddMetaData
#' @export AddMetaData
#'
#' @aliases SeuratAccess
#'
#' @concept seurat
#'
#' @examples
#' cluster_letters <- LETTERS[Idents(object = pbmc_small)]
#' names(cluster_letters) <- colnames(x = pbmc_small)
#' pbmc_small <- AddMetaData(
#'   object = pbmc_small,
#'   metadata = cluster_letters,
#'   col.name = 'letter.idents'
#' )
#' head(x = pbmc_small[[]])
#'
AddMetaData <- function(object, metadata, col.name = NULL) {
  UseMethod(generic = 'AddMetaData', object = object)
}

#' Coerce to a \code{Graph} Object
#'
#' Convert a \code{\link[base]{matrix}} (or \code{\link[Matrix]{Matrix}}) to
#' a \code{\link{Graph}} object
#'
#' @template param-dots-ignored
#' @param x The matrix to convert
#'
#' @return A \code{\link{Graph}} object
#'
#' @rdname as.Graph
#' @export as.Graph
#'
#' @family graph
#'
as.Graph <- function(x, ...) {
  UseMethod(generic = "as.Graph", object = x)
}

#' Convert Segmentation Layers
#'
#' @inheritParams CreateCentroids
#' @template param-dots-method
#' @param x An object
#'
#' @return \code{as.Centroids}: A
#' \code{\link[SeuratObject:Centroids-class]{Centroids}} object
#'
#' @export
#'
#' @concept spatial
#'
as.Centroids <- function(x, nsides = NULL, radius = NULL, theta = NULL, ...) {
  UseMethod(generic = "as.Centroids", object = x)
}

#' Coerce to a \code{Neighbor} Object
#'
#' Convert objects to \code{\link{Neighbor}} objects
#'
#' @template param-dots-method
#' @param x An object to convert to \code{\link{Neighbor}}
#'
#' @return A \code{\link{Neighbor}} object
#'
#' @rdname as.Neighbor
#' @export as.Neighbor
#'
#' @concept neighbor
#'
as.Neighbor <- function(x, ...) {
  UseMethod(generic = 'as.Neighbor', object = x)
}

#' @return \code{as.Segmentation}: A
#' \code{\link[SeuratObject:Segmentation-class]{Segmentation}} object
#'
#' @rdname as.Centroids
#' @export
#'
as.Segmentation <- function(x, ...) {
  UseMethod(generic = 'as.Segmentation', object = x)
}

#' Coerce to a \code{Seurat} Object
#'
#' Convert objects to Seurat objects
#'
#' @template param-dots-method
#' @param x An object to convert to class \code{Seurat}
#'
#' @return A \code{\link{Seurat}} object generated from \code{x}
#'
#' @rdname as.Seurat
#' @export as.Seurat
#'
#' @concept seurat
#'
as.Seurat <- function(x, ...) {
  UseMethod(generic = 'as.Seurat', object = x)
}

#' Cast to Sparse
#'
#' Convert dense objects to sparse representations
#'
#' @template param-dots-method
#' @param x An object
#'
#' @return A sparse representation of the input data
#'
#' @rdname as.sparse
#' @export as.sparse
#'
#' @concept utils
#'
as.sparse <- function(x, ...) {
  UseMethod(generic = 'as.sparse', object = x)
}

#' Query Specific Object Types
#'
#' List the names of \code{\link{Assay}}, \code{\link{DimReduc}},
#' \code{\link{Graph}}, \code{\link{Neighbor}} objects
#'
#' @template param-dots-ignored
#' @param object A \code{\link{Seurat}} object
#' @param slot Name of component object to return
#'
#' @return If \code{slot} is \code{NULL}, the names of all component objects
#' in this \code{Seurat} object. Otherwise, the specific object specified
#'
#' @rdname ObjectAccess
#'
#' @export
#'
#' @concept data-access
#'
#' @examples
#' Assays(pbmc_small)
#'
Assays <- function(object, ...) {
  UseMethod(generic = "Assays", object = object)
}

#' Get, Set, and Query Segmentation Boundaries
#'
#' @template param-dots-method
#' @param object An object
#'
#' @name Boundaries
#' @return \code{Boundaries}: The names of all segmentation boundaries present
#' within \code{object}
#'
#' @rdname Boundaries
#' @export
#'
#' @concept spatial
#'
Boundaries <- function(object, ...) {
  UseMethod(generic = 'Boundaries', object = object)
}

#' Cast Assay Layers
#'
#' Cast layers in v5 assays to other classes
#'
#' @param object An object
#' @param to Either a class name or a function that takes a layer and returns
#' the same layer as a new class
#' @param ... If \code{to} is a function, arguments passed to \code{to}
#'
#' @return \code{object} with the layers cast to class specified by \code{to}
#'
#' @export
#'
#' @concept assay5
#'
CastAssay <- function(object, to, ...) {
  UseMethod(generic = 'CastAssay', object = object)
}

#' Cell and Feature Names
#'
#' Get the cell and feature names of an object
#'
#' @template param-dots-method
#' @param x An object
#'
#' @return \code{Cell}: A vector of cell names
#'
#' @rdname Cells
#' @export Cells
#'
#' @concept data-access
#' @family dimnames
#'
#' @examples
#' Cells(x = pbmc_small)
#'
Cells <- function(x, ...) {
  UseMethod(generic = 'Cells', object = x)
}

#' Check Matrix Validity
#'
#' @template param-dots-method
#' @param object A matrix
#' @param checks Type of checks to perform, choose one or more from:
#' \itemize{
#'  \item \dQuote{\code{infinite}}: Emit a warning if any value is infinite
#'  \item \dQuote{\code{logical}}: Emit a warning if any value is a logical
#'  \item \dQuote{\code{integer}}: Emit a warning if any value is \emph{not}
#'   an integer
#'  \item \dQuote{\code{na}}: Emit a warning if any value is an \code{NA}
#'   or \code{NaN}
#' }
#'
#' @return Emits warnings for each test and invisibly returns \code{NULL}
#'
#' @name CheckMatrix
#' @rdname CheckMatrix
#'
#' @keywords internal
#'
#' @export
#'
#' @concept utils
#'
CheckMatrix <- function(object, checks, ...) {
  UseMethod(generic = 'CheckMatrix', object = object)
}

#' Get SeuratCommands
#'
#' Pull information on previously run commands in the Seurat object.
#'
#' @template param-dots-method
#' @param object An object
#'
#' @return Either a SeuratCommand object or the requested parameter value
#'
#' @rdname Command
#' @export Command
#'
#' @concept data-access
#'
Command <- function(object, ...) {
  UseMethod(generic = 'Command', object = object)
}

#' Create a \code{\link[SeuratObject:Centroids-class]{Centroids}} Objects
#'
#' @param coords The coordinates of cell/spot centroids
#' @param nsides The number of sides to represent cells/spots; pass
#' \code{\link[base]{Inf}} to plot as circles
#' @param radius Radius of shapes when plotting
#' @param theta Angle to adjust shapes when plotting
#'
#' @return A \code{\link[SeuratObject:Centroids-class]{Centroids}} object
#'
#' @export
#'
#' @concept spatial
#'
CreateCentroids <- function(coords, nsides, radius, theta) {
  UseMethod(generic = 'CreateCentroids', object = coords)
}

#' Create Spatial Coordinates
#'
#' @template param-dots-method
#' @param coords Spatial coordinates
#'
#' @return A \code{\link{FOV}} object
#'
#' @export
#'
#' @concept spatial
#'
#' @seealso \code{\link{FOV-class}}
#'
CreateFOV <- function(coords, ...) {
  UseMethod(generic = 'CreateFOV', object = coords)
}

#' Create a \code{\link{Molecules}} Object
#'
#' @template param-dots-method
#' @param coords Spatial coordinates for molecules; should be a data frame
#' with three columns:
#' \itemize{
#'  \item \dQuote{\code{x}}: x-coordinates for each molecule
#'  \item \dQuote{\code{y}}: y-coordinates for each molecule
#'  \item \dQuote{\code{gene}}: gene name for each molecule
#' }
#'
#' @return A \code{\link{Molecules}} object
#'
#' @export
#'
#' @concept spatial
#'
CreateMolecules <- function(coords, ...) {
  UseMethod(generic = 'CreateMolecules', object = coords)
}

#' Create a \code{\link[SeuratObject:Segmentation-class]{Segmentation}} Objects
#'
#' @param coords The coordinates of cell segmentations
#' @param compact Logical indicating whether or not the object should only store segmentation data
#' in the \code{sf.data} slot; see \link{Segmentation-class} for details.
#'
#' @return A \code{\link[SeuratObject:Segmentation-class]{Segmentation}} object
#'
#' @export
#'
#' @concept spatial
#'
CreateSegmentation <- function(coords, compact = FALSE) {
  UseMethod(generic = 'CreateSegmentation', object = coords)
}

#' Create a \code{Seurat} object
#'
#' Create a \code{Seurat} object from raw data
#'
#' @inheritParams CreateAssayObject
#' @template param-dots-method
#' @param counts Either a \code{\link[base]{matrix}}-like object with
#' unnormalized data with cells as columns and features as rows or an
#' \code{\link{Assay}}-derived object
#' @param project \link{Project} name for the \code{Seurat} object
#' @param assay Name of the initial assay
#' @param names.field For the initial identity class for each cell, choose this
#' field from the cell's name. E.g. If your cells are named as
#' BARCODE_CLUSTER_CELLTYPE in the input matrix, set \code{names.field} to 3 to
#' set the initial identities to CELLTYPE.
#' @param names.delim For the initial identity class for each cell, choose this
#' delimiter from the cell's column name. E.g. If your cells are named as
#' BARCODE-CLUSTER-CELLTYPE, set this to \dQuote{-} to separate the cell name
#' into its component parts for picking the relevant field.
#' @param meta.data Additional cell-level metadata to add to the Seurat object.
#' Should be a \code{\link[base]{data.frame}} where the rows are cell names and
#' the columns are additional metadata fields. Row names in the metadata need
#' to match the column names of the counts matrix.
#'
#' @note In previous versions (<3.0), this function also accepted a parameter to
#' set the expression threshold for a \sQuote{detected} feature (gene). This
#' functionality has been removed to simplify the initialization
#' process/assumptions. If you would still like to impose this threshold for
#' your particular dataset, simply filter the input expression matrix before
#' calling this function.
#'
#' @return A \code{\link{Seurat}} object
#'
#' @rdname CreateSeuratObject
#' @export
#'
#' @concept seurat
#'
#' @examples
#' \dontrun{
#' pbmc_raw <- read.table(
#'   file = system.file('extdata', 'pbmc_raw.txt', package = 'Seurat'),
#'   as.is = TRUE
#' )
#' pbmc_small <- CreateSeuratObject(counts = pbmc_raw)
#' pbmc_small
#' }
#'
CreateSeuratObject <- function(
  counts,
  assay = 'RNA',
  names.field = 1,
  names.delim = '_',
  meta.data = NULL,
  project = 'CreateSeuratObject',
  ...
) {
  UseMethod(generic = 'CreateSeuratObject', object = counts)
}

#' Crop Coordinates
#'
#' @template param-dots-method
#' @param object An object
#' @param x,y Range to crop x/y limits to; if \code{NULL}, uses full range of
#' \code{x}/\code{y}
#' @param coords Coordinate system to execute crop; choose from:
#' \itemize{
#'  \item \dQuote{\code{plot}}: Coordinates as shown when plotting
#'  \item \dQuote{\code{tissue}}: Coordinates from
#'   \code{\link{GetTissueCoordinates}}
#' }
#'
#' @return \code{object} cropped to the region specified by \code{x}
#' and \code{y}
#'
#' @export
#'
#' @concept spatial
#'
Crop <- function(object, x = NULL, y = NULL, coords = c('plot', 'tissue'), ...) {
  UseMethod(generic = 'Crop', object = object)
}

#' Default Assay
#'
#' Get and set the default assay
#'
#' @template param-dots-method
#' @param object An object
#'
#' @return \code{DefaultAssay}: The name of the default assay
#'
#' @rdname DefaultAssay
#' @export DefaultAssay
#'
#' @concept data-access
#'
DefaultAssay <- function(object, ...) {
  UseMethod(generic = 'DefaultAssay', object = object)
}

#' @param value Name of assay to set as default
#'
#' @return \code{DefaultAssay<-}: An object with the default assay updated
#'
#' @rdname DefaultAssay
#' @export DefaultAssay<-
#'
"DefaultAssay<-" <- function(object, ..., value) {
  UseMethod(generic = 'DefaultAssay<-', object = object)
}

#' @return \code{DefaultBoundary}: The name of the default
#' segmentation boundary
#'
#' @rdname Boundaries
#'
#' @export
#'
DefaultBoundary <- function(object) {
  UseMethod(generic = 'DefaultBoundary', object = object)
}

#' @param value The name of a segmentation boundary to set as default
#'
#' @return \code{DefaultBoundary<-}: \code{object} with the default
#' segmentation boundary set to \code{value}
#'
#' @rdname Boundaries
#'
#' @export
#'
"DefaultBoundary<-" <- function(object, ..., value) {
  UseMethod(generic = 'DefaultBoundary<-', object = object)
}

#' Set Default Dimensionality Reduction
#' Set the default dimensionality reduction for a specific assay in a Seurat object.
#' @param object A Seurat object.
#' @param value Character string specifying the name of the dimensionality reduction to set as default.
#' Set to NULL to clear the default for the current assay.
#' @return A Seurat object with the default dimensionality reduction updated.
#' @details
#' This function stores the default dimensionality reduction on a per-assay basis.
#' This function only needs to be run if users want to override the default DimReduc selection logic that
#' Seurat employs.
#' When `DefaultDimReduc()` is called, it will return the assay-specific default if one has been set.
#' If none was explicitly set then default Seurat logic will be used to select default.
#'
#' @rdname DefaultDimReduc
#' @export DefaultDimReduc<-
#'
"DefaultDimReduc<-" <- function(object, ..., value) {
  UseMethod(generic = 'DefaultDimReduc<-', object = object)
}


#' Get and Set the Default FOV
#'
#' @template param-dots-method
#' @param object A \code{\link{Seurat}} Object
#'
#' @return \code{DefaultFOV}: The name of the default \code{\link{FOV}}
#'
#' @name DefaultFOV
#' @rdname DefaultFOV
#'
#' @export
#'
#' @concept spatial
#'
DefaultFOV <- function(object, ...) {
  UseMethod(generic = 'DefaultFOV', object = object)
}

#' @param value The name of the \code{\link{FOV}} to set as the default
#'
#' @return \code{DefaultFOV<-}: \code{object} with the default FOV set
#' to \code{value}
#'
#' @rdname DefaultFOV
#'
#' @export
#'
"DefaultFOV<-" <- function(object, ..., value) {
  UseMethod(generic = 'DefaultFOV<-', object = object)
}

#' Default Layer
#'
#' Get and set the default layer
#'
#' @template param-dots-method
#' @param object An object
#'
#' @return \code{DefaultLayer}: The name of the default layer
#'
#' @rdname DefaultLayer
#' @export DefaultLayer
#'
#' @concept assay5
#'
DefaultLayer <- function(object, ...) {
  UseMethod(generic = 'DefaultLayer', object = object)
}

#' @param value Name of layer to set as default
#'
#' @return \code{DefaultLayer<-}: An object with the default layer updated
#'
#' @rdname DefaultLayer
#' @export DefaultLayer<-
#'
"DefaultLayer<-" <- function(object, ..., value) {
  UseMethod(generic = 'DefaultLayer<-', object = object)
}

#' Get the Neighbor nearest neighbors distance matrix
#'
#' @template param-dots-method
#' @param object An object
#'
#' @return The distance matrix
#'
#' @rdname Distances
#' @export Distances
#'
#' @concept data-access
#'
Distances <- function(object, ...) {
  UseMethod(generic = 'Distances', object = object)
}

#' Get Cell Embeddings
#'
#' @template param-dots-method
#' @param object An object
#'
#' @return The embeddings matrix
#'
#' @rdname Embeddings
#' @export Embeddings
#'
#' @concept data-access
#'
Embeddings <- function(object, ...) {
  UseMethod(generic = 'Embeddings', object = object)
}

#' Access cellular data
#'
#' Retrieves data (feature expression, PCA scores, metrics, etc.) for a set
#' of cells in a Seurat object
#'
#' @template param-dots-method
#' @param object An object
#'
#' @export FetchData
#'
#' @concept data-access
#'
FetchData <- function(object, ...) {
  UseMethod(generic = 'FetchData', object = object)
}

#' @return \code{Features}: A vector of feature names
#'
#' @rdname Cells
#' @export Features
#'
Features <- function(x, ...) {
  UseMethod(generic = 'Features', object = x)
}

#' Get and Set Assay Data
#'
#' General accessor and setter functions for \code{\link{Assay}} objects.
#' \code{GetAssayData} can be used to pull information from any of the
#' expression matrices (eg. \dQuote{counts}, \dQuote{data}, or
#' \dQuote{scale.data}). \code{SetAssayData} can be used to replace one of these
#' expression matrices
#'
#' @template param-dots-method
#' @param object An object
#' @param layer Name of layer to get or set
#' @param slot \Sexpr[stage=build,results=rd]{lifecycle::badge("deprecated")} Specific assay data to get or set
#'
#' @return \code{GetAssayData}: returns the specified assay data
#'
#' @template lifecycle-superseded
#'
#' @section Lifecycle:
#'
#' \code{GetAssayData} and \code{SetAssayData} have been superseded. To fetch
#' expression matrices, use \code{\link{LayerData}}; to set expression data,
#' use \code{\link{LayerData<-}}
#'
#' @name AssayData
#' @rdname AssayData
#' @export GetAssayData
#'
#' @order 1
#'
#' @concept data-access
#'
GetAssayData <- function(object, ...) {
  UseMethod(generic = 'GetAssayData', object = object)
}

#' Get image data
#'
#' @template param-dots-method
#' @param object An object
#' @param mode How to return the image; should accept one of \dQuote{grob},
#' \dQuote{raster}, \dQuote{plotly}, or \dQuote{raw}
#'
#' @return Image data, varying depending on the value of \code{mode}:
#' \describe{
#'  \item{\dQuote{grob}}{
#'   An object representing image data inheriting from \code{grob} objects
#'   (eg. \code{rastergrob})
#'  }
#'  \item{\dQuote{raster}}{An object of class \code{raster}}
#'  \item{\dQuote{plotly}}{
#'   A list with image data suitable for Plotly rendering, see
#'   \code{\link[plotly:layout]{plotly::layout}} for more details
#'  }
#'  \item{\dQuote{raw}}{The raw image data as stored in the object}
#' }
#'
#' @seealso \code{\link[plotly]{layout}}
#'
#' @rdname GetImage
#' @export GetImage
#'
#' @concept data-access
#'
GetImage <- function(object, mode = c('grob', 'raster', 'plotly', 'raw'), ...) {
  mode <- mode[1]
  mode <- match.arg(arg = mode)
  UseMethod(generic = 'GetImage', object = object)
}

#' Get tissue coordinates
#' 
#' Retrieve tissue coordinates from spatial objects. 
#'
#' Spatial classes may have specific implementations; 
#' refer to those documentation pages for more details.
#'
#' @template param-dots-method
#' @param object An object
#'
#' @return A data frame with tissue coordinates
#'
#' @rdname GetTissueCoordinates
#' @export GetTissueCoordinates
#'
#' @concept data-access
#'
GetTissueCoordinates <- function(object, ...) {
  UseMethod(generic = 'GetTissueCoordinates', object = object)
}

#' Highly Variable Features
#'
#' Get and set variable feature information for an \code{\link{Assay}} object.
#' \code{HVFInfo} and \code{VariableFeatures} utilize generally variable
#' features, while \code{SVFInfo} and \code{SpatiallyVariableFeatures} are
#' restricted to spatially variable features
#'
#' @template param-dots-method
#' @param object An object
#' @param method Which method to pull. For \code{HVFInfo} and
#' \code{VariableFeatures}, choose one from one of the
#' following:
#' \itemize{
#'  \item \dQuote{vst}
#'  \item \dQuote{sctransform} or \dQuote{sct}
#'  \item \dQuote{mean.var.plot}, \dQuote{dispersion}, \dQuote{mvp}, or
#'   \dQuote{disp}
#' }
#' For \code{SVFInfo} and \code{SpatiallyVariableFeatures}, choose from:
#' \itemize{
#'  \item \dQuote{markvariogram}
#'  \item \dQuote{moransi}
#' }
#' @param status Add variable status to the resulting data frame
#' @param selection.method \Sexpr[stage=build,results=rd]{lifecycle::badge("deprecated")}
#'
#' @return \code{HVFInfo}: A data frame with feature means, dispersion, and
#' scaled dispersion
#'
#' @rdname VariableFeatures
#' @export HVFInfo
#'
#' @order 1
#'
#' @concept data-access
#'
HVFInfo <- function(object, method, status = FALSE, ...) {
  UseMethod(generic = 'HVFInfo', object = object)
}

#' Get, set, and manipulate an object's identity classes
#'
#' @param x,object An object
#' @param ... Arguments passed to other methods; for \code{RenameIdents}: named
#' arguments as \code{old.ident = new.ident}; for \code{ReorderIdent}: arguments
#' passed on to \code{\link{FetchData}}
#'
#' @return \code{Idents}: The cell identities
#'
#' @rdname Idents
#' @export Idents
#'
#' @concept seurat
#'
#' @examples
#' # Get cell identity classes
#' Idents(pbmc_small)
#'
Idents <- function(object, ... ) {
  UseMethod(generic = 'Idents', object = object)
}

#' @param value The name of the identities to pull from object metadata or the
#' identities themselves
#'
#' @return \code{Idents<-}: \code{object} with the cell identities changed
#'
#' @rdname Idents
#' @export Idents<-
#'
#' @examples
#' # Set cell identity classes
#' # Can be used to set identities for specific cells to a new level
#' Idents(pbmc_small, cells = 1:4) <- 'a'
#' head(Idents(pbmc_small))
#'
#' # Can also set idents from a value in object metadata
#' colnames(pbmc_small[[]])
#' Idents(pbmc_small) <- 'RNA_snn_res.1'
#' levels(pbmc_small)
#'
"Idents<-" <- function(object, ..., value) {
  UseMethod(generic = 'Idents<-', object = object)
}

#' Get Neighbor algorithm index
#'
#' @template param-dots-method
#' @param object An object
#'
#' @return Returns the value in the alg.idx slot of the Neighbor object
#'
#' @rdname NNIndex
#' @export Index
#'
#' @concept data-access
#'
Index <- function(object, ...) {
  UseMethod(generic = "Index", object = object)
}

#' @param value The index to store
#'
#' @return \code{Idents<-}: A Neighbor object with the index stored
#'
#' @rdname NNIndex
#' @export Index<-
#'
"Index<-" <- function(object, ..., value) {
  UseMethod(generic = 'Index<-', object = object)
}

#' Get Neighbor nearest neighbor index matrices
#'
#' @template param-dots-method
#' @param object An object
#'
#' @return A matrix with the nearest neighbor indices
#'
#' @rdname Indices
#' @export Indices
#'
#' @concept data-access
#'
Indices <- function(object, ...) {
  UseMethod(generic = "Indices", object = object)
}

#' Is an object global/persistent?
#'
#' Typically, when removing \code{Assay} objects from an \code{Seurat} object,
#' all associated objects (eg. \code{DimReduc}, \code{Graph}, and
#' \code{SeuratCommand} objects)
#' are removed as well. If an associated object is marked as global/persistent,
#' the associated object will remain even if its original assay was deleted
#'
#' @template param-dots-method
#' @param object An object
#'
#' @return \code{TRUE} if the object is global/persistent otherwise \code{FALSE}
#'
#' @rdname IsGlobal
#' @export IsGlobal
#'
#' @concept data-access
#'
#' @examples
#' IsGlobal(pbmc_small[['pca']])
#'
IsGlobal <- function(object, ...) {
  UseMethod(generic = 'IsGlobal', object = object)
}

#' Check if a matrix is empty
#'
#' Takes a matrix and asks if it's empty (either 0x0 or 1x1 with a value of NA)
#'
#' @param x A matrix
#'
#' @return Whether or not \code{x} is empty
#'
#' @rdname IsMatrixEmpty
#' @export IsMatrixEmpty
#'
#' @concept utils
#'
#' @seealso \code{\link{EmptyMatrix}()}
#'
#' @examples
#' IsMatrixEmpty(new("matrix"))
#' IsMatrixEmpty(matrix())
#' IsMatrixEmpty(matrix(1:3))
#'
IsMatrixEmpty <- function(x) {
  UseMethod(generic = 'IsMatrixEmpty', object = x)
}

#' Split and Join Layers Together
#'
#' @param object An object
#' @template param-dots-method
#'
#' @return \code{object} with the layers specified joined
#'
#' @rdname SplitLayers
#' @export JoinLayers
#'
#' @concept assay5
#'
JoinLayers <- function(object, ...) {
  UseMethod(generic = 'JoinLayers', object = object)
}

#' Get and set JackStraw information
#'
#' @template param-dots-method
#' @param object An object
#'
#' @return \code{JS}: either a \code{\link{JackStrawData}} object or the
#' specified jackstraw data
#'
#' @rdname JS
#' @export JS
#'
#' @concept jackstraw
#'
JS <- function(object, ...) {
  UseMethod(generic = 'JS', object = object)
}

#' @param value JackStraw information
#'
#' @return \code{JS<-}: \code{object} with the update jackstraw information
#'
#' @rdname JS
#' @export JS<-
#'
"JS<-" <- function(object, ..., value) {
  UseMethod(generic = 'JS<-', object = object)
}

#' Get and set object keys
#'
#' @template param-dots-method
#' @param object An object
#'
#' @return \code{Key}: the object key
#'
#' @rdname Key
#' @export Key
#'
#' @concept data-access
#'
Key <- function(object, ...) {
  UseMethod(generic = 'Key', object = object)
}

#' @return \code{Keys}: a named vector of keys of sub-objects
#'
#' @rdname Key
#' @export
#'
Keys <- function(object, ...) {
  UseMethod(generic = 'Keys', object = object)
}

#' @param value Key value
#'
#' @return \code{Key<-}: \code{object} with an updated key
#'
#' @rdname Key
#' @export Key<-
#'
#' @concept data-access
#'
"Key<-" <- function(object, ..., value) {
  UseMethod(generic = 'Key<-', object = object)
}

#' Query and Manipulate Assay Layers
#'
#' @template param-dots-method
#' @param object An object
#' @param layer Name of layer to fetch or set
#' @param slot \Sexpr[stage=build,results=rd]{lifecycle::badge("deprecated")}
#'
#' @return \code{LayerData}: the layer data for \code{layer} from \code{object}
#'
#' @rdname Layers
#' @export LayerData
#'
LayerData <- function(object, layer, ...) {
  UseMethod(generic = 'LayerData', object = object)
}

#' @param value New two-dimensional data to be added as a layer
#'
#' @return \code{Layer<-}: \code{object} with \code{value} added as a layer
#' named \code{layer}
#'
#' @rdname Layers
#' @export LayerData<-
#'
"LayerData<-" <- function(object, layer, ..., value) {
  UseMethod(generic = 'LayerData<-', object = object)
}

#' @return \code{Layers}: the names of the layers present in \code{object}
#'
#' @rdname Layers
#' @export Layers
#'
#' @concept data-access
#'
Layers <- function(object, ...) {
  UseMethod(generic = 'Layers', object = object)
}

#' Get and set feature loadings
#'
#' @template param-dots-method
#' @param object An object
#'
#' @return \code{Loadings}: the feature loadings for \code{object}
#'
#' @rdname Loadings
#' @export Loadings
#'
#' @concept data-access
#'
Loadings <- function(object, ...) {
  UseMethod(generic = 'Loadings', object = object)
}

#' @param value Feature loadings to add
#'
#' @return \code{Loadings<-}: \code{object} with the updated loadings
#'
#' @rdname Loadings
#' @export Loadings<-
#'
"Loadings<-" <- function(object, ..., value) {
  UseMethod(generic = 'Loadings<-', object = object)
}

#' Match Cells
#'
#' @param new A vector of new cells
#' @param orig A vector of existing cells
#' @param ordered Sort the result to the same order as \code{orig}
#'
#' @return A numeric vector with new cells in order of the original cells; if
#' no match can be found, returns \code{NULL}
#'
#' @keywords internal
#'
#' @export
#'
#' @concept utils
#'
MatchCells <- function(new, orig, ordered = FALSE) {
  if (!is.character(x = orig)) {
    stop("'orig' must be a character vector")
  }
  UseMethod(generic = 'MatchCells', object = new)
}

#' Get and set miscellaneous data
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Miscellaneous data
#'
#' @rdname Misc
#' @export Misc
#'
#' @concept data-access
#'
Misc <- function(object, ...) {
  UseMethod(generic = 'Misc', object = object)
}

#' @param value Data to add
#'
#' @return An object with miscellaneous data added
#'
#' @rdname Misc
#' @export Misc<-
#'
"Misc<-" <- function(object, ..., value) {
  UseMethod(generic = 'Misc<-', object = object)
}

#' @return \code{Molecules}: The names of all molecule sets present within
#' \code{object}
#'
#' @rdname Boundaries
#' @export
#'
Molecules <- function(object, ...) {
  UseMethod(generic = 'Molecules', object = object)
}

#' Overlay \code{Spatial} Objects Over One Another
#'
#' Create an overlay of some query spatial object (\code{x}) against some
#' target object (\code{y}). Basically, find all components of a query that
#' fall within the bounds of a target spatial region
#'
#' @template param-dots-ignored
#' @param x Query \code{Spatial} object
#' @param y Target \code{Spatial} object
#' @param invert Invert the overlay and return only the components of \code{x}
#' that fall \emph{outside} the bounds of \code{y}
#'
#' @return \code{x} with only the components that fall within the
#' bounds of \code{y}
#'
#' @templateVar pkg sf
#' @template note-reqdpkg
#'
#' @export
#'
#' @concept spatial
#'
setGeneric(
  name = 'Overlay',
  def = function(x, y, invert = FALSE, ...) {
    standardGeneric(f = 'Overlay')
  },
  signature = c('x', 'y')
)

#' Get and set project information
#'
#' @template param-dots-method
#' @param object An object
#'
#' @return Project information
#'
#' @rdname Project
#' @export Project
#'
#' @concept seurat
#'
Project <- function(object, ...) {
  UseMethod(generic = 'Project', object = object)
}

#' @param value Project information to set
#'
#' @return An object with project information added
#'
#' @rdname Project
#' @export Project<-
#'
"Project<-" <- function(object, ..., value) {
  UseMethod(generic = 'Project<-', object = object)
}

#' Get the spot radius from an image
#'
#' @param object An image object
#' @param ... Arguments passed to other methods
#'
#' @return The radius size
#'
#' @rdname Radius
#' @export Radius
#'
#' @concept spatialimage
#'
Radius <- function(object, ...) {
  UseMethod(generic = 'Radius', object = object)
}

#' Rename cells
#'
#' Change the cell names in all the different parts of an object. Can be useful
#' before combining multiple objects.
#'
#' @template param-dots-method
#' @param object An object
#'
#' @return An object with new cell names
#'
#' @rdname RenameCells
#' @export RenameCells
#'
#' @concept seurat
#'
RenameCells <- function(object, ...) {
  UseMethod(generic = 'RenameCells', object = object)
}

#' @return \code{RenameIdents}: An object with selected identity classes renamed
#'
#' @rdname Idents
#' @export RenameIdents
#' @aliases RenameIdent
#'
#' @examples
#' # Rename cell identity classes
#' # Can provide an arbitrary amount of idents to rename
#' levels(pbmc_small)
#' pbmc_small <- RenameIdents(pbmc_small, '0' = 'A', '2' = 'C')
#' levels(pbmc_small)
#'
RenameIdents <- function(object, ...) {
  UseMethod(generic = 'RenameIdents', object = object)
}

#' @param var Feature or variable to order on
#'
#' @return \code{ReorderIdent}: An object with
#'
#' @rdname Idents
#' @export ReorderIdent
#' @aliases ReorderIdent
#'
#' @examples
#' \dontrun{
#' head(Idents(pbmc_small))
#' pbmc_small <- ReorderIdent(pbmc_small, var = 'PC_1')
#' head(Idents(pbmc_small))
#' }
#'
ReorderIdent <- function(object, var, ...) {
  UseMethod(generic = 'ReorderIdent', object = object)
}

#' S4/List Conversion
#'
#' Convert S4 objects to lists and vice versa. Useful for declassing an S4
#' object while keeping track of it's class using attributes (see section
#' \strong{S4 Class Definition Attributes} below for more details). Both
#' \code{ListToS4} and \code{S4ToList} are recursive functions, affecting
#' all lists/S4 objects contained as sub-lists/sub-objects
#'
#' @param object An S4 object
#' @param x A list with an S4 class definition
#' (\dQuote{\code{classDef}}) attribute
#'
#' @return \code{S4ToList}: A list with an S4 class definition attribute
#'
#' @section S4 Class Definition Attributes:
#' S4 classes are scoped to the package and class name. In order to properly
#' track which class a list is generated from in order to build a new one,
#' these function use an \code{\link[base:attr]{attribute}} to denote the
#' class name and package of origin. This attribute is stored as 1-length
#' character vector named \dQuote{\code{classDef}} and takes the form
#' of \dQuote{\code{package:class}}
#'
#' @name s4list
#' @rdname s4list
#'
#' @keywords internal
#'
#' @export
#'
#' @concept utils
#' @family s4list
#'
#' @examples
#' # Turn an S4 object into a list
#' pbmc.list <- S4ToList(pbmc_small)
#' class(pbmc.list)
#' attributes(pbmc.list)
#'
#' str(pbmc.list$reductions)
#'
S4ToList <- function(object) {
  if (!(isS4(object) || is_bare_list(x = object))) {
    return(object)
  }
  UseMethod(generic = 'S4ToList', object = object)
}

#' @param new.data New assay data to add
#'
#' @return \code{SetAssayData}: \code{object} with the assay data set
#'
#' @rdname AssayData
#' @export SetAssayData
#'
#' @order 2
#'
SetAssayData <- function(object, layer, new.data, slot = deprecated(), ...) {
  UseMethod(generic = 'SetAssayData', object = object)
}

#' @return \code{SetIdent}: An object with new identity classes set
#'
#' @rdname Idents
#' @export SetIdent
#'
#' @examples
#' # Set cell identity classes using SetIdent
#' cells.use <- WhichCells(pbmc_small, idents = '1')
#' pbmc_small <- SetIdent(pbmc_small, cells = cells.use, value = 'B')
#'
SetIdent <- function(object, ...) {
  UseMethod(generic = 'SetIdent', object = object)
}

#' Simplify Geometry
#'
#' @param coords ...
#'
#' @return A simplified version of \code{coords}
#'
#' @export
#'
#' @concept spatial
#'
Simplify <- function(coords, tol, topologyPreserve = TRUE) {
  UseMethod(generic = 'Simplify', object = coords)
}

#' @return \code{SpatiallyVariableFeatures}: a character vector of the spatially
#' variable features
#'
#' @rdname VariableFeatures
#' @export SpatiallyVariableFeatures
#'
#' @order 5
#'
SpatiallyVariableFeatures <- function(object, method, ...) {
  UseMethod(generic = 'SpatiallyVariableFeatures', object = object)
}

# @rdname SplitLayers
# @export SplitLayers
#
# @order 1
#
SplitLayers <- function(object, ...) {
  UseMethod(generic = 'SplitLayers', object = object)
}

#' @return \code{StashIdent}: An object with the identities stashed
#'
#' @rdname Idents
#' @export StashIdent
#'
#' @examples
#' head(pbmc_small[[]])
#' pbmc_small <- StashIdent(pbmc_small, save.name = 'idents')
#' head(pbmc_small[[]])
#'
StashIdent <- function(object, save.name, ...) {
  UseMethod(generic = 'StashIdent', object = object)
}

#' Get the standard deviations for an object
#'
#' @template param-dots-method
#' @param object An object
#'
#' @return The standard deviations
#'
#' @rdname Stdev
#' @export Stdev
#'
#' @concept data-access
#'
Stdev <- function(object, ...) {
  UseMethod(generic = 'Stdev', object = object)
}

#' Stitch Matrices Together
#'
#' @template param-dots-method
#' @param x A matrix
#' @param y One or more matrices of the same class or coercible to the
#' same class as \code{x}
#' @param rowmap,colmap \code{\link{LogMap}s} describing the row and cell
#' membership of each matrix; the \code{LogMap} entries are assumed to be in
#' the order of \code{c(x, y)}
#'
#' @return A single matrix of type \code{class(x)} consisting of all values
#' in component matrices
#'
#' @export
#'
#' @concept utils
#'
StitchMatrix <- function(x, y, rowmap, colmap, ...) {
  if (!inherits(x = rowmap, what = 'LogMap')) {
    abort(message = "'rowmap' must be a 'LogMap'")
  } else if (!inherits(x = colmap, what = 'LogMap')) {
    abort(message = "'colmap' must be a 'LogMap'")
  }
  UseMethod(generic = 'StitchMatrix', object = x)
}

#' @return \code{SVFInfo}: a data frame with the spatially variable features
#'
#' @rdname VariableFeatures
#' @export SVFInfo
#'
#' @order 4
#'
SVFInfo <- function(object, method, status, ...) {
  UseMethod(generic = 'SVFInfo', object = object)
}

#' Get the offset angle
#'
#' @param object An object
#'
#' @rdname Theta
#' @export
#'
#' @concept spatial
#'
Theta <- function(object) {
  UseMethod(generic = 'Theta', object = object)
}

#' Get and Set Additional Tool Data
#'
#' Use \code{Tool} to get tool data. If no additional arguments are provided,
#' will return a vector with the names of tools in the object.
#'
#' @template param-dots-method
#' @param object An object
#'
#' @return If no additional arguments, returns the names of the tools in the
#' object; otherwise returns the data placed by the tool requested
#'
#'@note For developers: set tool data using \code{Tool<-}. \code{Tool<-} will
#'automatically set the name of the tool to the function that called
#'\code{Tool<-}, so each function gets one entry in the tools list and cannot
#'overwrite another function's entry. The automatic naming will also remove any
#'method identifiers (eg. \code{RunPCA.Seurat} will become \code{RunPCA});
#'please plan accordingly
#'
#' @rdname Tool
#' @export Tool
#'
#' @aliases Tools
#'
#' @concept data-access
#'
#' @examples
#' # Example function that adds unstructured data to tools
#' MyTool <- function(object) {
#'   sample.tool.output <- matrix(rnorm(n = 16), nrow = 4)
#'   # Note: `Tool<-` must be called from within a function
#'   # and the name of the tool will be generated from the function name
#'   Tool(object) <- sample.tool.output
#'   return(object)
#' }
#'
#' # Run our tool
#' set.seed(42L)
#' pbmc_small <- MyTool(pbmc_small)
#'
#' # Get a list of tools run
#' Tool(pbmc_small)
#'
#' # Access specific tool data
#' Tool(pbmc_small, slot = "MyTool")
#'
Tool <- function(object, ...) {
  UseMethod(generic = 'Tool', object = object)
}

#' @param value Information to be added to tool list
#'
#' @rdname Tool
#' @export Tool<-
#'
"Tool<-" <- function(object, ..., value) {
  UseMethod(generic = 'Tool<-', object = object)
}

#' @return \code{VariableFeatures}: a vector of the variable features
#'
#' @rdname VariableFeatures
#' @export VariableFeatures
#'
#' @order 2
#'
VariableFeatures <- function(object, method = NULL, ...) {
  UseMethod(generic = 'VariableFeatures', object = object)
}

#' @param value A character vector of variable features
#'
#' @order 3
#'
#' @rdname VariableFeatures
#' @export VariableFeatures<-
#'
"VariableFeatures<-" <- function(object, ..., value) {
  UseMethod(generic = 'VariableFeatures<-', object = object)
}

#' Get Version Information
#'
#' @template param-dots-method
#' @param object An object
#'
#' @rdname Version
#' @export Version
#'
#' @concept data-access
#'
#' @examples
#' Version(pbmc_small)
#'
Version <- function(object, ...) {
  UseMethod(generic = "Version", object = object)
}

#' Identify cells matching certain criteria
#'
#' Returns a list of cells that match a particular set of criteria such as
#' identity class, high/low values for particular PCs, etc.
#'
#' @template param-dots-method
#' @param object An object
#'
#' @return A vector of cell names
#'
#' @rdname WhichCells
#' @export WhichCells
#'
#' @concept data-access
#'
#' @seealso \code{\link{FetchData}}
#'
#' @examples
#' WhichCells(pbmc_small, idents = 2)
#' WhichCells(pbmc_small, expression = MS4A1 > 3)
#' levels(pbmc_small)
#' WhichCells(pbmc_small, idents = c(1, 2), invert = TRUE)
#'
WhichCells <- function(object, ...) {
  UseMethod(generic = 'WhichCells', object = object)
}
