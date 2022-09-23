#' @slot layers A named list containing expression matrices; each matrix should
#' be a two-dimensional object containing some subset of cells and features
#' defined in the \code{cells} and \code{features} slots. Cell and feature
#' membership is recorded in the \code{cells} and \code{features} slots,
#' respectively
#' @slot cells A \link[LogMap]{logical mapping} of cell names and layer
#' membership; this map contains all the possible cells that this assay can
#' contain. New layers must have some subset of cells present in this map
#' @slot features A \link[LogMap]{logical mapping} of feature names and layer
#' membership; this map contains all the possible features that this assay can
#' contain. New layers must have some subset of features present in this map
#' @slot default A one-length integer with the end index of the
#' \link[DefaultLayer]{default layer}; the default layer be all layers up to
#' and including the layer at index \code{default}
#' @slot assay.orig Original assay that this assay is based off of;
#' used to track assay provenance
#' @slot meta.data A \link[base:data.frame]{data frame} with feature-level
#' meta data; should have the same number of rows as \code{features}
