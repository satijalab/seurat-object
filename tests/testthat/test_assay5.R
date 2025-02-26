#' Returns a random counts matrix.
get_test_assay <- function(ncells, nfeatures, assay_version) {
  # Use the `assay_version` param to choose the correct assay builder.
  create_assay <- switch(assay_version,
                         v3 = CreateAssayObject,
                         v5 = CreateAssay5Object,
                         stop("`assay_version` should be one of 'v3', 'v5'")
  )
  
  # Populate a `nfeatures` x `ncells` matrix with zeros.
  counts <- matrix(0, ncol = ncells, nrow = nfeatures)
  
  # Assign column and row names to the matrix to label cells and genes.
  colnames(counts) <- paste0("cell", seq(ncol(counts)))
  row.names(counts) <- paste0("gene", seq(nrow(counts)))
  
  # Convert `counts` to a `dgCMatrix`.
  counts_sparse <- as.sparse(counts)
  # Build an assay of the specified type.
  assay <- create_assay(counts_sparse)
  
  return(assay)
}

#' Mocks "highly variable feature" annotations and adds them to the 
#' feature-level metadata of `assay`.
add_hvf_info <- function(assay, nfeatures, method_name, layer_name) {
  variable_features <- sample(
    rownames(assay), 
    size = nfeatures, 
    replace = FALSE
  )

  all_features <- rownames(assay)
  constant_features <- setdiff(all_features, variable_features)
  
  # Add a column to the assay's feature-level metadata indicating which
  # features are variable.
  is_variable <- all_features %in% variable_features
  names(is_variable) <- all_features
  variable_column <- paste("vf", method_name, layer_name, "variable", sep = "_")
  assay[[variable_column]] <- is_variable
  
  # Add a column to the assay's feature-level metadata indicating the order
  # (i.e. "rank") that each variable feature was selected in.
  variable_rank <- seq_along(variable_features)
  names(variable_rank) <- variable_features
  constant_rank <- seq_along(constant_features) + length(variable_features)
  names(constant_rank) <- constant_features
  feature_rank <- c(variable_rank, constant_rank)
  rank_column <- paste("vf", method_name, layer_name, "rank", sep = "_")
  assay[[rank_column]] <- feature_rank

  # Add a column to the assay's feature-level metadata containing a random
  # values in reverse order from `variable_rank` (i.e. the lowest rank
  # has the highest value).
  feature_values <- runif(length(all_features))
  feature_values <- sort(feature_values, decreasing = TRUE)
  names(feature_values) <- names(feature_rank)
  value_column <- paste("vf", method_name, layer_name, "value", sep = "_")
  assay[[value_column]] <- feature_values
  
  return(assay)
}
