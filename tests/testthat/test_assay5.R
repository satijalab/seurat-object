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
