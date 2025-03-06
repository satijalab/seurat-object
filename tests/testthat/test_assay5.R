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
add_hvf_info <- function(
  assay, 
  nfeatures = NULL,
  features = NULL,
  method_name, 
  layer_name
) {
  if (is.null(nfeatures) & is.null(features)) {
    # Ensure that one of `nfeatures` or `features` is set.
    stop("One of `nfeatures` or `features` must be provided.")
  } else if (!is.null(nfeatures) & !is.null(features)) {
    # Ensure that only one of `nfeatures` or `features` is set.
    stop("Only one of `nfeatures` or `feature` may be provided.")
  } else if (!is.null(nfeatures)) {
    # If `nfeatures` is provided, randomly sample from the `assay`'s features.
    variable_features <- sample(
      rownames(assay),
      size = nfeatures,
      replace = FALSE
    )
  } else {
    # If `features` was provided, use it.
    variable_features <- features
  }

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

context("HVFInfo")

test_that("`HVFInfo.Assay5` works with a single set of metadata", {
  # Populate an assay with random values for testing.
  assay <- get_test_assay(
    ncells = 10,
    nfeatures = 10,
    assay_version = "v5"
  )
  # Add similarly random HVF metadata to `assay`.
  assay <- add_hvf_info(
    assay,
    nfeatures = 10,
    method_name = "vst",
    layer_name = "counts"
  )

  # Extract the expected HVFInfo and rename the columns.
  info_columns <- c(
    "vf_vst_counts_variable",
    "vf_vst_counts_rank",
    "vf_vst_counts_value"
  )
  expected_info <- assay[[]][, info_columns]
  colnames(expected_info) <- c("variable", "rank", "value")

  # Check the base case where `method` and `layer` are both set and valid.
  result <- HVFInfo(assay, method = "vst", layer = "counts")
  expect_identical(result, expected_info["value"])

  # Check the same case with all relevant HVF columns returned.
  result <- HVFInfo(assay, method = "vst", layer = "counts", status = TRUE)
  expect_identical(result, expected_info)

  # Check that `layer` can be omitted.
  result <- HVFInfo(assay, method = "vst")
  expect_identical(result, expected_info["value"])

  # Check that `layer` can be `NULL`.
  result <- HVFInfo(assay, method = "vst", layer = NULL)
  expect_identical(result, expected_info["value"])

  # Check that `layer` can be `NA`.
  result <- HVFInfo(assay, method = "vst", layer = NA)
  expect_identical(result, expected_info["value"])

  # Check that the `method` parameter can be omitted.
  result <- HVFInfo(assay)
  expect_identical(result, expected_info["value"])

  # Check that `NULL` is returned if `method` does not point HVF metadata.
  result <- expect_warning(
    HVFInfo(assay, method = "not-a-method"),
    paste(
      "Unable to find highly variable feature information for",
      "method='not-a-method' and layer='NA'."
    )
  )
  expect_null(result)
  result <- expect_warning(
    HVFInfo(assay, method = "not-a-method", layer = "counts"),
    paste(
      "Unable to find highly variable feature information for",
      "method='not-a-method' and layer='counts'."
    )
  )
  expect_null(result)
})

test_that("`HVFInfo.Assay5` works with multiple methods run on different layers", {
  # Populate an assay with random values for testing.
  assay <- get_test_assay(
    ncells = 10,
    nfeatures = 10,
    assay_version = "v5"
  )
  # Add similarly random HVF metadata to `assay`.
  assay <- add_hvf_info(
    assay,
    nfeatures = 5,
    method_name = "vst",
    layer_name = "counts"
  )
  # Add a second "data" layer to the assay by duplicating "count".
  LayerData(assay, layer = "data") <- LayerData(assay, layer = "counts")
  # Add a second set of HVF columns to `assay`.
  assay <- add_hvf_info(
    assay,
    nfeatures = 5,
    method_name = "mvp",
    layer_name = "data"
  )

  # Extract the first set of HVF info and rename the columns.
  vst_columns <- c(
    "vf_vst_counts_variable",
    "vf_vst_counts_rank",
    "vf_vst_counts_value"
  )
  vst_info <- assay[[]][, vst_columns]
  colnames(vst_info) <- c("variable", "rank", "value")
  # Extract the first set of HVF info and rename the columns.
  mvp_columns <- c(
    "vf_mvp_data_variable",
    "vf_mvp_data_rank",
    "vf_mvp_data_value"
  )
  mvp_info <- assay[[]][, mvp_columns]
  colnames(mvp_info) <- c("variable", "rank", "value")

  # Check the base case where `method` and `layer` are both set and valid.
  result <- HVFInfo(assay, method = "vst", layer = "counts")
  expect_identical(result, vst_info["value"])
  result <- HVFInfo(assay, method = "mvp", layer = "data")
  expect_identical(result, mvp_info["value"])

  # Check that `layer` can be omitted. In this case, `layer` will default to
  # `NA` which will be in turn be interpreted as `Layers(assay)` (i.e. all layers).
  result <- HVFInfo(assay, method = "vst")
  expect_identical(result, vst_info["value"])
  result <- HVFInfo(assay, method = "mvp")
  expect_identical(result, mvp_info["value"])

  # Check that `layer` can be `NULL`. In this case, `layer` will be interpreted
  # as `DefaultLayer(assay)`. Thus, we are expected `layer` to resolve to "counts".
  result <- HVFInfo(assay, method = "vst", layer = NULL)
  expect_identical(result, vst_info["value"])
  result <- expect_warning(
    HVFInfo(assay, method = "mvp", layer = NULL),
    paste(
      "Unable to find highly variable feature information for",
      "method='mvp' and layer='NULL'."
    )
  )
  expect_null(result)

  # Check that `layer` can be `NA`. In this case, `layer` will be interpreted
  # as `Layers(assay)` (i.e. all layers).
  result <- HVFInfo(assay, method = "vst", layer = NA)
  expect_identical(result, vst_info["value"])
  result <- HVFInfo(assay, method = "mvp", layer = NA)
  expect_identical(result, mvp_info["value"])
})

test_that("`HVFInfo.Assay5` works with a single method run on multiple layers", {
  # Populate an assay with random values for testing.
  assay <- get_test_assay(
    ncells = 10,
    nfeatures = 10,
    assay_version = "v5"
  )
  # Split the "counts" layer in half across it's columns and drop the original layer.
  LayerData(assay, layer = "counts.1") <- LayerData(assay, layer = "counts")[, 1:5]
  LayerData(assay, layer = "counts.2") <- LayerData(assay, layer = "counts")[, 5:10]
  # Since it's first, "counts.1" would be chosen as the default layer when
  # "counts" is dropped but we'll do it explicitly (also avoids a warning).
  DefaultLayer(assay) <- "counts.1"
  LayerData(assay, layer = "counts") <- NULL

  # Add random HVF metadata for each layer. By adding "counts.2" before "counts.1"
  # we introduce discrepancy between the behavior of `layer = NULL` and `layer = NA`.
  assay <- add_hvf_info(
    assay,
    nfeatures = 5,
    method_name = "vst",
    layer_name = "counts.2"
  )
  assay <- add_hvf_info(
    assay,
    nfeatures = 5,
    method_name = "vst",
    layer_name = "counts.1"
  )

  # Extract the first set of HVF info and rename the columns.
  vst.1_columns <- c(
    "vf_vst_counts.1_variable",
    "vf_vst_counts.1_rank",
    "vf_vst_counts.1_value"
  )
  vst.1_info <- assay[[]][, vst.1_columns]
  colnames(vst.1_info) <- c("variable", "rank", "value")
  # Extract the second set of HVF info and rename the columns.
  vst.2_columns <- c(
    "vf_vst_counts.2_variable",
    "vf_vst_counts.2_rank",
    "vf_vst_counts.2_value"
  )
  vst.2_info <- assay[[]][, vst.2_columns]
  colnames(vst.2_info) <- c("variable", "rank", "value")

  # Check the base case where `method` and `layer` are both set and valid.
  result <- HVFInfo(assay, method = "vst", layer = "counts.1")
  expect_identical(result, vst.1_info["value"])
  result <- HVFInfo(assay, method = "vst", layer = "counts.2")
  expect_identical(result, vst.2_info["value"])

  # Check that `layer` can be omitted. In this case, `layer` will default to
  # `NULL` which will be in turn be interpreted as `Layers(assay)`
  # (i.e. all layers). Thus, we are expected `layer` to resolve to "counts.2".
  result <- HVFInfo(assay, method = "vst")
  expect_identical(result, vst.2_info["value"])

  # Check that `layer` can be `NULL`. In this case, `layer` will be interpreted
  # as `DefaultLayer(assay)`. Thus, we are expected `layer` to resolve to "counts".
  result <- HVFInfo(assay, method = "vst", layer = NULL)
  expect_identical(result, vst.1_info["value"])

  # Check that `layer` can be `NA`. In this case, `layer` will be interpreted
  # as `Layers(assay)` (i.e. all layers) and the first one associated with
  # `method` will be returned.
  result <- HVFInfo(assay, method = "vst", layer = NA)
  expect_identical(result, vst.2_info["value"])
})

test_that("`HVFInfo.Assay5` works with multiple methods run on the same layer", {
  # Populate an assay with random values for testing.
  assay <- get_test_assay(
    ncells = 10,
    nfeatures = 10,
    assay_version = "v5"
  )
  # Add similarly random HVF metadata to `assay`.
  assay <- add_hvf_info(
    assay,
    nfeatures = 5,
    method_name = "vst",
    layer_name = "counts"
  )
  # Add a second set of HVF columns to `assay`.
  assay <- add_hvf_info(
    assay,
    nfeatures = 5,
    method_name = "mvp",
    layer_name = "counts"
  )

  # Extract the first set of HVF info and rename the columns.
  vst_columns <- c(
    "vf_vst_counts_variable",
    "vf_vst_counts_rank",
    "vf_vst_counts_value"
  )
  vst_info <- assay[[]][, vst_columns]
  colnames(vst_info) <- c("variable", "rank", "value")
  # Extract the first set of HVF info and rename the columns.
  mvp_columns <- c(
    "vf_mvp_counts_variable",
    "vf_mvp_counts_rank",
    "vf_mvp_counts_value"
  )
  mvp_info <- assay[[]][, mvp_columns]
  colnames(mvp_info) <- c("variable", "rank", "value")

  # Check the base case where `method` and `layer` are both set and valid.
  result <- HVFInfo(assay, method = "vst", layer = "counts")
  expect_identical(result, vst_info["value"])
  result <- HVFInfo(assay, method = "mvp", layer = "counts")
  expect_identical(result, mvp_info["value"])

  # Check the same case with all relevant HVF columns returned.
  result <- HVFInfo(assay, method = "vst", layer = "counts", status = TRUE)
  expect_identical(result, vst_info)
  result <- HVFInfo(assay, method = "mvp", layer = "counts", status = TRUE)
  expect_identical(result, mvp_info)

  # Check that `layer` can be omitted. In this case, `layer` will default to
  # `NULL` which will be in turn be interpreted as `DefaultLayer(assay)`.
  # Thus, we are expected `layer` to resolve to "counts".
  result <- HVFInfo(assay, method = "vst")
  expect_identical(result, vst_info["value"])
  result <- HVFInfo(assay, method = "mvp")
  expect_identical(result, mvp_info["value"])

  # Check that `layer` can be `NULL`. In this case, `layer` will be interpreted
  # as `DefaultLayer(assay)`. Thus, we are expected `layer` to resolve to "counts".
  result <- HVFInfo(assay, method = "vst", layer = NULL)
  expect_identical(result, vst_info["value"])
  result <- HVFInfo(assay, method = "mvp", layer = NULL)
  expect_identical(result, mvp_info["value"])

  # Check that `layer` can be `NA`. In this case, `layer` will be interpreted
  # as `Layers(assay)` (i.e. all layers).
  result <- HVFInfo(assay, method = "vst", layer = NA)
  expect_identical(result, vst_info["value"])
  result <- HVFInfo(assay, method = "mvp", layer = NA)
  expect_identical(result, mvp_info["value"])
})

context("VariableFeatures")

test_that("`VariableFeatures.Assay5` getter/setter works as expected", {
  # Populate an assay with random values for testing.
  assay <- get_test_assay(
    ncells = 10,
    nfeatures = 10,
    assay_version = "v5"
  )

  # Set initial variable features (including some missing values).
  test_features <- c(gene1 = "gene1", gene3 = "gene3", gene5 = "gene5")
  assay[["var.features"]] <- test_features
  expected_features <- c("gene1", "gene3", "gene5")
  expect_identical(VariableFeatures(assay), expected_features)
  expect_identical(VariableFeatures(assay, nfeatures = 2), expected_features[1:2])

  # Set the ranking for the variable features.
  test_rank <- c(gene1 = 3, gene3 = 1, gene5 = 2)
  assay[["var.features.rank"]] <- test_rank
  expected_features <- c("gene3", "gene5", "gene1")
  expect_identical(VariableFeatures(assay), expected_features)
  expect_identical(VariableFeatures(assay, nfeatures = 2), expected_features[1:2])

  # Use the VariableFeatures setter to overwrite the current variable features.
  expected_features <- c("gene10", "gene2", "gene8")
  VariableFeatures(assay) <- expected_features
  expect_identical(VariableFeatures(assay), expected_features)
})
