#' Returns a random counts matrix.
get_random_counts <- function(ncells, nfeatures) {
  # Instantiate an empty nfeatures by ncells matrix.
  counts <- matrix(NA, ncol = ncells, nrow = nfeatures)
  
  # Populate the matrix one row at a time by generating a negative 
  # binomial distribution using a random `mu` between 0.1 and 10.
  for (i in 1:nfeatures) {
    mu <- runif(1, min = 0.1, max = 10)
    counts[i, ] <- values <- rnbinom(ncells, mu = mu, size = 1)
  }

  # Assign column and row names to the matrix to label cells and genes.
  colnames(counts) <- paste0("cell", seq(ncol(counts)))
  row.names(counts) <- paste0("gene", seq(nrow(counts)))

  # Convert `counts` to a `dgCMatrix`.
  counts_sparse <- as.sparse(counts)

  return(counts_sparse)
}

context("merge")

test_that("`merge` works with multi-assay inputs", {
  # Run checks against `Assay` and `Assay5` instances.
  for (assay_version in c("v5")) {
    create_assay <- switch(assay_version,
      v3 = CreateAssayObject,
      v5 = CreateAssay5Object,
      stop("`assay_version` should be one of 'v3', 'v5'")
    )
    
    # Populate a `Seurat` instance with two assays, "LEFT" and "RIGHT".
    counts_LEFT_A <- get_random_counts(ncells = 10, nfeatures = 10)
    input_A <- CreateSeuratObject(
      create_assay(counts_LEFT_A), 
      assay = "LEFT", 
      project = "A"
    )
    counts_RIGHT_A <- get_random_counts(ncells = 10, nfeatures = 10)
    input_A[["RIGHT"]] <- create_assay(counts_RIGHT_A)
    
    # Populate a second `Seurat` instance with a single assay, "LEFT".
    counts_LEFT_B <- get_random_counts(ncells = 5, nfeatures = 5)
    input_B <- CreateSeuratObject(
      create_assay(counts_LEFT_B), 
      assay = "LEFT", 
      project = "B"
    )
    
    # Populate a third `Seurat` instance with two assays, "LEFT" and "RIGHT".
    counts_LEFT_C <- get_random_counts(ncells = 5, nfeatures = 5)
    input_C <- CreateSeuratObject(
      create_assay(counts_LEFT_C), 
      assay = "LEFT", 
      project = "C"
    )
    counts_RIGHT_C <- get_random_counts(ncells = 5, nfeatures = 5)
    input_C[["RIGHT"]] <- create_assay(counts_RIGHT_C)

    # Merge the three `Seurat` instances.
    result <- expect_warning(
      merge(input_A, c(input_B, input_C)),
      paste(
        "Some cell names are duplicated across objects provided.",
        "Renaming to enforce unique cell names."
      )
    )

    # If we have an `Assay` input, split it into multiple layers
    if (assay_version == "v3") {
      result <- expect_warning(split(result, f = Idents(result)))
    }
    
    # Check that the result layers are named according to their project values.
    expected_layers <- c("counts.A", "counts.B", "counts.C")
    expect_identical(Layers(result[["LEFT"]]), expected_layers)
    expected_layers <- c("counts.A", "counts.C")
    expect_identical(Layers(result[["RIGHT"]]), expected_layers)
    # Check that the "LEFT" assay contains the expected cell names.
    expected_cells <- c(
      paste0("cell", seq(10), "_1"),
      paste0("cell", seq(5), "_2"),
      paste0("cell", seq(5), "_3")
    )
    expect_identical(Cells(result[["LEFT"]]), expected_cells)
    # Check that the "RIGHT" assay contains the expected cell names.
    expected_cells <- c(
      paste0("cell", seq(10), "_1"),
      paste0("cell", seq(5), "_3")
    )
    expect_identical(Cells(result[["RIGHT"]]), expected_cells)
    # Check that the values for "counts.A" in the "LEFT" assay are preserved.
    expected_counts <- counts_LEFT_A
    colnames(expected_counts) <- paste0("cell", seq(10), "_1")
    expect_identical(
      LayerData(result, assay = "LEFT", layer = "counts.A"),
      expected_counts
    )
    # Check that the values for "counts.B" in the "LEFT" assay are preserved.
    expected_counts <- counts_LEFT_B
    colnames(expected_counts) <- paste0("cell", seq(5), "_2")
    observed_counts <- LayerData(result, assay = "LEFT", layer = "counts.B")
    expect_identical(
      observed_counts,
      expected_counts,
      info = dim(expected_counts)
    )
    # Check that the values for "counts.C" in the "LEFT" assay are preserved.
    expected_counts <- counts_LEFT_C
    colnames(expected_counts) <- paste0("cell", seq(5), "_3")
    expect_identical(
      LayerData(result, assay = "LEFT", layer = "counts.C"),
      expected_counts,
      info = assay_version
    )
    # Check that the values for "counts.A" in the "RIGHT" assay are preserved.
    expected_counts <- counts_RIGHT_A
    colnames(expected_counts) <- paste0("cell", seq(10), "_1")
    expect_identical(
      LayerData(result, assay = "RIGHT", layer = "counts.A"),
      expected_counts
    )
    # Check that the values for "counts.C" in the "RIGHT" assay are preserved.
    expected_counts <- counts_RIGHT_C
    colnames(expected_counts) <- paste0("cell", seq(5), "_3")
    expect_identical(
      LayerData(result, assay = "RIGHT", layer = "counts.C"),
      expected_counts,
      info = assay_version
    )
  }
})

context("RenameCells")

test_that("`RenameCells` works with multi-assay inputs", {
  counts <- get_random_counts(ncells = 10, nfeatures = 10)

  # Run checks against `Assay` and `Assay5` instances.
  for (assay_version in c("v3", "v5")) {
    create_assay <- switch(assay_version,
      v3 = CreateAssayObject,
      v5 = CreateAssay5Object,
      stop("`assay_version` should be one of 'v3', 'v5'")
    )

    counts_left <- counts[, 1:5]
    assay_left <- create_assay(counts_left)
    left <- CreateSeuratObject(assay_left)

    counts_right <- counts[, 6:10]
    assay_right <- create_assay(counts_right)
    right <- CreateSeuratObject(assay_right, assay = "right")

    test_case <- merge(left, right)
    old_names <- colnames(test_case)

    # Rename all the cells in the input using `add.cell.id`.
    new_names <- paste0("new_cell", 1:10)
    result <- RenameCells(test_case, add.cell.id = "new")
    expect_identical(colnames(result), new_names)

    # Rename all the cells in the input using `new.names`.
    new_names <- paste0("new-cell", 1:10)
    result <- RenameCells(test_case, new.names = new_names)
    expect_identical(colnames(result), new_names)

    # Rename just the cells in the input's default assay.
    new_names <- paste0("new-cell", 1:5)
    expected_names <- c(new_names, tail(old_names, n = 5))
    result <- RenameCells(test_case, new.names = new_names)
    expect_identical(colnames(result), expected_names)

    # Rename a subset of cells across both assays.
    new_names <- paste0("new-cell", 2:9)
    expected_names <- c(head(old_names, n = 1), new_names, tail(old_names, n = 1))
    # Requires a named vector to be passed to `new.names`.
    expect_error(RenameCells(test_case, new.names = new_names))
    names(new_names) <- old_names[2:9]
    result <- RenameCells(test_case, new.names = new_names)
    expect_identical(colnames(result), expected_names)

    # All cells provided in a named vector mapping must exist.
    new_names <- paste0("new-cell", 1:10)
    names(new_names) <- c("not-a-cell", tail(old_names, n = 1))
    expect_error(RenameCells(test_case, new.names = new_names))
  }
})
