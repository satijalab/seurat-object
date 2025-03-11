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
