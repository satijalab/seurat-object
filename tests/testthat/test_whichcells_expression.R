set.seed(42)

# An expression naming nothing the object has left FetchData() with an empty
# request, and it reported "None of the requested variables were found: " with
# nothing after the colon
data("pbmc_small", package = "SeuratObject", envir = environment())

test_that("a name the object does not have is reported", {
  expect_error(
    WhichCells(pbmc_small, expression = nope > 1),
    "None of the names in the expression are in the object: 'nope'"
  )
  expect_error(
    WhichCells(pbmc_small, expression = nope > 1),
    "features, meta data columns, or start with the key"
  )
})

test_that("expressions that do name something still work", {
  cells <- WhichCells(pbmc_small, expression = nCount_RNA > 100)
  expect_true(length(cells) > 0)
  expect_true(all(cells %in% colnames(pbmc_small)))
  expect_true(all(pbmc_small$nCount_RNA[cells] > 100))

  # a value from the calling environment alongside a column
  threshold <- 100
  expect_identical(WhichCells(pbmc_small, expression = nCount_RNA > threshold), cells)

  # a feature
  expect_true(length(WhichCells(pbmc_small, expression = MS4A1 > 1)) > 0)
})
