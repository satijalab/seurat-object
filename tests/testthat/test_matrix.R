# Tests for matrix checks

data("pbmc_small")

pbmc_small <- suppressWarnings(suppressMessages(UpdateSeuratObject(pbmc_small)))

counts <- GetAssayData(pbmc_small, assay = "RNA", layer = "counts")

test_that("CheckMatrix works for valid matrix", {
  expect_null(CheckMatrix(counts))
})

test_that("CheckMatrix warns for non-integer counts", {
  counts[1,1] <- 1/2
  expect_warning(CheckMatrix(counts))
})

test_that("CheckMatrix warns for NA values", {
  counts[1,1] <- NA
  expect_warning(CheckMatrix(counts))
})

test_that("CheckMatrix warns for NaN values", {
  counts[1,1] <- NaN
  expect_warning(CheckMatrix(counts))
})

test_that("CheckMatrix warns for logical values", {
  counts <- as(object = counts, Class = "lsparseMatrix")
  expect_warning(CheckMatrix(counts))
})
