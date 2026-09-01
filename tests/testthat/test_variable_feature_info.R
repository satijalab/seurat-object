set.seed(42)

# HVFInfo(status = TRUE) and SVFInfo(status = TRUE) filled their status columns
# with object[[name]], which returns a one-column data frame, so the result had
# data frames nested inside it and anything sorting on them failed with
# "cannot xtfrm data frames"
build_assay <- function(nfeat = 20L, ncell = 10L) {
  counts <- matrix(
    rpois(nfeat * ncell, lambda = 3),
    nrow = nfeat,
    dimnames = list(paste0("g", seq_len(nfeat)), paste0("c", seq_len(ncell)))
  )
  CreateAssayObject(counts = as.sparse(counts))
}

test_that("HVFInfo status columns are plain vectors", {
  assay <- build_assay()
  features <- rownames(assay)
  assay[["vst.mean"]] <- seq_along(features)
  assay[["vst.variance"]] <- rev(seq_along(features))
  assay[["vst.variance.standardized"]] <- seq_along(features) / 2
  assay[["vst.variable"]] <- features %in% features[1:5]

  info <- HVFInfo(assay, method = "vst", status = TRUE)
  expect_true(is.logical(info$variable))
  expect_false(is.data.frame(info$variable))
  expect_identical(sum(info$variable), 5L)
  expect_no_error(info[order(info$variance), ])
})

test_that("SVFInfo status columns are plain vectors and can be sorted", {
  assay <- build_assay()
  features <- rownames(assay)
  assay[["moransi.spatially.variable"]] <- features %in% features[1:4]
  assay[["moransi.spatially.variable.rank"]] <- seq_along(features)
  assay[["MoransI_observed"]] <- runif(length(features))
  assay[["MoransI_p.value"]] <- runif(length(features))

  info <- SVFInfo(assay, method = "moransi", status = TRUE)
  expect_true(is.logical(info$variable))
  expect_true(is.numeric(info$rank))
  expect_no_error(info[order(info$rank), ])
  expect_identical(SpatiallyVariableFeatures(assay, method = "moransi"), features[1:4])
})

test_that("variable features are still selected from the status column", {
  assay <- build_assay()
  features <- rownames(assay)
  assay[["vst.mean"]] <- seq_along(features)
  assay[["vst.variance"]] <- rev(seq_along(features))
  assay[["vst.variance.standardized"]] <- seq_along(features) / 2
  assay[["vst.variable"]] <- features %in% features[3:6]
  expect_identical(VariableFeatures(assay, method = "vst"), features[3:6])
})

# Asking for information that was never computed reached the meta data with a
# column name it does not have: "undefined columns selected"
test_that("information that was never computed is reported as such", {
  assay <- build_assay()
  expect_error(
    SVFInfo(assay, method = "markvariogram", status = TRUE),
    "Unable to find spatially variable feature information"
  )
  expect_error(
    SVFInfo(assay, method = "markvariogram", status = TRUE),
    "Run FindSpatiallyVariableFeatures"
  )
  expect_error(
    SpatiallyVariableFeatures(assay, method = "moransi"),
    "Unable to find spatially variable feature information"
  )
})

test_that("half-written information is reported as such", {
  assay <- build_assay()
  features <- rownames(assay)
  # the ranks stored without the status column, as an interrupted run leaves it
  assay[["moransi.spatially.variable.rank"]] <- seq_along(features)
  assay[["MoransI_observed"]] <- runif(length(features))
  expect_error(
    SVFInfo(assay, method = "moransi", status = TRUE),
    "no 'moransi.spatially.variable' column"
  )
})
