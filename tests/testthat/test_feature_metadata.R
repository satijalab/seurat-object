set.seed(42)

# An unnamed vector was named by its own values and then intersected with the
# feature names, so anything that was not a list of features was rejected with
# "No feature overlap between new meta data and assay"
build_assay <- function(nfeat = 20L, ncell = 10L) {
  counts <- matrix(
    rpois(nfeat * ncell, lambda = 3),
    nrow = nfeat,
    dimnames = list(paste0("g", seq_len(nfeat)), paste0("c", seq_len(ncell)))
  )
  CreateAssay5Object(counts = as.sparse(counts))
}

test_that("an unnamed vector is taken in feature order", {
  assay <- build_assay()
  assay <- AddMetaData(assay, metadata = seq_len(nrow(assay)), col.name = "score")
  expect_identical(assay[[]]$score, seq_len(nrow(assay)))
  expect_identical(rownames(assay[[]]), Features(assay))

  assay[["type"]] <- rep(c("a", "b"), each = 10)
  expect_identical(assay[[]]$type, rep(c("a", "b"), each = 10))
})

test_that("a named vector and a data frame line up by name", {
  assay <- build_assay()
  scores <- setNames(seq_len(nrow(assay)), rev(Features(assay)))
  assay[["score"]] <- scores
  expect_identical(assay[[]][rev(Features(assay)), "score"], seq_len(nrow(assay)))

  frame <- data.frame(other = seq_len(nrow(assay)), row.names = rev(Features(assay)))
  assay <- AddMetaData(assay, metadata = frame)
  expect_identical(assay[[]][rev(Features(assay)), "other"], seq_len(nrow(assay)))
})

test_that("variable features still round-trip", {
  assay <- build_assay()
  VariableFeatures(assay) <- c("g5", "g2", "g9")
  expect_identical(VariableFeatures(assay), c("g5", "g2", "g9"))

  # every feature variable, which is as long as the assay
  VariableFeatures(assay) <- rev(Features(assay))
  expect_identical(VariableFeatures(assay), rev(Features(assay)))
})

test_that("a scalar is recycled and a wrong length is explained", {
  assay <- build_assay()
  assay[["flag"]] <- TRUE
  expect_true(all(assay[[]]$flag))
  expect_equal(nrow(assay[[]]), nrow(assay))

  expect_error(assay[["short"]] <- 1:5, "length 5 with no names")
})
