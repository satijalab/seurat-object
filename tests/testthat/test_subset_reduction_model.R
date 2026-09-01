set.seed(42)

# A UMAP built with return.model = TRUE keeps the embedding it was fit on in
# misc$model. Subsetting left that untouched, so the model and the reduction
# described different sets of cells, and a query projected onto the reduction
# was placed against the wrong coordinates
build_reduction <- function(ncell = 40L, ndim = 3L) {
  cells <- paste0("c", seq_len(ncell))
  embeddings <- matrix(
    rnorm(ncell * ndim),
    nrow = ncell,
    dimnames = list(cells, paste0("umap_", seq_len(ndim)))
  )
  reduction <- CreateDimReducObject(embeddings = embeddings, key = "umap_", assay = "RNA")
  # a model of the shape uwot stores
  Misc(reduction, slot = "model") <- list(
    embedding = embeddings,
    n_neighbors = 5L,
    metric = list(euclidean = list())
  )
  reduction
}

test_that("subsetting a reduction subsets the model it carries", {
  reduction <- build_reduction()
  kept <- Cells(reduction)[seq(from = 11, to = 40)]
  subset.reduction <- subset(reduction, cells = kept)

  expect_identical(Cells(subset.reduction), kept)
  expect_identical(rownames(Misc(subset.reduction, slot = "model")$embedding), kept)
  # and the coordinates themselves are the ones the reduction has
  expect_identical(
    Misc(subset.reduction, slot = "model")$embedding,
    Embeddings(subset.reduction)
  )
})

test_that("the rest of the model is left alone", {
  reduction <- build_reduction()
  subset.reduction <- subset(reduction, cells = Cells(reduction)[1:20])
  model <- Misc(subset.reduction, slot = "model")
  expect_identical(model$n_neighbors, 5L)
  expect_identical(model$metric, list(euclidean = list()))
})

test_that("a reduction with no model is unchanged", {
  reduction <- build_reduction()
  Misc(reduction, slot = "model") <- NULL
  subset.reduction <- subset(reduction, cells = Cells(reduction)[1:20])
  expect_null(Misc(subset.reduction, slot = "model"))
  expect_identical(nrow(Embeddings(subset.reduction)), 20L)
})
