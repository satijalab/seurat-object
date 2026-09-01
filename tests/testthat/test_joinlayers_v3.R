# A v3 assay holds one counts, data and scale.data matrix rather than layers, so
# there is nothing to join. JoinLayers() had no method for it, and merging SCT
# objects leaves an SCT assay that is still v3, so JoinLayers() on the result
# came back with "no applicable method ... class SCTAssay" instead
test_that("joining a v3 assay returns it unchanged", {
  object <- pbmc_small
  expect_identical(JoinLayers(object[["RNA"]]), object[["RNA"]])
  expect_identical(Layers(JoinLayers(object)[["RNA"]]), Layers(object[["RNA"]]))
})

test_that("joining an object whose default assay is v3 works", {
  object <- pbmc_small
  expect_no_error(JoinLayers(object))
  expect_identical(LayerData(JoinLayers(object), layer = "counts"),
                   LayerData(object, layer = "counts"))
})

test_that("a v5 assay in the same object still joins", {
  object <- pbmc_small
  object[["v5"]] <- CreateAssay5Object(counts = LayerData(object, layer = "counts"))
  object[["v5"]] <- split(object[["v5"]], f = rep(c("a", "b"), length.out = ncol(object)))
  expect_length(Layers(object[["v5"]]), 2L)
  joined <- JoinLayers(object, assay = "v5")
  expect_identical(Layers(joined[["v5"]]), "counts")
  # and the v3 assay alongside it is untouched
  expect_identical(Layers(joined[["RNA"]]), Layers(object[["RNA"]]))
})
