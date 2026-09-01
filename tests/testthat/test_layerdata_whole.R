set.seed(42)

build <- function(nfeat = 40, ncell = 30) {
  counts <- matrix(rpois(nfeat * ncell, lambda = 2), nrow = nfeat, ncol = ncell)
  dimnames(counts) <- list(
    paste0("g", seq_len(nfeat)),
    paste0("c", seq_len(ncell))
  )
  scaled <- matrix(rnorm(nfeat * ncell), nrow = nfeat, dimnames = dimnames(counts))
  assay <- CreateAssay5Object(counts = as.sparse(counts))
  suppressWarnings(LayerData(assay, layer = "scale.data") <- scaled)
  list(assay = assay, counts = counts, scaled = scaled)
}

test_that("a whole layer comes back unchanged", {
  # fetching a layer whole used to subset it with every row and column index,
  # copying the matrix for nothing
  b <- build()
  expect_equal(as.matrix(LayerData(b$assay, layer = "counts")), b$counts)
  expect_equal(as.matrix(LayerData(b$assay, layer = "scale.data")), b$scaled)
})

test_that("dimnames and class are preserved", {
  b <- build()
  got <- LayerData(b$assay, layer = "counts")
  expect_identical(dimnames(got), dimnames(b$counts))
  expect_s4_class(got, "dgCMatrix")
  expect_true(is.matrix(LayerData(b$assay, layer = "scale.data")))
})

test_that("feature and cell subsets are unaffected", {
  b <- build()
  expect_equal(
    as.matrix(LayerData(b$assay, layer = "counts", features = paste0("g", 5:10))),
    b$counts[paste0("g", 5:10), ]
  )
  expect_equal(
    as.matrix(LayerData(b$assay, layer = "counts", cells = paste0("c", 3:8))),
    b$counts[, paste0("c", 3:8)]
  )
  expect_equal(
    as.matrix(LayerData(
      b$assay, layer = "counts",
      features = paste0("g", 1:4), cells = paste0("c", 1:4)
    )),
    b$counts[1:4, 1:4]
  )
})

test_that("split layers come back whole and correct", {
  b <- build()
  assay <- split(b$assay, f = rep(c("x", "y"), length.out = ncol(b$assay)))
  for (layer in Layers(assay)) {
    got <- LayerData(assay, layer = layer)
    src <- if (startsWith(layer, "scale.data")) b$scaled else b$counts
    expect_equal(
      as.matrix(got),
      src[rownames(got), colnames(got)],
      info = layer
    )
  }
})

test_that("the returned layer is independent of the assay", {
  # the whole-layer path returns the stored matrix rather than a copy, so check
  # that modifying the result does not reach back into the assay
  b <- build()
  got <- LayerData(b$assay, layer = "scale.data")
  got[1, 1] <- 12345
  expect_equal(LayerData(b$assay, layer = "scale.data")[1, 1], b$scaled[1, 1])
})
