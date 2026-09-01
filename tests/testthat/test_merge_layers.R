set.seed(42)

make_layer <- function(genes, cells, seed) {
  set.seed(seed)
  m <- Matrix::rsparsematrix(length(genes), length(cells), density = 0.4)
  m <- abs(m) * 10
  m@x <- round(m@x) + 1
  dimnames(m) <- list(genes, cells)
  as(m, "dgCMatrix")
}

test_that("layers with different features merge on the union of features", {
  a <- make_layer(paste0("gene", 1:20), paste0("a", 1:15), seed = 1)
  b <- make_layer(paste0("gene", 11:30), paste0("b", 1:12), seed = 2)
  merged <- RowMergeSparseMatrices(mat1 = a, mat2 = b)

  expect_s4_class(merged, "dgCMatrix")
  expect_identical(rownames(merged), union(rownames(a), rownames(b)))
  expect_identical(colnames(merged), c(colnames(a), colnames(b)))
  # values carried over unchanged, and features absent from a layer are zero
  expect_equal(
    as.matrix(merged[rownames(a), colnames(a)]),
    as.matrix(a)
  )
  expect_equal(
    as.matrix(merged[rownames(b), colnames(b)]),
    as.matrix(b)
  )
  only.in.b <- setdiff(rownames(b), rownames(a))
  expect_true(all(as.matrix(merged[only.in.b, colnames(a)]) == 0))
  # total non-zeros are preserved
  expect_equal(length(merged@x), length(a@x) + length(b@x))
})

test_that("more than two layers merge correctly", {
  a <- make_layer(paste0("gene", 1:20), paste0("a", 1:10), seed = 3)
  b <- make_layer(paste0("gene", 11:30), paste0("b", 1:10), seed = 4)
  d <- make_layer(paste0("gene", 25:40), paste0("d", 1:10), seed = 5)
  merged <- RowMergeSparseMatrices(mat1 = a, mat2 = list(b, d))

  expect_identical(
    rownames(merged),
    Reduce(union, list(rownames(a), rownames(b), rownames(d)))
  )
  expect_identical(colnames(merged), c(colnames(a), colnames(b), colnames(d)))
  expect_equal(length(merged@x), sum(length(a@x), length(b@x), length(d@x)))
  for (layer in list(a, b, d)) {
    expect_equal(
      as.matrix(merged[rownames(layer), colnames(layer)]),
      as.matrix(layer)
    )
  }
})

test_that("layers sharing all features take the cbind path and still agree", {
  genes <- paste0("gene", 1:20)
  a <- make_layer(genes, paste0("a", 1:10), seed = 6)
  b <- make_layer(genes, paste0("b", 1:10), seed = 7)
  merged <- RowMergeSparseMatrices(mat1 = a, mat2 = b)
  expect_identical(rownames(merged), genes)
  expect_equal(as.matrix(merged), as.matrix(cbind(a, b)))
})
