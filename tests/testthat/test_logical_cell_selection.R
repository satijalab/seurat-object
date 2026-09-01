set.seed(11)
counts <- as.sparse(matrix(
  rpois(20 * 12, lambda = 3),
  nrow = 20,
  dimnames = list(paste0("g", seq_len(20)), paste0("c", seq_len(12)))
))
object <- suppressWarnings(CreateSeuratObject(counts))
mask <- c(rep(TRUE, 6), rep(FALSE, 6))

# `object[, mask]` has always worked. Everything taking a `cells` argument
# resolved positions to names but passed a mask through untouched, so the mask
# was compared against cell names, matched nothing, and the selection silently
# became every cell or none of them
test_that("cells can be selected with a logical mask", {
  expect_identical(WhichCells(object, cells = mask), colnames(object)[mask])
  expect_identical(colnames(subset(object, cells = mask)), colnames(object)[mask])
  expect_identical(colnames(subset(object[["RNA"]], cells = mask)), colnames(object)[mask])
  # the same selection through each of the three accepted forms
  expect_identical(
    colnames(subset(object, cells = mask)),
    colnames(subset(object, cells = which(mask)))
  )
  expect_identical(
    colnames(subset(object, cells = mask)),
    colnames(subset(object, cells = colnames(object)[mask]))
  )
})

test_that("features can be selected with a logical mask", {
  fmask <- c(rep(TRUE, 5), rep(FALSE, 15))
  expect_identical(rownames(subset(object[["RNA"]], features = fmask)), rownames(object)[fmask])
})

test_that("a mask of the wrong length is refused", {
  expect_error(subset(object, cells = c(TRUE, FALSE)), "one value per cell")
  expect_error(subset(object[["RNA"]], cells = c(TRUE, FALSE)), "one value per cell")
})

test_that("dimensional reductions take a mask too", {
  object[["pca"]] <- CreateDimReducObject(
    embeddings = matrix(
      rnorm(ncol(object) * 3), ncol = 3,
      dimnames = list(colnames(object), paste0("pc_", 1:3))
    ),
    key = "pc_",
    assay = "RNA"
  )
  expect_identical(Cells(subset(object[["pca"]], cells = mask)), colnames(object)[mask])
})
