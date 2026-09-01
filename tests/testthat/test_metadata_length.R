set.seed(42)

# An unnamed vector of the wrong length was recycled to the number of cells,
# silently, so cells were given values belonging to other cells
build_object <- function(ncell = 20L, nfeat = 30L) {
  counts <- matrix(
    rpois(nfeat * ncell, lambda = 3),
    nrow = nfeat,
    dimnames = list(paste0("g", seq_len(nfeat)), paste0("c", seq_len(ncell)))
  )
  suppressWarnings(CreateSeuratObject(counts = as.sparse(counts)))
}

test_that("a vector of the wrong length is refused", {
  object <- build_object()
  expect_error(object$three <- 1:3, "Cannot add 3 values as meta data for 20 cells")
  expect_error(object$five <- 1:5, "name the values by cell")
  # a factor takes the same path
  expect_error(object$fac <- factor(rep("a", 5)), "Cannot add 5 values")
  expect_error(
    AddMetaData(object, metadata = 1:5, col.name = "five"),
    "Cannot add 5 values"
  )
})

test_that("one value per cell, and a single value, still work", {
  object <- build_object()
  object$all <- seq_len(ncol(object))
  expect_identical(object$all, setNames(seq_len(ncol(object)), colnames(object)))

  object$one <- 1
  expect_true(all(object$one == 1))
  expect_length(object$one, ncol(object))

  object$flag <- factor(rep("a", ncol(object)))
  expect_true(all(object$flag == "a"))
})

test_that("named values are still matched by name", {
  object <- build_object()
  object$partial <- setNames(1:5, colnames(object)[1:5])
  expect_identical(unname(object$partial[1:5]), 1:5)
  expect_true(all(is.na(object$partial[6:ncol(object)])))
})
