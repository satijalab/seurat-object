set.seed(42)

make_object <- function(ncell = 10, nfeat = 6) {
  m <- matrix(rpois(nfeat * ncell, lambda = 2), nrow = nfeat, ncol = ncell)
  dimnames(m) <- list(
    paste0("gene", seq_len(nfeat)),
    paste0("c", seq_len(ncell))
  )
  suppressWarnings(CreateSeuratObject(counts = as.sparse(m)))
}

test_that("NA identities do not leak into the other groups", {
  object <- make_object()
  Idents(object) <- c(rep("A", 4), rep("B", 4), NA, NA)
  # comparing an NA identity with a level gives NA, which indexed an NA into
  # every group and inflated their lengths
  groups <- suppressWarnings(CellsByIdentities(object))
  expect_setequal(names(groups), c("A", "B", "NA"))
  expect_identical(groups[["A"]], paste0("c", 1:4))
  expect_identical(groups[["B"]], paste0("c", 5:8))
  expect_identical(groups[["NA"]], paste0("c", 9:10))
  expect_false(any(vapply(groups, anyNA, logical(1))))
})

test_that("identities without NA are unaffected", {
  object <- make_object()
  Idents(object) <- c(rep("A", 5), rep("B", 5))
  groups <- suppressWarnings(CellsByIdentities(object))
  expect_setequal(names(groups), c("A", "B"))
  expect_identical(groups[["A"]], paste0("c", 1:5))
  expect_identical(groups[["B"]], paste0("c", 6:10))
})

test_that("every cell is accounted for exactly once", {
  object <- make_object()
  Idents(object) <- c("A", NA, "B", "A", NA, "B", "A", "B", NA, "A")
  groups <- suppressWarnings(CellsByIdentities(object))
  expect_setequal(unlist(groups, use.names = FALSE), colnames(object))
  expect_equal(sum(lengths(groups)), ncol(object))
})

test_that("restricting to a subset of cells still excludes NA", {
  object <- make_object()
  Idents(object) <- c(rep("A", 4), rep("B", 4), NA, NA)
  groups <- suppressWarnings(
    CellsByIdentities(object, cells = paste0("c", c(1:2, 9:10)))
  )
  expect_identical(groups[["A"]], paste0("c", 1:2))
  expect_false(anyNA(groups[["A"]]))
  expect_identical(groups[["NA"]], paste0("c", 9:10))
})
