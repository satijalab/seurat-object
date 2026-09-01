set.seed(42)

make_counts <- function(nfeat = 20, ncell = 10) {
  m <- matrix(rpois(nfeat * ncell, lambda = 2), nrow = nfeat, ncol = ncell)
  dimnames(m) <- list(paste0("g", seq_len(nfeat)), paste0("c", seq_len(ncell)))
  as.sparse(m)
}

test_that("duplicate features are reported as features", {
  counts <- make_counts()
  rownames(counts) <- c(paste0("g", 1:18), "g1", "g2")
  # previously: invalid class "LogMap" object: Duplicate rownames not allowed
  expect_error(CreateSeuratObject(counts = counts), "Duplicate features")
  expect_error(CreateSeuratObject(counts = counts), "g1")
  expect_error(CreateSeuratObject(counts = counts), "g2")
  expect_error(CreateSeuratObject(counts = counts), "make.unique")
})

test_that("duplicate cell names are reported as cell names", {
  counts <- make_counts()
  colnames(counts) <- c(paste0("c", 1:9), "c1")
  expect_error(CreateSeuratObject(counts = counts), "Duplicate cell names")
  expect_error(CreateSeuratObject(counts = counts), "c1")
})

test_that("a single duplicate reads in the singular", {
  counts <- make_counts()
  rownames(counts) <- c(paste0("g", 1:19), "g1")
  expect_error(CreateSeuratObject(counts = counts), "1 is repeated")
})

test_that("many duplicates are summarised", {
  counts <- make_counts(nfeat = 20)
  rownames(counts) <- c(paste0("g", 1:12), paste0("g", 1:8))
  message <- tryCatch(CreateSeuratObject(counts = counts), error = conditionMessage)
  expect_match(message, "8 are repeated")
  expect_match(message, "3 more")
})

test_that("make.unique resolves it", {
  counts <- make_counts()
  rownames(counts) <- make.unique(c(paste0("g", 1:18), "g1", "g2"))
  object <- suppressWarnings(CreateSeuratObject(counts = counts))
  expect_equal(nrow(object), 20L)
  expect_false(anyDuplicated(rownames(object)) > 0)
})

test_that("unique names are unaffected", {
  counts <- make_counts()
  object <- suppressWarnings(CreateSeuratObject(counts = counts))
  expect_equal(dim(object), c(20L, 10L))
  expect_identical(rownames(object), paste0("g", 1:20))
})

test_that("LogMap itself reports duplicates", {
  expect_error(LogMap(c("a", "b", "a")), "Duplicate names")
  expect_error(LogMap(c("a", "b", "a"), what = "features"), "Duplicate features")
  expect_s4_class(LogMap(c("a", "b")), "LogMap")
})
