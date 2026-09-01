set.seed(42)

make_object <- function(tag, nfeat = 6, ncell = 4) {
  m <- matrix(rpois(nfeat * ncell, lambda = 3), nrow = nfeat, ncol = ncell)
  dimnames(m) <- list(paste0("g", seq_len(nfeat)), paste0(tag, seq_len(ncell)))
  suppressWarnings(CreateSeuratObject(counts = as.sparse(m)))
}

test_that("collapse = TRUE reports that it is not implemented", {
  # the documentation described collapse as if it worked; this pins the two
  # together so they cannot drift apart again
  expect_error(
    merge(make_object("a"), make_object("b"), collapse = TRUE),
    "not yet supported"
  )
})

test_that("collapse = FALSE keeps the layers separate", {
  merged <- merge(make_object("a"), make_object("b"), collapse = FALSE)
  expect_setequal(Layers(merged), c("counts.1", "counts.2"))
  expect_equal(ncol(merged), 8L)
})

test_that("the default matches collapse = FALSE", {
  expect_setequal(
    Layers(merge(make_object("a"), make_object("b"))),
    Layers(merge(make_object("a"), make_object("b"), collapse = FALSE))
  )
})
