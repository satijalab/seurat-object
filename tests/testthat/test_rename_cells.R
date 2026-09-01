set.seed(42)

make_counts <- function(ncell, nfeat = 10) {
  m <- matrix(rpois(nfeat * ncell, lambda = 3), nrow = nfeat, ncol = ncell)
  dimnames(m) <- list(
    paste0("gene", seq_len(nfeat)),
    paste0("old", seq_len(ncell))
  )
  as.sparse(m)
}

test_that("a v5 assay renames from a plain vector of names", {
  # RenameCells() documents new.names as "vector of new cell names"; the v5
  # method treated it as a lookup keyed by the old names, so an unnamed vector
  # produced NA for every cell
  for (n in c(1L, 2L, 5L)) {
    assay <- CreateAssay5Object(counts = make_counts(n))
    renamed <- RenameCells(assay, new.names = paste0("new", seq_len(n)))
    expect_identical(Cells(renamed), paste0("new", seq_len(n)), info = n)
    expect_false(anyNA(Cells(renamed)))
    expect_identical(
      colnames(LayerData(renamed, layer = "counts")),
      paste0("new", seq_len(n))
    )
  }
})

test_that("a v5 assay still accepts a named lookup", {
  assay <- CreateAssay5Object(counts = make_counts(4))
  lookup <- setNames(paste0("new", 1:4), Cells(assay))
  expect_identical(Cells(RenameCells(assay, new.names = lookup)), paste0("new", 1:4))
})

test_that("a v5 assay rejects the wrong number of names", {
  assay <- CreateAssay5Object(counts = make_counts(3))
  expect_error(
    RenameCells(assay, new.names = c("a", "b")),
    "one entry per cell"
  )
})

test_that("a v3 assay with a single cell is renamed", {
  # the layer loop skipped anything with <= 1 column, so a one-cell assay was
  # silently left alone
  assay <- suppressWarnings(CreateAssayObject(counts = make_counts(1)))
  expect_identical(Cells(assay), "old1")
  renamed <- RenameCells(assay, new.names = "new1")
  expect_identical(Cells(renamed), "new1")
  expect_identical(colnames(GetAssayData(renamed, layer = "counts")), "new1")
})

test_that("a v3 assay with several cells is unaffected", {
  assay <- suppressWarnings(CreateAssayObject(counts = make_counts(4)))
  renamed <- RenameCells(assay, new.names = paste0("new", 1:4))
  expect_identical(Cells(renamed), paste0("new", 1:4))
})

test_that("renaming preserves the data", {
  counts <- make_counts(3)
  for (assay in list(CreateAssay5Object(counts = counts),
                     suppressWarnings(CreateAssayObject(counts = counts)))) {
    renamed <- RenameCells(assay, new.names = paste0("new", 1:3))
    expect_equal(
      unname(as.matrix(GetAssayData(renamed, layer = "counts"))),
      unname(as.matrix(counts))
    )
  }
})
