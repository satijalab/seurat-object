set.seed(42)

# Renaming to names that are not unique left the frames with duplicate row
# names, and the error, "duplicate 'row.names' are not allowed", named neither
# the cells nor the renaming
build_object <- function(ncell = 20L, nfeat = 30L) {
  counts <- matrix(
    rpois(nfeat * ncell, lambda = 3),
    nrow = nfeat,
    dimnames = list(paste0("g", seq_len(nfeat)), paste0("c", seq_len(ncell)))
  )
  suppressWarnings(CreateSeuratObject(counts = as.sparse(counts)))
}

test_that("renaming every cell to the same name is refused", {
  object <- build_object()
  expect_error(
    RenameCells(object, new.names = rep("dup", ncol(object))),
    "Renaming would give the same name to more than one cell"
  )
  expect_error(
    RenameCells(object, new.names = rep("dup", ncol(object))),
    "'dup'"
  )
})

test_that("a single collision is refused and named", {
  object <- build_object()
  names <- colnames(object)
  names[1] <- names[2]
  expect_error(
    RenameCells(object, new.names = names),
    sprintf("1 name is repeated: '%s'", names[2])
  )
})

test_that("unique names still rename", {
  object <- build_object()
  renamed <- RenameCells(object, new.names = paste0(colnames(object), "_x"))
  expect_identical(colnames(renamed), paste0(colnames(object), "_x"))
  expect_identical(rownames(renamed[[]]), colnames(renamed))

  prefixed <- RenameCells(object, add.cell.id = "sample")
  expect_identical(colnames(prefixed), paste0("sample_", colnames(object)))
})
