set.seed(42)

# A layer joined from several on-disk matrices records all of their paths in one
# comma-separated entry. SaveSeuratRds() passed that entry to .FileMove() whole,
# which then reported "Can't find path:" for a path that is really two
test_that("an object with joined on-disk layers can be saved and read back", {
  skip_if_not_installed("BPCells")
  counts <- as.sparse(matrix(
    rpois(200 * 120, lambda = 2),
    nrow = 200,
    dimnames = list(paste0("g", seq_len(200)), paste0("c", seq_len(120)))
  ))
  store <- tempfile("bpcells")
  on.exit(unlink(store, recursive = TRUE), add = TRUE)
  suppressWarnings(BPCells::write_matrix_dir(counts, store))

  object <- suppressWarnings(CreateSeuratObject(BPCells::open_matrix_dir(store)))
  object <- split(object, f = rep(c("a", "b"), length.out = ncol(object)))
  object[["RNA"]] <- JoinLayers(object[["RNA"]])
  # the joined layer points at more than one store
  expect_true(grepl(",", .FilePath(LayerData(object, layer = "counts")), fixed = TRUE))

  rds <- file.path(tempfile("saved"), "object.rds")
  dir.create(dirname(rds), recursive = TRUE)
  on.exit(unlink(dirname(rds), recursive = TRUE), add = TRUE)
  expect_no_error(SaveSeuratRds(object, rds))

  restored <- LoadSeuratRds(rds)
  layer <- LayerData(restored, layer = "counts")
  expect_identical(dim(layer), dim(counts))
  expect_equal(
    as.matrix(layer[, colnames(counts)]),
    as.matrix(counts)
  )
})
