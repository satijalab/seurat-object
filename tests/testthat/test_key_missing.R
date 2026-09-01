set.seed(42)

# Key() vapplies over the keyed members with FUN.VALUE = character(1), so a
# member carrying no key at all returned length 0 and the whole call failed
# with "values must be length 1, but FUN(X[[2]]) result is length 0", which
# took down FetchData and everything that plots
build_object <- function() {
  counts <- matrix(
    rpois(200, lambda = 3),
    nrow = 20,
    dimnames = list(paste0("g", seq_len(20)), paste0("c", seq_len(10)))
  )
  object <- suppressWarnings(CreateSeuratObject(counts = as.sparse(counts)))
  object[["other"]] <- CreateAssayObject(counts = counts)
  object
}

drop_key <- function(object, assay) {
  slot(slot(object, name = "assays")[[assay]], name = "key") <- character(0)
  object
}

test_that("an assay with no key does not break Key()", {
  object <- drop_key(build_object(), "other")
  keys <- Key(object)
  expect_length(keys, 3L)
  expect_identical(unname(keys[["other"]]), NA_character_)
  expect_identical(unname(keys[["RNA"]]), "rna_")
  expect_identical(Keys(object), keys)
})

test_that("an assay with no key does not break fetching data", {
  object <- drop_key(build_object(), "other")
  fetched <- FetchData(object, vars = "g1")
  expect_identical(rownames(fetched), colnames(object))
  expect_identical(ncol(fetched), 1L)
})

test_that("keys are unchanged when every object has one", {
  object <- build_object()
  keys <- Key(object)
  expect_false(anyNA(keys))
  expect_identical(unname(keys[["other"]]), "other_")
})
