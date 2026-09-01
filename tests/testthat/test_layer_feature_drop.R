set.seed(42)

make_assay <- function(nfeat = 4, ncell = 3) {
  m <- matrix(1, nrow = nfeat, ncol = ncell)
  dimnames(m) <- list(paste0("g", seq_len(nfeat)), paste0("c", seq_len(ncell)))
  CreateAssay5Object(counts = as.sparse(m))
}

replace_counts <- function(assay, value) {
  warnings <- character()
  result <- withCallingHandlers(
    {
      x <- assay
      LayerData(x, layer = "counts") <- value
      x
    },
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(assay = result, warnings = warnings)
}

test_that("replacing a layer names the features whose data is discarded", {
  assay <- make_assay()
  renamed <- as.matrix(LayerData(assay, layer = "counts"))
  rownames(renamed)[1] <- "RENAMED"
  # the old warning said only that the features differed, so the fact that a
  # feature's data had been dropped had to be inferred from a row count
  out <- replace_counts(assay, as.sparse(renamed))
  expect_length(out$warnings, 1L)
  expect_match(out$warnings, "RENAMED")
  expect_match(out$warnings, "discarded")
})

test_that("several dropped features are counted and truncated", {
  assay <- make_assay(nfeat = 12)
  renamed <- as.matrix(LayerData(assay, layer = "counts"))
  rownames(renamed)[1:8] <- paste0("NEW", 1:8)
  out <- replace_counts(assay, as.sparse(renamed))
  expect_match(out$warnings, "8 features are not in the assay")
  expect_match(out$warnings, "3 more")
})

test_that("a single dropped feature reads in the singular", {
  assay <- make_assay()
  renamed <- as.matrix(LayerData(assay, layer = "counts"))
  rownames(renamed)[1] <- "RENAMED"
  out <- replace_counts(assay, as.sparse(renamed))
  expect_match(out$warnings, "1 feature is not in the assay")
})

test_that("replacing with the same features warns about nothing", {
  assay <- make_assay()
  same <- LayerData(assay, layer = "counts")
  out <- replace_counts(assay, same)
  expect_length(out$warnings, 0L)
  expect_equal(
    as.matrix(LayerData(out$assay, layer = "counts")),
    as.matrix(same)
  )
})

test_that("replacing with a subset of features still warns, without claiming a discard", {
  assay <- make_assay()
  subset <- LayerData(assay, layer = "counts")[1:2, , drop = FALSE]
  out <- replace_counts(assay, subset)
  expect_length(out$warnings, 1L)
  expect_match(out$warnings, "Different features")
  expect_false(any(grepl("discarded", out$warnings)))
})
