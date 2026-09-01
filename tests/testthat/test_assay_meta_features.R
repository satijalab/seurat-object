set.seed(42)

make_assay5 <- function(nfeat = 40, ncell = 25) {
  counts <- matrix(rpois(nfeat * ncell, lambda = 2), nrow = nfeat, ncol = ncell)
  dimnames(counts) <- list(
    paste0("gene", seq_len(nfeat)),
    paste0("cell", seq_len(ncell))
  )
  assay <- CreateAssay5Object(counts = as.sparse(counts))
  Key(assay) <- "rna_"
  list(assay = assay, counts = counts)
}

test_that("assigning an empty feature meta.data frame is a no-op", {
  built <- make_assay5()
  assay <- as(object = built$assay, Class = "Assay")
  expect_equal(ncol(assay[[]]), 0L)
  # A data frame with no columns is what an assay that has never had feature
  # metadata computed hands over; 1:ncol() walked indices 1 and 0 for it
  expect_no_error(assay[[]] <- data.frame(row.names = rownames(assay)))
  expect_equal(ncol(assay[[]]), 0L)
  expect_identical(rownames(assay[[]]), rownames(assay))
})

test_that("a v5 assay with no feature metadata coerces to a v3 assay", {
  built <- make_assay5()
  expect_equal(ncol(built$assay[[]]), 0L)
  converted <- expect_no_error(suppressWarnings(as(object = built$assay, Class = "Assay")))
  expect_s4_class(converted, "Assay")
  expect_equal(
    as.matrix(GetAssayData(converted, layer = "counts")),
    built$counts
  )
  expect_identical(rownames(converted), rownames(built$assay))
  expect_identical(colnames(converted), colnames(built$assay))
})

test_that("feature metadata is carried across when it is present", {
  built <- make_assay5()
  assay <- built$assay
  score <- setNames(seq_len(nrow(assay)) / 10, rownames(assay))
  flag <- setNames(rep(c(TRUE, FALSE), length.out = nrow(assay)), rownames(assay))
  assay[["score"]] <- score
  assay[["flag"]] <- flag
  converted <- suppressWarnings(as(object = assay, Class = "Assay"))
  expect_true(all(c("score", "flag") %in% colnames(converted[[]])))
  expect_equal(converted[["score"]][, 1], unname(score))
  expect_equal(converted[["flag"]][, 1], unname(flag))
})

test_that("assigning a non-empty frame still replaces the metadata", {
  built <- make_assay5()
  assay <- suppressWarnings(as(object = built$assay, Class = "Assay"))
  df <- data.frame(
    score = seq_len(nrow(assay)),
    row.names = rownames(assay)
  )
  assay[[]] <- df
  expect_identical(colnames(assay[[]]), "score")
  expect_equal(assay[["score"]][, 1], df$score)
})
