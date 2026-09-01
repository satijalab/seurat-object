# `Idents(obj) <- "group"` where no such meta data column exists recycles the
# string over every cell. That looks like a successful grouping until the
# levels are inspected, and it destroys the existing identities.

test_that("a string that is not a meta data column warns", {
  object <- pbmc_small
  expect_warning(Idents(object) <- "group", "not a column")
  expect_identical(levels(object), "group")
  expect_true(all(as.character(Idents(object)) == "group"))
})

test_that("the warning names the value and the way to silence it", {
  object <- pbmc_small
  expect_warning(Idents(object) <- "group", "group")
  expect_warning(Idents(object) <- "group", "I\\(")
})

test_that("an existing meta data column is used, silently", {
  object <- pbmc_small
  expect_no_warning(Idents(object) <- "letter.idents")
  expect_setequal(levels(object), c("A", "B"))
})

test_that("I() sets a literal identity without warning", {
  object <- pbmc_small
  expect_no_warning(Idents(object) <- I("group"))
  expect_identical(levels(object), "group")
})

test_that("I() wins over a column of the same name", {
  object <- pbmc_small
  expect_true("letter.idents" %in% colnames(object[[]]))
  expect_no_warning(Idents(object) <- I("letter.idents"))
  # the literal, not the column's A/B values
  expect_identical(levels(object), "letter.idents")
})

test_that("vectors of identities are unaffected", {
  object <- pbmc_small
  value <- rep(c("x", "y"), length.out = ncol(object))
  expect_no_warning(Idents(object) <- value)
  expect_setequal(levels(object), c("x", "y"))

  numeric.object <- pbmc_small
  expect_no_warning(Idents(numeric.object) <- seq_len(ncol(numeric.object)))
})

test_that("a single-cell object can still take a scalar identity", {
  object <- pbmc_small[, 1]
  expect_no_warning(Idents(object) <- "only")
  expect_identical(levels(object), "only")
})
