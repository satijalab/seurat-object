set.seed(42)

make_coords <- function(ncell = 3, npoint = 6) {
  do.call(rbind, lapply(seq_len(ncell), function(i) {
    data.frame(
      x = runif(npoint) * 10 + i * 20,
      y = runif(npoint) * 10,
      cell = paste0("cell", i)
    )
  }))
}

test_that("well-formed segmentations are unaffected", {
  coords <- make_coords()
  segmentation <- expect_no_warning(CreateSegmentation(coords))
  expect_s4_class(segmentation, "Segmentation")
  expect_setequal(Cells(segmentation), paste0("cell", 1:3))
})

test_that("NA coordinates name the cells responsible", {
  coords <- make_coords()
  coords$x[4] <- NA          # cell1
  coords$y[13] <- NA         # cell3
  # sp reports only "NA values in coordinates", which is unactionable when a
  # handful of rows in a large segmentation file are at fault
  expect_error(CreateSegmentation(coords), "cell1")
  expect_error(CreateSegmentation(coords), "cell3")
  expect_error(CreateSegmentation(coords), "2 cells")
  expect_error(CreateSegmentation(coords), "Remove or repair")
})

test_that("many bad cells are summarised rather than listed in full", {
  coords <- make_coords(ncell = 12)
  coords$x[coords$cell %in% paste0("cell", 1:8)] <- NA
  message <- tryCatch(CreateSegmentation(coords), error = conditionMessage)
  expect_match(message, "8 cells")
  expect_match(message, "3 more")
})

test_that("cells with too few boundary points are named once", {
  coords <- rbind(
    make_coords(),
    data.frame(x = c(1, 2), y = c(1, 2), cell = "tiny")
  )
  # sp warns once per polygon without saying which one
  warnings <- character()
  segmentation <- withCallingHandlers(
    CreateSegmentation(coords),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(warnings, 1L)
  expect_match(warnings, "tiny")
  expect_match(warnings, "fewer than 4")
  # the cell is still kept
  expect_true("tiny" %in% Cells(segmentation))
})

test_that("a single bad cell is reported in the singular", {
  coords <- make_coords()
  coords$x[1] <- NA
  expect_error(CreateSegmentation(coords), "1 cell:")
})
