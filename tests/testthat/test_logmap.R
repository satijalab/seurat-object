# Tests for the LogMap class

test_that("`labels` generic works as expected for `LogMap` instances", {
  # Instantiate and populate a LogMap instance for testing.
  map <- LogMap(paste0("value_", 1:6))
  map[["label_a"]] <- c(1, 3)
  map[["label_b"]] <- c(2, 4)
  map[["label_c"]] <- c(2, 4, 6)
  map[["label_d"]] <- c(2, 4, 6)
  map[["label_e"]] <- 2
  map[["label_f"]] <- 4

  # Labels can be fetched for specified values.
  values <- c("value_1", "value_3")
  result_key <- c(value_1 = "label_a", value_3 = "label_a")
  expect_identical(result_key, labels(map, values = values))

  # For values with multiple labels, the first label is returned by default.
  values <- c("value_2", "value_4")
  result_key <- c(value_2 = "label_b", value_4 = "label_b")
  result <- labels(map, values = values)
  expect_identical(result_key, result)

  # The last value for each label can also be fetched.
  values <- c("value_1", "value_2", "value_4")
  result_key <- c(value_1 = "label_a", value_2 = "label_e", value_4 = "label_f")
  result <- labels(map, values = values, select = "last")
  expect_identical(result_key, result)

  # It is also possible to fetch the label that is shared by the most values
  # in the requested set. If multiple labels are equally common, the first
  # is returned.
  values <- c("value_2", "value_4", "value_6")
  result_key <- c(value_2 = "label_c", value_4 = "label_c", value_6 = "label_c")
  result <- labels(map, values = values, select = "common")
  expect_identical(result_key, result)

  # Label resolution is based on the column order of the underlying matrix
  label_order <- c(
    "label_a",
    "label_e",
    "label_b",
    "label_c",
    "label_e",
    "label_f"
  )
  map <- map[, label_order, drop = FALSE]
  values <- c("value_2", "value_4")
  result_key <- c(value_2 = "label_e", value_4 = "label_b")
  result <- labels(map, values = values)
  expect_identical(result_key, result)

  # The output order is taken from the `values` parameter.
  values <- c("value_3", "value_1")
  result_key <- c(value_3 = "label_a", value_1 = "label_a")
  result <- labels(map, values = values)
  expect_identical(result_key, result)

  # Values without labels are excluded in the return.
  values <- c("value_1", "value_5", "value_6")
  result_key <- c(value_1 = "label_a", value_6 = "label_c")
  result <- labels(map, values = values)
  expect_identical(result_key, result)

  # If no labels are found, an error is raised.
  expect_no_error(labels(map, "value_5"))
})
