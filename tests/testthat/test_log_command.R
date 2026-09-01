set.seed(42)

make_object <- function(nfeat = 6, ncell = 5) {
  m <- matrix(rpois(nfeat * ncell, lambda = 3), nrow = nfeat, ncol = ncell)
  dimnames(m) <- list(paste0("g", seq_len(nfeat)), paste0("c", seq_len(ncell)))
  suppressWarnings(CreateSeuratObject(counts = as.sparse(m)))
}

logged <- function(object) names(slot(object, "commands"))

MyStep <- function(object, alpha = 1) {
  suppressWarnings(LogSeuratCommand(object))
}

WrappedStep <- function(object, alpha = 1) {
  withCallingHandlers(
    suppressWarnings(LogSeuratCommand(object)),
    warning = function(w) invokeRestart("muffleWarning")
  )
}

test_that("a direct call is named after its function", {
  expect_identical(logged(MyStep(make_object())), "MyStep")
})

test_that("a call wrapped in withCallingHandlers is named after the function", {
  # the wrapper frame sits between the function and LogSeuratCommand, and was
  # being recorded as the command
  expect_identical(logged(WrappedStep(make_object())), "WrappedStep")
})

test_that("a call made through do.call is named after the function", {
  object <- make_object()
  expect_identical(logged(do.call(MyStep, list(object))), "MyStep")
  # and when do.call is not itself at the bottom of the stack
  nested <- function(object) do.call(MyStep, list(object))
  expect_identical(logged(nested(object)), "MyStep")
})

test_that("other evaluation wrappers do not become the command name", {
  object <- make_object()
  expect_identical(logged(suppressWarnings(MyStep(object))), "MyStep")
  expect_identical(logged(suppressMessages(MyStep(object))), "MyStep")
  expect_identical(logged(tryCatch(MyStep(object), error = function(e) e)), "MyStep")
  expect_identical(logged(eval(quote(MyStep(object)))), "MyStep")
})

test_that("the recorded arguments survive the change", {
  object <- suppressWarnings(MyStep(make_object(), alpha = 7))
  command <- slot(object, "commands")[["MyStep"]]
  expect_identical(slot(command, "name"), "MyStep")
  expect_equal(slot(command, "params")$alpha, 7)
})
