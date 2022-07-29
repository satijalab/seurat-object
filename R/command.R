#' @include zzz.R
#' @include generics.R
#' @importFrom methods new slot slot<-
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The \code{SeuratCommand} Class
#'
#' The \code{SeuratCommand} is used for logging commands that are run
#' on a \code{Seurat} object; it stores parameters and timestamps
#'
#' @slot name Command name
#' @slot time.stamp Timestamp of when command was tun
#' @slot assay.used Optional name of assay used to generate
#' \code{SeuratCommand} object
#' @slot call.string String of the command call
#' @slot params List of parameters used in the command call
#'
#' @name SeuratCommand-class
#' @rdname SeuratCommand-class
#' @exportClass SeuratCommand
#'
#' @family command
#'
#' @aliases SeuratCommand
#'
setClass(
  Class = 'SeuratCommand',
  slots = c(
    name = 'character',
    time.stamp = 'POSIXct',
    assay.used = 'OptionalCharacter',
    call.string = 'character',
    params = 'ANY'
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Log a command
#'
#' Logs command run, storing the name, timestamp, and argument list. Stores in
#' the Seurat object
#'
#' @param object Name of Seurat object
#' @param return.command Return a \code{\link{SeuratCommand}} object instead
#'
#' @return If \code{return.command}, returns a \code{\link{SeuratCommand}}
#' object; otherwise, returns the Seurat object with command stored
#'
#' @export
#'
#' @family command
#'
#' @seealso \code{\link{Command}}
#'
LogSeuratCommand <- function(object, return.command = FALSE) {
  time.stamp <- Sys.time()
  object <- UpdateSlots(object = object)
  #capture function name
  which.frame <- sys.nframe() - 1
  if (which.frame < 1) {
    stop("'LogSeuratCommand' cannot be called at the top level", call. = FALSE)
  }
  if (as.character(x = sys.calls()[[1]])[1] == "do.call") {
    call.string <- deparse(expr = sys.calls()[[1]])
    command.name <- as.character(x = sys.calls()[[1]])[2]
  } else {
    command.name <- as.character(x = deparse(expr = sys.calls()[[which.frame]]))
    command.name <- gsub(
      pattern = "\\.Seurat",
      replacement = "",
      x = command.name
    )
    call.string <- command.name
    command.name <- ExtractField(
      string = command.name,
      field = 1,
      delim = "\\("
    )
  }
  #capture function arguments
  argnames <- names(x = formals(fun = sys.function(which = sys.parent(n = 1))))
  argnames <- grep(
    pattern = "object",
    x = argnames,
    invert = TRUE,
    value = TRUE
  )
  argnames <- grep(
    pattern = "anchorset",
    x = argnames,
    invert = TRUE,
    value = TRUE
  )
  argnames <- grep(
    pattern = "\\.\\.\\.",
    x = argnames,
    invert = TRUE,
    value = TRUE
  )
  params <- list()
  p.env <- parent.frame(n = 1)
  argnames <- intersect(x = argnames, y = ls(name = p.env))
  # fill in params list
  for (arg in argnames) {
    param_value <- get(x = arg, envir = p.env)
    if (inherits(x = param_value, what = 'Seurat')) {
      next
    }
    #TODO Institute some check of object size?
    params[[arg]] <- param_value
  }
  # check if function works on the Assay and/or the DimReduc Level
  assay <- params[["assay"]]
  reduction <- params[["reduction"]]
  # Get assay used for command
  cmd.assay <- assay %||% (reduction %iff% if (inherits(x = reduction, what = 'DimReduc')) {
    DefaultAssay(object = reduction)
  } else if (reduction %in% Reductions(object = object)) {
    DefaultAssay(object = object[[reduction]])
  })
  if (inherits(x = reduction, what = 'DimReduc')) {
    reduction <- 'DimReduc'
  }
  # rename function name to include Assay/DimReduc info
  if (length(x = assay) == 1) {
    command.name <- paste(command.name, assay, reduction, sep = '.')
  }
  command.name <- sub(
    pattern = "[\\.]+$",
    replacement = "",
    x = command.name,
    perl = TRUE
  )
  command.name <- sub(pattern = "\\.\\.", replacement = "\\.", x = command.name, perl = TRUE)
  # store results
  seurat.command <- new(
    Class = 'SeuratCommand',
    name = command.name,
    params = params,
    time.stamp = time.stamp,
    call.string = call.string,
    assay.used = cmd.assay
  )
  if (isTRUE(x = return.command)) {
    return(seurat.command)
  }
  object[[command.name]] <- seurat.command
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname DefaultAssay
#' @export
#' @method DefaultAssay SeuratCommand
#'
DefaultAssay.SeuratCommand <- function(object, ...) {
  object <- UpdateSlots(object = object)
  return(slot(object = object, name = 'assay.used'))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @inherit .DollarNames.Seurat title
#'
#' @description Autocompletion for \code{$} access on a
#' \code{\link{SeuratCommand}} object
#'
#' @inheritParams utils::.DollarNames
#' @param x A \code{\link{SeuratCommand}} object
#'
#' @return The parameter name matches for \code{pattern}
#'
#' @importFrom utils .DollarNames
#'
#' @method .DollarNames SeuratCommand
#' @export
#'
#' @family command
#'
".DollarNames.SeuratCommand" <- function(x, pattern = '') {
  return(.DollarNames(x = slot(object = x, name = "params"), pattern = pattern))
}

#' Command Log Parameter Access
#'
#' Pull parameter values from a \code{\link{SeuratCommand}} object
#'
#' @inheritParams .DollarNames.SeuratCommand
#' @param i A parameter name
#'
#' @return The value for parameter \code{i}
#'
#' @method $ SeuratCommand
#' @export
#'
#' @family command
#'
#' @examples
#' cmd <- pbmc_small[["NormalizeData.RNA"]]
#' cmd$normalization.method
#'
"$.SeuratCommand" <- function(x, i, ...) {
  params <- slot(object = x, name = "params")
  return(params[[i]])
}

#' Command Log Data Access
#'
#' Access data from a \code{SeuratCommand} object
#'
#' @inheritParams .DollarNames.SeuratCommand
#' @param i The name of a command log slot
#'
#' @return \code{[}: Slot \code{i} from \code{x}
#'
#' @method [ SeuratCommand
#' @export
#'
#' @family command
#'
#' @examples
#' cmd <- pbmc_small[["NormalizeData.RNA"]]
#' cmd["call.string"]
#'
"[.SeuratCommand" <- function(x, i, ...) {
  i <- arg_match(arg = i, values = slotNames(x = x))
  return(slot(object = x, name = i))
}

#' Coerce a \code{SeuratCommand} to a list
#'
#' @inheritParams .DollarNames.SeuratCommand
#' @param complete Include slots besides just parameters
#' (eg. call string, name, timestamp)
#'
#' @return A list with the parameters and, if \code{complete = TRUE},
#' the call string, name, and timestamp
#'
#' @method as.list SeuratCommand
#' @export
#'
#' @family command
#'
#' @examples
#' cmd <- pbmc_small[["NormalizeData.RNA"]]
#' as.list(cmd)
#' as.list(cmd, complete = TRUE)
#'
as.list.SeuratCommand <- function(x, complete = FALSE, ...) {
  CheckDots(...)
  cmd <- slot(object = x, name = 'params')
  if (isTRUE(x = complete)) {
    cmd <- append(
      x = cmd,
      values = sapply(
        X = setdiff(x = slotNames(x = x), y = 'params'),
        FUN = slot,
        object = x,
        simplify = FALSE,
        USE.NAMES = TRUE
      ),
      after = 0
    )
  }
  for (i in seq_along(along.with = cmd)) {
    if (is.character(x = cmd[[i]])) {
      cmd[[i]] <- paste(trimws(x = cmd[[i]]), collapse = ' ')
    }
  }
  return(cmd)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Command Log Overview
#'
#' Overview of a \code{\link{SeuratCommand}} object
#'
#' @template return-show
#'
#' @keywords internal
#'
#' @concept command
#'
#' @examples
#' cmd <- pbmc_small[["NormalizeData.RNA"]]
#' cmd
#'
setMethod(
  f = 'show',
  signature = 'SeuratCommand',
  definition = function(object) {
    params <- slot(object = object, name = "params")
    params <- params[sapply(X = params, FUN = class) != "function"]
    cat(
      "Command: ", slot(object = object, name = "call.string"), '\n',
      "Time: ", as.character(slot(object = object, name = "time.stamp")), '\n',
      sep = ""
    )
    for (p in seq_len(length.out = length(x = params))) {
      cat(
        names(params[p]), ":", params[[p]], "\n"
      )
    }
  }
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
