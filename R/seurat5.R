#' @include zzz.R
#' @include assay5.R
#' @include logmap.R
#'
#' @importFrom methods setGeneric
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setClass(
  Class = 'Seurat5',
  slots = c(
    assays = 'list',
    reductions = 'list',
    idents = 'factor',
    cells = 'LogMap',
    meta.data = 'data.frame',
    project = 'character',
    version = 'package_version'
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @method Assays Seurat5
#' @export
#'
Assays.Seurat5 <- function(object, ...) {
  return(names(x = slot(object = object, name = 'assays')))
}

#' @method CastAssay Seurat5
#' @export
#'
CastAssay.Seurat5 <- function(object, to, assay = NULL, layers = NULL, ...) {
  assay <- assay[1L] %||% DefaultAssay(object = object)
  assay <- match.arg(arg = assay, choices = Assays(object = object))
  object[[assay]] <- CastAssay(
    object = object[[assay]],
    to = to,
    layers = layers,
    ...
  )
  validObject(object = object)
  return(object)
}

#' @method Cells Seurat5
#' @export
#'
Cells.Seurat5 <- function(x, assay = NA, ...) {
  assay <- assay[1L] %||% DefaultAssay(object = x)
  if (is.na(x = assay)) {
    return(rownames(x = slot(object = x, name = 'cells')))
  }
  assay <- match.arg(arg = assay, choices = Assays(object = x))
  return(slot(object = x, name = 'cells')[[assay]])
}

#' @method CreateSeurat5Object default
#' @export
#'
CreateSeurat5Object.default <- function(
  counts,
  assay = 'RNA',
  names.field = 1L,
  names.delim = '_',
  meta.data = NULL,
  min.cells = 0,
  min.features = 0,
  project = 'SeuratProject',
  ...
) {
  assay.data <- CreateAssay5Object(
    counts = counts,
    min.cells = min.cells,
    min.features = min.features,
   ...
  )
  return(CreateSeurat5Object(
    counts = assay.data,
    assay = assay,
    names.field = names.field,
    names.delim = names.delim,
    meta.data = meta.data,
    project = project
  ))
}

#' @method CreateSeurat5Object StdAssay
#' @export
#'
CreateSeurat5Object.StdAssay <- function(
  counts,
  assay = 'RNA',
  names.field = 1L,
  names.delim = '_',
  meta.data = NULL,
  project = 'SeuratProject',
  ...
) {
  cells <- LogMap(y = colnames(x = counts))
  cells[[assay]] <- colnames(x = counts)
  if (!isTRUE(x = nzchar(x = Key(object = counts)))) {
    Key(object = counts) <- Key(object = tolower(x = assay), quiet = TRUE)
  }
  assay.list <- list(counts)
  names(x = assay.list) <- assay
  idents <- factor(x = rep_len(x = project, length.out = ncol(x = counts)))
  names(x = idents) <- colnames(x = counts)
  op <- options(Seurat.object.validate = FALSE)
  on.exit(expr = options(op), add = TRUE)
  object <- suppressWarnings(expr = new(
    Class = 'Seurat5',
    assays = assay.list,
    reductions = list(),
    idents = idents,
    cells = cells,
    meta.data = EmptyDF(n = nrow(x = cells)),
    project = project,
    version = utils::packageVersion(pkg = "SeuratObject")
  ))
  options(op)
  # TODO: Calculate nCount and nFeature
  n.calc <- CalcN5(object = counts)
  if (!is.null(x = n.calc)) {
    names(x = n.calc) <- paste(names(x = n.calc), assay, sep = '_')
  }
  validObject(object = object)
  return(object)
}

#' @method DefaultAssay Seurat5
#' @export
#'
DefaultAssay.Seurat5 <- function(object, ...) {
  return(Assays(object = object)[1L])
}

#' @method DefaultAssay<- Seurat5
#' @export
#'
"DefaultAssay<-.Seurat5" <- function(object, ..., value) {
  value <- value[1L]
  assays <- Assays(object = object)
  value <- match.arg(arg = value, choices = assays)
  idx <- which(x = assays == value)
  slot(object = object, name = 'assays') <- c(
    slot(object = object, name = 'assays')[idx],
    slot(object = object, name = 'assays')[-idx]
  )
  return(object)
}

#' @method Features Seurat5
#' @export
#'
Features.Seurat5 <- function(x, assay = NULL, ...) {
  assay <- assay[1L] %||% DefaultAssay(object = x)
  assay <- match.arg(arg = assay, choices = Assays(object = x))
  return(Features(x = x[[assay]], ...))
}

#' @method FetchData Seurat5
#' @export
#'
FetchData.Seurat5 <- function(
  object,
  vars,
  cells = NULL,
  layer = NULL,
  slot = deprecated(),
  ...
) {
  if (is_present(arg = slot)) {
    deprecate_soft(
      when = '4.9.0',
      what = 'FetchData(slot = )',
      with = 'FetchData(layer = )'
    )
    layer <- layer %||% slot
  }
  # Find cells to use
  cells <- cells %||% Cells(x = object, assay = NULL)
  if (is.numeric(x = cells)) {
    cells <- Cells(x = object, assay = NULL)[cells]
  }
  if (is.null(x = vars)) {
    df <- EmptyDF(n = length(x = cells))
    rownames(x = df) <- cells
    return(df)
  }
  # Get a list of all objects to search through and their keys
  object.keys <- Key(object = object)
  # Find all vars that are keyed
  keyed.vars <- sapply(
    X = object.keys,
    FUN = function(key) {
      if (length(x = key) == 0 || !nzchar(x = key)) {
        return(character(length = 0L))
      }
      return(grep(pattern = paste0('^', key), x = vars, value = TRUE))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  keyed.vars <- Filter(f = length, x = keyed.vars)
  # return(keyed.vars)
  # TODO: Handle FOVs
  # Find all keyed. vars
  data.fetched <- lapply(
    X = names(x = keyed.vars),
    FUN = function(x) {
      data.return <- if (x == 'meta.data') {
        md <- gsub(pattern = '^md_', replacement = '', x = keyed.vars[[x]])
        df <- object[[md]][cells, , drop = FALSE]
        names(x = df) <- paste0('md_', names(x = df))
        df
      } else {
        FetchData(
          object = object[[x]],
          vars = keyed.vars[[x]],
          cells = cells,
          layer = layer,
          ...
        )
      }
      return(as.list(x = data.return))
    }
  )
  data.fetched <- unlist(x = data.fetched, recursive = FALSE)
  # data.fetched <- do.call(what = 'cbind', args = data.fetched)
  # Pull vars from object metadata
  meta.vars <- intersect(x = vars, y = names(x = object[[]]))
  meta.vars <- setdiff(x = meta.vars, y = names(x = data.fetched))
  meta.default <- intersect(x = meta.vars, y = rownames(x = object))
  if (length(x = meta.default)) {
    warning(
      "The following variables were found in both object metadata and the default assay: ",
      paste0(meta.default, collapse = ', '),
      "\nReturning metadata; if you want the feature, please use the assay's key (eg. ",
      paste0(Key(object = object)[DefaultAssay(object = object)], meta.default[1L]),
      ")",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  # Pull vars from the default assay
  default.assay <- DefaultAssay(object = object)
  default.vars <- intersect(x = vars, y = rownames(x = object))
  if (length(x = default.vars)) {
    data.fetched[default.vars] <- as.list(x = FetchData(
      object = object[[default.assay]],
      vars = default.vars,
      cells = cells,
      layer = layer,
      ...
    ))
  }
  # Pull identities
  if ('ident' %in% vars && !'ident' %in% colnames(x = object[[]])) {
    data.fetched[['ident']] <- Idents(object = object)[cells]
  }
  # Try to find ambiguous vars
  missing.vars <- setdiff(x = vars, y = names(x = data.fetched))
  # Warn about missing vars
  m2 <- if (length(x = missing.vars) > 10) {
    paste0(' (10 out of ', length(x = missing.vars), ' shown)')
  } else {
    ''
  }
  if (length(x = missing.vars) == length(x = vars)) {
    stop("None of the requested variables were found", call. = FALSE)
  } else if (length(x = missing.vars)) {
    warning(
      "The following requested variables were not found",
      m2,
      ": ",
      paste(head(x = missing.vars, n = 10L), collapse = ', '),
      call. = FALSE,
      immediate. = TRUE
    )
  }
  data.fetched <- as.data.frame(x = data.fetched, row.names = cells)
  return(data.fetched)
}

#' @method HVFInfo Seurat5
#' @export
#'
HVFInfo.Seurat5 <- function(
  object,
  assay = NULL,
  method = NULL,
  status = FALSE,
  layer = NULL,
  strip = TRUE,
  ...
) {
  assay <- assay[1L] %||% DefaultAssay(object = object)
  assay <- match.arg(arg = assay, choices = Assays(object = object))
  return(HVFInfo(
    object = object[[assay]],
    method = method,
    status = status,
    layer = layer,
    strip = strip,
    ...
  ))
}

#' @method Idents Seurat5
#' @export
#'
Idents.Seurat5 <- function(object, ...) {
  return(slot(object = object, name = 'idents'))
}

#' @method Idents<- Seurat5
#' @export
#'
"Idents<-.Seurat5" <- function(object, cells = NULL, drop = FALSE, ..., value) {
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  cells <- intersect(x = cells, y = colnames(x = object))
  cells <- match(x = cells, table = colnames(x = object))
  if (length(x = cells) == 0) {
    warning("Cannot find cells provided", call. = FALSE, immediate. = TRUE)
    return(object)
  }
  idents.new <- if (length(x = value) == 1 && value %in% colnames(x = object[[]])) {
    unlist(x = object[[value]], use.names = FALSE)[cells]
  } else {
    if (is.list(x = value)) {
      value <- unlist(x = value, use.names = FALSE)
    }
    rep_len(x = value, length.out = length(x = cells))
  }
  new.levels <- if (is.factor(x = idents.new)) {
    levels(x = idents.new)
  } else {
    unique(x = idents.new)
  }
  old.levels <- levels(x = object)
  levels <- c(new.levels, old.levels)
  idents.new <- as.vector(x = idents.new)
  idents <- as.vector(x = Idents(object = object))
  idents[cells] <- idents.new
  idents[is.na(x = idents)] <- 'NA'
  levels <- intersect(x = levels, y = unique(x = idents))
  names(x = idents) <- colnames(x = object)
  missing.cells <- which(x = is.na(x = names(x = idents)))
  if (length(x = missing.cells) > 0) {
    idents <- idents[-missing.cells]
  }
  idents <- factor(x = idents, levels = levels)
  slot(object = object, name = 'idents') <- idents
  if (isTRUE(x = drop)) {
    object <- droplevels(x = object)
  }
  validObject(object = object)
  return(object)
}

#' @method Key Seurat5
#' @export
#'
Key.Seurat5 <- function(object, ...) {
  return(c(
    sapply(
      X = .FilterObjects(
        object = object,
        classes.keep = c('KeyMixin')
      ),
      FUN = function(x) {
        return(Key(object = object[[x]]))
      }
    ),
    meta.data = 'md_'
  ))
}

#' @method Project Seurat5
#' @export
#'
Project.Seurat5 <- function(object, ...) {
  return(slot(object = object, name = 'project'))
}

#' @method VariableFeatures Seurat5
#' @export
#'
VariableFeatures.Seurat5 <- function(
  object,
  assay = NULL,
  method = NULL,
  layer = NULL,
  ...
) {
  assay <- assay[1L] %||% DefaultAssay(object = object)
  assay <- match.arg(arg = assay, choices = Assays(object = object))
  return(VariableFeatures(
    object = object[[assay]],
    method = method,
    layer = layer,
    ...
  ))
}

#' @method VariableFeatures<- Seurat5
#' @export
#'
"VariableFeatures<-.Seurat5" <- function(object, assay = NULL, ..., value) {
  assay <- assay[1L] %||% DefaultAssay(object = object)
  assay <- match.arg(arg = assay, choices = Assays(object = object))
  VariableFeatures(object = object[[assay]], ...) <- value
  return(object)
}


# LayerIntersectFeatures <- function(
#   object,
#   assay = NULL,
#   layers
#   ) {
#   assay <- assay %||% DefaultAssay(object = object)
#   feature.df <- object[[assay]]@features@.Data[,layers]
#   features.inte <- rownames(feature.df)[rowSums2(feature.df)  == length(layers)]
#   return(features.inte)
# }




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom utils .DollarNames
#'
#' @method .DollarNames Seurat5
#' @export
#'
.DollarNames.Seurat5 <- .DollarNames.Seurat

#' @method $ Seurat5
#' @export
#'
"$.Seurat5" <- `$.Seurat`

#' @method $<- Seurat5
#' @export
#'
"$<-.Seurat5" <- `$<-.Seurat`

#' @method [[ Seurat5
#' @export
#'
"[[.Seurat5" <- function(x, i, ..., drop = FALSE, na.rm = FALSE) {
  # Pull all meta data
  if (missing(x = i)) {
    df <- slot(object = x, name = 'meta.data')
    row.names(x = df) <- colnames(x = x)
    return(df)
  }
  # Correct invalid `i`
  meta.cols <- names(x = slot(object = x, name = 'meta.data'))
  if (is.numeric(x = i)) {
    stopifnot(any(i <= length(x = meta.cols)))
    i <- i[i <= length(x = meta.cols)]
    i <- meta.cols[as.integer(x = i)]
  }
  stopifnot(is.character(x = i))
  slot.use <- .FindObject(object = x, name = i)
  # Pull cell-level meta data
  if (is.null(x = slot.use)) {
    ic <- intersect(x = i, y = meta.cols)
    if (!length(x = ic)) {
      stop("Cannot find '", i, "' in this Seurat object", call. = FALSE)
    }
    data.return <- slot(object = x, name = 'meta.data')[, i, drop = FALSE, ...]
    row.names(x = data.return) <- colnames(x = x)
    if (isTRUE(x = na.rm)) {
      idx.na <- apply(X = is.na(x = data.return),MARGIN = 1L, FUN = all)
      data.return <- data.return[!idx.na, , drop = FALSE]
    } else {
      idx.na <- rep_len(x = FALSE, length.out = ncol(x = x))
    }
    if (isTRUE(x = drop)) {
      data.return <- unlist(x = data.return, use.names = FALSE)
      names(x = data.return) <- rep.int(
        x = colnames(x = x)[!idx.na],
        times = length(x = i)
      )
    }
    return(data.return)
  }
  # Pull a sub object
  data.return <- slot(object = x, name = slot.use)[[i]]
  return(data.return)
}

#' @method dim Seurat5
#' @export
#'
dim.Seurat5 <- function(x) {
  return(c(
    tryCatch(
      expr = nrow(x = x[[DefaultAssay(object = x)]]),
      error = function(...) {
        return(0L)
      }
    ),
    nrow(x = slot(object = x, name = 'cells'))
  ))
}

#' @method dimnames Seurat5
#' @export
#'
dimnames.Seurat5 <- function(x) {
  return(list(
    DefaultAssay(object = x) %iff% rownames(x = x[[DefaultAssay(object = x)]]),
    Cells(x = x, assay = NA)
  ))
}

#' @importFrom rlang is_bare_list
#'
#' @method dimnames<- Seurat5
#' @export
#'
"dimnames<-.Seurat5" <- function(x, value) {
  op <- options(Seurat.object.validate = FALSE)
  on.exit(expr = options(op), add = TRUE)
  msg <- "Invalid 'dimnames given for a Seurat object"
  if (!is_bare_list(x = value, n = 2L)) {
    stop(msg, call. = FALSE)
  } else if (!all(sapply(X = value, FUN = length) == dim(x = x))) {
    stop(msg, call. = FALSE)
  }
  value <- lapply(X = value, FUN = as.character)
  onames <- dimnames(x = x)
  # Rename cells at the Seurat level
  rownames(x = slot(object = x, name = 'cells')) <- value[[2L]]
  # Rename features/cells at the Assay level
  for (assay in Assays(object = x)) {
    anames <- dimnames(x = x[[assay]])
    afeatures <- MatchCells(
      new = onames[[1L]],
      orig = anames[[1L]],
      ordered = TRUE
    )
    if (!is.null(x = afeatures)) {
      idx <- MatchCells(new = anames[[1L]], orig = onames[[1L]])
      anames[[1L]][idx] <- value[[1L]][afeatures]
    }
    acells <- MatchCells(new = onames[[2L]], orig = anames[[2L]])
    anames[[2L]] <- value[[2L]][acells]
    suppressWarnings(expr = dimnames(x = x[[assay]]) <- anames)
  }
  # Rename idents
  names(x = slot(object = x, name = 'idents')) <- value[[2L]]
  # TODO: Rename features/cells at the DimReduc level
  # Validate and return the Seurat object
  options(op)
  validObject(object = x)
  return(x)
}

#' @method droplevels Seurat5
#' @export
#'
droplevels.Seurat5 <- function(x, ...) {
  Idents(object = x) <- droplevels(x = Idents(object = x), ...)
  return(x)
}

#' @importFrom utils head
#'
#' @method head Seurat5
#' @export
#'
head.Seurat5 <- .head

#' @method levels Seurat5
#' @export
#'
levels.Seurat5 <- levels.Seurat

#' @method levels<- Seurat5
#' @export
#'
"levels<-.Seurat5" <- `levels<-.Seurat`

#' @method merge Seurat5
#' @export
#'
merge.Seurat5 <- function(
  x,
  y,
  # labels = NULL,
  add.cell.ids = NULL,
  collapse = FALSE,
  project = 'SeuratProject',
  ...
) {
  objects <- c(x, y)
  names(x = objects) <- NULL
  projects <- vapply(
    X = objects,
    FUN = Project,
    FUN.VALUE = character(length = 1L)
  )
  # Check cell names
  if (isTRUE(x = is.na(x = add.cell.ids))) {
    add.cell.ids <- as.character(x = seq_along(along.with = objects))
  } else if (isTRUE(x = add.cell.ids)) {
    add.cell.ids <- projects
  }
  if (!is.null(x = add.cell.ids)) {
    if (length(x = add.cell.ids) != length(x = objects)) {
      stop("add.cell.ids length", call. = FALSE)
    } else if (anyDuplicated(x = add.cell.ids)) {
      stop("add.cell.ids duplicate", call. = FALSE)
    }
    for (i in seq_along(along.with = add.cell.ids)) {
      colnames(x = objects[[i]]) <- paste(
        colnames(x = objects[[i]]),
        add.cell.ids[i],
        sep = '_'
      )
    }
  }
  objects <- CheckDuplicateCellNames(object.list = objects)
  cmat <- LogMap(y = unlist(x = lapply(X = objects, FUN = colnames)))
  # Merge identities
  idents.all <- unlist(x = lapply(X = objects, FUN = slot, name = 'idents'))
  idents.all <- idents.all[rownames(x = cmat)]
  # Initialize the combined object
  op <- options(Seurat.object.validate = FALSE)
  on.exit(expr = options(op), add = TRUE)
  obj.combined <- suppressWarnings(expr = new(
    Class = 'Seurat5',
    assays = list(),
    reductions = list(),
    idents = idents.all,
    cells = cmat,
    meta.data = EmptyDF(n = nrow(x = cmat)),
    project = project,
    version = packageVersion(pkg = "SeuratObject")
  ))
  options(op)
  # TODO: Merge assays
  assays.all <- c(
    DefaultAssay(object = x),
    setdiff(
      x = unique(x = unlist(x = lapply(X = objects, FUN = Assays))),
      y = DefaultAssay(object = x)
    )
  )
  for (assay in assays.all) {
    assay.objs <- which(x = vapply(
      X = lapply(X = objects, FUN = names),
      FUN = '%in%',
      FUN.VALUE = logical(length = 1L),
      x = assay
    ))
    if (length(x = assay.objs) == 1L) {
      obj.combined[[assay]] <- objects[[assay.objs]][[assay]]
      next
    }
    idx.x <- assay.objs[[1L]]
    idx.y <- setdiff(x = assay.objs, y = idx.x)
    obj.combined[[assay]] <- merge(
      x = objects[[idx.x]][[assay]],
      y = lapply(X = objects[idx.y], FUN = '[[', assay),
      labels = projects,
      add.cell.ids = NULL,
      collapse = collapse,
      ...
    )
  }
  # Merge meta data
  for (i in seq_along(along.with = objects)) {
    md <- objects[[i]][[]]
    if (!ncol(x = md)) {
      next
    }
    obj.combined[[]] <- md
  }
  # Validate and return the merged object
  validObject(object = obj.combined)
  return(obj.combined)
}

#' @method names Seurat5
#' @export
#'
names.Seurat5 <- function(x) {
  return(.FilterObjects(
    object = x,
    classes.keep = c('StdAssay','DimReduc')
  ))
}

#' @importFrom utils tail
#'
#' @method tail Seurat5
#' @export
#'
tail.Seurat5 <- .tail

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Check for Duplicate Names
#'
#' @inheritParams base::paste
#' @param names A list of character vectors
#' @param stop Throw an error if any duplicate names are found
#'
#' @return \code{names} with all values enforced to be unique across all entries
#'
#' @importFrom rlang is_bare_list
#'
#' @keywords internal
#'
#' @noRd
#'
.CheckDuplicateNames <- function(names, stop = FALSE, sep = '_') {
  if (!is_bare_list(x = names)) {
    stop("'names' must be a list", call. = FALSE)
  } else if (length(x = names) == 1L) {
    return(names)
  }
  if (length(x = Reduce(f = intersect, x = names))) {
    if (isTRUE(x = stop)) {
      stop("Duplicate names provided", call. = FALSE)
    }
    warning(
      "Duplicate names provided, adjusting to enfoce unique names",
      call. = FALSE,
      immediate. = TRUE
    )
    for (i in seq_along(along.with = names)) {
      names[[i]] <- paste(names, i, sep = '_')
    }
  }
  return(names)
}

#' @importFrom rlang caller_env
#'
.DuplicateError <- function(name, cls, error = TRUE) {
  letter <- tolower(x = substr(x = cls[1L], start = 1L, stop = 1L))
  article <- ifelse(test = letter %in% .Vowels(), yes = 'an', no = 'a')
  msg <- paste0("'", name[1L], "' already taken for ", paste(article, cls))
  if (isTRUE(x = error)) {
    abort(message = msg, call = caller_env())
    stop(msg, call. = FALSE)
  }
  return(msg)
}

#' @importFrom methods slotNames
#' @importFrom rlang is_bare_list
#'
.ObjectNames <- function(object) {
  slots <- setdiff(x = slotNames(x = object), y = c('misc', 'tools'))
  slots <- Filter(
    f = function(x) {
      return(is_bare_list(x = slot(object = object, name = x)))
    },
    x = slots
  )
  objects <- lapply(
    X = slots,
    FUN = function(x) {
      obj.names <- names(x = slot(object = object, name = x))
      cls <- switch(
        EXPR = x,
        'assays' = c('StdAssay', 'Assay'),
        'reductions' = 'DimReduc',
        'commands' = 'SeuratCommand',
        x
      )
      objs <- rep_len(x = cls, length.out = length(x = obj.names))
      names(x = objs) <- obj.names
      return(objs)
    }
  )
  return(unlist(x = objects))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Add cell-level meta data
setMethod(
  f = '[[<-',
  signature = c(
    x = 'Seurat5',
    i = 'character',
    j = 'missing',
    value = 'data.frame'
  ),
  definition = function(x, i, ..., value) {
    # Check that the `i` we're adding are present in the data frame
    if (!is.null(x = names(x = value))) {
      i <- match.arg(arg = i, choices = names(x = value), several.ok = TRUE)
    } else if (length(x = i) != ncol(x = value)) {
      stop(
        "Cannot assign ",
        length(x = i),
        " names to ",
        ncol(x = value),
        " bits of meta data",
        call. = FALSE
      )
    }
    # Handle meta data for different cells
    names.intersect <- intersect(x = row.names(x = value), y = colnames(x = x))
    if (length(x = names.intersect)) {
      value <- value[names.intersect, , drop = FALSE]
      if (!nrow(x = value)) {
        stop(
          "None of the cells provided are in this Seurat object",
          call. = FALSE
        )
      }
    } else if (nrow(x = value) == ncol(x = x)) {
      # When no cell names are provided in value, assume it's in cell order
      row.names(x = value) <- colnames(x = x)
    } else {
      # Throw an error when no cell names provided and cannot assume cell order
      stop(
        "Cannot add more or less meta data without cell names",
        call. = FALSE
      )
    }
    # Add the cell-level meta data using the `value = vector` method
    for (n in i) {
      v <- value[[n]]
      names(x = v) <- row.names(x = value)
      x[[n]] <- v
    }
    return(x)
  }
)

# Add cell-level meta data
setMethod(
  f = '[[<-',
  signature = c(
    x = 'Seurat5',
    i = 'missing',
    j = 'missing',
    value = 'data.frame'
  ),
  definition = function(x, ..., value) {
    # Allow removing all meta data
    if (IsMatrixEmpty(x = value)) {
      x[[names(x = x[[]])]] <- NULL
      return(x)
    }
    # If no `i` provided, use the column names from value
    x[[names(x = value)]] <- value
    return(x)
  }
)

# Add dimensional reductions
setMethod(
  f = '[[<-',
  signature = c(
    x = 'Seurat5',
    i = 'character',
    j = 'missing',
    value = 'DimReduc'
  ),
  definition = function(x, i, ..., value) {
    i <- make.names(names = i)
    # Check for duplicate names
    if (i %in% names(x = x)) {
      # Checks for overwriting DimReducs
      if (inherits(x = x[[i]], what = 'DimReduc')) {
        ''
      } else {
        .DuplicateError(name = i, cls = class(x = x[[i]]))
      }
    }
    # Check keys
    Key(object = value) <- .CheckKey(
      key = Key(object = value),
      existing = Key(object = x),
      name = i
    )
    # TODO: Check cells
    # Add the DimReduc
    slot(object = x, name = 'reductions')[[i]] <- value
    slot(object = x, name = 'reductions') <- Filter(
      f = Negate(f = is.null),
      x = slot(object = x, name = 'reductions')
    )
    return(x)
  }
)

# Add cell-level meta data
#' @importFrom methods selectMethod
#'
setMethod(
  f = '[[<-',
  signature = c(x = 'Seurat5', i = 'character', j = 'missing', value = 'factor'),
  definition = function(x, i, ..., value) {
    # Reuse the `value = vector` method
    f <- slot(
      object = selectMethod(
        f = '[[<-',
        signature = c(
          x = 'Seurat5',
          i = 'character',
          j = 'missing',
          value = 'vector'
        )
      ),
      name = '.Data'
    )
    return(f(x = x, i = i, value = value))
  }
)

# Remove objects and cell-level meta data
setMethod(
  f = '[[<-',
  signature = c(x = 'Seurat5', i = 'character', j = 'missing', value = 'NULL'),
  definition = function(x, i, ..., value) {
    # Allow removing multiple objects or bits of cell-level meta data at once
    for (name in i) {
      # Determine the slot to use
      # If no subobject found, check cell-level meta data
      slot.use <- .FindObject(object = x, name = name) %||% 'meta.data'
      switch(
        EXPR = slot.use,
        'meta.data' = {
          # If we can't find the cell-level meta data, throw a warning and move
          # to the next name
          if (!name %in% colnames(x = x[[]])) {
            warning(
              "Cannot find cell-level meta data named ",
              name,
              call. = FALSE,
              immediate. = TRUE
            )
            next
          }
          # Remove the column of meta data
          slot(object = x, name = 'meta.data')[, name] <- value
        },
        'assays' = {
          # Cannot remove the default assay
          if (isTRUE(x = name == DefaultAssay(object = x))) {
            stop("Cannot delete default assay", call. = FALSE)
          }
          # Remove the assay
          slot(object = x, name = slot.use)[[i]] <- value
          # Remove the assay entry from the LogMap
          slot(object = x, name = 'cells') <- droplevels(x = slot(
            object = x,
            name = 'cells'
          ))
        },
        # Remove other subobjects
        slot(object = x, name = slot.use)[[name]] <- value
      )
    }
    validObject(object = x)
    return(x)
  }
)

# Add assays
setMethod(
  f = '[[<-',
  signature = c(
    x = 'Seurat5',
    i = 'character',
    j = 'missing',
    value = 'StdAssay'
  ),
  definition = function(x, i, ..., value) {
    i <- make.names(names = i)
    # Checks for if the assay or name already exists
    if (i %in% names(x = x)) {
      # Check for duplicate names
      if (!inherits(x = x[[i]], what = 'StdAssay')) {
        .DuplicateError(name = i, cls = class(x = x[[i]]))
      }
      if (!identical(x = dim(x = value), y = dim(x = x[[i]]))) {
        warning("different cells/features from existing", call. = FALSE, immediate. = TRUE)
      }
      # if (!all(colnames(x = value) == Cells(x = x, assay = i))) {
      #   stop("different cells", call. = FALSE)
      # }
    }
    if (!all(colnames(x = value) %in% colnames(x = x))) {
        stop("new cells", call. = FALSE)
    }
    # TODO: enable reordering cells in assay
    if (is.unsorted(x = MatchCells(new = colnames(x = value), orig = colnames(x = x), ordered = TRUE))) {
      stop("unorderd cells", call. = FALSE)
    }
    # Check keys
    Key(object = value) <- .CheckKey(
      key = Key(object = value),
      existing = Key(object = x),
      name = i
    )
    # Add the assay
    slot(object = x, name = 'assays')[[i]] <- value
    slot(object = x, name = 'assays') <- Filter(
      f = Negate(f = is.null),
      x = slot(object = x, name = 'assays')
    )
    # TODO: Update the cells LogMap
    slot(object = x, name = 'cells')[[i]] <- colnames(x = value)
    # Return the Seurat object
    validObject(object = x)
    return(x)
  }
)

# Add multiple objects or cell-level meta data
setMethod(
  f = '[[<-',
  signature = c(x = 'Seurat5', i = 'character', j = 'missing', value = 'vector'),
  definition = function(x, i, ..., value) {
    # Add multiple objects
    if (length(x = i) > 1L) {
      value <- rep_len(x = value, length.out = length(x = i))
      for (idx in seq_along(along.with = i)) {
        x[[i[idx]]] <- value[[idx]]
      }
      return(x)
    }
    # Add a column of cell-level meta data
    if (is.null(x = names(x = value))) {
      # Handle cases where new meta data is unnamed
      value <- rep_len(x = value, length.out = ncol(x = x))
      names(x = value) <- colnames(x = x)
    } else {
      # Check cell names for new objects
      names.intersect <- intersect(x = names(x = value), y = colnames(x = x))
      if (!length(x = names.intersect)) {
        stop(
          "No cell overlap between new meta data and Seurat object",
          call. = FALSE
        )
      }
      value <- value[names.intersect]
    }
    df <- EmptyDF(n = ncol(x = x))
    row.names(x = df) <- colnames(x = x)
    df[[i]] <- if (i %in% names(x = x[[]])) {
      x[[i]]
    } else {
      NA
    }
    df[names(x = value), i] <- value
    slot(object = x, name = 'meta.data')[, i] <- df[[i]]
    validObject(object = x)
    return(x)
  }
)

setMethod(
  f = 'show',
  signature = c(object = 'Seurat5'),
  definition = function(object) {
    nfeatures <- sum(vapply(
      X = Assays(object = object),
      FUN = function(x) {
        return(nrow(x = object[[x]]))
      },
      FUN.VALUE = numeric(length = 1L),
      USE.NAMES = FALSE
    ))
    nassays <- length(x = Assays(object = object))
    cat(
      "A Seurat (v5) object for the",
      Project(object = object),
      "project\n"
    )
    cat(
      nfeatures,
      'features across',
      ncol(x = object),
      'samples within',
      nassays,
      ifelse(test = nassays == 1, yes = 'assay', no = 'assays'),
      '\n'
    )
    cat(
      'Active assay:',
      DefaultAssay(object = object),
      paste0('(', nrow(x = object), ' features)')
      # paste0(
      #   "(",
      #   nrow(x = object),
      #   ' features, ',
      #   length(x = VariableFeatures(object = object)),
      #   ' variable features)'
      # )
    )
    nother <- nassays - 1L
    if (as.logical(x = nother)) {
      cat(
        '\n',
        nother,
        'other',
        ifelse(test = nother == 1, yes = 'assay', no = 'assays'),
        'present:',
        strwrap(x = paste(
          setdiff(
            x = Assays(object = object),
            y = DefaultAssay(object = object)
          ),
          collapse = ', '
        ))
      )
    }
    return(invisible(x = NULL))
  }
)

setValidity(
  Class = 'Seurat5',
  method = function(object) {
    if (isFALSE(x = getOption(x = "Seurat.object.validate", default = TRUE))) {
      warning("Not validating Seurat objects", call. = FALSE, immediate. = TRUE)
      return(TRUE)
    }
    valid <- NULL
    # TODO: Check assays
    if (!IsNamedList(x = slot(object = object, name = 'assays'))) {
      valid <- c(valid, "'assays' must be a named list")
    } else {
      for (assay in Assays(object = object)) {
        if (!inherits(x = object[[assay]], what = 'StdAssay')) {
          valid <- c(
            valid,
            "Objects in the 'assays' list must be 'Assay' objects"
          )
          break
        }
        if (!assay %in% colnames(x = slot(object = object, name = 'cells'))) {
          valid <- c(
            valid,
            "All assays must have an entry in 'cells'"
          )
        } else if (any(colnames(x = object[[assay]]) != Cells(x = object, assay = assay))) {
          valid <- c(
            valid,
            paste0(
              "All assays must have the same cells as listed in 'cells' (offending: ",
              assay,
              ")"
            )
          )
        }
        if (!isTRUE(x = nzchar(x = Key(object = object[[assay]])))) {
          valid <- c(valid, "All assays must have a key")
        }
      }
    }
    # TODO: Check reductions
    # TODO: Check metadata
    # TODO: Check idents
    idents <- Idents(object = object)
    if (length(x = idents) != ncol(x = object)) {
      valid <- c(valid, "'idents' must have the same length as number of cells")
    } else if (is.null(x = names(x = idents))) {
      valid <- c(valid, "'idents' must be named")
    } else if (!all(names(x = idents) == colnames(x = object))) {
      valid <- c(valid, "'idents' must be named with the cells")
    }
    # TODO: Check project
    proj <- Project(object = object)
    if (length(x = proj) != 1L) {
      valid <- c(valid, "'project' must be a 1-length character")
    } else if (is.na(x = proj)) {
      valid <- c(valid, "'project' cannot be NA")
    } else if (!nzchar(x = proj)) {
      valid <- c(valid, "'project' cannot be an empty character")
    }
    # TODO: Check version
    if (length(x = slot(object = object, name = 'version')) > 1) {
      valid <- c(valid, "Only one version is allowed")
    }
    return(valid %||% TRUE)
  }
)

