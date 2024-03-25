#' A small example version of the PBMC dataset
#'
#' A subsetted version of 10X Genomics' 3k PBMC dataset
#'
#' @format A Seurat object with the following slots filled
#' \describe{
#'   \item{assays}{
#'   \itemize{Currently only contains one assay ("RNA" - scRNA-seq expression data)
#'   \item{counts - Raw expression data}
#'   \item{data - Normalized expression data}
#'   \item{scale.data - Scaled expression data}
#'   \item{var.features - names of the current features selected as variable}
#'   \item{meta.features - Assay level metadata such as mean and variance}
#'    }}
#'   \item{meta.data}{Cell level metadata}
#'   \item{active.assay}{Current default assay}
#'   \item{active.ident}{Current default idents}
#'   \item{graphs}{Neighbor graphs computed, currently stores the SNN}
#'   \item{reductions}{Dimensional reductions: currently PCA and tSNE}
#'   \item{version}{Seurat version used to create the object}
#'   \item{commands}{Command history}
#' }
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k}
#'
"pbmc_small"

.LoadPkgData <- function(ds, ..., mode = c('load', 'resave'), env = NULL) {
  resave_data_others <- function(src) {
    env <- new.env()
    sys.source(file = src, envir = env, chdir = TRUE)
    return(env)
  }
  data.dir <- system.file('data', package = 'SeuratObject', mustWork = TRUE)
  datasets <- list.files(path = data.dir, pattern = "\\.R$")
  if (!length(x = datasets)) {
    rlang::warn(
      message = "Loading datasets by function works only under `devtools::load_all()`"
    )
    return(invisible(x = NULL))
  }
  ds <- match.arg(arg = ds, choices = datasets)
  mode <- match.arg(arg = mode)
  ds.env <- resave_data_others(src = file.path(data.dir, ds))
  if (mode == 'resave') {
    save(
      list = ls(envir = ds.env, all.names = TRUE),
      file = file.path(
        data.dir,
        sub(pattern = '\\.R$', replacement = '.rda', x = ds)
      ),
      compress = TRUE,
      compression_level = 9L,
      envir = ds.env
    )
    return(invisible(x = NULL))
  }
  ds.env <- as.list(x = ds.env)
  if (is.environment(x = env)) {
    for (i in names(x = ds.env)) {
      env[[i]] <- ds.env[[i]]
    }
    return(invisible(x = ds.env))
  }
  if (length(x = ds.env) == 1L) {
    return(ds.env[[1L]])
  }
  return(ds.env)
}

.PBMCsmall <- \(mode = 'load', env = .GlobalEnv) .LoadPkgData(
  ds = 'pbmc_small',
  mode = mode,
  env = env
)
