
pbmc_small <- local({
  callcheck <- "resave_data_others" %in% unlist(x = lapply(
    X = sys.calls(),
    FUN = as.character
  ))
  if (!isTRUE(x = callcheck)) {
    return(NULL)
  }
  # Check required packages
  pkgcheck <- requireNamespace("rprojroot", quietly = TRUE) &&
    requireNamespace("Seurat", quietly = TRUE) &&
    utils::packageVersion(pkg = "Seurat") >= "5.0.0"
  if (!isTRUE(x = pkgcheck)) {
    return(NULL)
  }

  library(SeuratObject)
  op <- options(Seurat.object.assay.version = "v3")

  # Find the raw inputs
  root <- rprojroot::find_package_root_file()
  raw <- file.path(root, "inst", "extdata", "raw", "pbmc_small")
  filecheck <-  dir.exists(raw) &&
    all(file.exists(file.path(raw, c("counts.mtx", "features.txt", "cells.txt"))))
  if (!isTRUE(x = filecheck)) {
    return(NULL)
  }

  # Read in the raw data
  mat <- methods::as(
    object = Matrix::readMM(file = file.path(raw, "counts.mtx")),
    Class = "CsparseMatrix"
  )
  dimnames(x = mat) <- list(
    readLines(con = file.path(raw, "features.txt")),
    readLines(con = file.path(raw, "cells.txt"))
  )

  # Construct the `Seurat` object
  pbmc_small <- CreateSeuratObject(counts = mat, project = "pbmc_small")
  if (!inherits(x = pbmc_small[["RNA"]], what = "Assay")) {
    return(NULL)
  }

  # Process the object
  pbmc_small <- Seurat::NormalizeData(
    object = pbmc_small,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  ) |>
    Seurat::FindVariableFeatures(
      selection.method = "vst",
      loess.span = 0.3,
      clip.max = "auto",
      num.bin = 20L,
      binning.method = "equal_width",
      nfeatures = 20L
    ) |>
    Seurat::ScaleData(
      model.use = "linear",
      use.umi = FALSE,
      do.scale = TRUE,
      do.center = TRUE,
      scale.max = 10L,
      block.size = 1000L,
      min.cells.to.block = 80L
    ) |>
    Seurat::RunPCA(
      npcs = 20L,
      rev.pca = FALSE,
      weight.by.var = TRUE,
      reduction.name = "pca",
      reduction.key = Key("PC", quiet = TRUE),
      seed.use = 42L
    ) |>
    Seurat::JackStraw(dims = 10L, num.replicate = 10L) |>
    Seurat::ScoreJackStraw(dims = 1:5, score.thresh = 1e-5) |>
    Seurat::FindNeighbors(
      dims = 1:10,
      k.param = 30L,
      prune.SNN = 1/15,
      nn.eps = 0L
    ) |>
    Seurat::FindClusters(resolution = c(0.8, 1)) |>
    Seurat::RunTSNE(
      dims = 1:5,
      perplexity = 5L,
      seed.use = 1L,
      check_duplicates = FALSE,
      reduction.name = "tsne",
      reduction.key = Key("tSNE", quiet = TRUE)
    ) |>
    Seurat::RunUMAP(
      dims = 1:5,
      reduction.name = "umap",
      reduction.key = Key("UMAP", quiet = TRUE)
    ) |> Seurat::ProjectDim()

  # Add some meta data
  set.seed(seed = 42L)
  pbmc_small[[]] <- data.frame(
    letter.idents = factor(x = sample(
      x = c("A", "B"),
      size = ncol(x = pbmc_small),
      replace = TRUE
    )),
    groups = sample(
      x = c("g1", "g2"),
      size = ncol(x = pbmc_small),
      replace = TRUE
    ),
    row.names = colnames(x = pbmc_small)
  )

  # Add a v5 assay
  rna5 <- methods::as(object = pbmc_small[["RNA"]], Class = "Assay5")
  Key(rna5) <- Key("RNA5", quiet = TRUE)
  pbmc_small[["RNA5"]] <- rna5

  # Return
  pbmc_small
})
