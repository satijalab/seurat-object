.SetSeuratCompat <- local({
  seurat.version <- NULL
  function(...) {
    current <- .RoundVersion(current = packageVersion(pkg = 'Seurat'))
    seurat.version <<- paste(current, collapse = '.')
    return(invisible(x = NULL))
  }
})

.GetSeuratCompat <- local(
  envir = environment(fun = .SetSeuratCompat),
  function() {
    if (is.null(x = seurat.version) && isNamespaceLoaded(name = 'Seurat')) {
      .SetSeuratCompat()
    }
    return(seurat.version %||% '5.0.0')
  }
)

.SeuratCompatMessage <- local(
  envir = environment(fun = .SetSeuratCompat),
  function(...) {
    seurat <- .GetSeuratCompat()
    if (!is.null(x = seurat) && seurat < '5.0.0') {
      options(
        Seurat.object.assay.brackets = 'v3',
        Seurat.object.assay.version = 'v3'
      )
      version <- paste0('v', substr(x = seurat, start = 1L, stop = 1L))
      packageStartupMessage(paste(
        strwrap(x = paste(
          "Seurat",
          version,
          "was just loaded with SeuratObject v5;",
          "disabling v5 assays and validation routines,",
          "and ensuring assays work in strict v3/v4 compatibility mode"
        )),
        collapse = '\n'
      ))
    }
    return(invisible(x = NULL))
  }
)
