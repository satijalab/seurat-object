.SetSeuratCompat <- local({
  seurat.version <- NULL
  function(pkgname, pkgpath) {
    current <- .RoundVersion(current = packageVersion(pkg = pkgname))
    if (pkgname == 'Signac') {
      if (is.null(x = seurat.version)) {
        seurat.version <<- ifelse(
          test = paste(current, collapse = '.') >= '1.12.9000',
          yes = '5.0.0',
          no = '4.4.0'
        )
      }
      return(invisible(x = NULL))
    }
    seurat.version <<- paste(current, collapse = '.')
    if (!is.null(x = seurat.version) && seurat.version < '5.0.0') {
      options(
        Seurat.object.assay.brackets = 'v3',
        Seurat.object.assay.version = 'v3'
      )
    }
    return(invisible(x = NULL))
  }
})

.GetSeuratCompat <- local(
  envir = environment(fun = .SetSeuratCompat),
  function(pkgname = 'Seurat') {
    if (is.null(x = seurat.version) && isNamespaceLoaded(name = pkgname)) {
      .SetSeuratCompat(pkgname = pkgname)
    }
    return(seurat.version %||% '5.0.0')
  }
)

.SeuratCompatMessage <- local(
  envir = environment(fun = .SetSeuratCompat),
  function(pkgname, pkgpath) {
    seurat <- .GetSeuratCompat()
    if (!is.null(x = seurat) && seurat < '5.0.0') {
      options(
        Seurat.object.assay.brackets = 'v3',
        Seurat.object.assay.version = 'v3'
      )
      version <- paste0('v', substr(x = seurat, start = 1L, stop = 1L))
      packageStartupMessage(paste(
        strwrap(x = paste(
          pkgname,
          switch(
            EXPR = pkgname,
            Seurat = version,
            "built for for SeuratObject v4"
          ),
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
