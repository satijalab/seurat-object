# SeuratObject v5.3.0

## Test environments
* local Ubuntu 20.04.6 install (R 4.3.2)
* local macOS 15.6 (R 4.5.1)
* win-builder
  * R-oldrelease (R 4.4.3)
  * R-release (R 4.5.2)
  * R-devel (R 4.6.0)

## R CMD check results

There were no ERRORs or WARNINGs.

There were two NOTEs given by win-builder (oldrel):
- - -
```
* checking CRAN incoming feasibility ... [21s] NOTE
Maintainer: 'Rahul Satija <seurat@nygenome.org>'

Suggests or Enhances not in mainstream repositories:
  BPCells
Availability using Additional_repositories specification:
  BPCells   yes   https://bnprks.r-universe.dev
```

Rahul Satija is the maintainer and the email is correct. 

BPCells is hosted on R-universe and used conditionally in SeuratObject.

```
* checking package dependencies ... NOTE
Package suggested but not available for checking: 'BPCells'
```

BPCells is hosted on R-universe and used conditionally in SeuratObject.

## Downstream dependencies

We checked 63 reverse dependencies (25 from CRAN + 38 from Bioconductor), comparing R CMD check results across CRAN and dev versions of this package.

Reverse dependency checks gave ERRORs in 5 packages: Platypus (examples), SCpubr (tests), Seurat (examples, tests), XYomics (vignettes), scpoisson (examples, vignettes). All of these ERRORs are regarding the defunctation of the `slot` argument in `GetAssayData` and `SetAssayData`, which was deprecated in our package three releases ago (v5.0.0, 2023-10-26) and announced then as slated for eventual defunctation.

Besides this point, this update does not affect the functionality of any of the dependencies.