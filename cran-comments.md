# SeuratObject v5.0.1

## Test environments
* local ubuntu 20.04 install, R 4.3.2
* local macOS 14.1, R 4.4.0
* win-builder (oldrelease, release, devel)
* mac-builder (devel)

We were unable to test on r-release on mac-builder because the portal seemed to point to the wrong version.

## R CMD check results

There were no ERRORs or WARNINGs

There were two NOTEs

> * checking CRAN incoming feasibility ... NOTE
> Maintainer: 'Rahul Satija <seurat@nygenome.org>'
> 
> New maintainer:
>   Rahul Satija <seurat@nygenome.org>
> Old maintainer(s):
>   Paul Hoffman <seurat@nygenome.org>
> Suggests or Enhances not in mainstream repositories:
>   BPCells
> Availability using Additional_repositories specification:
>   BPCells   yes  https://bnprks.r-universe.dev

The new maintainer is Rahul Satija, the email address has remained the same.

> * checking package dependencies ... NOTE
> Package suggested but not available for checking: 'BPCells'

BPCells is hosted on R-universe and used conditionally in SeuratObject.

## Downstream dependencies
There are 2 packages that depend on SeuratObject: Seurat, and tidyseurat; this update does not impact their functionality.

There are 9 packages that import Seurat: APackOfTheClones, bbknnR, CAMML, Platypus, scAnnotate, scaper, scCustomize, scpoisson, and Signac; this update does not impact their functionality.

There are 10 packages that suggest Seurat: cellpypes, CytoSimplex, RESET, rliger, scOntoMatch, SCpubr, singleCellHaystack, SpaTopic, SPECK, and VAM; this update does not impact their functionality.
