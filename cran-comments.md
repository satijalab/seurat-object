# SeuratObject v5.2.0

## Test environments
* local ubuntu 20.04 install (R 4.3.2)
* local macOS 15.6 (R 4.4.0, R 4.5.1)
* win-builder (oldrelease, release, devel)
* mac-builder (release, devel)

## R CMD check results

There were no ERRORs or WARNINGs.

There were three NOTEs given by win-builder (oldrel):
- - -
> * checking CRAN incoming feasibility ... [21s] NOTE
Maintainer: 'Rahul Satija <seurat@nygenome.org>'

The maintainer remains Rahul Satija and the email is correct.
 
> Suggests or Enhances not in mainstream repositories:
>   BPCells
> Availability using Additional_repositories specification:
>   BPCells   yes   https://bnprks.r-universe.dev

> * checking package dependencies ... NOTE
> Package suggested but not available for checking: 'BPCells'

BPCells is hosted on R-universe and used conditionally in SeuratObject.

> * checking DESCRIPTION meta-information ... NOTE
> Author field differs from that derived from Authors@R

There seems to be a slight difference in ORCID formatting between Author and Authors@R; the information in both is the same.
## Downstream dependencies
There are 2 packages that depend on SeuratObject: Seurat and tidyseurat. This update does not impact their functionality.

There are 11 packages that import SeuratObject: APackOfTheClones, bbknnR, CAMML, Platypus, RepeatedHighDim, scAnnotate, scaper, scCustomize, scpoisson, SeuratExplorer, and Signac. This update does not impact their functionality.

There are 14 packages that suggest SeuratObject: cellpypes, ClustAssess, CytoSimplex, mxfda, RESET, rliger, scOntoMatch, SCpubr, SingleCellComplexHeatmap, singleCellHaystack, spatialGE, SpaTopic, SPECK, and VAM. This update does not impact their functionality.
