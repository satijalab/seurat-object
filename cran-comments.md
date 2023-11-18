# SeuratObject v5.0.1

## Test environments
* local ubuntu 22.04 install, R 4.2.3
* win-builder (oldrelease, release)

We were unable to test on r-devel on win-builder due to insufficient Matrix version

## R CMD check results

There were no ERRORs or WARNINGs

There were two NOTEs

> * checking CRAN incoming feasibility ... NOTE
> Suggests or Enhances not in mainstream repositories:
>   BPCells
> Availability using Additional_repositories specification:
>   BPCells   yes  https://bnprks.r-universe.dev

> * checking package dependencies ... NOTE
> Package suggested but not available for checking: 'BPCells'

BPCells is hosted on R-universe and used conditionally in SeuratObject

## Downstream dependencies

The following reverse dependency is impacted by this release of SeuratObject:

- tidyseurat
  - test failures: SeuratObject changes the order of the results, but not the actual values. The authors of tidyseurat are aware of this, but the functionality of tidyseurat is not impacted
