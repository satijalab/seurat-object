# SeuratObject v5.4.0

## Test environments
* local Ubuntu 20.04.6 LTS (R-4.3.2)
* local Ubuntu 24.04.1 LTS (R-4.4.2)
* local macOS Sequoia 15.6.1 (R-4.5.1)
* [win-builder](https://win-builder.r-project.org/)
  * R-oldrelease (R-4.4.3 patched)
  * R-release (R-4.5.2 patched)
  * R-devel (R-4.6.0)
* [macos-builder](https://mac.r-project.org/macbuilder/submit.html)
  * R-devel (R-4.6.0)

## R CMD check results

There were no ERRORs or WARNINGs.

There were three NOTEs given by win-builder (oldrel):
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

```
* checking DESCRIPTION meta-information ... NOTE
Author field differs from that derived from Authors@R
```
There seems to be a slight difference in ORCID formatting; the information in both is the same.

## Downstream dependencies

We checked 29 reverse dependencies (12 from CRAN + 17 from Bioconductor), comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 1 package (XYomics) due to a failure during package installation. CRAN checks for XYomics as of 9 April confirm this is unrelated to changes here.

This update does not affect the functionality of any dependencies.