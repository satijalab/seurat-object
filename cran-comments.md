# SeuratObject v5.0.0

## Test environments
* local ubuntu 22.04 install, R 4.2.3
* win-builder (oldrelease, release, devel)

## R CMD check results

There were no ERRORs, WARNINGs, or NOTEs

## Downstream dependencies

The following reverse dependencies are impacted by this release of SeuratObject:

- Seurat:
  - new test failures and warnings: Seurat's tests expect errors that are now handled by SeuratObject. Warnings occur due to use of deprecated, but still accepted arguments. The authors of Seurat are aware of these changes, but the functionality of Seurat is not impacted
  - S3 generic/method consistency: Seurat defines a method for a generic defined in SeuratObject. In the latest update, SeuratObject changes one of the parameters in the method signature, but still accepts the old arguments. The functionality of Seurat is not impacted by this update

- Signac
  - new test failures: SeuratObject changes the order of the results, but not the actual values. The authors of Signac are aware of this, but the functionality of Signac is not impacted

- tidyseurat
  - new test failures: SeuratObject changes the order of the results, but not the actual values. The authors of tidyseurat are aware of this, but the functionality of tidyseurat is not impacted
