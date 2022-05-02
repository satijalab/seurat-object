# SeuratObject 4.1.0
## Changed
- Allow `UpdateSeuratObject` to work when `data` is `NULL` (#38)
- Fix superclass issue with R-devel 4.3.x (#42)

## Added
- New `FOV`, `Segmentations`, `Centroids`, and `Molecules` classes for imaging-based spatial datasets

# SeuratObject 4.0.4
## Changed
- `CreateSeuratObject.Assay` sets Assay key when not present (#29)
- Ignore warnings when creating an `Assay` from a data frame (#32)

## Added
- New `CheckMatrix` generic for validating expression matrices

# SeuratObject 4.0.3

## Changed
- Export utility functions (#22)
- Bug fix in names with `Key.Seurat` (#26)
- Improved duplicate key checking and resolution

# SeuratObject 4.0.2

## Changed
- Provide default option for `Seurat.checkdots` option if option is not set (#16)

# SeuratObject 4.0.1

## Added
- `head` and `tail` methods for `Seurat` and `Assay` objects (#5)
- New utility functions (#6):
  - `AttachDeps` to attach required imported dependencies on package attachment
  - `IsMatrixEmpty` to test if a matrix is empty or not

## Changed
- Allow super classes to replace child classes (#1). For example, allows `Assay`
  objects to replace `Seurat::SCTAssay` or `Signac::ChromatinAssay` objects of
  the same name
- Better support for creating sparse matrices from `data.table`/`tibble`
  objects (#4)
- Improved error messages for clashing object names (#7)
- Allow returning a `NULL` if a subset results in zero cells (#9)

## Removed
- SCT-specific code (#2)

# SeuratObject 4.0.0

- Initial release of SeuratObject
