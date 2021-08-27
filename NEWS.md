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
