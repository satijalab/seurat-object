# Unreleased

## Changes:
- Properly re-export `%||%` from rlang (#178)

# SeuratObject 5.0.1

## Changes:
- Update internal calls to `GetAssayData()` to use `layer` instead of `slot` (#160)
- Update Matrix version to 1.6-2 (#164)
- Change layer-saving in `SaveSeuratRds()` to move all layers instead of just those in `tempdir()` (#169)
- Update internal calls to `SetAssayData()` to use `layer` instead of `slot` (#171)
- Replace internal calls of `FilterObjects()` to `.FilterObjects()` (#171)

# SeuratObject 5.0.0
## Added
- New `Assay5` class with support for layers; layers provide support for:
  - arbitrary expression matrix names and number
  - arbitrary expression matrix shape
  - disk-backed expression matrices
- New `$` method for `Assay` and `Assay5` objects to pull expression matrices, replacing informal usage of `@`
- New `LayerData()` and `LayerData()<-` functions to replace `GetAssayData()` and `SetAssayData()`, respectively
- Support for renaming cells and features with `dimnames()<-` (changing feature names does not apply to v3 `Assay` objects)
- New `SaveSeuratRds()` and `LoadSeuratRds()` to save and load `Seurat` objects with disk-backed layers
- New `droplevels.LogMap()` to drop unused entries from a `LogMap`
- New ability to split (`split()`) and rejoin layers (`JoinLayers()`) within `Assay` and `Assay5` objects based on grouping factor

## Changes
- `slot` argument deprecated in all contexts; where applicable, replaced with `layer` argument
- `[` for `Assay` and `Assay5` objects take a layer name to pull an expression matrix
  - option `Seurat.object.assay.brackets` allows restoring v3/v4 behavior of subsetting the main expression matrix (eg. `data`)
- Stricter object validation routines at all levels
- `PackageCheck()` deprecated in favor of `rlang::check_installed()`
- `AttachDeps()` deprecated in favor of using the `Depends` field of `DESCRIPTION`
- Subobjects within a `Seurat` object may have subsets of cells present at the object level
- Begun replacement of `stop()` and `warning()` with `rlang::abort()` and `rlang::warn()` for easier debugging
- Expanded validation and utility of `KeyMixin` objects

## Removed
- Unused object constructors (eg. `Assay()`, `Seurat()`)

# SeuratObject 4.1.4
## Changes
- Fixes for `CellsByIdentities` (#80)
- Remove {rgeos} from Suggests and replace with {sf} due to {rgeos} package retirement
- New check for potential binary breaks between dependencies and SeuratObject

# SeuratObject 4.1.3
## Changes
- Move {rgeos} to Suggests; segmentation simplification now requires {rgeos} to be installed manually
- Move {sp} to Depends

## Added
- Add keys to `Assays` and `DimReducs` in `UpdateSeuratObject` when missing

# SeuratObject 4.1.2
## Changed
- Bump required Matrix version to >= 1.5.0

# SeuratObject 4.1.1
## Changed
- Update sparse matrix coersions due to Matrix deprecations

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
