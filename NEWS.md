# SeuratObject 5.3.0

## Changes:
 - Update `subset.Seurat` to no longer call `droplevels` on the input's cell-level `meta.data` slot; update `subset.Assay` to no longer call `droplevels` on the input's feature-level `meta.features` slot; update `subset.StdAssay` to no longer call `droplevels` on the input's feature-level `meta.data` slot (reverting #251, see discussion in #247)
 - Add `drop = FALSE` when retrieving data in `LayerData.Assay` to preserve dimensions when subsetting to the requested cells (#261)
 - Update `sf.data` slot in the `Segmentation` class to store a `data.frame` in accordance with changes in loading Visium segmentations (#267)
 - Add `compact` slot to the `Segmentation` class to denote whether the object only stores segmentation information in the `sf.data` slot, or also the `sp`-inherited slots (#267)
 - Update `CreateSegmentation`, `CreateSegmentation.data.frame`, `CreateSegmentation.sf`, and methods interacting with `Segmentation` objects (`Cells.Segmentation`, `RenameCells.Segmentation`, `Crop.Segmentation`, `subset.Segmentation`, `[[<-`, `coordinates`, `Overlay`, `show`, `setValidity`) to account for the addition of `compact` and the update to `sf.data`
   - For a detailed description of changes, see [#267](https://github.com/satijalab/seurat-object/pull/267)
 - Add `coords_x_orientation` slot to the `FOV` class to hold the orientation of the x-axis for spatial data (to mark whether the coordinate system of a spatial object has been updated)
 - Update `UpdateSeuratObject` to check whether an object requires an update to the coordinate system based on the existence and value of the `coords_x_orientation` slot (currently only relevant for Visium objects)
 - Add method `safeValidityCheck` to show an improved error message suggesting to run `UpdateSeuratObject` when a "slots in class definition but not in object" error is thrown by R's internal object validation; use in `FOV`-interacting methods

# SeuratObject 5.2.0

## Changes:
- Add `sf.data` slot to the `Segmentation` class to store an [`sf`](https://r-spatial.github.io/sf/articles/sf1.html) object (#258)
  - `sf.data` will represent segmentation boundaries for a given image inside the `images` slot of a Seurat object
  - Add `CreateSegmentation.sf`, `[[<-`, `setValidity` for interacting with `Segmentation` objects 
  - Update `RenameCells.Segmentation`, `subset.Segmentation`, `[`, `UpdateSeuratObject`
- Add optional `misc` slot to `SpatialImage` to store additional info associated with an object in a list (#258)

# SeuratObject 5.1.0

## Changes:
- Update `subset.Seurat` to call `droplevels` on the input's cell-level `meta.data` slot; update `subset.Assay` to call `droplevels` on the input's feature-level `meta.features` slot; update `subset.StdAssay` to call `droplevels` on the input's feature-level `meta.data` slot (#251)
- Update `UpdateSeuratObject` to call `droplevels` on the input's cell-level `meta.data` slot (@samuel-marsh, #247)
- Drop `Seurat` from `Enhances`; update `.IsFutureSeurat` to avoid calling `requireNamespace('Seurat', ...)` (#250)
- Update the `VariableFeatures.StdAssay` setter to apply a speedup (#240)
- Add `SVFInfo.Assay5` & `SpatiallyVariableFeatures.Assay5` (#242)
- Fix bug in `UpdateSeuratObject` (@neanderthalensis, #210)
- Fix bug in `WhichCells.Seurat` (@maxim-h, #219)
- Fix bug in `SpatiallyVariableFeatures.Assay` (#242)
- Fix bug in `merge.Seurat` (#246)
- Fix bug in `VariableFeatures.StdAssay` (#245)
- Fix bug in `HVFInfo.StdAssay` (#244)
- Fix bug in `RenameCells.Seurat` (#237)
- Fix bug in `subset.StdAssay` (#214)

# SeuratObject 5.0.2

## Changes:
- Properly re-export `%||%` from rlang (#178)
- Class key-based warnings (#180)
- Require R 4.1 (#180)
- Fix errors in `UpdateSeuratObject` (@ddiez, #182)
- Add `...` to call signature for `Radius` generic (#190)
- Fix bug in `PolyVtx` (#194)
- Fix bug in feature-level subsetting (#200)
- Update `UpdateSeuratObject` to run without `Seurat` installed (#199)
- Add warning in `Layers.Assay()` when the search returns no results (@maxim-h, #189)
- Fix bug in `subset` to allow empty images to be dropped (#204)

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
