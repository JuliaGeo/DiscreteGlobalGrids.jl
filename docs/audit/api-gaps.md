# API gaps for later implementation

This page describes four operations that need API design or implementation.
Each entry starts with a task a user wants to complete.
**API-001, API-003, API-004, and API-005 are open. API-002 is withdrawn.**
This documentation rewrite does not implement these operations.

Code blocks labeled **Existing syntax** describe current calls.
`DGG` abbreviates `DiscreteGlobalGrids`; input variables stand for objects the user already has.
Blocks labeled **Proposed behavior** describe a missing operation; they are not runnable Julia or accepted API names.

Evidence comes from the historical [core sweep](sweep-core.md) and [workflow sweep](sweep-workflows.md).
Immediate documentation corrections are tracked in [findings.md](findings.md).
Source line numbers in earlier reports are snapshots; locate the symbols before implementation.
A workaround does not close an entry.
Close an entry after its public behavior is decided, implemented, documented, and tested, or record a decision to withdraw it.

Priorities describe impact: P0 means a silent mismatch with the intended behavior, P1 means a major missing operation, and P2 means inconvenient usage.
They do not authorize implementation during this rewrite.

## API-001: Select cells whose centroids lie in a target

**Task:** Calculate a summary for a region using cells whose centers lie inside its outline.
A centroid is the representative center point of a cell.

**Current behavior:** `Intersects` selects cells that touch or overlap the region.
`Within` selects cells wholly inside it.
Both test the cell geometry, so neither expresses the center rule.

**Wanted behavior:** Select cells by testing their centroids against the region.
Make this rule available through both `query` and the `Cells` dimension.
The selector name and syntax still need a design decision.

**Existing syntax**, with `grid`, `field`, and `region` supplied by the user:

```julia
DGG.query(grid, DGG.Intersects(region))       # Cell touches or overlaps region.
field[DGG.Cells(DGG.Within(region))]         # Whole cell lies inside region.
```

**Proposed behavior, not implemented; pseudocode:**

```text
select cells from grid whose centroids lie in region
select the same cells from field through its Cells dimension
```

**Consequence:** Border cells can change the summary.
Users currently need their own centroid filtering to apply this rule.

**Design and acceptance:** Decide whether a centroid on the boundary counts, and specify its coordinate frame.
Test interior, boundary, hole, and exterior cases.
Check that `query` and dimension selection agree, and that results differ from `Intersects` and `Within` where expected.

**Tracking:** P1; SC-006; DOC-011; [UC-012, UC-018, UC-070](use-cases.md).
Evidence: `src/engine/query.jl` (`_matches`) tests cell geometries; `docs/src/tutorials/zonal.jl` describes the missing center rule.
**Status:** Open; design and implementation deferred.

## API-002: Construct an owning region from ordinary cell IDs

**Withdrawn decision:** Do not add an ordinary-user convenience constructor for `PartialGrid`.
Users should not need to construct `PartialGrid` directly.
Keep its constructor requirements in the advanced documentation.

This ID remains for links from earlier reports. It is not an active API gap.
Historical tracking: SC-010; DOC-006; [UC-007, UC-014, UC-017](use-cases.md).

## API-003: Extract polygons consistently from public cell collections

**Task:** Obtain a cell polygon to pass to code that accepts GeoInterface geometry.

**Current behavior:** `cell_boundary(grid, cell_id)` returns boundary vertices.
`cell_polygon(grid, cell_id)` already wraps those vertices as a GeoInterface polygon with a closed ring.
However, its docstring labels it internal.
Users lack a clear public promise for obtaining a polygon across the supported cell representations.

**Wanted behavior:** Make `cell_polygon` a public operation with a consistent GeoInterface polygon result.
Decide which cell representations it accepts and document those inputs explicitly.
Keep `cell_boundary` as the operation that returns vertices.

**Existing syntax**, with a matching `grid` and typed `cell_id` supplied by the user:

```julia
vertices = DGG.cell_boundary(grid, cell_id)  # Boundary vertices.
polygon = DGG.cell_polygon(grid, cell_id)   # Works, but documented as internal.
```

The existing polygon uses unit-sphere `(x, y, z)` coordinates.
The design must make the coordinate frame explicit for each supported input.

**Consequence:** Users who need a polygon must rely on an internal-documented operation or build their own wrapper.
A documented public contract would remove that uncertainty.

**Design and acceptance:** List the supported cell representations and required grid context.
Test GeoInterface polygon behavior, ring closure, and coordinates for each representation.
Broad input coverage does not imply that a collection becomes one polygon.
This entry does not propose expanding the plural `cell_polygons` helper; its collection behavior is a separate question.

**Tracking:** P1; SC-012; DOC-012; [UC-009, UC-014, UC-022, UC-061, UC-069](use-cases.md).
Evidence: `src/interface/grid.jl` (`cell_boundary`, `cell_polygon`) defines the vertex/polygon distinction and labels the polygon wrapper internal.
The historical finding concerned collection documentation; this entry now records the public single-cell polygon requirement.
**Status:** Open; public status and supported input coverage need design and implementation review.

## API-004: Compare target and candidate resolution over the same region

**Task:** Choose a grid level that matches a global longitude/latitude raster near Greenland.
The raster's cells have different physical sizes at different latitudes.
The user needs a match near Greenland, where the analysis takes place.

**Current behavior:** `levelfor(sys, target; over=region)` samples candidate levels of `sys` within `region`.
It measures a spatial `target` over its full extent.
Thus, the two sides of the comparison can describe different places.

**Wanted behavior:** Use `region` when measuring both the candidate grid and the spatial target.
The call syntax already exists; the missing feature is regional measurement of the target.
A numeric target in metres has no spatial variation and needs no change.

**Existing syntax**, with `sys` supplied by the user and `global_raster` standing for an existing global longitude/latitude raster:

```julia
import Extents
greenland = Extents.Extent(X = (-74.0, -10.0), Y = (59.0, 84.0))
level = DGG.levelfor(sys, global_raster; over = greenland)
# Today: candidates near Greenland versus the raster's full extent.
# Wanted: candidates near Greenland versus raster cells in the same region.
```

This example uses a bounding box around Greenland. It does not create or load a raster.

**Consequence:** The selected level can match the raster's overall cell size instead of its cell size near Greenland.
The difference matters for spatial targets whose cell sizes vary by location.

**Design and acceptance:** Define regional sampling for grid, raster, and regridding-space targets.
Decide what happens when the target and region do not overlap.
Test a target whose regional and full-extent median areas differ.
Preserve the default without `over` and the behavior for numeric metre targets.

**Tracking:** P0; SC-015; DOC-013; [UC-003, UC-004, UC-074](use-cases.md).
Evidence: `src/sizing.jl` (`levelfor`, `_targetarea`) passes `over` to candidate sampling but not target measurement.
**Status:** Open; design and implementation deferred.

## API-005: Request radius-k or exact-ring neighborhoods in global sweeps

**Task:** Calculate a value for every cell using neighbors up to three adjacency steps away, or only those exactly three steps away.
An adjacency step moves from a cell to one of its immediate neighbors.

**Current behavior:** Direct `neighbors(grid, cell_id, k)` queries support wider neighborhoods.
Direct `ring(grid, cell_id, k)` queries select the exact ring.
However, `mapneighbors`, `foreachneighbors`, and `mapneighbors!` process only the immediate neighbors of each center cell.
They have no radius or ring selector.

**Wanted behavior:** Let each sweep request either a disk through radius `k` or the exact ring at `k`.
Here, a disk means all neighbors one through `k` steps away; an exact ring means only neighbors `k` steps away.
The callback receives the center separately.
The selector names and syntax remain undecided.

**Existing syntax**, with a `grid`, typed `cell_id`, and cell collection `cells` supplied by the user:

```julia
nearby = DGG.neighbors(grid, cell_id, 3)  # Distances 1, 2, and 3; no center.
outer = DGG.ring(grid, cell_id, 3)        # Distance 3 only; no center.
DGG.mapneighbors((center, nearby) -> length(nearby), cells)
# The sweep above counts only immediate neighbors within cells.
```

**Proposed behavior, not implemented; pseudocode:**

```text
mapneighbors: visit each center with all neighbors at distances 1 through 3
mapneighbors: visit each center with only neighbors at distance 3
```

**Consequence:** A wider data halo does not make the sweep visit more neighbors.
A halo is extra surrounding data loaded for a chunk.
Repeated one-step averaging also gives a different calculation from applying one callback to a radius-three neighborhood.
Users currently need custom loops for the latter task.

**Design and acceptance:** Define zero-radius behavior, center handling, order, and duplicate exclusion.
Test pentagons, grid seams, subset clipping, and all callback forms.
Check that eager and chunked sweeps agree.
Reject chunk plans whose halos cannot cover the requested radius.

**Tracking:** P1; SW-001, SW-002; DOC-014; [UC-025, UC-026, UC-027, UC-030, UC-052, UC-054](use-cases.md).
Evidence: `src/engine/neighborhood.jl` requests `neighbors(..., 1)`; `src/chunks.jl` delegates to that sweep.
Direct radius and ring queries are defined in `src/interface/grid.jl`.
**Status:** Open; design and implementation deferred.

## Future scope candidates, not confirmed API gaps

These tasks need decisions about what the package should provide.
They are not promised APIs and have no new API IDs.
Documentation should describe current support without presenting custom user code as a built-in feature.

| User task | Current support and decision needed |
| --- | --- |
| [UC-069](use-cases.md): Export to common file formats | The package writes its own Zarr layouts and exposes GeoInterface geometry. The audit found no direct writer or tested handoff for GeoJSON, Shapefile, GeoTIFF, NetCDF, FITS, or standard MOC. Decide which formats the package should support. See [UCDOC-018](use-case-findings.md#ucdoc-018--document-tested-exchange-format-handoffs). |
| [UC-070](use-cases.md): Calculate with neighbors within a distance in metres | Users can combine cell selection, centroid distances, and weights. The package promises no physical-radius sweep, Gaussian filter, or physical Laplacian. Decide whether to add these operations. API-005 instead measures distance in adjacency steps. See [UCDOC-028](use-case-findings.md#ucdoc-028--add-a-physical-radius-neighborhood-task-route). |
| [UC-071](use-cases.md): Build a solver that refines cells during a simulation | Mixed-level membership and adjacency exist. A complete solver also needs face orientation, fluxes, value transfer, refinement rules, time stepping, and mixed-level halos. Decide which parts belong in this package. See UCDOC-019 in [use-case-findings.md](use-case-findings.md). |
| [UC-072](use-cases.md): Run calculations on GPUs or clusters | Current documentation covers CPU threads and logical partitions. It promises no GPU kernels, automatic differentiation, or automatic cluster scheduling. Decide which execution features belong here. See UCDOC-020 in [use-case-findings.md](use-case-findings.md). |
| [UC-073](use-cases.md): Append to remote data or coordinate remote writes | Storage supports remote reads and limited writes. It promises no general transaction, retry, conflicting-write, or append protocol. Decide whether this package or the storage backend should provide those guarantees. See UCDOC-021 in [use-case-findings.md](use-case-findings.md). |
| [UC-074](use-cases.md): Use other planetary shapes and projections | Custom authalic ellipsoids and explicit projected-raster transformations exist. The package promises no general planar grid or planetary exchange format. Document supported cases before expanding scope. See UCDOC-003 in [use-case-findings.md](use-case-findings.md). |

The ledger now has four active API gaps and one withdrawn decision. These scope candidates add no API commitments.
