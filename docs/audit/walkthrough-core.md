# Newcomer walkthrough: UC-001–024

## Method

This walkthrough starts at each entry in [the use-case inventory](use-cases.md).
It follows links and rendered public docstrings before it reads implementation code.
Each question below has a separate result.

`Static` means that this pass inspected text, source, or tests without executing it.
`Executed` means that this pass ran a small example in the current package environment.
The execution used Julia 1.13.0 and the named `walkthrough_core` daemon.
This pass did not download data, build graphics, or test installation from the network.

The result labels have these meanings:

- **Discoverable:** The public documentation answers the question from the listed start path.
- **DOC GAP:** The implementation or tests answer the question, but the public path does not.
- **API GAP:** The natural operation is absent. The API gap section records it separately.
- **Runtime failure:** A documented operation failed during this pass.
- **Untested:** This pass did not execute the operation.

## A. Start a project and choose a grid

### UC-001 — Install the package and run a first cell query

**Discovery path:** `docs/src/index.md` → repository `README.md` →
`docs/src/tutorials/choosing_a_grid.jl` → `Project.toml`.

- **Which Julia version and installation source work?** The index and README require Julia 1.11 or later. Both install from the GitHub repository. `Project.toml` confirms `julia = "1.11"`. The URL has no tag or revision, so it does not define a reproducible package version. **DOC GAP UCOR-001.**
- **Are local workspace dependencies required?** The public pages do not state this directly. The root project declares all core dependencies. The workspace lists test, docs, benchmark, and library projects. A normal package user should not need those workspace projects. Source inspection was required for this answer. **DOC GAP UCOR-001.**
- **Which imports belong in a minimal example?** The README quick start uses only `import DiscreteGlobalGrids as DGG`. It then constructs a grid and calls `cellat`. **Discoverable.**
- **Which optional packages activate later features?** The index has one table for Zarr, Makie, Metis, KaHyPar, and Scotch. It separately names the visualization package. This validates the correction tracked by [DOC-018](findings.md#doc-018-provide-one-optional-integration-discovery-table). **Discoverable.**
- **Discovered question: does the documented install work in a clean environment?** This pass could not answer it without a network download. [DOC-020](findings.md#doc-020-add-reliable-local-documentation-validation) also remains open. **Untested.**

**Evidence:** Static: `docs/src/index.md:74-88,131-145`, `README.md:24-68`, and `Project.toml:6-53`. Executed: the already active project loaded the package, constructed `HEALPixSystem()` level 2, and located Zürich as `LevelIndex(2, 10)`. That run does not validate installation.

**Outcome:** **DOC GAP** and **Untested**. The first query works in the existing workspace. The clean installation path is not reproducible or tested here.

### UC-002 — Choose a system for geometry, topology, and compatibility

**Discovery path:** `docs/src/tutorials/choosing_a_grid.jl` → `docs/src/all_dggs.md` →
`systems` docstring → `docs/src/internals/system-capabilities.md`.

- **Which systems exist?** The chooser lists IGeo7, H3, A5, HEALPix, ISEA4R, and S2. It lists CopernicusDEM separately. **Discoverable.**
- **Why does `systems()` omit CopernicusDEM?** The chooser and `systems` docstring explain that its latitude-dependent raster lattice does not fit uniform-size registry sweeps. **Discoverable.**
- **Which cells have equal areas?** The gallery and capability table label HEALPix, A5, and ISEA4R as equal-area. The chooser instead calls A5 “nearly equal.” The A5 system docstring resolves the distinction: A5 is equal-area on its ellipsoid, while unit-sphere `cell_area` varies by about 1%. The overview tables do not state the measurement surface and therefore give conflicting guidance. **DOC GAP UCOR-012.**
- **Where do exceptional cell shapes occur?** The chooser identifies the twelve IGeo7 and H3 pentagons. The capability table records exceptional ISEA4R degrees and A5 connectivity. **Discoverable.**
- **Which identifiers match external datasets?** The chooser states H3 compatibility and HEALPix nested/ring conversion. It disclaims native S2 and external ISEA4R compatibility. **Discoverable.**
- **Discovered question: are equal level numbers comparable?** The chooser says no and directs the reader to `cellsize` and `levelfor`. **Discoverable.**

**Evidence:** Static: `docs/src/tutorials/choosing_a_grid.jl:20-102,131-153`, `docs/src/all_dggs.md:54-87`, `src/DiscreteGlobalGrids.jl:243-267`, `src/systems/A5/system.jl:4-14`, and `docs/src/internals/system-capabilities.md:17-49,113-123`. Executed: A5 level 2 had a maximum-to-minimum unit-sphere area ratio of `1.0128715`; HEALPix level 2 had ratio `1.0`. The system list and compatibility guidance validate most of [DOC-002](findings.md#doc-002-provide-one-complete-system-choice-overview), but the area label does not.

**Outcome:** **DOC GAP**. System discovery and compatibility are clear. The equal-area decision is inconsistent for A5.

### UC-003 — Choose a level from a physical resolution

**Discovery path:** `docs/src/tutorials/choosing_a_grid.jl` →
`docs/src/api/grid-interface.md` → rendered `cellsize` and `levelfor` docstrings.

- **Does size mean edge length, diameter, or equivalent-area width?** It is the side of a square whose area equals the sampled median cell area. **Discoverable.**
- **What are its units and sampling assumptions?** The result is metres on the WGS84 authalic-radius sphere by default. The default median uses at most 256 sampled cells. **Discoverable.**
- **How does rounding work?** `levelfor` compares area ratios. It selects the geometrically nearer adjacent level. **Discoverable.**
- **What happens beyond supported levels?** It selects the nearest endpoint of `levels(sys)`. A nonpositive size raises `ArgumentError`. **Discoverable.**
- **How much storage will the result need?** `ncells` gives a cell count, but the task path does not connect that count to array element size, extra dimensions, or cell-axis overhead. **DOC GAP UCOR-002.**
- **Discovered question: does `over` measure both datasets over the same region?** No. It restricts candidate-grid sampling only. [DOC-013](findings.md#doc-013-correct-the-current-levelfor-area-of-interest-contract) now states this limitation. The natural symmetric operation remains [API-004](api-gaps.md#api-004-compare-target-and-candidate-resolution-over-the-same-region). **API GAP.**

**Evidence:** Static: `src/sizing.jl:11-90`, `docs/src/api/grid-interface.md:53-63`, and `test/interface/sizing.jl:31-69`. Executed: HEALPix level 2 returned `1.629905451188185e6` metres; a 100 km request selected level 6.

**Outcome:** **DOC GAP** and **API GAP**. Size selection is discoverable. Storage estimation and symmetric regional matching are not complete workflows.

### UC-004 — Match the resolution of an existing raster or grid

**Discovery path:** `docs/src/tutorials/choosing_a_grid.jl` →
`docs/src/tutorials/between_grids.jl` → `docs/src/api/grid-interface.md` →
`cellsize` and `levelfor` docstrings → `src/sizing.jl`.

- **Which raster dimensions and coordinate systems are accepted?** The docstring says a `DimensionalData.AbstractDimArray` measured through X/Y lookups. It does not state required lookup sampling, longitude/latitude units, projected-coordinate handling, or failure cases. Source reading was required. **DOC GAP UCOR-003.**
- **How is nonuniform resolution summarized?** The implementation takes a median of sampled source cell areas. The docstring gives the sample count but not enough raster-axis requirements to reproduce the sample. **DOC GAP UCOR-003.**
- **Does latitude change the estimate?** Yes. The sizing tests show the change for CopernicusDEM longitude/latitude boxes. The public docstring says location can matter, but the task path has no raster example that demonstrates the effect. **DOC GAP UCOR-003.**
- **Is automatic target selection reproducible?** It is deterministic for the same axes and `samples` value. The source uses a fixed irrational stride. The public docs do not state that sampling rule. **DOC GAP UCOR-003.**
- **Discovered question: can `over` compare the same region on both sides?** No. The documented asymmetry links to [API-004](api-gaps.md#api-004-compare-target-and-candidate-resolution-over-the-same-region). **API GAP.**

**Evidence:** Static: `src/sizing.jl:54-105,115-145`, `test/interface/sizing.jl:42-69`, and `docs/src/api/grid-interface.md:53-63`. The existing user edit in `docs/src/tutorials/between_grids.jl` was only read. No raster example was executed.

**Outcome:** **DOC GAP**, **API GAP**, and **Untested**.

### UC-005 — Align spherical grids with geodetic Earth coordinates

**Discovery path:** `docs/src/tutorials/choosing_a_grid.jl` →
`docs/src/api/grid-interface.md` → `AuthalicSystem` docstring → `docs/src/abstractions.md`.

- **Are input latitudes geodetic or authalic?** Plain spherical systems use their sphere's latitude. `AuthalicSystem` accepts geodetic latitude and maps it to the authalic sphere. **Discoverable.**
- **Which ellipsoid and radius apply?** The wrapper defaults to WGS84. `cellsize` uses the WGS84 authalic radius unless `radius` is supplied. **Discoverable.**
- **Which methods transform coordinates?** The chooser lists `cellat`, `cell_boundary`, `cell_centroid`, `cell_area`, and `cell_cap`. **Discoverable.**
- **Why is wrapping A5 rejected?** A5 already converts its geometry to geodetic latitude. A wrapper would convert twice. **Discoverable.**
- **Does the wrapper change identifiers or neighbors?** No. The chooser states that it preserves IDs, hierarchy, order, and neighbor topology. **Discoverable.**
- **Discovered question: does the system name prove a producer's latitude convention?** No. The chooser explicitly tells users to check producer metadata. **Discoverable.**

**Evidence:** Static: `docs/src/tutorials/choosing_a_grid.jl:167-224`, `docs/src/api/grid-interface.md:128-140`, and `test/fallbacks/authalic.jl`. This validates the task path tracked by [DOC-002](findings.md#doc-002-provide-one-complete-system-choice-overview) and [DOC-007](findings.md#doc-007-show-geographic-coordinates-beside-spherical-geometry-apis). The authalic round trip was not rerun.

**Outcome:** **Discoverable** by static inspection.

## B. Identify cells and inspect geometry

### UC-006 — Locate a point and read its cell

**Discovery path:** `docs/src/api/grid-interface.md` → rendered `cellat` and
`localindex` docstrings → `docs/src/tutorials/zonal.jl`.

- **Is a pair longitude/latitude or latitude/longitude?** It is `(longitude, latitude)`. **Discoverable.**
- **Are angles degrees or radians?** Numeric longitude/latitude calls use degrees. A `UnitSphericalPoint` supplies a point on the unit sphere. **Discoverable.**
- **When must I use `UnitSphericalPoint`?** Use it when the caller already has unit-sphere Cartesian coordinates or needs to avoid geographic conversion. The examples show geographic conversion, but do not state this choice directly. **Discoverable after the docstring.**
- **What happens on seams, poles, or exact boundaries?** The `cellat` docstring promises an incident cell and deterministic platform-local ties. It warns that exact ties can differ across platforms. It does not give each shipped system's tie rule, although it says each system documents one. **DOC GAP UCOR-004.**
- **What happens outside a region?** `cellat` and `localindex` return `nothing`. **Discoverable.**
- **Discovered question: can a tie return `nothing` on a complete grid?** The docstring says no. **Discoverable.**

**Evidence:** Static: `src/interface/grid.jl:288-305`, `docs/src/api/grid-interface.md:18-38`, and `test/systems/crosssystem/cell_vector.jl:188-210`. Executed: a point outside a 16-cell subtree returned `nothing`; Zürich returned a level-2 HEALPix ID.

**Outcome:** **DOC GAP**. Ordinary and absent-region point location worked. Per-system exact-boundary rules remain missing.

### UC-007 — Keep cell identities separate from array positions

**Discovery path:** `docs/src/abstractions.md` → `docs/src/api/grid-interface.md` →
`docs/src/api/neighbor-fields.md` → public identity docstrings → collection tests.

- **What does an `Int` mean in each call?** A bare integer is a local position unless an API explicitly asks for a global index. The neighbor-field table distinguishes `Local()`, `Global()`, and typed ID requests. **Discoverable.**
- **Does a typed identity encode its level and system?** Every typed ID encodes its level. The system can be implicit in a concrete type, but `LevelIndex` is shared by several systems. The public task path does not warn that a `LevelIndex` alone does not identify its system. **DOC GAP UCOR-005.**
- **Is a missing lookup `nothing` or an exception?** `localindex` and `globalindex` return `nothing`. A dimensional `At` selector raises `SelectorError`. The selector difference is not stated beside the identity introduction. **DOC GAP UCOR-009.**
- **Can an indexed handle move to another subset?** The `cellid` docstring says to extract the ID and resolve it in the other collection. A stale subset offset must not move directly. **Discoverable.**
- **Discovered question: can local and global positions differ?** Yes. The grid-interface docstring explains this explicitly. **Discoverable.**

**Evidence:** Static: `docs/src/abstractions.md:20-28`, `src/interface/grid.jl:88-230`, `docs/src/api/neighbor-fields.md:10-36`, and `test/systems/crosssystem/dimensionaldata.jl:166-179`. Executed: in IDs `[1,4,8]`, the second cell had local position 2 and complete-level position 5; an absent ID returned `nothing`.

**Outcome:** **DOC GAP**. The main identity distinction is clear. Shared ID types and selector errors need one user-level explanation. [API-002](api-gaps.md#api-002-construct-an-owning-region-from-ordinary-cell-ids) remains relevant to safe movement into owned regions.

### UC-008 — Import and export native cell identifiers

**Discovery path:** `docs/src/api/grid-interface.md` → rendered `rawid`, `reindex`,
`cellindextype`, and `cellindextypes` docstrings → system module docstrings →
`docs/src/all_dggs.md`.

- **Which identifier types and representations are accepted?** `cellindextypes(sys)` is the authority. `reindex` converts only between listed schemes. **Discoverable.**
- **Are IDs zero-based, sparse, hexadecimal, or ordinal?** The system docstrings define H3 and Z7 `UInt64` encodings and HEALPix nested/ring forms. `rawid` warns that the integer is an encoding, not an array position. **Discoverable.**
- **How are invalid IDs rejected?** Unsupported schemes raise `ArgumentError`. Invalid cells resolve to `nothing` at the index gate where applicable. **Discoverable from public docstrings.**
- **Can external S2 or ISEA4R identifiers be used directly?** No compatibility is claimed. S2 uses package ordinal IDs. ISEA4R needs a fixture-derived permutation before interchange. **Discoverable.**
- **Discovered question: is `rawid` sufficient to reconstruct an ID?** No. It can lose the level and type context. The docstring recommends typed IDs. **Discoverable.**

**Evidence:** Static: `src/interface/grid.jl:88-164`, `docs/src/tutorials/choosing_a_grid.jl:34-43`, and `docs/src/internals/system-capabilities.md:113-123`. This validates the identifier part of [DOC-002](findings.md#doc-002-provide-one-complete-system-choice-overview). Codec fixtures were not rerun.

**Outcome:** **Discoverable** and **Untested** in this pass.

### UC-009 — Measure and export a cell footprint

**Discovery path:** `docs/src/api/grid-interface.md` → rendered geometry docstrings →
`docs/src/abstractions.md` → `docs/src/api/boundaries.md`.

- **What units and coordinate types are returned?** Boundaries and centroids use unit-sphere points. Areas use steradians. Extents use longitude/latitude degrees. **Discoverable.**
- **Are boundaries closed?** They are implicitly closed. The first point is not repeated. `cell_polygon` closes the ring. **Discoverable.**
- **Are edges geodesic or sampled curves?** Consumers join returned vertices with great-circle arcs. Systems must densify boundaries curved in their native chart. **Discoverable from the rendered `cell_boundary` docstring.**
- **Is an extent exact or conservative?** It is conservative for pole and antimeridian cells. Such cells use the full longitude span. **Discoverable.**
- **How do I obtain GeoInterface-compatible geometry?** `getcell` returns a unit-sphere GeoInterface polygon for a local index. `cell_polygon` is documented as internal. Geographic conversion uses `GeometryOps.GeographicFromUnitSphere`. **Discoverable.**
- **Discovered question: can I extract polygons from any public collection?** No. `cell_polygons` supports `MultiOrderCellSet` only. This is [API-003](api-gaps.md#api-003-extract-polygons-consistently-from-public-cell-collections), with the current limit stated by [DOC-012](findings.md#doc-012-narrow-the-collection-polygon-claim-to-supported-inputs). **API GAP.**

**Evidence:** Static: `src/interface/grid.jl:52-87,232-285`, `docs/src/api/grid-interface.md:18-51`, and `docs/src/api/boundaries.md:20-29`. Executed: a HEALPix boundary had 32 points and `first(boundary) != last(boundary)`, consistent with implicit closure. The area was `0.06544984694978735` steradians.

**Outcome:** **Discoverable** for one cell; **API GAP** for general collection extraction.

### UC-010 — Use a tree to accelerate spatial work

**Discovery path:** `docs/src/api/grid-interface.md` → `docs/src/architecture.md` →
`docs/src/tutorials/multiorder.jl` → public `treeify` and cursor docstrings →
`src/engine/cursor.jl`.

- **What does a tree node represent?** A node represents a hierarchy cell and the stored leaf indices beneath it. Sparse subsets can narrow that set. **Discoverable only after internal material.**
- **Is there a synthetic root?** Yes. `node_cell(root)` returns `nothing`, and its extent covers the sphere. The optional multi-order traversal shows this. **Discoverable.**
- **Are caps exact or inflated?** HEALPix, S2, and ISEA4R have exact subtree caps. IGeo7, H3, and A5 use inflated caps. **Discoverable after the capability reference.**
- **Which tree interfaces are external?** `treeify` and `getcell` extend ConservativeRegridding.Trees. Node traversal uses `GeometryOps.SpatialTreeInterface`. The optional example imports the latter. **Discoverable.**
- **When does a query reuse a tree?** The public path does not explain object lifetime or whether callers can supply a prepared tree. Source inspection shows `query` owns traversal and has no tree argument. **DOC GAP UCOR-006.**
- **Discovered question: are node indices local or global?** Cursor nodes own ascending indices in the original grid's local index space. This required source reading. **DOC GAP UCOR-006.**

**Evidence:** Static: `docs/src/api/grid-interface.md:103-112`, `docs/src/tutorials/multiorder.jl:315-372`, `docs/src/internals/system-capabilities.md:37-45`, and `src/engine/cursor.jl:107-121,267-288,449-461`. This task still crosses the boundary tracked by [DOC-003](findings.md#doc-003-separate-task-reference-from-implementation-contracts).

**Outcome:** **DOC GAP** and **Untested**. The optional traversal is usable, but reuse and index-space guidance require source reading.

## C. Build and select regional datasets

### UC-011 — Select cells intersecting a polygon, box, or spherical cap

**Discovery path:** `docs/src/api/selecting-cells.md` → `docs/src/tutorials/zonal.jl` →
`docs/src/tutorials/healpix_astronomy.jl` → rendered query docstrings.

- **Which geometry types and coordinate conventions work?** Query accepts GeoInterface geometries, longitude/latitude extents in degrees, and unit-sphere spherical caps. **Discoverable.**
- **Which predicate is the default?** `covering` means intersection coverage. `query` requires an explicit predicate. **Discoverable.**
- **Does coverage include boundary spill?** Yes. `Intersects` and `Covering` include cells that meet the boundary. The zonal tutorial calls this a rim. **Discoverable.**
- **What are the return type and index space?** Fixed-level `query` returns sorted typed IDs. `covering_indices` returns local collection positions. **Discoverable.**
- **Can the result stay compressed?** `covering` on a `CellVector` returns a compressed subset. A query ID vector can be converted to `CellVector` after satisfying its order contract. **Discoverable.**
- **Discovered question: are cap predicates as broad as polygon predicates?** No. Caps support only `Intersects`, `Disjoint`, and `Within`. **Discoverable.**

**Evidence:** Static: `docs/src/api/selecting-cells.md:7-22,58-96`, `src/engine/query.jl:480-520`, and both tutorials. Executed: an extent around Zürich returned a sorted `Vector{LevelIndex}` with two cells and valid local indices.

**Outcome:** **Discoverable**. The revised result table validates [DOC-001](findings.md#doc-001-state-the-exact-query-result-containers), and predicate direction validates [DOC-004](findings.md#doc-004-make-predicate-direction-and-target-support-explicit).

### UC-012 — Select by an exact spatial relationship

**Discovery path:** `docs/src/api/selecting-cells.md` → rendered `query` docstring →
`docs/src/tutorials/zonal.jl`.

- **Which object is the left operand?** The cell is the left operand. Read `Predicate(target)` as “cell relation target.” **Discoverable.**
- **Which relationships include boundaries?** The predicate table distinguishes `Within`/`Contains` from `CoveredBy`/`Covers`. **Discoverable.**
- **Which results can be empty by dimensionality?** The table defines relations but does not teach the dimensional constraints beyond `Overlaps`. DE9IM users can infer them; newcomers cannot. **DOC GAP UCOR-007.**
- **Does each predicate test the cell polygon or a representative point?** Predicates test the footprint. There is no centroid selector. **Discoverable.**
- **Discovered question: is every exported predicate implemented?** `Crosses` is exported but query rejects it. The table and error contract say so. **Discoverable.**

**Evidence:** Static: `docs/src/api/selecting-cells.md:67-96` and `src/engine/query.jl:376-416,480-498`. The absent centroid operation remains [API-001](api-gaps.md#api-001-select-cells-whose-centroids-lie-in-a-target), while [DOC-011](findings.md#doc-011-identify-the-absent-centroid-selection-operation) now states the limit.

**Outcome:** **DOC GAP** for newcomer dimensional interpretation and **API GAP** for centroid selection. Predicate direction is discoverable.

### UC-013 — Query difficult global geometry

**Discovery path:** `docs/src/api/selecting-cells.md` → query and geometry docstrings →
query tests → `src/engine/query.jl`.

- **What happens across the antimeridian?** The public selection page does not state seam semantics. Multi-order tests cover a seam-crossing polygon. Source and tests were required. **DOC GAP UCOR-007.**
- **Are polar caps supported?** Spherical caps are supported. The docs do not collect pole behavior or cap radius limits in a user guide. **DOC GAP UCOR-007.**
- **Are holes and multipolygons supported?** Cross-system multi-order tests cover both. The public selection page only says “GeoInterface geometry” and gives no supported-shape matrix. **DOC GAP UCOR-007.**
- **Are empty geometries or large spherical polygons supported?** Tests cover a wide polygon. The public path does not state empty-geometry behavior or winding/complement rules. **DOC GAP UCOR-007.**
- **Must a projected polygon be transformed first?** The query path expects geographic degree coordinates or unit-sphere geometry. It does not give a preparation instruction for projected GeoInterface input. **DOC GAP UCOR-007.**
- **What input validity is required?** The public path does not state ring closure, winding, self-intersection, or validity requirements. **DOC GAP UCOR-007.**

**Evidence:** Static: `docs/src/api/selecting-cells.md:7-14,67-96`, `src/engine/query.jl:266-416`, and `test/systems/crosssystem/multiorder_polygons.jl:334-408`. No difficult geometry was rerun.

**Outcome:** **DOC GAP** and **Untested**. Tests show broader support than newcomers can discover.

### UC-014 — Construct a region from an explicit cell list

**Discovery path:** `docs/src/abstractions.md` → `docs/src/api/selecting-cells.md` →
rendered collection docstrings → collection tests.

- **Which container should I choose?** The conversion table distinguishes `PartialGrid`, `CellVector`, `CellLookup`, and `MultiOrderCellSet`. **Discoverable.**
- **Must cells be sorted, unique, and at one level?** Yes for `PartialGrid(sys,l,ids)` and `CellVector(sys,l,ids)`. IDs must be strictly ascending, unique, and at `l`. **Discoverable.**
- **Are duplicates preserved or removed?** Constructors reject duplicates because strict ascending order is required. `vcat` can preserve duplicates, but it is not a constructor normalization step. **Discoverable after public docstrings.**
- **Can a region be empty?** The implementation accepts a correctly typed empty vector. This is not stated in the task path. **DOC GAP UCOR-008.**
- **Does construction reorder associated values?** No. The guide warns that conversions do not sort an unrelated value array. **Discoverable.**
- **Discovered question: does the region own the input ID vector?** `PartialGrid` retains it by reference. Mutation can invalidate the grid. **Discoverable.**

**Evidence:** Static: `docs/src/api/selecting-cells.md:24-48`, `src/engine/partial_grid.jl:68-90`, and `test/systems/crosssystem/cell_vector.jl`. Executed: reversed IDs raised `ArgumentError: ids must be strictly ascending` for both `CellVector` and `PartialGrid`.

**Outcome:** **DOC GAP** for empty input. Safe normalization and ownership remain [API-002](api-gaps.md#api-002-construct-an-owning-region-from-ordinary-cell-ids). General polygon extraction remains [API-003](api-gaps.md#api-003-extract-polygons-consistently-from-public-cell-collections).

### UC-015 — Attach a cell axis to a vector or multidimensional cube

**Discovery path:** `docs/src/tutorials/store_io.jl` →
`docs/src/tutorials/healpix_astronomy.jl` → `docs/src/api/selecting-cells.md` →
rendered `Cells` and `CellLookup` docstrings.

- **How do cell IDs align with values?** Axis position `i` names the value at position `i`. The tutorials use values that encode positions and compare them after IO. **Discoverable.**
- **Can `Cells` occupy any dimension?** The implementation and tests support a non-leading cell dimension. The task path only demonstrates a leading or sole cell dimension. **DOC GAP UCOR-008.**
- **How are time, bands, names, metadata, and missing values preserved?** Regridding and storage pages discuss parts of this contract. The selection entry does not give one preservation statement for ordinary cube construction and slicing. **DOC GAP UCOR-008.**
- **When can I use a plain vector instead?** The HEALPix tutorial shows `query` plus `globalindex` for a complete-level vector. The chooser table distinguishes `CellVector` from `CellLookup`. **Discoverable.**
- **Discovered question: what happens when axis and data lengths differ?** Construction raises `DimensionMismatch`. This is tested but not stated in the task path. **DOC GAP UCOR-008.**

**Evidence:** Static: `docs/src/tutorials/store_io.jl:18-53`, `docs/src/tutorials/healpix_astronomy.jl:176-184`, and `test/systems/crosssystem/dimensionaldata.jl:441`. Executed: a `2 × 3` array accepted `Cells` as its second dimension and preserved the values and lookup.

**Outcome:** **DOC GAP**. Basic alignment is discoverable, but general cube construction needs one explicit contract. [DOC-015](findings.md#doc-015-compare-all-dimensional-array-pass-modes-beside-the-sweep) helps later stencil work, not this construction gap.

### UC-016 — Select and mask a cell-indexed array

**Discovery path:** `docs/src/tutorials/zonal.jl` →
`docs/src/tutorials/healpix_astronomy.jl` → `docs/src/tutorials/store_io.jl` →
`docs/src/api/selecting-cells.md` → dimensional selection tests.

- **How do I avoid the two different `Contains` bindings?** Use `DD.Contains` for the dimensional point or ID selector. Use `DGG.Contains` for the DE9IM cell-footprint predicate. The zonal tutorial states this. **Discoverable.**
- **Does selection preserve a cell lookup?** The storage and HEALPix tutorials show that boolean and region selections retain a `Cells` lookup. **Discoverable.**
- **Does a point select one cell?** `DD.Contains((lon,lat))` selects the one stored cell holding the point. **Discoverable.**
- **What happens when a requested cell is absent?** `At` and point `Contains` raise `SelectorError`; geometry selection can return an empty labeled array. This distinction appears in tests, not the task path. **DOC GAP UCOR-009.**
- **Does a view remain lazy?** The guides describe lazy store values, but they do not state whether each selector returns a view, a lazy wrapper, or an allocated result. **DOC GAP UCOR-009.**
- **Discovered question: what does an empty boolean mask return?** It returns a zero-length array with a zero-length cell lookup. This is tested but not documented. **DOC GAP UCOR-009.**

**Evidence:** Static: `docs/src/tutorials/zonal.jl:127-188`, `docs/src/tutorials/store_io.jl:55-79`, and `test/systems/crosssystem/dimensionaldata.jl:140-269`. Executed: a missing `At` raised `SelectorError`; a false mask returned length zero; a partial mask retained a two-cell `CellLookup`.

**Outcome:** **DOC GAP**. Common selectors work, but cardinality, failure, and laziness need one public contract.

### UC-017 — Combine, intersect, and compare regions

**Discovery path:** `docs/src/api/selecting-cells.md` → rendered `CellVector` docstring →
`docs/src/abstractions.md` → region-algebra tests → `src/engine/region_algebra.jl`.

- **Which inputs must share a system and level?** `intersect` and containment require compatible one-level collections. Mismatched levels or systems raise `ArgumentError`. The user guide does not state this. **DOC GAP UCOR-010.**
- **Does concatenation differ from union?** Yes. `vcat` preserves repeated cells and order. `union` returns unique canonical membership. The public region-algebra section renders only `expand` and `compact`. Source reading was required. **DOC GAP UCOR-010.**
- **What order does the result use?** Set results use canonical order. Concatenation follows input order. **DOC GAP UCOR-010.**
- **Are values combined automatically?** No. These operations act on cell collections. The guide says data aggregation is separate only indirectly. **DOC GAP UCOR-010.**
- **Which set operations are absent?** The task path has no supported-operation list. **DOC GAP UCOR-010.**
- **Discovered question: do incompatible collections return empty results?** No. They raise an error to avoid a misleading empty set. Tests state this, but the guide does not. **DOC GAP UCOR-010.**

**Evidence:** Static: `docs/src/api/selecting-cells.md:128-137`, `src/engine/region_algebra.jl`, and `test/systems/crosssystem/cell_vector.jl:243-277`. Executed: union removed overlap; intersection returned two cells; `vcat` retained both repeated overlap cells.

**Outcome:** **DOC GAP**. [DOC-006](findings.md#doc-006-explain-collection-choices-conversions-and-ownership) improved conversion guidance, but region algebra still requires source reading. [API-002](api-gaps.md#api-002-construct-an-owning-region-from-ordinary-cell-ids) remains open for safe construction.

### UC-018 — Compute zonal statistics from a global field

**Discovery path:** `docs/src/tutorials/zonal.jl` →
`docs/src/api/regridding-methods.md` → `docs/src/tutorials/choosing_a_grid.jl`.

- **Is the statistic based on selected cells or polygon overlap fractions?** The tutorial computes means over selected cells. It warns that this is not an exact polygon statistic. It does not give an overlap-weighted zonal recipe. **DOC GAP UCOR-011.**
- **Should I weight by cell area?** The chooser says equal-area grids permit an ordinary mean and other grids need `cell_area` weights. The zonal tutorial does not apply that rule to its cross-system example. **DOC GAP UCOR-011.**
- **How do coastal missing data affect the result?** The example uses `skipmissing`, but it does not define the denominator or distinguish missing source coverage from selected geometry. **DOC GAP UCOR-011.**
- **How do overlapping zones affect the result?** The tutorial does not discuss independent versus exclusive zones or double counting. **DOC GAP UCOR-011.**
- **Can the same recipe run on another grid?** Yes. The tutorial repeats the coverage mean on IGeo7. **Discoverable.**
- **Discovered question: can I use a center-in-zone rule?** No natural selector exists. This is [API-001](api-gaps.md#api-001-select-cells-whose-centroids-lie-in-a-target). **API GAP.**

**Evidence:** Static: `docs/src/tutorials/zonal.jl:75-175`, `docs/src/tutorials/choosing_a_grid.jl:110-166`, and `docs/src/api/regridding-methods.md`. No downloaded raster or overlap calculation ran.

**Outcome:** **DOC GAP**, **API GAP**, and **Untested**. The selection mean is clear. A complete zonal-statistics decision path is absent.

## D. Use hierarchy and mixed resolutions

### UC-019 — Traverse parents, children, and descendants

**Discovery path:** `docs/src/api/grid-interface.md` → hierarchy docstrings →
`docs/src/tutorials/multiorder.jl` → system capability reference.

- **Is the requested level absolute or relative?** `ancestor` and `descendants` take absolute level numbers. `parent` and `children` move one level. **Discoverable.**
- **Does descendant order match storage order?** It does when `has_sorted_subtrees` holds. `descendant_range` then gives complete-level positions. **Discoverable.**
- **Which systems provide a contiguous range?** All registered systems except A5. **Discoverable.**
- **Do children geometrically tile their parent?** Only HEALPix, S2, and ISEA4R guarantee congruent refinement. **Discoverable.**
- **What happens at level limits?** The hierarchy docstrings constrain ancestor levels and reject children beyond `maxlevel`. **Discoverable.**
- **Discovered question: does `descendants` promise a mutable vector?** No. The revised docs describe an iterable collection without promising allocation. This validates [DOC-009](findings.md#doc-009-describe-descendants-without-promising-allocation-or-mutability). **Discoverable.**

**Evidence:** Static: `docs/src/api/grid-interface.md:65-101`, rendered hierarchy docstrings, and `docs/src/internals/system-capabilities.md:40-45`. Executed: a HEALPix level-2 cell had four level-3 children and sixteen level-4 descendants. Their complete-grid positions matched `descendant_range`.

**Outcome:** **Discoverable**. The usage path validates [DOC-008](findings.md#doc-008-give-hierarchy-operations-a-runnable-usage-path) and [DOC-010](findings.md#doc-010-use-explicit-numeric-level-direction).

### UC-020 — Turn one ancestor into a regional fine grid

**Discovery path:** `docs/src/api/boundaries.md` → `subtree` docstring →
`docs/src/api/grid-interface.md` → `docs/src/tutorials/store_io.jl` → collection guide.

- **Does construction enumerate descendants?** `subtree` returns a `PartialGrid` that can retain a root and compressed membership. It does not need a materialized ID vector on sorted-subtree systems. **Discoverable from docstrings and linked internals.**
- **What does preserving the subtree root enable?** It enables subtree-aware traversal and boundary algorithms. **Discoverable after the boundary guide.**
- **Are returned indices local or complete-level positions?** `cellindex` uses local positions. `globalindex` returns complete-level positions. **Discoverable.**
- **Can the region carry field values?** Convert it to `CellLookup` and wrap it in `Cells`. The storage tutorial demonstrates the axis workflow. **Discoverable.**
- **Discovered question: what is the storage cost?** `CellVector` can store descendant ranges. A5 can require explicit selected indices. The linked collection contracts state this. **Discoverable after the internal reference.**

**Evidence:** Static: `docs/src/api/boundaries.md:96-103`, `docs/src/api/selecting-cells.md:24-48`, and `docs/src/internals/collection-contracts.md:10-74`. Executed: a HEALPix level-2 root produced a 16-cell level-4 subtree and retained `(root_id, root_level)`.

**Outcome:** **Discoverable**. No runtime failure occurred.

### UC-021 — Expand or compact a region across levels

**Discovery path:** `docs/src/api/selecting-cells.md` → rendered `expand` and `compact`
docstrings → `docs/src/tutorials/multiorder.jl` → region-algebra tests.

- **Is level movement geometric or hierarchical?** It follows hierarchy membership. It does not recompute geometric overlap. **Discoverable.**
- **Does compaction require complete sibling groups?** Yes. It replaces complete sibling groups with their parent. **Discoverable.**
- **Which output can contain mixed levels?** `compact` returns `MultiOrderCellSet`. `expand` returns a one-level `CellVector`. **Discoverable.**
- **Can expansion recover the original fine-level membership?** Yes when expansion returns to the original level. This is a membership round trip, not a polygon-union promise. **Discoverable.**
- **Discovered question: can expansion move to a shallower level?** No. It raises `ArgumentError` if the target is shallower than a stored member. **Discoverable from the docstring and tests.**

**Evidence:** Static: `docs/src/api/selecting-cells.md:128-137`, `docs/src/tutorials/multiorder.jl:267-283`, and `test/systems/crosssystem/region_algebra.jl:37-108`. Executed: four HEALPix siblings compacted to their level-2 parent and expanded back to the same level-3 membership.

**Outcome:** **Discoverable**. This validates the relevant parts of [DOC-005](findings.md#doc-005-separate-multi-order-guarantees-from-refinement-budgets) and [DOC-006](findings.md#doc-006-explain-collection-choices-conversions-and-ownership).

### UC-022 — Cover a geometry with cells at several levels

**Discovery path:** `docs/src/tutorials/multiorder.jl` →
`docs/src/api/selecting-cells.md` → multi-order docstrings → system capability reference.

- **How do the two multi-order types differ?** `MultiOrderCoverage` is a query request. `MultiOrderCellSet` is the returned mixed-level collection. **Discoverable.**
- **Which levels can appear?** Fixed-level mode uses cells from system roots through the requested finest level. Budget mode uses roots through the deepest reached level or `maxlevel`. **Discoverable.**
- **How is interior coverage separated from the boundary?** `iscontained(set,i)` marks members proven `Within`; `coarsest_contained` finds the shallowest proven member. Maximum-depth cells can remain unproven. **Discoverable.**
- **Do footprints overlap on noncongruent systems?** The tutorial explains that IGeo7 children can overlap a coarser neighbor at a level change. It distinguishes polygon geometry from represented descendants. **Discoverable.**
- **How do I plot or regrid the result?** The tutorial passes the set to `dggpoly!` and to `regrid(...; to=set)`. **Discoverable.**
- **Discovered question: can every set use `level_ranges`?** No. A5 lacks sorted subtrees. `expand` still works by enumerating descendants. **Discoverable.**

**Evidence:** Static: `docs/src/tutorials/multiorder.jl:30-89,210-303`, `src/engine/multiorder.jl:5-111`, and `src/engine/multiorder_set.jl:23-179`. Executed: a 20-degree HEALPix cap produced 69 members across levels 2 through 4 and 69 polygons.

**Outcome:** **Discoverable**. General collection polygon extraction remains [API-003](api-gaps.md#api-003-extract-polygons-consistently-from-public-cell-collections), but this supported set type is complete.

### UC-023 — Cover a region within a maximum cell budget

**Discovery path:** `docs/src/tutorials/multiorder.jl` → multi-order constructor and
query docstrings → `docs/src/internals/system-capabilities.md` → budget tests.

- **Is the budget strict?** It is strict after the seed. If the initial coarsest seed already exceeds the budget, that seed is returned over budget. **Discoverable.**
- **What happens when root coverage already exceeds it?** The result keeps the seed and exceeds `maxcells`. **Discoverable.**
- **How does refinement prioritize cells?** It refines the coarsest crossing cells first. Ties follow curve order. **Discoverable from the tutorial for the main rule; tie order needs the implementation notes.**
- **How do noncongruent children affect coverage?** IGeo7, H3, and A5 have no budget-mode polygon or leaf-coverage guarantee. **Discoverable.**
- **What does a larger budget improve?** It permits more refinement and usually yields finer boundary representation. The docs correctly avoid a monotonic accuracy guarantee on noncongruent systems. **Discoverable.**
- **Discovered question: can I cap traversal depth?** Yes. `maxlevel` limits budget mode. **Discoverable.**

**Evidence:** Static: `docs/src/tutorials/multiorder.jl:36-47,166-208`, `src/engine/multiorder.jl:79-111`, `src/engine/multiorder_budget.jl:1-43`, and `test/systems/crosssystem/multiorder_budget.jl`. Executed: a HEALPix cap with `maxcells=15,maxlevel=6` returned 15 members across levels 2 through 4.

**Outcome:** **Discoverable**. This validates the budget correction in [DOC-005](findings.md#doc-005-separate-multi-order-guarantees-from-refinement-budgets).

### UC-024 — Test membership and find neighbors across levels

**Discovery path:** `docs/src/api/selecting-cells.md` →
`docs/src/api/neighbors.md` → `member_neighbors`, `iscontained`, and
`coarsest_contained` docstrings → system capability reference.

- **Does membership mean exact presence or ancestor coverage?** Iteration shows stored members. Expansion and `CellLookup` use ancestor coverage at the reference level. `iscontained` does not test membership; it reports proof against the query target by member position. **Discoverable.**
- **Can one coarse member touch several fine members?** Yes. `member_neighbors` returns all adjacent members across the level boundary. **Discoverable.**
- **When is adjacency geometric?** It is geometric for congruent HEALPix, S2, and ISEA4R. It follows the hierarchy relation for noncongruent IGeo7, H3, and A5. **Discoverable.**
- **Why is there no single-level halo for a mixed-level set?** `halo`, `interior`, and ordinary adjacency require one level. A mixed-level set has no single level at which to define that ring. Use `member_neighbors` for member adjacency. **Discoverable after the neighbor page and architecture link.**
- **Discovered question: what happens for a cell outside the set?** `member_neighbors` raises `ArgumentError`. **Discoverable from the rendered docstring.**
- **Discovered question: what does `iscontained` take?** It takes a stored-member position, not a cell ID. The public signature states this, but the similar name makes it easy to mistake for a membership test. **Discoverable, with terminology risk.**

**Evidence:** Static: `docs/src/api/neighbors.md:132-140`, `src/engine/multiorder_set.jl:23-73`, `docs/src/internals/system-capabilities.md:108-111`, and `test/systems/crosssystem/stencils.jl:652-808`. Executed: a compact set containing one level-2 HEALPix member returned no in-set member neighbors, and its containment flag was false because the set came from algebra rather than a coverage proof.

**Outcome:** **Discoverable**. The public docs answer the semantics, though `iscontained` is a proof flag rather than a general membership function.

## New documentation gaps

These identifiers belong to this walkthrough. They are findings for later consolidation.
This pass does not edit [findings.md](findings.md).

| ID | Use cases | Gap | Evidence |
| --- | --- | --- | --- |
| UCOR-001 | UC-001 | State whether the package has a registered release or pin a repository revision. State that normal users do not need workspace projects. Add a clean-install check. | `docs/src/index.md:74-88`; `README.md:24-31`; `Project.toml:6-7` |
| UCOR-002 | UC-003 | Connect `ncells` to a practical storage estimate for element type and extra dimensions. | `docs/src/tutorials/choosing_a_grid.jl:77-109`; `docs/src/api/grid-interface.md:53-63` |
| UCOR-003 | UC-004 | Define accepted raster axis sampling, coordinate units, projection limits, sampling reproducibility, and latitude-dependent examples for `cellsize` and `levelfor`. | `src/sizing.jl:54-105`; `test/interface/sizing.jl:42-69` |
| UCOR-004 | UC-006 | Document the exact-boundary tie rule for each shipped system, or remove the public promise that every system documents it. | `src/interface/grid.jl:288-305` |
| UCOR-005 | UC-007 | Explain that all typed IDs encode level, while shared types such as `LevelIndex` still need system context. | `src/interface/grid.jl:88-106`; `docs/src/architecture.md:47-68` |
| UCOR-006 | UC-010 | Add a task-level tree contract for node index space, prepared-tree reuse, and query ownership. | `docs/src/api/grid-interface.md:103-112`; `src/engine/cursor.jl:107-121,449-461` |
| UCOR-007 | UC-012, UC-013 | Add a supported global-geometry matrix. Cover dimensional predicate outcomes, seam and pole rules, holes, multipart and empty geometries, winding, validity, and projected input preparation. | `docs/src/api/selecting-cells.md:67-96`; `test/systems/crosssystem/multiorder_polygons.jl:334-408` |
| UCOR-008 | UC-014, UC-015 | Document empty collection construction, non-leading `Cells` dimensions, dimension-length errors, and preservation of names, metadata, missing values, and extra dimensions. | `docs/src/api/selecting-cells.md:24-48,139-150`; `test/systems/crosssystem/dimensionaldata.jl:441` |
| UCOR-009 | UC-007, UC-016 | State selector cardinality, absent-cell errors, empty-mask results, lookup preservation, and whether results allocate or stay lazy. | `test/systems/crosssystem/dimensionaldata.jl:140-269`; `docs/src/tutorials/store_io.jl:55-79` |
| UCOR-010 | UC-017 | Document `union`, `intersect`, `vcat`, `issubset`, compatibility errors, order, duplicates, and separation from value aggregation. | `docs/src/api/selecting-cells.md:128-137`; `src/engine/region_algebra.jl`; `test/systems/crosssystem/cell_vector.jl:243-277` |
| UCOR-011 | UC-018 | Add a zonal-statistics decision path for cell selection, area weights, overlap fractions, missing coverage, and overlapping zones. | `docs/src/tutorials/zonal.jl:75-175`; `docs/src/tutorials/choosing_a_grid.jl:110-166` |
| UCOR-012 | UC-002 | Qualify A5 area claims by measurement surface. The overview tables say equal-area, the chooser says nearly equal, and unit-sphere `cell_area` varies about 1%. | `docs/src/all_dggs.md:70-79`; `docs/src/internals/system-capabilities.md:17-24`; `docs/src/tutorials/choosing_a_grid.jl:43-50,131-153`; `src/systems/A5/system.jl:4-14` |

## API gaps confirmed by this walkthrough

This section keeps API outcomes separate from documentation gaps. It does not edit
[api-gaps.md](api-gaps.md), and it does not count a custom loop as completion.

| Existing ID | Use cases in this pass | Walkthrough result |
| --- | --- | --- |
| [API-001](api-gaps.md#api-001-select-cells-whose-centroids-lie-in-a-target) | UC-012, UC-018 | Still open. The docs now state that no centroid-in-target selector exists. |
| [API-002](api-gaps.md#api-002-construct-an-owning-region-from-ordinary-cell-ids) | UC-007, UC-014, UC-017 | Still open. Current constructors require normalized IDs, and `PartialGrid` aliases its input. |
| [API-003](api-gaps.md#api-003-extract-polygons-consistently-from-public-cell-collections) | UC-009, UC-014, UC-022 | Still open. `cell_polygons` accepts `MultiOrderCellSet`, including UC-022, but not the other common collections. |
| [API-004](api-gaps.md#api-004-compare-target-and-candidate-resolution-over-the-same-region) | UC-003, UC-004 | Still open. `over` restricts candidate sampling only. |

No new API gap was confirmed in UC-001–024.

## Outcome count

Counts use the primary outcome of each use case. A use case can also carry an API gap or an untested execution path.

| Primary outcome | Count | Use cases |
| --- | ---: | --- |
| Discoverable | 9 | UC-005, UC-008, UC-011, UC-019, UC-020, UC-021, UC-022, UC-023, UC-024 |
| DOC GAP | 15 | UC-001, UC-002, UC-003, UC-004, UC-006, UC-007, UC-009, UC-010, UC-012, UC-013, UC-014, UC-015, UC-016, UC-017, UC-018 |
| Runtime failure | 0 | None |

Four use cases also expose existing API gaps: UC-003/004, UC-007/014/017,
UC-009/014/022, and UC-012/018. Installation, codec fixtures, difficult global
geometry, raster matching, authalic round trips, tree traversal, and zonal overlap
remain untested in this pass. All executed lightweight examples completed without a
package runtime failure. Two first-attempt audit snippets had Julia parse or caller
construction errors; they were corrected before evidence was recorded and were not
failures of documented package operations.
