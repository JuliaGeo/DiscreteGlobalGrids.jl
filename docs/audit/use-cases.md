# Newcomer use-case inventory

This inventory defines tasks for a documentation audit of DiscreteGlobalGrids.jl.
It records questions that a newcomer must answer before a task can succeed.
It does not certify the current documentation or report completed execution tests.

## How to use this inventory

- Keep each `UC-###` identifier when assigning work or recording findings.
- Begin at the listed documentation entry. Follow its links before reading implementation code.
- Record each unanswered question, misleading answer, missing prerequisite, and failed example against its use-case identifier.
- Record the exact input, observed result, expected result, and source evidence for a failure.
- Run a small synthetic example before using downloaded data or a global high-resolution grid.
- A successful outcome includes correct cell identity, value alignment, units, and boundary semantics.

**Status vocabulary:** `Supported` means a local API exists. `Composed` means callers combine APIs or external packages.
`Uncertain` means the audit must establish the supported boundary. `Unsupported` means the inspected surface states a limitation or lacks the requested operation.
Missing natural API support belongs in `docs/audit/api-gaps.md` for later implementation.
Do not substitute a custom loop or workaround and mark the requested API workflow complete.
These labels describe the task, not the quality of its documentation.

Documentation paths below are relative to `docs/src/` unless they begin with `src/`, `test/`, or `lib/`.
Tutorial `.jl` files are the editable Literate sources. Their published pages use `.md`.
Questions separated by semicolons are separate audit checks.

## A. Start a project and choose a grid

### UC-001 — Install the package and run a first cell query

- **Status / entry points:** Supported; package loading, `systems`, `levelgrid`, `cellat`.
- **Questions:** Which Julia version and installation source work? Are local workspace dependencies required? Which imports belong in a minimal example? Which optional packages activate later features?
- **Success:** A clean environment constructs a small grid and locates a known point without an undeclared dependency.
- **Start reading:** `index.md`, repository `README.md`, `tutorials/choosing_a_grid.jl`.
- **Validate:** Run the published setup in a temporary environment. Compare package metadata with documented compatibility requirements.

### UC-002 — Choose a system for geometry, topology, and compatibility

- **Status / entry points:** Supported; `systems`, seven system constructors, capability traits.
- **Questions:** Which systems exist? Why does `systems()` omit CopernicusDEM? Which cells have equal areas? Where do exceptional cell shapes occur? Which identifiers match external datasets?
- **Success:** A reader chooses a system from explicit requirements and understands its main exceptions.
- **Start reading:** `tutorials/choosing_a_grid.jl`, `all_dggs.md`, `systems` docstring.
- **Validate:** Compare the overview with constructors, system traits, and system tests. Check all systems, including S2 and CopernicusDEM.

### UC-003 — Choose a level from a physical resolution

- **Status / entry points:** Supported; `cellsize`, `levelfor`, `levels`, `maxlevel`, `ncells`.
- **Questions:** Does size mean edge length, diameter, or equivalent-area width? What are its units and sampling assumptions? How does rounding work? What happens beyond supported levels? How much storage will the result need?
- **Success:** A requested size selects a valid level with an explained approximation and cell count.
- **Start reading:** `tutorials/choosing_a_grid.jl`, `api/grid-interface.md`.
- **Validate:** Check sizes between adjacent levels, endpoint levels, invalid sizes, and unequal-area systems against `test/interface/sizing.jl`.

### UC-004 — Match the resolution of an existing raster or grid

- **Status / entry points:** Supported; `levelfor`, `cellsize`, bare-system regridding targets.
- **Questions:** Which raster dimensions and coordinate systems are accepted? How is nonuniform resolution summarized? Does latitude change the estimate? Is automatic target selection reproducible?
- **Success:** A reader can explain the level selected from a source and override it when necessary.
- **Start reading:** `tutorials/choosing_a_grid.jl`, `tutorials/between_grids.jl`, `api/grid-interface.md`.
- **Validate:** Compare explicit and inferred levels for a small raster, a complete grid, and a regional source.

### UC-005 — Align spherical grids with geodetic Earth coordinates

- **Status / entry points:** Supported; `AuthalicSystem`, `AuthalicGrid`, `cellat`, geometry methods.
- **Questions:** Are input latitudes geodetic or authalic? Which ellipsoid and radius apply? Which methods transform coordinates? Why is wrapping A5 rejected? Does the wrapper change identifiers or neighbors?
- **Success:** Point location, plotted boundaries, and area calculations use one documented coordinate convention.
- **Start reading:** `tutorials/choosing_a_grid.jl`, `api/grid-interface.md`, `abstractions.md`.
- **Validate:** Use a non-equatorial point. Check coordinate round trips and unchanged identity against `test/fallbacks/authalic.jl`.

## B. Identify cells and inspect geometry

### UC-006 — Locate a point and read its cell

- **Status / entry points:** Supported; `cellat`, `cellindex`, `getcell`, `cell_centroid`.
- **Questions:** Is a pair longitude/latitude or latitude/longitude? Are angles degrees or radians? When must I use `UnitSphericalPoint`? What happens on seams, poles, exact boundaries, or outside a region?
- **Success:** The returned cell contains the intended point, and an absent regional cell has a documented result.
- **Start reading:** `api/grid-interface.md`, `tutorials/zonal.jl`.
- **Validate:** Check an ordinary point, an antimeridian point, a pole, and a point outside a subset.

### UC-007 — Keep cell identities separate from array positions

- **Status / entry points:** Supported; `AbstractCellIndex`, `LevelIndex`, `localindex`, `globalindex`, `cellindex`, `rawid`, `cellid`.
- **Questions:** What does an `Int` mean in each call? Does a typed identity encode its level and system? Is a missing lookup `nothing` or an exception? Can an indexed handle move to another subset?
- **Success:** Values remain attached to the same cells after selecting, sorting, and moving between collections.
- **Start reading:** `abstractions.md`, `api/grid-interface.md`, `api/neighbor-fields.md`.
- **Validate:** Use a noncontiguous subset whose local positions differ from global positions. Resolve a handle through `cellid` on another axis.

### UC-008 — Import and export native cell identifiers

- **Status / entry points:** Supported with system-specific limits; `rawid`, `reindex`, `cellindextype`, `cellindextypes`, H3/A5/Z7 codecs.
- **Questions:** Which identifier types and representations are accepted? Are IDs zero-based, sparse, hexadecimal, or ordinal? How are invalid IDs rejected? Can external S2 or ISEA4R identifiers be used directly?
- **Success:** A supported external ID round-trips without changing its geographic cell. Unsupported codecs are identified explicitly.
- **Start reading:** `api/grid-interface.md`, system module docstrings, `all_dggs.md`.
- **Validate:** Use committed system fixtures. Include invalid IDs; do not infer native S2 compatibility or external ISEA4R numbering.

### UC-009 — Measure and export a cell footprint

- **Status / entry points:** Supported; `cell_boundary`, `cell_centroid`, `cell_area`, `cell_extent`, `getcell`; qualified geometry helpers.
- **Questions:** What units and coordinate types are returned? Are boundaries closed? Are edges geodesic or sampled curves? Is an extent exact or conservative? How do I obtain GeoInterface-compatible geometry?
- **Success:** A polygon and its measurements have documented units and suitable coordinates for the next consumer.
- **Start reading:** `api/grid-interface.md`, `abstractions.md`.
- **Validate:** Compare area totals on a small complete level. Inspect curved edges, seam cells, and an authalic wrapper.

### UC-010 — Use a tree to accelerate spatial work

- **Status / entry points:** Supported; `treeify`, `node_extent`, `HierarchicalGridCursor`, qualified `node_cell`.
- **Questions:** What does a tree node represent? Is there a synthetic root? Are caps exact or inflated? Which tree interfaces are external? When does a query reuse a tree?
- **Success:** A traversal visits the intended cells and uses bounds without treating them as exact cell polygons.
- **Start reading:** `api/grid-interface.md`, `architecture.md`, `tutorials/multiorder.jl`.
- **Validate:** Compare tree-assisted selection with exhaustive selection on a tiny grid. Check root and leaf semantics.

## C. Build and select regional datasets

### UC-011 — Select cells intersecting a polygon, box, or spherical cap

- **Status / entry points:** Supported; `query`, `covering`, `covering_indices`, `Covering`, DE9IM predicates.
- **Questions:** Which geometry types and coordinate conventions work? Which predicate is the default? Does coverage include boundary spill? What are the return type and index space? Can the result stay compressed?
- **Success:** A reproducible selection uses an explicit geometric rule and indexes the correct values.
- **Start reading:** `api/selecting-cells.md`, `tutorials/zonal.jl`, `tutorials/healpix_astronomy.jl`.
- **Validate:** Compare a cap, an extent, and a polygon against a small exhaustive geometry check.

### UC-012 — Select by an exact spatial relationship

- **Status / entry points:** Supported with predicate-specific semantics; `Intersects`, `Within`, `Contains`, `Covers`, `CoveredBy`, `Touches`, `Crosses`, `Overlaps`, `Equals`, `Disjoint`.
- **Questions:** Which object is the left operand? Which relationships include boundaries? Which results can be empty by dimensionality? Does each predicate test the cell polygon or a representative point?
- **Success:** The reader can distinguish full containment, contact, and intersection without guessing from predicate names.
- **Start reading:** `api/selecting-cells.md`, `query` docstring.
- **Validate:** Use one geometry with interior, crossing, exterior, and boundary-touching cells. Check each documented predicate.

### UC-013 — Query difficult global geometry

- **Status / entry points:** Uncertain at individual geometry boundaries; `query`, `covering`, supported geometry adapters.
- **Questions:** What happens across the antimeridian? Are polar caps, holes, multipolygons, empty geometries, or large spherical polygons supported? Must a projected polygon be transformed first? What input validity is required?
- **Success:** Supported geometry cases return correct cells; unsupported cases have a clear preparation step or limitation.
- **Start reading:** `api/selecting-cells.md`, geometry/query docstrings.
- **Validate:** Use synthetic seam and polar geometries, one hole, and one multipart shape. Compare exact predicates with conservative coverage.

### UC-014 — Construct a region from an explicit cell list

- **Status / entry points:** Supported; `PartialGrid`, `CellVector`, `CellLookup`, `region`, `cellset`.
- **Questions:** Which container should I choose? Must cells be sorted, unique, and at one level? Are duplicates preserved or removed? Can a region be empty? Does construction reorder associated values?
- **Success:** The collection represents the intended membership with documented ordering and data alignment.
- **Start reading:** `abstractions.md`, `api/selecting-cells.md`.
- **Validate:** Try unsorted, repeated, empty, and mixed-level inputs. Compare accepted behavior with cell-vector and dimensional-data tests.

### UC-015 — Attach a cell axis to a vector or multidimensional cube

- **Status / entry points:** Supported; `Cells`, `CellLookup`, `DimArray`, Raster integration.
- **Questions:** How do cell IDs align with the values? Can `Cells` occupy any dimension? How are time, bands, names, metadata, and missing values preserved? When can I use a plain vector instead?
- **Success:** A cell-time cube supports spatial selection without mixing cell and time positions.
- **Start reading:** `tutorials/store_io.jl`, `tutorials/healpix_astronomy.jl`, `api/selecting-cells.md`.
- **Validate:** Use a non-leading cell dimension and values that encode their cell/time coordinates.

### UC-016 — Select and mask a cell-indexed array

- **Status / entry points:** Supported; `Cells(At(...))`, dimensional `Contains`, `Covering`, predicate selectors, boolean indexing.
- **Questions:** How do I avoid the two different `Contains` bindings? Does selection preserve a cell lookup? Does a point select one cell? What happens when a requested cell is absent? Does a view remain lazy?
- **Success:** Point, identity, geometry, and boolean selectors return correctly labeled data.
- **Start reading:** `tutorials/zonal.jl`, `tutorials/healpix_astronomy.jl`, `tutorials/store_io.jl`.
- **Validate:** Compare selector results with explicit cell positions. Check an empty mask and a missing member.

### UC-017 — Combine, intersect, and compare regions

- **Status / entry points:** Supported; `union`, `intersect`, `vcat`, `issubset`, `cellset`.
- **Questions:** Which inputs must share a system and level? Does concatenation differ from union? What order does the result use? Are values combined automatically? Which set operations are absent?
- **Success:** Membership algebra has predictable semantics, and data aggregation remains an explicit separate step.
- **Start reading:** `api/selecting-cells.md`, region-algebra docstrings.
- **Validate:** Compare with small explicit sets, including overlapping regions, incompatible levels, and empty collections.

### UC-018 — Compute zonal statistics from a global field

- **Status / entry points:** Composed; geometry selection, reductions, `cell_area`, regridding to geometry targets where supported.
- **Questions:** Is the statistic based on selected cells or polygon overlap fractions? Should I weight by cell area? How do coastal missing data and overlapping zones affect the result? Can the same recipe run on another grid?
- **Success:** The statistic states its coverage and weighting rule and matches a small independent calculation.
- **Start reading:** `tutorials/zonal.jl`, `api/regridding-methods.md`.
- **Validate:** Use a constant field and a two-valued field over a partial polygon. Compare selection means with overlap-weighted results.

## D. Use hierarchy and mixed resolutions

### UC-019 — Traverse parents, children, and descendants

- **Status / entry points:** Supported; `parent`, `children`, `ancestor`, `descendants`, `rootcells`, `descendant_range`.
- **Questions:** Is the requested level absolute or relative? Does descendant order match storage order? Which systems provide a contiguous range? Do children geometrically tile their parent? What happens at level limits?
- **Success:** The reader moves between levels without confusing hierarchical ancestry with geometric coverage.
- **Start reading:** `api/grid-interface.md`, `tutorials/multiorder.jl`.
- **Validate:** Check ancestor round trips and descendant counts on congruent and noncongruent systems, including A5's range limitation.

### UC-020 — Turn one ancestor into a regional fine grid

- **Status / entry points:** Supported; `subtree`, `region`, `CellVector`, `PartialGrid`.
- **Questions:** Does construction enumerate descendants? What does preserving the subtree root enable? Are returned indices local or complete-level positions? Can the region carry field values?
- **Success:** A subtree supports queries, data axes, and boundary operations with documented storage costs.
- **Start reading:** `api/boundaries.md`, `api/grid-interface.md`, `tutorials/store_io.jl`.
- **Validate:** Compare subtree membership with `descendants` at a shallow depth. Check root-preserving and explicit-list forms.

### UC-021 — Expand or compact a region across levels

- **Status / entry points:** Supported; `expand`, `compact`, `cellset`, `cellindices`.
- **Questions:** Is level movement geometric or hierarchical? Does compaction require complete sibling groups? Which output can contain mixed levels? Can expansion recover the original fine-level membership?
- **Success:** The reader understands what information compaction preserves and what it does not preserve geometrically.
- **Start reading:** `api/selecting-cells.md`, `tutorials/multiorder.jl`.
- **Validate:** Test full and partial sibling groups. Check compact/expand membership round trips on multiple systems.

### UC-022 — Cover a geometry with cells at several levels

- **Status / entry points:** Supported; `MultiOrderCoverage`, `MultiOrderCellSet`, `level_ranges`, `cell_polygons`.
- **Questions:** How do the two multi-order types differ? Which levels can appear? How is interior coverage separated from the boundary? Do footprints overlap on noncongruent systems? How do I plot or regrid the result?
- **Success:** A compact coverage has an explicit finest level and documented geometric interpretation.
- **Start reading:** `tutorials/multiorder.jl`, `api/selecting-cells.md`.
- **Validate:** Expand to the finest level and compare with reference coverage. Include HEALPix and IGeo7.

### UC-023 — Cover a region within a maximum cell budget

- **Status / entry points:** Supported; `MultiOrderCellSet(...; maxcells=...)`.
- **Questions:** Is the budget strict? What happens when root coverage already exceeds it? How does refinement prioritize cells? How do noncongruent children affect coverage? What does a larger budget improve?
- **Success:** The result respects the documented budget contract and exposes the coverage-versus-resolution tradeoff.
- **Start reading:** `tutorials/multiorder.jl`, multi-order constructor docstrings.
- **Validate:** Check tiny and moderate budgets against `test/systems/crosssystem/multiorder_budget.jl`.

### UC-024 — Test membership and find neighbors across levels

- **Status / entry points:** Supported with hierarchy-specific semantics; `iscontained`, `coarsest_contained`, `member_neighbors`.
- **Questions:** Does membership mean exact presence or ancestor coverage? Can one coarse member touch several fine members? When is adjacency geometric? Why is there no single-level halo for a mixed-level set?
- **Success:** Cross-level answers distinguish stored members from cells represented through ancestry.
- **Start reading:** `api/selecting-cells.md`, `api/neighbors.md`, `systems` docstring.
- **Validate:** Construct a coarse member beside refined members. Compare congruent and noncongruent systems.

## E. Compute with neighborhoods and boundaries

### UC-025 — Find first-order neighbors or an exact higher-order ring

- **Status / entry points:** Supported; `neighbors`, `ring`, `neighborcount`, `Vertex`, `Edge`.
- **Questions:** Does `neighbors(..., k)` include every distance through `k`? Is the center included? How does `ring` differ? What does order zero mean? Is distance topological or metric? How are duplicates handled?
- **Success:** The reader chooses a ring or a neighborhood with the intended membership.
- **Start reading:** `api/neighbors.md`, `tutorials/stencils.jl`.
- **Validate:** Compare explicit graph distances on a tiny complete grid. Include zero, one, two, and exceptional cells.

### UC-026 — Convolve a global dataset over second-order neighbors

- **Status / entry points:** API gap; `neighbors(..., 2)` and `ring(..., 2)` select cells, but `mapneighbors` supplies only one-ring neighborhoods.
- **Questions:** Does the task need the radius-two disk, the exact second ring, or two successive diffusion steps? Is the center included? How should weights depend on distance? How do pentagons and seams change kernel width? What natural radius/selection API is missing from the sweep?
- **Success:** The audit records the requested operation precisely and identifies the missing sweep support in `docs/audit/api-gaps.md`. A larger halo is not documented as a substitute for a larger callback neighborhood.
- **Start reading:** `tutorials/stencils.jl`, `api/neighbors.md`, sweep docstrings.
- **Validate:** Check signatures and callback membership against an explicit small-grid reference. Distinguish a radius-two average from repeated one-ring averaging; they generally use different weights.

### UC-027 — Apply a first-order stencil to a field

- **Status / entry points:** Supported; `mapneighbors`, `foreachneighbors`, `Values`, `NeighborSlices`.
- **Questions:** What arguments reach the callback? Where is the center value? What type does the callback return? Can the input be a plain vector? Does output preserve dimensions? What happens for an isolated cell?
- **Success:** A smoothing or difference kernel runs on a vector and a dimensional array with matching results.
- **Start reading:** `tutorials/stencils.jl`, `api/neighbors.md`.
- **Validate:** Compare a hand-computed cell, a constant field, and a one-cell region with `mapneighbors` tests.

### UC-028 — Read values, centroids, and identifiers in one stencil

- **Status / entry points:** Supported; `needs`, qualified `Cell`, `Index`, `Local`, `Global`, `Value`, `Centroid`, `cellfield`.
- **Questions:** Are callback records organized by field or by neighbor? Which indices address which arrays? How do several quantities stay aligned? Can a quantity be computed lazily? Can a partial field supply only needed cells?
- **Success:** A geometry-aware kernel reads consistent values and coordinates without relying on undocumented tuple positions.
- **Start reading:** `api/neighbor-fields.md`, `tutorials/stencils.jl`.
- **Validate:** Compare computed centroids with precomputed fields. Use distinct local/global indices and a partial value source.

### UC-029 — Use neighbor order for a directional kernel

- **Status / entry points:** Supported with system-specific constraints; `winding`, adjacency slots, `directioncode`, `StorageOrder`.
- **Questions:** Is order counterclockwise? Where does the first slot begin? Does clipping preserve slot numbers? Is a slot a compass direction? Which IGeo7 direction codes are portable to other cells or systems?
- **Success:** Directional output relies only on a documented ordering contract and handles missing directions explicitly.
- **Start reading:** `api/neighbors.md`, `api/neighbor-fields.md`, IGeo7 docstrings.
- **Validate:** Check clockwise geometry, exceptional cells, and clipped versus marked adjacency rows. Do not assume six fixed slots.

### UC-030 — Run a stencil on a region with holes

- **Status / entry points:** Supported; regional `neighbors`, `ring`, `mapneighbors`, `grow`.
- **Questions:** Are neighbors clipped to membership? Does a hole change graph distance or only remove results? How do I provide outside values? Does a cell outside the region produce an error? What boundary condition does my kernel imply?
- **Success:** The kernel has an explicit treatment of absent neighbors and does not confuse clipped system distance with paths around obstacles.
- **Start reading:** `tutorials/stencils.jl`, `api/boundaries.md`, `api/neighbors.md`.
- **Validate:** Remove one interior cell from a small region. Compare with complete-level neighbors filtered by membership.

### UC-031 — Cache adjacency for repeated graph calculations

- **Status / entry points:** Supported; `adjacency`, `AdjacencyTable`, `halocells`, `haloindices`.
- **Questions:** What is the row and column index space? How do `halo=0`, `halo=1`, and `halo=:mark` differ? What does zero mark? Why are wider halos rejected? Can the table be reused with another ordering?
- **Success:** A caller indexes a clipped or completed neighbor table without mixing halo positions and cell identities.
- **Start reading:** `api/neighbors.md`, `tutorials/stencils.jl`.
- **Validate:** Compare every row with direct neighbors. Check the `[region; halo]` buffer and missing-slot markers.

### UC-032 — Find the border, interior, and outside halo of a region

- **Status / entry points:** Supported; `border`, `interior`, `halo`, `halo_indices`.
- **Questions:** Which side of the boundary does each verb return? Are holes included in the outside halo? Is iteration lazy and sorted? Is `length` available? What changes across subtree and arbitrary-subset inputs?
- **Success:** Border and interior partition the region; its halo consists of external adjacent cells in the documented index space.
- **Start reading:** `api/boundaries.md`, `internals/boundary-engines.md`.
- **Validate:** Compare with explicit membership on a subtree, a region with a hole, an empty region, and a complete globe.

### UC-033 — Grow a region for a multi-step computation

- **Status / entry points:** Supported; `grow`, `halo`, region algebra.
- **Questions:** Does growth return the original region plus new cells? Which connectivity applies? Does repeated growth match a wider neighborhood? How do I map original output positions into the grown input?
- **Success:** A kernel receives enough surrounding data while output remains associated with the original owned cells.
- **Start reading:** `api/boundaries.md`, `api/selecting-cells.md`, `api/chunk-sweep.md`.
- **Validate:** Compare repeated one-step growth with an explicit breadth-first reference and check value remapping.

### UC-034 — Choose traversal order and handle callback failures

- **Status / entry points:** Supported; `StorageOrder`, sweep ordering/thread options, `NeighborCallbackError`.
- **Questions:** Which buffers can callbacks retain? Can callbacks mutate input or shared state? Which iteration order is guaranteed? How are exceptions attributed to a cell? How does threading affect determinism?
- **Success:** A documented callback contract supports correct serial and parallel execution, including useful error reporting.
- **Start reading:** `api/neighbor-fields.md`, sweep docstrings.
- **Validate:** Compare serial and threaded pure kernels. Trigger one controlled callback exception and inspect its location information.

### UC-035 — Compute a cost-distance field through neighboring cells

- **Status / entry points:** Composed; neighborhood APIs plus the tutorial's priority-queue algorithm.
- **Questions:** Which connectivity defines a route? How are edge distances and costs computed? Are missing cells barriers? Does the shortest path use the subset graph? Is the algorithm a package API or example code?
- **Success:** A seed produces a reproducible distance field with documented unreachable-cell behavior.
- **Start reading:** `tutorials/stencils.jl`.
- **Validate:** Use uniform costs, an obstacle, and a disconnected region. Compare with an independent tiny graph.

## F. Transfer data between spatial grids

### UC-036 — Regrid a raster onto a DGGS

- **Status / entry points:** Supported; `regrid`, a grid target, `DGGSpace`.
- **Questions:** Which raster types, dimensions, and coordinate systems work? How are input footprints defined? Does the target cover the globe or source region? Which result wrapper and cell axis are returned?
- **Success:** A small raster becomes a labeled cell field with the intended coverage and interpretation.
- **Start reading:** `tutorials/regridding.jl`, `api/regridding-methods.md`.
- **Validate:** Regrid a constant and a spatially varying synthetic raster. Check destination cells and metadata.

### UC-037 — Move a field between DGGS systems or levels

- **Status / entry points:** Supported; `regrid`, explicit `from`, grid/lookup/collection targets.
- **Questions:** When can the source space be inferred? When is `from` required? Can a regional source feed a complete destination? How do finer and coarser targets differ? Does an equal level number imply equal resolution?
- **Success:** Data moves between two systems with correctly specified spaces and comparable physical resolution.
- **Start reading:** `tutorials/between_grids.jl`, `src/regridding.jl` docstrings.
- **Validate:** Compare explicit and inferred source spaces. Use unequal system levels with similar cell sizes.

### UC-038 — Choose area conservation, interpolation, or containing-cell transfer

- **Status / entry points:** Supported; `Conservative`, `BarycentricPoint`, `NearestCell`, `DirectNearest`.
- **Questions:** Are values area averages, point samples, categories, or totals? Which methods preserve integrals? Does nearest mean nearest centroid or containing cell? When is a reusable operator useful? Which methods retain source values exactly?
- **Success:** The method matches the data meaning, and the example explains the result at coarse-to-fine resolution.
- **Start reading:** `api/regridding-methods.md`, `tutorials/between_grids.jl`.
- **Validate:** Test constant preservation, source-site reproduction, categorical transfer, and a conservative area integral.

### UC-039 — Control missing values and partial source coverage

- **Status / entry points:** Supported; `Weighted`, `Extensive`, `missingpolicy`, `missingval`, `regrid!`.
- **Questions:** Is a threshold a fraction of total weight? Is an uncovered cell different from a missing source value? Which sentinel appears in each output wrapper? Can the destination element type hold it? How do totals differ from averages?
- **Success:** Coastal or partially covered cells use an explicit policy and a representable missing sentinel.
- **Start reading:** `api/regridding-methods.md`.
- **Validate:** Use known overlap fractions and missing samples. Check `missing`, `NaN`, and a numeric raster sentinel.

### UC-040 — Reuse a regridding plan across time or variables

- **Status / entry points:** Supported; `plan_regrid`, `regrid!`.
- **Questions:** Which geometry, ordering, and dimensions must remain fixed? Can a plan be reused after changing a mask? How is destination storage allocated? Does each field need the same missing policy? Is concurrent application supported?
- **Success:** Repeated plan application matches separate regrids and preserves the intended time or variable axes.
- **Start reading:** `tutorials/between_grids.jl`, `tutorials/regridding.jl`.
- **Validate:** Reuse one plan for two distinct fields and compare each result with an independent call.

### UC-041 — Regrid onto a regional or mixed-level destination

- **Status / entry points:** Supported; `PartialGrid`, `CellVector`, `CellLookup`, `MultiOrderCellSet` targets.
- **Questions:** Which target types infer their system and level? Does output follow target order? How do multi-level footprints affect interpretation? Can an empty destination be used? Which collections qualify as data axes?
- **Success:** Only requested cells receive values, in a documented order, with mixed-level meaning stated explicitly.
- **Start reading:** `tutorials/multiorder.jl`, regridding target docstrings.
- **Validate:** Compare a small subset against corresponding full-target results. Include a mixed-level target.

### UC-042 — Regrid back to a raster or another external space

- **Status / entry points:** Composed with GlobalRegridding; DGGS source adaptation and external destination spaces.
- **Questions:** Which destination forms belong to GlobalRegridding? How are raster resolution, extent, projection, and dimension order specified? Can the DGGS source be regional or chunked? Which package documents the full call?
- **Success:** A DGGS field transfers to one documented external destination without inventing a target constructor.
- **Start reading:** `tutorials/regridding.jl`, `lib/GlobalRegridding/README.md`, regridding API docstrings.
- **Validate:** Run a small DGGS-to-raster transfer and check coordinate alignment and constant preservation.

### UC-043 — Regrid large or tiled inputs within a memory budget

- **Status / entry points:** Supported through GlobalRegridding integration; `DGGSpace`, `PerChunk`, `Spilled`, chunked plans.
- **Questions:** Which storage option limits weights versus source data? Where are temporary weights stored? Which chunking parameters belong to the source and target? How do I estimate peak memory? What cleanup and restart behavior is supported?
- **Success:** A small forced-chunk run matches an in-memory result and exposes the relevant resource controls.
- **Start reading:** `tutorials/out_of_core.jl`, `api/partitioning.md`, `lib/GlobalRegridding/README.md`.
- **Validate:** Force several chunks and a spill path. Compare values with an in-memory plan using acceptance-test patterns.

### UC-044 — Interpret interpolation near boundaries, poles, and degeneracies

- **Status / entry points:** Supported with documented method-specific rules; point methods and dual-cell machinery.
- **Questions:** What happens outside a sample hull? How are missing corners, irregular polygons, pole cells, and regional rims handled? Which fallback preserves a value? Can interpolation overshoot source values?
- **Success:** The caller can predict edge behavior and knows which invariants remain valid.
- **Start reading:** `api/regridding-methods.md`, `test/systems/crosssystem/regrid_dual.jl`.
- **Validate:** Check exact sample sites, a regional rim, an exceptional cell, and a polar destination against point-method tests.

## G. Read, write, and inspect stored cubes

### UC-045 — Write and reopen a cell-indexed cube

- **Status / entry points:** Supported through the Zarr extension; `dggwrite`, `dggread`, `using Zarr`.
- **Questions:** What activates IO? Which input wrappers and dimensions work? How are variables, metadata, missing values, chunks, and cell axes stored? Does reading load data immediately? What happens if the path exists?
- **Success:** A cube round-trips with values, cell identities, dimensions, and essential metadata preserved.
- **Start reading:** `tutorials/store_io.jl`, `api/store-io.md`.
- **Validate:** Round-trip two variables and a time axis in a temporary directory. Check laziness separately from value equality.

### UC-046 — Select a small region from a large or remote store

- **Status / entry points:** Supported; `dggread`, `ChunkedCellLookup`, spatial selectors, `chunkof`, `chunkbounds`, `nchunks`.
- **Questions:** Which URL schemes work? What is lazy at open and selection time? Which cell-ID or data chunks must be fetched? Does the selected lookup stay chunked? How are absent cells handled?
- **Success:** A regional read returns correct labels and values without materializing the full source cube.
- **Start reading:** `tutorials/store_io.jl`, `api/store-io.md`.
- **Validate:** Use a local store with instrumented reads before testing a public URL. Compare selected values with an in-memory reference.

### UC-047 — Choose a cell-ID encoding and preserve interoperability

- **Status / entry points:** Supported; `DenseEncoding`, `RangesEncoding`, `ImplicitEncoding`, `encoding`, merge modes.
- **Questions:** Which encodings can be written versus only read? What does `:auto` choose? How do integer adjacency and valid-ID rank differ? Must IDs be sorted and unique? Which external readers understand rank-merged intervals?
- **Success:** Encoding reduces storage only where its identity and interoperability contracts hold.
- **Start reading:** `tutorials/store_io.jl`, `api/store-io.md`, encoding docstrings.
- **Validate:** Round-trip dense and range forms on a sparse-ID system. Check unsupported write requests and malformed intervals.

### UC-048 — Recognize an existing store and diagnose format errors

- **Status / entry points:** Supported; `describe_store`, `Detection`, conventions, `DGGSFormatError`.
- **Questions:** Which conventions are recognized? What metadata identifies system, level, encoding, and axis? How are conflicting or incomplete declarations reported? Can I inspect a store without reading data?
- **Success:** A reader can identify the supported convention or act on a specific format error.
- **Start reading:** `api/store-io.md`, store-description and convention docstrings.
- **Validate:** Use convention fixtures plus malformed metadata. Check that descriptions distinguish detection from successful decoding.

### UC-049 — Update or incrementally populate a store

- **Status / entry points:** Supported within writer contracts; `dggwrite!`, open `ZGroup` targets, incremental subzone writing.
- **Questions:** Which operation creates versus updates? Can I write a variable, a region, or a whole cube? Must the axis already exist? How are incompatible dimensions rejected? Are append, transactions, and concurrent writes supported?
- **Success:** A documented update changes the intended data region and preserves the store's axis contract.
- **Start reading:** `api/store-io.md`, `api/subzone-layout.md`, writer docstrings.
- **Validate:** Update a tiny existing store and reopen it. Check incompatible axes; classify append and transaction behavior separately.

### UC-050 — Store fine cells by ancestor and subzone

- **Status / entry points:** Supported; `SubzoneLayout`, `subzonestore`, qualified column/row and subzone helpers.
- **Questions:** Which systems and level pairs qualify? Why is storage two-dimensional? How are unequal descendant counts padded? Which coordinate maps a cell to a row and column? How does the layout reduce regional reads?
- **Success:** Fine-cell values round-trip through ancestor columns without treating padding as real cells.
- **Start reading:** `api/subzone-layout.md`, `tutorials/store_io.jl`.
- **Validate:** Include an exceptional ancestor with fewer descendants. Check all cell-to-row/column round trips and padding exclusion.

### UC-051 — Register a custom store convention, encoding, or grid reference

- **Status / entry points:** Supported extension surface; `register_convention!`, `register_encoding!`, `register_grid!`, qualified detect/decode/encode hooks.
- **Questions:** Which hooks and metadata must an extension provide? How are detection conflicts resolved? Does registration affect writes? What validation protects cell identity? Which names are public but qualified?
- **Success:** A minimal custom format extension has a reproducible decode/encode contract and explicit registration behavior.
- **Start reading:** `api/store-io.md`, `src/io/conventions.jl`, `src/io/encodings.jl`.
- **Validate:** Use an isolated registry fixture or documented restoration strategy. Round-trip a tiny custom example without polluting unrelated tests.

## H. Compute out of core and distribute chunks

### UC-052 — Apply a neighborhood kernel to a stored global dataset

- **Status / entry points:** Supported; `chunkplan`, `mapneighbors!`, `chunkcube`.
- **Questions:** How do I choose halo width? Does the plan follow existing storage chunks? What output shape and chunking are required? Are halo values read from the source? Is in-place aliasing safe?
- **Success:** A chunked run matches the in-memory kernel, including cells along chunk boundaries.
- **Start reading:** `tutorials/out_of_core.jl`, `api/chunk-sweep.md`.
- **Validate:** Force several small chunks. Compare every value with UC-027 and include a region with a hole.

### UC-053 — Write a custom chunk callback

- **Status / entry points:** Supported; `foreachchunk`, `MapChunk`, `ChunkCube`, `ownedindices`, `localindices`, `axisindices`, `chunkhalo`, `halowidth`.
- **Questions:** Which positions belong to the source axis, chunk buffer, or owned output? Does a buffer include halos first or last? Can callbacks retain buffers? How do extra dimensions and empty chunks behave?
- **Success:** A callback reads surrounding data and writes each owned output cell exactly once.
- **Start reading:** `api/chunk-sweep.md`, chunk type docstrings.
- **Validate:** Use values equal to source-axis positions. Verify ownership disjointness and full output coverage.

### UC-054 — Run wider or repeated stencils across chunk boundaries

- **Status / entry points:** API gap for direct wider sweeps; `chunkplan` can provide wider input halos, but `mapneighbors!` still supplies one-ring callback neighborhoods.
- **Questions:** What natural API should request a wider neighborhood? How much input halo does the requested operation need? How does repeated diffusion differ from one radius-two kernel? When do intermediate halo values require another pass or exchange?
- **Success:** The audit separates the input halo width from callback neighborhood width. It records missing wider-sweep support and any missing multi-pass orchestration in `docs/audit/api-gaps.md`.
- **Start reading:** `api/chunk-sweep.md`, `api/neighbors.md`, `tutorials/out_of_core.jl`.
- **Validate:** Inspect callback membership with `halo=2` on tiny chunks. Verify that documentation does not claim this changes the one-ring sweep into UC-026.

### UC-055 — Assign chunks to threads or distributed workers

- **Status / entry points:** Supported planning plus composed execution; `partitionproblem`, `partition`, plan slicing, `partindices`, `partchunks`, `partsources`, `partweights`.
- **Questions:** Does partitioning launch workers? What data and packages must each worker open? Which output writes are disjoint? Can plans be serialized? How do I avoid nested threading and duplicated source reads?
- **Success:** Every task runs once, and parallel output matches serial output with documented worker setup.
- **Start reading:** `api/partitioning.md`, `tutorials/out_of_core.jl`.
- **Validate:** Use two partitions and a small store. Follow `test/partitioning/distributed.jl` for process-level checks.

### UC-056 — Balance work and reduce shared source reads

- **Status / entry points:** Supported; `WeightedContiguous`, `MetisPartition`, `KaHyParPartition`, `ScotchPartition`, capacities and weights.
- **Questions:** What objective does each algorithm optimize? Which optional package enables it? What do source weights, work weights, capacities, and imbalance mean? Which seed controls reproducibility? How does edge limiting affect a backend?
- **Success:** The reader selects an available algorithm and can inspect both load balance and source replication.
- **Start reading:** `api/partitioning.md`.
- **Validate:** Use a known dependency graph and unequal worker capacities. Check missing-backend errors and deterministic seeded behavior where promised.

### UC-057 — Partition a chunked regridding plan

- **Status / entry points:** Supported; `partitionproblem` for GlobalRegridding plans/dependency graphs and `partition`.
- **Questions:** Which destination chunks form tasks? Which source chunks create dependencies? How are assigned chunks executed? Does partitioning preserve plan reuse and spill behavior? Who coordinates writing results?
- **Success:** Partitioned regridding produces the same destination as the original plan and exposes its source dependencies.
- **Start reading:** `api/partitioning.md`, `lib/GlobalRegridding/README.md`.
- **Validate:** Partition a tiny many-to-many chunk graph. Check exact task coverage and compare serial results.

## I. Complete domain workflows and display results

### UC-058 — Analyze elevation, roughness, slope, and flow directions

- **Status / entry points:** Composed; regridding, cell fields, neighborhood sweeps, tutorial hydrology algorithms.
- **Questions:** Are DEM values samples or area means? What distance and elevation units produce slope? How are pits, flats, missing cells, and pentagons handled? Which operations belong to this package versus companion hydrology tools?
- **Success:** A small DEM workflow has physically interpretable output and states its treatment of unresolved drainage cases.
- **Start reading:** `tutorials/hydrology.jl`, `tutorials/stencils.jl`, `api/neighbor-fields.md`.
- **Validate:** Use a synthetic plane or bowl before a downloaded DEM. Consult `test/integration/geomorphometry_synthetic.jl` and related integration fixtures.

### UC-059 — Use Copernicus DEM's native tile and pixel hierarchy

- **Status / entry points:** Supported specialized system; `CopernicusDEMSystem`, tiled cursor and source integration.
- **Questions:** Which product resolutions are accepted? What do its two levels represent? How do latitude bands change pixel dimensions and adjacency? Are sample posts or pixel footprints represented? Where are tile gaps, rims, and polar limits handled?
- **Success:** A small native tile source maps points and cells correctly before transfer to another grid.
- **Start reading:** CopernicusDEM module docstring, `tutorials/hydrology.jl`, `scripts/dagger_regrid/README.md`.
- **Validate:** Use local CopernicusDEM fixtures, including a latitude-band transition, a tile seam, and source-mode tests.

### UC-060 — Analyze a HEALPix sky map and perform a cone search

- **Status / entry points:** Composed; `HEALPixSystem`, `HEALPixRingIndex`, `reindex`, cell lookups, spherical-cap selection.
- **Questions:** How do `nside`, level, and pixel count relate? Is input ordering nested or ring? Does converting IDs also reorder values? Which celestial coordinate frame is assumed? Does this package read FITS itself?
- **Success:** Pixel values align with nested cells, and a cone selection has an explicit footprint rule and coordinate frame.
- **Start reading:** `tutorials/healpix_astronomy.jl`, HEALPix system docstrings.
- **Validate:** Use a vector whose values encode original pixel positions. Compare ring/nested mapping with fixture or Healpix.jl checks.

### UC-061 — Plot grids, regions, and cell values on a map or globe

- **Status / entry points:** Composed with visualization/Makie packages; `dggpoly`, Makie extension for cell collections.
- **Questions:** Which plotting package and backend must I load? Which plotting names belong to the companion package? How does value order match cells? How are missing colors, seams, projections, and outlines handled?
- **Success:** A complete grid and a selected field render with correct cell-value alignment in a documented plotting environment.
- **Start reading:** `all_dggs.md`, tutorials, `lib/DiscreteGlobalGridsVisualization/README.md`.
- **Validate:** Render a tiny categorical grid and subset. Use `test/plotting/runtests.jl` and companion tests for the applicable API.

### UC-062 — Display continuous surfaces or large interactive fields

- **Status / entry points:** Companion-package functionality; `dggsurface`, `dggresample`.
- **Questions:** Does display interpolation change the stored data? How are heights scaled? What limits interactive resolution? Which backends and GeoMakie axes work? Which plotting API is experimental?
- **Success:** A reader can choose cell polygons, an interpolated surface, or a zoom-dependent display with its limitations stated.
- **Start reading:** `lib/DiscreteGlobalGridsVisualization/README.md`.
- **Validate:** Run one small surface and one interactive resampling example in the documented companion environment.

## J. Extend the package and understand its boundaries

### UC-063 — Implement a new finite grid

- **Status / entry points:** Supported extension surface; `AbstractGrid`, base grid interface, generic fallbacks.
- **Questions:** Which methods are mandatory? What cell identity and geometry types are required? Which methods use local indices? Which queries and topology operations become available automatically? What costs do fallbacks introduce?
- **Success:** A minimal new grid passes the conformance suite and supports a documented query and geometry example.
- **Start reading:** `extending.md`, `api/grid-interface.md`, `lib/DiscreteGlobalGridsConformanceTesting/README.md`.
- **Validate:** Build a tiny fixture implementation and run `test_grid_interface`; compare the guide with the current required contracts.

### UC-064 — Implement a hierarchical or quad-face system

- **Status / entry points:** Supported extension surface; `AbstractHierarchicalGridSystem`, `AbstractQuadFaceGridSystem`, `HierarchicalLevelGrid`, hierarchy primitives and traits.
- **Questions:** Can a system use the generic level-grid type? Which primitive methods are mandatory? When may sorted-subtree, direct-location, or congruent-refinement traits be true? Which quad-face helpers are public contracts?
- **Success:** A small system satisfies hierarchy laws without making unsupported optimization claims.
- **Start reading:** `extending.md`, `architecture.md`, `api/grid-interface.md`.
- **Validate:** Run `test_hierarchical_system` and selected cross-system laws. Check trait declarations against actual behavior.

### UC-065 — Add an optimized boundary walker without changing semantics

- **Status / entry points:** Advanced extension surface; boundary engines, subtree walkers, cap bounds.
- **Questions:** Which hooks may a system override? What order, memory, and completeness contracts apply? When must a specialized path fall back? How are seams and exceptional cells validated? Which helpers are implementation details?
- **Success:** An optimization agrees with a small exhaustive boundary reference across its guards and exceptions.
- **Start reading:** `internals/boundary-engines.md`, `extending.md`.
- **Validate:** Use subtree border, interior, and halo suites at several depths. Include the specialization's fallback cases.

### UC-066 — Add a custom partitioning backend

- **Status / entry points:** Supported extension surface; `AbstractPartitioningAlgorithm`, `partitionlabels`, `PartitionProblem`.
- **Questions:** Which labels, weights, capacities, and task order must a backend honor? Who validates its output? How is optional loading arranged? What happens for empty work or more workers than tasks?
- **Success:** A minimal backend returns valid assignments and integrates with plan slicing and result accessors.
- **Start reading:** `api/partitioning.md`, `src/partitioning_backends.jl`.
- **Validate:** Reuse partitioning contract tests with a simple custom algorithm and invalid-output fixtures.

### UC-067 — Migrate older package examples

- **Status / entry points:** Supported migration where shims exist; `src/deprecated.jl`, deprecated `cellposition` and `globalindices` names.
- **Questions:** Which old names still work? What replaced them? Did index semantics change? Which examples refer to packages or plotting names that moved? Are deprecation warnings actionable?
- **Success:** An old minimal workflow uses current names without silently changing index spaces.
- **Start reading:** `index.md`, current API pages, `src/deprecated.jl`.
- **Validate:** Compare each shim with its documented replacement on a noncontiguous region or chunk plan.

### UC-068 — Run tests and build the documentation locally

- **Status / entry points:** Supported repository workflow; `test/runtests.jl`, `docs/make.jl`, project files.
- **Questions:** Which environment and local workspace dependencies are required? Which tests need optional packages or external data? Does a docs build execute tutorials? Which warnings fail the build? Does the build attempt deployment?
- **Success:** A contributor can run a focused check and distinguish a rendered page from a verified executable example.
- **Start reading:** repository `README.md`, project metadata, `docs/make.jl`, conformance package README.
- **Validate:** Read the actual runner configuration before executing it. Report exact commands, environments, exclusions, and warnings.

### UC-069 — Export standard vector, raster, astronomy, or DGGS exchange formats

- **Status / entry points:** Composed or unsupported by the core package; geometry adapters, regridding, external format packages.
- **Questions:** Is there a direct GeoJSON, Shapefile, GeoTIFF, NetCDF, FITS, or MOC writer? Which external package completes the workflow? Are mixed-level sets equivalent to a standard MOC? Which metadata and ID conversions are necessary?
- **Success:** Each requested format has a tested external handoff or an explicit statement that no direct writer exists.
- **Start reading:** `api/store-io.md`, `tutorials/healpix_astronomy.jl`, relevant geometry/regridding docs.
- **Validate:** Trace available methods before claiming support. Test one documented handoff; do not infer standards compliance from similar data structures.

### UC-070 — Apply arbitrary metric-radius or physical convolution kernels

- **Status / entry points:** Composed; spherical-cap queries, cell centroids, distance calculations, user-defined weights.
- **Questions:** Is a topological ring equivalent to a radius in kilometers? Does a cap select intersecting footprints or centers? How are area and distance weights normalized? Is there a built-in Gaussian kernel or physical Laplacian?
- **Success:** A distance-based filter states its selection and weighting rules instead of treating neighbor order as a metric distance.
- **Start reading:** `api/neighbors.md`, `api/selecting-cells.md`, `tutorials/stencils.jl`.
- **Validate:** Compare a small cap-based kernel with explicit spherical distances. Classify named kernels as user composition unless an API exists.

### UC-071 — Perform adaptive PDE, finite-volume, or mixed-level halo calculations

- **Status / entry points:** Unsupported as a complete solver; partial building blocks include geometry, adjacency, and `member_neighbors`.
- **Questions:** Are face fluxes, oriented edge geometry, time stepping, refinement criteria, or conservative restriction/prolongation provided? Can mixed-level sets have a halo? Do noncongruent refinements support the assumed mesh law?
- **Success:** Documentation separates available building blocks from solver functionality that callers must implement.
- **Start reading:** `api/neighbors.md`, `api/selecting-cells.md`, hierarchy trait docstrings.
- **Validate:** Check for actual public methods. Preserve the explicit absence of a single-level halo on `MultiOrderCellSet`.

### UC-072 — Use accelerators or automatic distributed execution

- **Status / entry points:** Uncertain or unsupported as a turnkey workflow; CPU sweeps and explicit chunk partitioning are visible.
- **Questions:** Are GPU arrays, automatic differentiation, device kernels, or automatic cluster scheduling supported? Does partitioning execute tasks? Which parts require scalar indexing or CPU-only native libraries?
- **Success:** A reader finds an explicit support boundary and a tested CPU or manually distributed route.
- **Start reading:** `api/chunk-sweep.md`, `api/partitioning.md`, package extension metadata.
- **Validate:** Require an implementation and a representative test before claiming accelerator or automatic scheduler support.

### UC-073 — Publish, append, or concurrently modify a remote store

- **Status / entry points:** Partly unsupported; documented remote reads and local/open-group writes do not establish remote publication or transaction support.
- **Questions:** Does writing a URL work? Who uploads a local directory? Are credentials managed here? Can writers append cells or coordinates? Are concurrent writes atomic, locked, or resumable?
- **Success:** The documented workflow identifies the external publication step and does not promise transaction behavior without evidence.
- **Start reading:** `tutorials/store_io.jl`, writer docstrings.
- **Validate:** Inspect extension methods and test coverage. Distinguish disjoint local chunk writes from remote transaction guarantees.

### UC-074 — Use arbitrary planetary bodies or projected coordinate systems

- **Status / entry points:** Partly supported through manifold wrappers; arbitrary system and projection support requires verification.
- **Questions:** Which radius or ellipsoid parameters are accepted? Which calculations remain unit-sphere quantities? Can source projections be transformed by the regridding stack? Is a planar grid accepted by the spherical contracts?
- **Success:** Units and coordinate transformations remain explicit, and the guide does not imply unrestricted CRS or planetary support.
- **Start reading:** `api/grid-interface.md`, authalic wrapper docstrings, `api/regridding-methods.md`.
- **Validate:** Test a documented nondefault ellipsoid or radius if available. Check unsupported projections and transformations against actual adapters.

### UC-075 — Estimate scalability before constructing a global high-resolution workflow

- **Status / entry points:** Supported planning from counts, compressed collections, iterators, chunks, and documented traits.
- **Questions:** Which objects are lazy versus materialized? When does a call allocate one entry per cell or edge? Which iterators lack exact lengths? Why can A5 take a slower boundary path? How do plan, halo, and data memory combine?
- **Success:** A reader can choose a small prototype and a chunked production route without accidentally enumerating an enormous global grid.
- **Start reading:** `systems` docstring, `api/boundaries.md`, `api/chunk-sweep.md`, `api/partitioning.md`.
- **Validate:** Inspect container sizes and allocation behavior at small levels. Use documented complexity bounds; do not extrapolate performance guarantees from one timing.

## Cross-cutting checks for every assigned workflow

These checks avoid repeating every system, shape, and storage combination as a separate use case.
An agent should choose the combinations that can expose a different contract.

| Axis | Required audit questions or representative cases |
| --- | --- |
| Prerequisites | Can a newcomer identify every import, optional extension, environment, and downloaded input? |
| Entry point | Does a page explain why this API solves the task before introducing its implementation? |
| Identity | Does every integer name a local position, complete-level position, chunk position, or raw identifier? |
| Coordinates | Are angle units, longitude/latitude order, manifold, area units, and distance units explicit? |
| Systems | Compare a hexagonal system and a congruent quadrilateral system. Add A5, S2, ISEA4R, or CopernicusDEM where their exceptions matter. |
| Geometry | Check a seam, pole, exceptional cell, region boundary, and hole when the algorithm can encounter them. |
| Collections | Compare a complete level, a noncontiguous subset, and a subtree. Use mixed levels only where explicitly accepted. |
| Data | Check a plain vector, a cell-indexed array, and extra dimensions when promised. Use values that expose ordering errors. |
| Empty and invalid inputs | Establish empty-region, absent-cell, invalid-level, invalid-ID, and incompatible-system behavior for the relevant entry point. |
| Missing data | Distinguish absent cells, missing samples, no overlap, and numeric missing sentinels. |
| Ordering and mutation | State sorting, duplicate handling, callback buffer lifetime, aliasing, and permitted mutation. |
| Resources | Identify lazy objects, full materialization, optional length, threading, temporary files, and chunk boundaries. |
| Interoperability | Distinguish a shared name from a compatible identifier scheme, coordinate frame, or exchange standard. |
| Outcome | End each runnable example with a useful result and one check that detects a plausible silent error. |

## Evidence map and validation notes

The initial sweep inspected the public/exported surface in `src/DiscreteGlobalGrids.jl`.
It also inspected the tutorial headings, API guides, core neighborhood contracts, extension surfaces, and test runners.
The inventory includes the visualization companion because current tutorials rely on its plotting API.
It includes GlobalRegridding because the package re-exports that package's regridding verbs and implements its space contract.

| Workflow group | Main implementation evidence | Existing validation entry points |
| --- | --- | --- |
| A–B: systems, sizing, identity, geometry | `src/interface/`, `src/fallbacks/`, `src/sizing.jl`, `src/systems/` | `test/interface/`, `test/fallbacks/`, `test/systems/*/runtests.jl` |
| C–D: regions, selection, hierarchy | `src/engine/query.jl`, `cell_vector.jl`, `region_algebra.jl`, `multiorder*.jl`, `src/dimensionaldata.jl` | `test/systems/crosssystem/cell_vector.jl`, `dimensionaldata.jl`, `region_algebra.jl`, `multiorder*.jl` |
| E: neighborhoods and boundaries | `src/engine/stencil.jl`, `neighborhood.jl`, `adjacency.jl`, `halo.jl`, `needs.jl`, `cellfield.jl` | `test/systems/crosssystem/stencils.jl`, `neighborhood.jl`, `mapneighbors.jl`, `needs.jl`, `subtree_iterators.jl`, `subtree_halos.jl` |
| F: regridding | `src/regridding.jl`, `src/dual_cells.jl`, `lib/GlobalRegridding/` | `test/systems/crosssystem/regrid*.jl`, `regridding_conservation.jl` |
| G: storage | `src/io/`, `ext/DiscreteGlobalGridsZarrExt/` | `test/io/` |
| H: chunk execution and assignment | `src/chunks.jl`, `src/partitioning*.jl`, optional partitioner extensions | `test/io/chunk_sweep.jl`, `test/partitioning/` |
| I: domain workflows and plotting | tutorials, `src/systems/CopernicusDEM/`, `ext/MakieExt/`, visualization companion | `test/integration/`, `test/systems/CopernicusDEM/`, `test/plotting/`, companion tests |
| J: extension contracts | `src/interface/`, fallback modules, registry and backend contracts | `lib/DiscreteGlobalGridsConformanceTesting/`, relevant extension tests |

No repository `AGENTS.md` file was found during this checkout's initial file sweep.
The session-provided instructions therefore remain the applicable audit instructions.
The inventory applies the simplified technical English skill requested for audit writing.

`test/runtests.jl` includes system, interface, IO, plotting, regridding, partitioning, and selected script tests.
Standalone integration files exist outside that runner; do not claim the main suite ran them automatically.
The conformance harness also has its own workspace package and self-tests.
Use the repository's actual project environment and applicable Julia execution skill before running Julia checks.

`docs/make.jl` generates tutorial Markdown with `Literate.markdown(...; execute=false)`.
Documenter can execute resulting example blocks, so distinguish Literate generation from the full documentation build.
The build uses `warnonly=true` and `checkdocs=:none`.
A zero exit status does not prove that all examples succeeded or all public docstrings appear on the site.
The script ends with `deploydocs`; inspect and control the build environment before treating it as a local rendering command.
Record warnings and verify executable examples separately from page generation.

## Suggested assignment order

1. Audit UC-001–024 to establish vocabulary, coordinates, identity, and selection.
2. Audit UC-025–044 for complete computational workflows, including the second-order global convolution case.
3. Audit UC-045–057 for persistence, chunk boundaries, and worker assignment.
4. Audit UC-058–075 for domain completion, extension contracts, and explicit support boundaries.
5. Re-run representative newcomer paths after documentation revisions. Record remaining gaps against the same identifiers.

This file is an inventory only. It does not modify tutorials, API behavior, or test results.
