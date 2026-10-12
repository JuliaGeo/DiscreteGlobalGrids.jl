# Core interface documentation sweep

This report records findings only. It does not change public documentation or source docstrings.

## Scope

The sweep covered grid construction, grid systems, identifiers, geometry, hierarchy,
spatial queries, regions, mixed-level sets, and level selection. It used local source
and tests as the authority.

The sweep excluded neighborhood sweeps, stencils, cell fields, store I/O, regridding,
chunk execution, and partitioning. Those areas belong to separate sweeps.

Severity has this meaning:

- **Critical:** The documentation states behavior that the API does not provide.
- **High:** A common task can produce a wrong choice or wrong interpretation.
- **Medium:** A common task requires source reading or avoidable trial and error.
- **Low:** The issue affects clarity but does not usually block a task.

Each finding is either a **DOC GAP** or an **API GAP**. An API gap needs a later
design decision. This report does not prescribe a temporary workaround.

## Findings



### SC-001 — The selection guide states the wrong return type

- **Kind:** DOC GAP
- **Severity:** Critical
- **Affected surface:** `query`, `CellVector`, `MultiOrderCellSet`
- **Evidence:** `docs/src/api/selecting-cells.md:12-16` says a one-level answer is a
`CellVector`. The implementation returns a sorted `Vector` of cell IDs at
`src/engine/query.jl:480-506`. Only a mixed-level query returns a
`MultiOrderCellSet`, at `src/engine/multiorder.jl:179-182`.
- **A newcomer needs:** The exact return container for each query form. They also
need to know when to convert a result to a `CellVector` or `PartialGrid`.
- **Current documentation:** The guide gives the wrong one-level container. The
reference signature gives the correct type, but the two statements conflict.
- **Proposed correction:** Correct the guide. Add a two-row table for fixed-level
and mixed-level queries. Show the next conversion for data indexing and region
operations.



### SC-002 — `systems()` claims to list shipped systems but omits Copernicus DEM

- **Kind:** DOC GAP
- **Severity:** High
- **Affected surface:** `systems`, `CopernicusDEMSystem`
- **Evidence:** The function docstring calls its result "The grid systems shipped
by this package" at `src/DiscreteGlobalGrids.jl:243-254`. Source comments state
that `CopernicusDEMSystem` is deliberately absent at
`src/DiscreteGlobalGrids.jl:239-242`. The package exports that system at
`src/DiscreteGlobalGrids.jl:497`.
- **A newcomer needs:** Whether `systems()` is exhaustive. They also need the route
to a system that does not appear in it.
- **Current documentation:** The reason for the omission exists only in a source
comment. The public docstring does not name the exception.
- **Proposed correction:** Define `systems()` as the uniform global-grid registry.
Name `CopernicusDEMSystem` as an exported exception. Link to its constructor and
explain why generic cross-system loops omit it.



### SC-003 — The grid-choice material omits two public systems from its decision table

- **Kind:** DOC GAP
- **Severity:** High
- **Affected surface:** `S2System`, `CopernicusDEMSystem`, system selection
- **Evidence:** The decision table lists five systems at
`docs/src/tutorials/choosing_a_grid.jl:39-50`. It omits S2 and Copernicus DEM.
The gallery includes S2 at `docs/src/all_dggs.md:24-27`, but it also omits
Copernicus DEM. The Copernicus constructor exists at
`src/systems/CopernicusDEM/bands.jl:7-26`.
- **A newcomer needs:** One complete list with intended use, cell geometry, area
behavior, coordinate frame, identifier compatibility, and major limitations.
- **Current documentation:** A reader must combine the tutorial, gallery, the very
long `systems()` docstring, and system source files. Copernicus DEM has no entry
in the user documentation.
- **Proposed correction:** Add S2 to the main choice table. Add a separate section
for the Copernicus DEM raster lattice. Keep one compact system matrix as the
source of truth.



### SC-004 — Public reference pages embed contributor contracts and internal hooks

- **Kind:** DOC GAP
- **Severity:** High
- **Affected surface:** Core grid reference and most core docstrings
- **Evidence:** The user-facing grid page includes `node_cell`, `node_indices`,
`one_ring`, `cap_inflation`, `authalic_stretch`, traversal engines,
`lattice_decode`, and tree hooks at `docs/src/api/grid-interface.md:59-141`.
The `node_extent` contract alone spans `src/interface/system.jl:146-190`.
`CellVector` mixes usage with storage representation and A5 implementation
details at `src/engine/cell_vector.jl:220-284`. The `systems()` docstring contains
engine algorithms at `src/DiscreteGlobalGrids.jl:272-357`.
- **A newcomer needs:** Short signatures, results, units, valid inputs, important
errors, and one normal example. A system author needs the full contracts.
- **Current documentation:** Both audiences receive the same text. Important task
details are difficult to find inside implementation rationale and performance
history.
- **Proposed correction:** Split usage reference from extension contracts. Keep
public docstrings short. Move invariants, traversal rules, and benchmarks to the
extension or internals pages. Link between the two layers.



### SC-005 — Query predicate direction is implicit

- **Kind:** DOC GAP
- **Severity:** High
- **Affected surface:** `query`, `Intersects`, `Contains`, `Within`, `Covers`,
`CoveredBy`, `Touches`, `Overlaps`, `Equals`, `Crosses`
- **Evidence:** The engine defines `Intersects(target)` as "cell intersects target"
and reverses asymmetric predicates at `src/engine/query.jl:376-388`. The query
docstring lists names and target types at `src/engine/query.jl:480-499`, but it
does not state the subject and object for each relation. The selection page lists
`Crosses` beside supported predicates at `docs/src/api/selecting-cells.md:41-59`.
The engine rejects `Crosses` at `src/engine/query.jl:390-395`. Cap targets support
only three predicates at `src/engine/query.jl:472-474`.
- **A newcomer needs:** A direct statement such as "return each cell for which
`cell RELATION target` is true." They also need a support matrix by target type.
- **Current documentation:** A reader can easily reverse `Contains` and `Within`.
The predicate index implies that every listed type works with `query`.
- **Proposed correction:** Add a small truth table with relation direction,
supported target types, and one example per asymmetric pair. Mark `Crosses` as
unsupported by this query engine.



### SC-006 — The zonal tutorial identifies a center-based rule that the API cannot express

- **Kind:** API GAP
- **Severity:** High
- **Affected surface:** Spatial selection by cell centroid
- **Evidence:** The tutorial compares three zonal rules at
`docs/src/tutorials/zonal.jl:127-138`. The center-in-zone row has no spelling.
The query engine tests DE9IM relations against the complete cell polygon at
`src/engine/query.jl:399-416`. No public predicate selects cells whose centroid
lies in the target.
- **A newcomer needs:** A direct selector for a standard zonal rule: keep a cell
when its representative point lies in a polygon.
- **Current API:** `Intersects` includes a boundary rim. `Within` requires the full
cell to fit. Neither operation means center-in-zone.
- **Proposed correction:** Design a named centroid selection predicate that works
with `query`, `predicate_indices`, and `Cells(...)`. Then give the missing table
entry a real spelling.



### SC-007 — Multi-order coverage guarantees contradict each other

- **Kind:** DOC GAP
- **Severity:** Critical
- **Affected surface:** `MultiOrderCoverage`, mixed-level `query`, `covering`
- **Evidence:** The type docstring promises leaf coverage without limiting that
promise to one mode at `src/engine/multiorder.jl:32-46`. The query docstring then
says every target point lies in an emitted cell at
`src/engine/multiorder.jl:144-152`. The same docstring reports misses for IGeo7,
H3, and A5 budget queries at `src/engine/multiorder.jl:154-166`. The tutorial
introduces a parent as standing for all target-level descendants at
`docs/src/tutorials/multiorder.jl:30-38`.
- **A newcomer needs:** Separate guarantees for `level` mode and `maxcells` mode.
They need to know which guarantees depend on congruent refinement.
- **Current documentation:** One paragraph promises coverage. A later table says
the implementation can miss the target. Both cannot be the public contract.
- **Proposed correction:** State the two modes separately. Give exact guarantees
only where the implementation and tests establish them. Label empirical miss
rates as measurements, not contracts. Repeat the distinction in `covering` and
the tutorial.



### SC-008 — `maxcells` is presented as a hard cap, but the seed can exceed it

- **Kind:** DOC GAP
- **Severity:** High
- **Affected surface:** `query(sys, MultiOrderCoverage(...); maxcells)`
- **Evidence:** The tutorial section is titled "Cap the cell count with `maxcells`"
at `docs/src/tutorials/multiorder.jl:157-170`. The function docstring states that
a seed larger than the budget returns over budget at
`src/engine/multiorder.jl:128-139`.
- **A newcomer needs:** Whether `length(result) <= maxcells` always holds.
- **Current documentation:** The exception appears deep in a long docstring. The
tutorial heading states the stronger rule.
- **Proposed correction:** Call `maxcells` a refinement budget near every first
example. State the seed exception in the first paragraph and in the tutorial.



### SC-009 — The region and collection types lack a task-oriented conversion map

- **Kind:** DOC GAP
- **Severity:** High
- **Affected surface:** `AbstractGrid`, `PartialGrid`, `CellVector`, `CellLookup`,
`MultiOrderCellSet`, `region`, `cellset`
- **Evidence:** The selection page lists type docstrings at
`docs/src/api/selecting-cells.md:18-30`. The conversion entry points are spread
across `src/engine/cell_vector.jl:220-265`,
`src/engine/cell_vector.jl:549-566`, and `src/dimensionaldata.jl:76-116`.
`cellset` has provenance behavior at `src/engine/cell_vector.jl:473-490`.
- **A newcomer needs:** Which type to use for geometry, cell IDs, a data axis, and
mixed-level storage. They also need the return type and cost of each conversion.
- **Current documentation:** Each type explains itself in depth. No page answers
"I have this type; which conversion gives the operation I need?"
- **Proposed correction:** Add a compact conversion table. Include the five main
types, `query` results, `region`, `CellVector`, `PartialGrid`, `CellLookup`,
`expand`, and `cellset`.



### SC-010 — `PartialGrid` has no safe convenience constructor for ordinary ID input

- **Kind:** API GAP
- **Severity:** Medium
- **Affected surface:** `PartialGrid(sys, level, ids)`
- **Evidence:** The constructor requires an exact canonical ID element type and
strict ascending order at `src/engine/partial_grid.jl:68-83`. It stores `ids` by
reference at `src/engine/partial_grid.jl:37-42`. The public selection page gives
no construction example at `docs/src/api/selecting-cells.md:18-30`.
- **A newcomer needs:** A constructor that accepts a normal collection of cell IDs
and produces a valid region without hidden aliasing requirements.
- **Current API:** The caller must prepare the exact type, order, and uniqueness.
Later mutation of the referenced vector can invalidate constructor checks.
- **Proposed correction:** Design an owning convenience constructor that normalizes
and validates IDs. Keep an explicit no-copy constructor for advanced callers who
need it. Document ownership in both signatures.



### SC-011 — Geographic geometry extraction is not shown in the core API

- **Kind:** DOC GAP
- **Severity:** Medium
- **Affected surface:** `cell_boundary`, `cell_centroid`, `UnitSphericalPoint`
- **Evidence:** The core contract returns unit-sphere points at
`src/interface/grid.jl:52-89`. The grid reference embeds those docstrings but
gives no geographic conversion example at `docs/src/api/grid-interface.md:17-34`.
Conversion appears later inside specialized tutorials, for example
`docs/src/tutorials/choosing_a_grid.jl:178-185`.
- **A newcomer needs:** How to obtain longitude and latitude from one centroid or
boundary. They also need the output units and longitude convention.
- **Current documentation:** The first API page exposes `(x, y, z)` geometry and
leaves the conversion spelling in unrelated tutorials.
- **Proposed correction:** Add one short conversion example beside the geometry
functions. State degrees and longitude range. Link to coordinate-frame guidance.



### SC-012 — The plural polygon API does not support the collections named by the guide

- **Kind:** API GAP
- **Severity:** High
- **Affected surface:** `cell_polygons`
- **Evidence:** The boundary guide tells readers to use `cell_polygons` "for a
collection" at `docs/src/api/boundaries.md:24-26`. The only public method accepts
`MultiOrderCellSet` at `src/engine/multiorder_set.jl:98-109`. There is no method
for a complete grid, `PartialGrid`, `CellVector`, or `CellLookup`.
- **A newcomer needs:** One plural geometry operation that works across the region
types used elsewhere in the package.
- **Current API:** The guide promises a general collection operation. Most public
collection types raise a `MethodError` for that call.
- **Proposed correction:** Define the intended common collection contract for
`cell_polygons`. Add methods for the supported region types, or narrow and rename
the operation if mixed-level sets are the only intended input.



### SC-013 — The hierarchy reference has no usage path

- **Kind:** DOC GAP
- **Severity:** Medium
- **Affected surface:** `parent`, `children`, `ancestor`, `descendants`,
`descendant_range`, `subtree`
- **Evidence:** The grid reference presents only an `@docs` list at
`docs/src/api/grid-interface.md:46-57`. `subtree` appears on a different page at
`docs/src/api/boundaries.md:95-102`. The source contracts define materially
different products at `src/interface/system.jl:474-517` and
`src/interface/system.jl:614-638`.
- **A newcomer needs:** A small example that starts with one cell and shows an
immediate relation, a target-level list, an index range, and a region.
- **Current documentation:** The reader must infer when to choose `descendants`,
`descendant_range`, or `subtree`. The A5 availability limit is easy to miss.
- **Proposed correction:** Add a hierarchy task table and one executable example.
Put `subtree` beside the other hierarchy operations. Show that
`descendant_range` is available only when `has_sorted_subtrees(sys)` is true.



### SC-014 — The `descendants` container and cost contract is inconsistent

- **Kind:** DOC GAP
- **Severity:** Medium
- **Affected surface:** `descendants`
- **Evidence:** The general docstring says the call materializes `O(subtree)` IDs at
`src/interface/system.jl:486-499`. Copernicus DEM returns a lazy read-only vector
and calls this an interface divergence at `src/systems/CopernicusDEM/system.jl:147-159`.
The general signature does not promise an owned `Vector`.
- **A newcomer needs:** The operations they may use on the result and whether they
can mutate it. They should not have to depend on one system's allocation model.
- **Current documentation:** The generic contract implies materialization. One
shipped system explicitly does something else.
- **Proposed correction:** Specify an `AbstractVector` or iterable return contract.
State ordering and mutability. Describe allocation as implementation-specific.
Remove the claim that the lazy implementation diverges if it satisfies the
revised contract.



### SC-015 — `levelfor(...; over=...)` does not restrict both sides as documented

- **Kind:** API GAP
- **Severity:** Critical
- **Affected surface:** `levelfor`
- **Evidence:** The docstring says `over` restricts both sides of the comparison at
`src/sizing.jl:54-71`. The implementation applies `over` only to candidate system
levels at `src/sizing.jl:73-86`. `_targetarea` receives no `over` argument and
measures a grid, raster, or regrid space globally at `src/sizing.jl:90-96`.
- **A newcomer needs:** A level choice that compares both datasets over the same
area when cell size varies by location.
- **Current API:** The documented call silently compares a regional system sample
with a global target sample. This can select a different level.
- **Proposed correction:** Extend target measurement to accept the same area of
interest, or narrow the public contract. Add a test with a location-dependent
target where regional and global medians differ.



### SC-016 — Level-direction wording is easy to reverse

- **Kind:** DOC GAP
- **Severity:** Low
- **Affected surface:** `expand`, `ancestor`, `descendants`, `subtree`
- **Evidence:** `expand` says "`l` above it throws" at
`src/engine/region_algebra.jl:57-68`, while the implementation accepts numeric
levels greater than or equal to the current level at
`src/engine/region_algebra.jl:77-85`. Other contracts use "deeper" and
"shallower" with explicit inequalities at `src/interface/system.jl:474-499`.
- **A newcomer needs:** One consistent rule: larger level numbers mean finer cells.
- **Current documentation:** Spatial words such as "above" require the reader to
infer tree orientation.
- **Proposed correction:** Use "coarser," "finer," and explicit inequalities in
every hierarchy docstring. Avoid "above" and "below" for level numbers.



## Coverage ledger

The following table makes the sweep boundary explicit.


| Area                              | Public source inspected                                                                                                                                                                                              | User documentation inspected                                                            | Result                                                           |
| --------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------- | ---------------------------------------------------------------- |
| Module surface and registry       | `src/DiscreteGlobalGrids.jl`                                                                                                                                                                                         | `docs/src/index.md`, `docs/src/all_dggs.md`                                             | SC-002, SC-003, SC-004                                           |
| Core types and index spaces       | `src/interface/types.jl`, `src/interface/grid.jl`                                                                                                                                                                    | `docs/src/abstractions.md`, `docs/src/api/grid-interface.md`                            | SC-004, SC-011, SC-016                                           |
| Hierarchical system contract      | `src/interface/system.jl`, `src/fallbacks/level_grid.jl`, `src/fallbacks/subtree.jl`, `src/fallbacks/subtree_iterators.jl`                                                                                           | `docs/src/api/grid-interface.md`, `docs/src/api/boundaries.md`, `docs/src/extending.md` | SC-013, SC-014, SC-016                                           |
| Geometry and location             | `src/fallbacks/geometry.jl`, `src/fallbacks/locate.jl`, `src/core/manifolds.jl`, `src/fallbacks/authalic_grid.jl`                                                                                                    | `docs/src/api/grid-interface.md`, `docs/src/tutorials/choosing_a_grid.jl`               | SC-011                                                           |
| Spatial queries                   | `src/engine/query.jl`                                                                                                                                                                                                | `docs/src/api/selecting-cells.md`, `docs/src/tutorials/zonal.jl`                        | SC-001, SC-005, SC-006                                           |
| Regions and collections           | `src/engine/partial_grid.jl`, `src/engine/cell_vector.jl`, selection-related parts of `src/dimensionaldata.jl`, `src/engine/region.jl`, `src/engine/region_algebra.jl`                                               | `docs/src/api/selecting-cells.md`, `docs/src/api/boundaries.md`                         | SC-009, SC-010, SC-012, SC-016                                   |
| Mixed-level coverage              | `src/engine/multiorder.jl`, `src/engine/multiorder_set.jl`, `src/engine/multiorder_budget.jl`                                                                                                                        | `docs/src/tutorials/multiorder.jl`, `docs/src/api/selecting-cells.md`                   | SC-007, SC-008, SC-009, SC-012                                   |
| Trees and cursors                 | `src/engine/index_tree.jl`, `src/engine/cursor.jl`, `src/engine/extent_memo.jl`, `src/cap_cached_tree.jl`                                                                                                            | `docs/src/api/grid-interface.md`, `docs/src/architecture.md`                            | SC-004                                                           |
| Cell size and level choice        | `src/sizing.jl`                                                                                                                                                                                                      | `docs/src/tutorials/choosing_a_grid.jl`, `docs/src/api/grid-interface.md`               | SC-003, SC-015                                                   |
| Shipped systems                   | Public system and identifier files for IGeo7, H3, HEALPix, A5, S2, ISEA4R, and Copernicus DEM                                                                                                                        | `docs/src/all_dggs.md`, `docs/src/tutorials/choosing_a_grid.jl`                         | SC-002, SC-003, SC-014                                           |
| Tests used as behavioral evidence | `test/interface`, system conformance tests, `test/systems/crosssystem/multiorder_polygons.jl`, `test/systems/crosssystem/region_algebra.jl`, and selection sections of `test/systems/crosssystem/dimensionaldata.jl` | Not applicable                                                                          | Confirmed hierarchy, region, selection, and mixed-level behavior |




## Deferred scope

The sweep did not assess these files for documentation quality:

- `src/engine/neighborhood.jl`, `src/engine/stencil.jl`, `src/engine/cellfield.jl`,
`src/engine/adjacency.jl`, and `src/engine/halo.jl`
- `src/io/**`, `src/regridding.jl`, `src/chunks.jl`, `src/partitioning*.jl`
- Their matching API pages and tutorials

The sweep read small references to those APIs only when a core docstring depended
on them. It did not create findings for those areas.