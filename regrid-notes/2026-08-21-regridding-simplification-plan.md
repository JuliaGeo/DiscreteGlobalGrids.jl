# Regridding simplification plan

Date: 2026-08-21

## Objective

Reduce the number of overlapping abstractions in `GlobalRegridding` while
retaining the scalability and performance work already present. Dependency
discovery should reuse the hierarchical indexes already provided by GeometryOps
and DiscreteGlobalGrids, be planned once, and then be shared by lazy execution,
scheduling, caching, and prefetch. Cell-tree fallbacks should likewise reuse one
GeometryOps extent protocol and packed spatial index rather than maintaining
package-local trees for each space type.

The intended division of work is:

```text
Raster chunks              Generic chunks               DGGS chunks
DiskArrays ranges +        GeometryOps packed      treeify(original grid),
CR quadtree cursor         FlexibleRTree           stop at chunklevel
              \                 |                 /
               +---------- chunk query ----------+
                          |
                          v
          materialized ChunkDependencyGraph
                          |
               +----------+----------+
               |          |          |
               v          v          v
        lazy execution  scheduling  source-cache policy
               |
               v
        on-demand weights and data
```

Planning is therefore eager about geometry and chunk dependencies, but remains
lazy about weight construction and source data.

## Scope

This plan covers:

- thin adapters onto existing hierarchical chunk spatial indexes;
- simplification of `RasterGrid`'s coordinate-transform and cell-tree layers,
  including task-local Proj transformations;
- a public GeometryOps spherical-extent protocol that ConservativeRegridding
  can use with cap and Cartesian-box tree nodes;
- one packed fallback cell tree for arbitrary/scattered cell subsets;
- consolidation and documentation of the actual cross-package `RegridSpace`
  extension surface;
- construction and ownership of the chunk dependency graph;
- graph-backed lazy execution;
- convergence of the two current dependency-discovery paths;
- one final weight-block construction and reference-weight representation for
  eager and chunked plans;
- removal of DGG-specific output wrapping and target-conversion indirection;
- consolidation of overlapping conservative-geometry caches; and
- simplification of private lazy source-residency and DiskArrays metadata state.

The nearest/bilinear/point-method interface is explicitly deferred. It will be
redesigned separately rather than constrained by this refactor.

`DirectPlan` and `ChunkedPlan` will remain separate. Their repeated fields are
small, while their construction, storage, and execution lifetimes are genuinely
different. `ShapedRegridArray` also remains unless a later array-design change
can preserve lazy multidimensional reads without moving equivalent complexity
into `LazyRegridArray`.

## Phase status and order

Phases 0 and 1 are complete. The current execution point is Phase 1A. Finish
the spatial/geometry layer through Phases 1A, 1B, and 1C before beginning graph
ownership, so Phase 2 benchmarks and APIs exercise the representations intended
to survive:

```text
0 baseline [complete] -> 1 native chunk indexes [complete]
                       -> 1A raster geometry/CRS
                       -> 1B spherical extents/fallback cell tree
                       -> 1C RegridSpace extension surface
                       -> 2 graph ownership -> 3 graph executor -> 4 retirement
                       -> 5 weights -> 6 DGG/API -> 7 geometry cache -> 8 lazy state
```

## Current state

### Lazy dependency discovery is repeated

After Phase 1, `LazyRegridArray` stores the source space's native chunk index
plus source and destination cap vectors. Each requested destination tile calls
`_connectedsource!`, which runs a fresh native-index query to rediscover its
source chunks. Separately, production scheduling materializes
`ChunkDependencyGraph` with the latitude-sorted cap join.

The two broad phases are both conservative but no longer necessarily identical.
The native raster and DGG hierarchies may safely remove relations caused only by
overlap between unused portions of covering caps; the first indexed graph
prototype removed 322 such apparent edges, subject to Phase 2's geometric
no-miss proof. More importantly, the two paths can express different policy: the
production driver's lon/lat-box `refine` was applied to the scheduling graph
only, so the refcount cache released source chunks the independent executor
still read — 8 uncredited demands over 11 chunks
(`regrid-notes/2026-08-21-dag-driver.md:81-119`;
`scripts/copdem_production.jl:640-668`). `refinegraph = false` has been the
default since (`copdem_production.jl:81`). The defect is that `refine` has two
possible owners; Phase 2 gives it one.

### Phase 1 starting state: chunk adapters ignored existing hierarchy

The following was the starting state that Phase 1 removed. Here, "flat"
describes the topology that was exposed to regridding, not the spatial
infrastructure available in its dependencies.

- `DGGChunkTree` was one leaf whose entries were every `(chunk number, cap)` pair.
- `RasterGrid` constructed a `RasterFlatTree` containing every chunk cap, also as
  one leaf.

This duplicated facilities that already existed:

- `RasterGrid` already has an implicit cell-lattice hierarchy and already reads
  its spatial storage ranges from `DiskArrays.eachchunk`. A lightweight
  `ConservativeRegridding.Trees.AbstractCurvilinearGrid` view can carry its
  arbitrary task-local transform into the existing `TopDownQuadtreeCursor`;
  traversal can stop as soon as a node is owned by one DiskArrays chunk. The
  lookup supplies coordinates only. It must not become a second source of
  raster storage chunking.
- GeometryOps `FlexibleRTrees.RTree` is a packed hierarchical spatial index,
  implements `SpatialTreeInterface`, accepts precomputed XYZ extents for
  arbitrary payloads, and preserves original payload indices through packing.
  It is the fallback for generic spaces without an implicit structured index,
  not the structured-raster implementation.
- A DGG chunk is already an ancestor cell. The original complete or
  `PartialGrid` can be traversed through `treeify` and
  `HierarchicalGridCursor`, stopping at `DGGSpace.chunklevel`. Those frontier
  nodes already have subtree-covering extents and correspond in canonical order
  to the compact chunk numbers `1:nchunks(space)`.
- GeometryOps implements `Extents.intersects` between a unit-spherical cap and
  a Cartesian XYZ box in both argument orders. Its contract is conservative:
  it may retain false positives but does not introduce false negatives.

Those adapters prevented hierarchical pruning and led
`chunk_dependency_graph` to maintain a separate latitude-sorted cap join.
Regridding should adapt the existing indexes instead of implementing another
raster tree, another DGG tree, or another generic cap tree.

The generic `chunkextents` path also traversed the flat wrapper and copied its
leaf caps back into a vector. `DGGSpace` avoids that copy with a specialization,
and Phase 1 removed its redundant one-node chunk tree. `RasterFlatTree` remains
only for scattered cell subsets and the legacy `chunktree` compatibility path,
which Phases 1B and 4 address separately.

### Weight construction has two result paths

The public method seam appends entries to `WeightCOO`. Conservative assembly,
however, already produces a sparse CSC matrix. Chunked construction copies that
matrix into COO and sparsifies it again, while eager construction has a separate
`wholeblock(::Conservative)` fast path solely to adopt the CSC directly.

The final reference weights are also represented more than once. `WeightCOO`
always allocates a denominator vector and separately records `hasdenom`;
`WeightBlock` represents the same distinction with `denom === nothing`; and
`CachedBlock` stores another `ref` vector which is a copy of the denominator for
Conservative blocks. Eager application recomputes that reference on every plan
application. Phase 5 makes the final block own it once.

### Spherical extent operations are split across packages

GeometryOps already implements conservative cap–Cartesian-box intersection in
both argument orders and private cap–cap intersection and cap merging.
GlobalRegridding nevertheless owns its own cap-to-XYZ bound, dilated cap
predicate, cap merging and cell-cap construction. DiscreteGlobalGrids wraps the
same private GeometryOps cap operations, while ConservativeRegridding selects a
private cap predicate for every spherical search and therefore cannot accept a
mixed cap/box tree even though GeometryOps can compare those extents.

The missing abstraction is not another regridding tree. It is one public
GeometryOps extent protocol that all three packages can call. Phase 1B promotes
the required cap primitives and uses them to replace the package-local generic
and scattered cell trees with one packed `FlexibleRTree` adapter.

### The actual space extension interface is scattered

`GlobalRegridding.jl` currently names five unexported functions as the external
extension surface, but `DGGSpace` also extends `chunkindex`,
`candidatechunks!`, and `chunkranges`. A method specialized by a different
package is not private in the compatibility sense, even when it remains
unexported. The declarations and contracts are spread across `spaces.jl`,
`conservative.jl`, `discovery.jl`, `lazy.jl`, `executor.jl`, and `api.jl`.
Phase 1C gathers and accurately documents that surface without adding an
interface wrapper or trait hierarchy.

### Conservative destination geometry is cached through several wrappers

`TileCells`, `CachedCellTree`, and destination `CellMemo` overlap in their
responsibility for reusing destination polygons. Raster and DGG cap decorators
cache a different kind of geometry and have measured performance justification;
they should not be collapsed blindly with polygon caching.

### Phase 0 diagnosis: the recorded test baselines used different trees

The local `Manifest.toml` was untracked (`.gitignore:2`) and stale in both git
sources: it paired GeometryOps `repo-rev = 36c853e0` with tree `c30ec910`, the
tree of `02750768`, which predates `intersection_area`; and
ConservativeRegridding `repo-rev = 66ed54c` with tree `873cc732`, the tree of
`6a4b997`, which predates both the windowed assembly and the cached dual DFS.
The revs were bumped without re-resolving. Consequences:

- On this machine `isdefined(GeometryOps, :intersection_area)` is `false`, so
  every Conservative weight block fails at
  `lib/GlobalRegridding/src/intersection_area.jl:31`; no installed GeometryOps
  tree carries the symbol.
- Every suite result recorded "at the new CR pin" since `5c3960b` — including
  2,295 / 0 / 1 broken at `regrid-notes/2026-08-21-merge-train-2.md:110-115` —
  ran against `6a4b997`'s code. `1183e6f`'s 2,202 / 1 broken was a clean-depot
  run at its own pins (CR `6a4b997`, GeometryOps `36c853e0`).
- No note records a 496-pass / 22-error run; that figure has no source.

Upstream state as of 2026-08-22: GeometryOps `36c853e0` is an orphan — the
branch `claude/perf-ladder-predicates` was rebased onto main after #473
(`intersection_area`) merged as `6bdb983e`; its tip `32c60581` is main plus
exactly the two perf-predicate commits, PR #476 open and green.
ConservativeRegridding `66ed54c` is the reachable tip of
`claude/cached-dual-dfs` and the only revision carrying `Trees.split_weight`,
the sparse-CSC `intersection_areas`, `_sparse_from_chunks`, and
`cached_dual_depth_first_search`; CR main (`89ec5ca5`) lacks `split_weight`
and cannot precompile GlobalRegridding. CR `66ed54c` defines methods on
`SpatialTreeInterface.node_extent_is_expensive`, which requires GeometryOps
at or after PR #474; both candidate GeometryOps pins satisfy that.

This package-environment problem was resolved before Phase 1. Commits `f245a0b`
and `3939dd8` also give the subpackages explicit test environments and correct
the conformance test's root source. The Phase 1 verification below is therefore
against the resolved sources, not the stale trees described here.

## Phase 0: reconcile and record the baseline — complete 2026-08-22

1. Repin GeometryOps from `36c853e0` to
   `32c60581afa09f19aeaefee26446d95693ec52c4` in the five `[sources]` blocks
   (`Project.toml:33`, `docs/Project.toml:34`, `test/Project.toml:29`,
   `lib/GlobalRegridding/Project.toml:24`,
   `lib/DiscreteGlobalGridsConformanceTesting/Project.toml:17`; `benchmark/`
   inherits through the workspace). Keep ConservativeRegridding at `66ed54c`
   and GeometryOpsCore at `rev = "main"`. Pin GeometryOps main `6bdb983e`
   instead only to deliberately drop the measured −15 % predicate work
   (`7627ee3`) until #476 merges.
2. `rm Manifest.toml` — mandatory, see above — then `Pkg.instantiate()`.
   Expected trees: GeometryOps `63274824…`, ConservativeRegridding
   `9469fdc5…`. Verified in a clean depot on 2026-08-22: both resolve and
   precompile, and `isdefined(GeometryOps, :intersection_area)`,
   `isdefined(ConservativeRegridding.Trees, :split_weight)` and
   `isdefined(ConservativeRegridding, :cached_dual_depth_first_search)` all
   hold.
3. Record both tree SHAs in the commit message; the Manifest is untracked, so
   the message is the only record (`7627ee3` is the precedent).
4. Run the complete `GlobalRegridding` tests and the DGG cross-system
   regridding and acceptance tests. Reference counts from `1183e6f`:
   GlobalRegridding 2,202 / 1 broken; `crosssystem/regrid.jl` 116;
   `regridding_conservation.jl` 76 / 12 broken; `systems/CopernicusDEM`
   16,162 / 3 broken; `regrid_acceptance.jl` 22. The new counts are the first
   at the real `66ed54c` tree and replace every earlier figure.
5. Create the baseline artifacts later phases gate on: graph construction
   time and bytes for the Copernicus DEM to IGeo7 graph (bytes has never been
   recorded); a test for `_TILE_CELL_CACHE_MAX`; a committed production
   Conservative benchmark; and a committed re-measurement of the chunked
   weight-block round trip's peak RSS (Phase 5).

Exit criterion: reachable pins, a re-resolved Manifest, a green regridding
test gate, and the four artifacts, all recorded before any simplification
commit. The pin commit is separate from every structural change.

## Phase 1: reuse the existing chunk indexes — complete 2026-08-22

Introduce one private chunk-candidate interface for spatial discovery. Phase 1
uses it for lazy single-tile queries; Phase 2 will move whole-graph construction
onto it only after the production performance and geometric-miss gates pass. It
should express the operation the regridder needs, not require every index to
expose the same internal extent type:

```julia
index = chunkindex(src_space)
candidatechunks!(out, index, dstcap; radius)
```

The current `chunktree` contract assumes that internal and leaf extents are all
`SphericalCap`s. That remains appropriate for cell trees, but it should not force
a packed XYZ-box R-tree to impersonate a cap tree. Chunk indexing may be a
private query abstraction over conservative extents.

### Structured raster spaces: DiskArrays ranges over the existing quadtree

Raster storage chunking comes from `DiskArrays`, not from the coordinate
lookup. `RasterGrid(A)` already records the spatial ranges returned by
`DiskArrays.eachchunk(A)` as `xchunks` and `ychunks`; explicit `chunks = ...`
on the dimension-only constructor is the test/manual equivalent. Preserve that
single ownership rule.

Expose the raster geometry through a lightweight
`ConservativeRegridding.Trees.AbstractCurvilinearGrid` view with:

- the existing X/Y cell-edge vectors;
- the raster's actual native-to-unit-sphere transform, including closures and
  other task-local callables;
- cell numbering in the source array's fastest-dimension order; and
- range extents delegated to the existing conservative raster rectangle-cap
  calculation.

Use ConservativeRegridding's existing `TopDownQuadtreeCursor` over that view.
For a cap query, descend the implicit cell lattice and stop as soon as the
cursor's cell rectangle lies wholly inside one `(xchunk, ychunk)` pair. Emit
that pair's DiskArrays chunk number. If a quadtree node crosses a storage chunk
boundary, keep descending; at the bottom, map the touched cells back to their
owning chunks and deduplicate. This makes index construction independent of the
number of chunks and does not materialize one extent per tiny chunk.

The view is an adapter to an existing quadtree, not a second raster tree. It may
later replace `RasterCellTree` as well, but Phase 1 changes chunk discovery only
so conservative cell-tree behavior and its memoization remain independently
testable.

### Generic spaces: `FlexibleRTrees`

For a `RegridSpace` without a structured chunk index, use GeometryOps' existing
packed R-tree rather than adding a new generic hierarchical cap tree.

1. Materialize the source chunk caps once.
2. Convert each cap to an outward-rounded Cartesian XYZ bounding box with one
   small, tested `cap_xyz_extent` helper. GeometryOps does not currently define
   `Extents.extent(::SphericalCap)`.
3. Construct `RTree(HPR(), caps; extents = boxes)`. Caps remain the payload, and
   the R-tree's original-index mapping keeps leaf hits equal to chunk numbers.
4. Query the tree with a destination cap whose angular radius is enlarged by
   `support_radius(method, src_space)`. Use GeometryOps' conservative mixed
   cap/box intersection during descent.
5. Apply the existing exact cap-cap support predicate to leaf candidates and
   sort the resulting chunk numbers.

The cap/box predicate is defined in both argument orders and documented to
permit false positives but no false negatives (GeometryOps
`src/utils/UnitSpherical/cap.jl:129-142`); it dispatches on the exact
`(:X, :Y, :Z)` extent keys. A query radius at or above `pi` is already a
whole-sphere query (`chord = 2 sin(min(r, pi) / 2)`, `cap.jl:137`). Only the
empty collection needs a fallback: `RTree` rejects empty input
(`FlexibleRTrees/types.jl:86`) and handles a singleton.

API details the helper must respect:

- `extents` is a concrete `Vector{<:Extents.Extent}`, one per item, and the
  tree takes ownership of it (`types.jl:79, 91-93`);
- a bare cap is not a predicate (`SpatialTreeInterface.jl:67-70`); query with
  `depth_first_search(Base.Fix1(Extents.intersects, cap), tree)`;
- `GeometryOps.query` is ambiguous — `SpatialTreeInterface` and
  `FlexibleRTrees` export different `query`s — so qualify it;
- `Extents.extent(::SphericalCap)` returns `nothing` silently through Extents'
  universal fallback rather than erroring; `cap_xyz_extent`'s tests should pin
  that nothing reaches for it;
- ConservativeRegridding's spherical broad phase is hard-wired to cap–cap
  `UnitSpherical._intersects` (`src/regridder/intersection_areas.jl:259-263`),
  so the box R-tree is only ever the private chunk index, never a cell tree.

`RasterFlatTree` remains only where it is genuinely a cell-tree fallback for
arbitrary scattered cell indices. Its chunk-lattice use is removed; raster
chunk discovery consumes only the DiskArrays-derived ranges.

### DGG spaces: expose the existing tree's chunk frontier

For a normally chunked `DGGSpace`, traverse the original grid tree and stop at
the chunk level:

```julia
cursor = treeify(space.grid)
# Descend through the existing cursor only until cursor.level == space.chunklevel.
# Emit that frontier node as one regridding chunk instead of descending to cells.
```

This distinction matters. `node_extent(system, cell)` is explicitly the cap of
the entire subtree and must cover every descendant's geometry through the
system's maximum level. `HierarchicalGridCursor` uses that covering extent above
the original grid's leaf level. It switches to the private tight `cell_cap` only
at literal cell leaves, where there are no stored descendants left to cover.

Constructing a new `PartialGrid` at `chunklevel` would incorrectly reclassify
chunk ancestors as literal leaves and select the tight boundary cap. Traversing
`treeify(space.grid)` preserves their role as internal ancestor nodes, so no
extent-correcting cursor wrapper is needed for ordinary hierarchical systems.

At the chunk frontier:

- the ancestor IDs occur in the same canonical order as `space.chunkids`;
- each cursor's `node_indices` equals the corresponding `space.ranges[chunk]`;
- compact chunk numbers can be assigned by traversal order, by searching
  `space.chunkids`, or via `chunkat(first(node_indices(cursor)))`; and
- sparse `PartialGrid` nodes may retain the cursor's existing sound
  selected-descendant tightening.

Complete grids, rooted partial grids, and unrooted sparse partial grids already
follow this path through their existing cursor implementations. If
`chunklevel == level(space.grid)`, each chunk is one literal cell and its tight
leaf cap is correct. If a rooted grid starts deeper than `chunklevel`, or the
system does not provide sorted subtrees, `DGGSpace` already produces one chunk;
use the trivial single-chunk path. Empty grids retain their existing full-sphere
fallback.

Two cases the cursor path does not cover as written:

- `treeify(space.grid)` is not always a `HierarchicalGridCursor`. CopernicusDEM
  returns its existing `MemoBlockCursor`, whose rectangular spatial hierarchy
  is essential: a sparse production grid has 26,475 selected level-0 tiles, so
  a synthetic hierarchical root would scan that flat fanout for every query.
  For its level-0 frontier, traverse `treeify(levelgrid(sys, chunklevel))` and
  filter leaf IDs through the sorted `space.chunkids`. The block cursor's leaf
  pad has the same descendant-geometry covering role as `node_extent`, and no
  source pixels or second tree are materialized. Higher-level rooted cases can
  use the ordinary cursor path. This is the same rule as for rasters: adapt the
  system's existing spatial lattice rather than forcing one universal cursor.
- The single-chunk case is the common one for rooted subtrees, not the rare
  one: the automatic `chunklevel` is chosen from cell counts independent of the
  root, so `DGGSpace(subtree(IGeo7, L2, 4))` gives `chunklevel = 0` and
  `nchunks = 1`. Dispatch on `nchunks(space) == 1` and `_ischunked`, not on
  `chunklevel` or `isempty(chunkids)` (`_wholechunk` stores `ID[]` with
  `ID = Any` when the grid has no system, `src/regridding.jl:50-56`).

The sparse tightening at `src/engine/cursor.jl:232-257` is sound on this path
because the stored ids are grid-level cells with nothing below them. It would
have been unsound on a `PartialGrid` built at `chunklevel`, whose stored ids
are ancestors and whose `cell_cap` under-covers overhang (IGeo7 L0: `cell_cap`
radius 0.652 against `node_extent` 0.783) — which is why that construction is
rejected above.

After cutover, remove `DGGChunkTree`; do not add a second DGG hierarchy or a
chunk-specific extent adapter.

### Required laws

- Every chunk number appears exactly once in an index.
- Every chunk-frontier extent covers the actual geometry of every stored cell in
  that chunk, including descendant overhang in aperture-7 systems.
- Every internal extent covers the actual geometry of all descendant chunks. It
  need not contain their deliberately over-covering cap objects.
- Broad-phase queries have no false negatives, including polar, antimeridian,
  tiny-cap, and whole-sphere cases.
- Generic R-tree exact-filtered candidates equal brute-force
  `DilatedIntersects` results. Raster and DGG hierarchy queries may safely omit
  pairs caused only by overlap between two chunks' unused cap over-coverage,
  but may omit no geometrically possible contribution.
- Raster results are independent of DiskArrays chunk sizes and quadtree split
  shape. Generic results are independent of R-tree packing and leaf capacity.
  All results are independent of thread count.

Exit criterion: raster queries use DiskArrays chunk ranges over the existing CR
quadtree, generic queries use `FlexibleRTrees`, DGG queries expose the
chunk-level frontier of `treeify(space.grid)`, and all satisfy the covering and
no-false-negative laws without new general-purpose tree implementations.

Completed in Phase 1. The production graph builder deliberately remains on its
latitude-sorted cap join; the first indexed-graph measurement below did not
justify combining that Phase 2 change with the hierarchy cutover.

## Phase 1A: simplify `RasterGrid` geometry and CRS transformations

`RasterGrid` currently mixes four concerns in one file:

1. extracting dimensions and DiskArrays chunk ranges;
2. converting native raster coordinates onto the unit sphere;
3. constructing spherical cells and caps; and
4. implementing two separate spatial trees over the same raster lattice.

That has produced both exact duplication and generic-looking APIs whose
contracts are actually geographic or unit-spherical. Resolve this before the
Phase 2 raster graph benchmarks, so those measurements exercise the intended
long-term geometry layer.

### Be precise about the supported geometry

Global regridding is currently unit-spherical: `RasterGrid` returns
`UnitSphericalPoint`s, its manifold is a unit `Spherical`, and chunk discovery
uses `SphericalCap`s. It is therefore not yet manifold- or output-CRS-generic.

The useful generic boundary is narrower and should be stated explicitly:

- the raster's **native CRS** may be geographic or projected;
- a forward transform maps native `(x, y)` coordinates to the unit sphere;
- its inverse maps a unit-sphere point back to native coordinates; and
- all cell polygons, spatial indexes, support radii, and intersections remain
  unit-spherical.

Rename ambiguous storage and helpers accordingly. In particular, prefer
`native_to_unit_sphere` and `unit_sphere_to_native` over `transform` and
`inverse`. A helper such as `chartspacing` that returns radians of spherical arc
must say so in its name and documentation. Geographic-only types and helpers
must also say `Geographic` rather than appearing to be valid for every chart.

Full manifold-generic raster regridding would require replacing cap-based
discovery and spherical intersection assumptions throughout the package. It is
not part of this refactor and should not be implied by the `RasterGrid` API.

### Reuse GeometryOps coordinate transformations

Delete the local `LonLatToSphere` and `SphereToLonLat` transformations. They
duplicate GeometryOps' existing:

```julia
GeometryOps.UnitSpherical.UnitSphereFromGeographic()
GeometryOps.UnitSpherical.GeographicFromUnitSphere()
```

Those types already implement the `CoordinateTransformations.Transformation`
interface and `inv`, and are already used elsewhere in this repository.
Specialize the raster chart metadata directly on the GeometryOps types:

- geographic X period and latitude bounds;
- spherical step bounds in radians; and
- optional geographic edge tables.

Internally normalize forward transforms to the CoordinateTransformations
single-coordinate calling convention, `transform((x, y))`. Preserve the
currently documented two-argument closure form with one construction-time
adapter, not conditional dispatch at every vertex. This allows plain callables,
GeometryOps transformations, composed transformations, and Proj-backed
transformations to pass through the same raster geometry code.

The default remains geographic longitude/latitude in degrees for compatibility,
but documentation must state that assumption. Projected X/Y dimensions must
carry an explicit CRS or explicit native-to-sphere transformation; do not infer
that arbitrary coordinates named `X` and `Y` are degrees.

### Use Proj.jl's API for task-local transformations

PROJ transformation objects and contexts may not be used concurrently. A Julia
task may migrate between OS threads, so thread-indexed clones are the wrong
ownership unit. Give each Julia `Task` its own context and its own transformations,
cached in task-local storage and reused for that task's geometry work.

Use only Proj.jl's exposed Julia API. Do not add a direct `ccall`, do not call
`PROJ_jll`, and do not reproduce a binding in GlobalRegridding. Proj.jl 1.9
provides the required operations:

```julia
ctx = Proj.proj_context_clone()

forward = Proj.Transformation(
    Proj.proj_clone(template.pj, ctx),
    template.direction,
)

reverse = Proj.Transformation(
    Proj.proj_clone(template.pj, ctx),
    inv(template.direction),
)
```

This exact sequence was validated locally against Proj.jl 1.9: the clone and
template returned identical coordinates. `proj_context_clone()` is preferred
to a fresh context because it inherits the configured default context. Cloning
the pipeline also preserves its selected operation, axis normalization, area of
use, and grid configuration instead of reconstructing it from incomplete CRS
metadata.

Implement this in a Proj extension. Core GlobalRegridding exposes a small
task-local transform-pair preparation hook whose default returns immutable or
ordinary Julia callables unchanged. The Proj specialization:

1. obtains the calling task's cache;
2. creates one cloned context on the first use of a template;
3. clones distinct forward and reverse `Proj.Transformation`s into it;
4. returns the cached pair on later calls from the same task; and
5. owns cleanup as one idempotent resource.

A `Proj.Transformation` owns and finalizes its cloned transformation pointer;
the context remains caller-owned. Cleanup order is therefore mandatory:

1. finalize or otherwise release both cloned transformations;
2. call `Proj.proj_context_destroy(ctx)`; and
3. mark the owner closed so repeated finalization is harmless.

Construction failures must release any successfully created clone before
destroying the context. No cloned transformation or context may be stored on the
shared `RasterGrid`; it stores only the immutable/native chart description or
template. No lock is needed because the cache and every clone belong to one
task. Network and search-path policy are configured before cloning and remain
outside the raster hot path.

For projected rasters, the conceptual forward pipeline is:

```text
native CRS --Proj--> geographic longitude/latitude --GeometryOps--> unit sphere
```

and the inverse reverses those steps. Proj constructors already accept strings
and `GeoFormatTypes` CRS values. Use `always_xy = true` when constructing raster
pipelines so their coordinate order agrees with the raster's X/Y axes. A later
multi-package Rasters/Proj extension may obtain the native CRS from raster
metadata, but the core constructor must continue to accept an explicit
transformation without depending on Proj or Rasters.

### Remove duplicate raster lattice and cap implementations

`RasterCellTree` and `RasterGridView` describe the same curvilinear lattice.
The latter already implements ConservativeRegridding's
`AbstractCurvilinearGrid` interface and is used by the Phase 1 chunk quadtree.
Make it the single structured raster geometry adapter:

- whole-space and rectangular `celltree` requests return restricted
  `TopDownQuadtreeCursor`s over `RasterGridView`;
- global cell numbering continues to follow storage dimension order;
- scattered cell positions may retain `RasterFlatTree` until Phase 1B replaces
  the generic/scattered fallback; and
- remove `RasterCellTree` once its leaf sizing, split-weight, and cell-access
  contracts are covered by the CR cursor.

`raster_tree_memo.jl` exists only around `RasterCellTree` and overlaps CR's
cached dual-depth-first traversal. Delete it if the existing raster benchmark
confirms no material regression. If a cache is still earned, adapt one small
cache to the CR cursor rather than retaining a second tree implementation.

The spherical cap machinery needs the same separation:

- the general curvilinear path should delegate range extents to
  ConservativeRegridding's existing `cell_range_extent`, which visits the
  transformed perimeter vertices through `getvertex`;
- the current finite `_boxcap` sampling must not be presented as conservative
  for an arbitrary projected transform;
- retain and clearly name the geographic wide-longitude/polar specialization
  only where it proves a tighter sound cap or a material performance win; and
- keep geographic sine/cosine edge tables only as a measured optimization,
  renamed so they do not masquerade as a generic chart facility.

Do not remove the wide-box fallback merely for line-count reduction. The generic
CR cap builder must first pass full-longitude, polar, antimeridian, projected,
and large-index-rectangle coverage laws. If those laws expose a generally useful
gap, fix it in ConservativeRegridding or GeometryOps rather than adding another
general cap builder here.

### Required tests and exit criterion

Add tests that establish:

- the default geographic transform is the GeometryOps type and retains current
  cell rings, lookups, and interpolation behavior;
- the compatibility adapter accepts existing two-argument closures while the
  internal contract is one-coordinate;
- a projected native CRS round-trips through the unit sphere within the
  projection's accuracy;
- concurrent Julia tasks obtain distinct Proj contexts and transformation
  clones, reuse them within a task, and return identical numerical results;
- task cleanup releases transformations before their context, including partial
  construction failure;
- no GlobalRegridding source contains a direct `ccall` to PROJ;
- general curvilinear and geographic-specialized node extents cover every
  boundary point owned by their leaves; and
- eager/lazy results, zero-source-read construction, allocation, and the
  small-spatial-chunk traversal benchmark do not regress materially.

Exit criterion: GeometryOps owns geographic/unit-sphere conversion, Proj.jl owns
all PROJ calls, every Proj object is task-local with explicit lifetime ordering,
`RasterGrid` names its unit-sphere contract honestly, and one CR curvilinear-grid
adapter replaces the duplicate structured raster trees without weakening extent
coverage.

## Phase 1B: converge spherical extents and fallback cell trees

Phase 1 reused native hierarchy for **chunk** discovery. Conservative cell-pair
discovery still has two package-local fallbacks: `CellCapTree` for a generic
restricted cell set and `RasterFlatTree` for scattered raster positions. Both
materialize per-cell caps and then implement another spatial tree. GeometryOps'
packed `FlexibleRTrees.RTree` already provides the required hierarchy and stable
original payload numbering, but ConservativeRegridding's spherical path
currently insists that both sides expose cap extents.

Resolve the protocol at its owner before replacing the trees.

### Promote the GeometryOps cap operations actually used across packages

Add public GeometryOps operations for:

- conservative cap–cap intersection through `Extents.intersects`, delegating to
  the existing numerically robust unit-spherical predicate;
- the outward-rounded Cartesian `(:X, :Y, :Z)` extent of a `SphericalCap`;
- cap containment and cap merging, which DiscreteGlobalGrids and
  GlobalRegridding currently expose or reproduce around private GeometryOps
  methods; and
- cap dilation by an angular radius if retaining that named operation is clearer
  than constructing a new cap at every caller.

Adjudicate whether the Cartesian bound should be
`Extents.extent(::SphericalCap)` or an explicitly named unit-spherical helper.
Do not overload `Extents.extent` if that would imply that a cap is ordinary
Cartesian geometry rather than a spherical query object. Whichever spelling is
chosen becomes the one implementation of the outward-rounding rules currently
in `cap_xyz_extent`.

The exact cap-lens area calculation used only by lazy wave costing may also move
to GeometryOps if it is a generally useful unit-spherical operation. It does
not block the tree convergence; retaining one task-specific cost estimator is
preferable to publishing an API solely to remove lines here.

After the public operations exist:

- `DilatedIntersects` delegates to a dilated cap plus the public robust cap–cap
  predicate rather than recomputing spherical distance;
- DiscreteGlobalGrids' `intersects_cap`, `cap_contains`, and `merge_caps`
  wrappers either become thin public aliases with no private calls or disappear
  where callers can use GeometryOps directly; and
- GlobalRegridding removes its duplicate `_mergecap` implementation and uses
  the GeometryOps cap builder for fallback leaf geometry.

### Let ConservativeRegridding compare mixed extent types

ConservativeRegridding currently chooses
`GeometryOps.UnitSpherical._intersects` for every spherical manifold and
`Extents.intersects` otherwise (`intersection_areas.jl:270-274`). Once
GeometryOps defines cap–cap intersection in the extent protocol, use
`Extents.intersects` uniformly. The dual search can then compare:

- cap against cap;
- cap against Cartesian XYZ box; and
- Cartesian XYZ box against Cartesian XYZ box.

This is an extent-protocol change, not a change to intersection geometry. Cell
polygons remain unit-spherical and the clipping operator is unchanged. Add the
mixed-extent combinations to ConservativeRegridding's dual-search tests before
using them here.

The multithreaded frontier's `pair_weight` has a cap-specific overlap estimate.
A packed box tree would take its generic subtree-size fallback today. Either add
a sound cheap mixed cap/box work estimate or establish by benchmark that the
generic estimate does not materially harm task balance. Do not trade the custom
tree's source complexity for an unmeasured production scheduling regression.

### Replace the generic and scattered cell-tree implementations

Introduce one small private adapter over `FlexibleRTrees.RTree` containing:

- the owning `RegridSpace` for `getcell` and manifold access;
- the packed tree;
- stable original cell positions as leaf payloads; and
- outward-rounded Cartesian extents derived by GeometryOps.

Use it for arbitrary cell-index subsets for which a space cannot return a
restricted native cursor. Then:

- remove `CellCapTree`;
- remove `RasterFlatTree`'s scattered-cell role after Phase 1A has moved
  rectangular raster ranges to `TopDownQuadtreeCursor`;
- let arbitrary non-windowable DGG ranges use the same fallback rather than a
  DGG-specific tree; and
- retain `RasterFlatTree` only as the temporary legacy `chunktree` adapter until
  Phase 4 removes that surface.

Do not replace the native structured cursors. Raster rectangles remain CR
quadtree cursors, DGG ranges with native subtree/window support remain DGG
cursors, and the packed tree is only the unstructured fallback.

### Required laws and exit criterion

- Public cap–cap results are identical to GeometryOps' existing robust private
  predicate, including nearly coincident centres and radii within a few ULP of
  contact.
- Cap–box and box–cap comparisons retain their no-false-negative contract.
- Cartesian cap extents contain extrema from tiny, polar, antimeridian,
  hemisphere, and whole-sphere caps with outward rounding.
- Packed fallback leaf positions equal the requested original cell positions,
  independent of packing and node capacity.
- Conservative blocks built from native cap cursors and packed box fallbacks
  are numerically identical.
- Serial and threaded candidate pairs and weights remain deterministic.
- The committed Conservative production benchmark shows no material regression
  in candidate count, task balance, construction time, or peak memory.

Exit criterion: GeometryOps owns the reusable spherical extent primitives,
ConservativeRegridding accepts mixed spatial extent types through one public
predicate, and GlobalRegridding owns no general-purpose fallback cell tree
other than a thin adapter over `FlexibleRTrees`.

## Phase 1C: consolidate the `RegridSpace` extension surface

Gather the cross-package declarations and contracts in `spaces.jl`, grouped by
responsibility rather than by the file that happens to consume them:

```text
cell geometry:       ncells, getcell, celltree, subtree
chunk geometry:      nchunks, cellindices, chunkextent,
                     chunkindex, candidatechunks!
array storage:       chunkranges
point/chart methods: cellat, cellcentroid, hascellchart
output/resolution:   destinationdims, dimsource, _asspace
```

Mark load-bearing qualified hooks `public` even when they remain unexported.
In particular, `chunkindex`, `candidatechunks!`, and `chunkranges` are external
interfaces because DiscreteGlobalGrids specializes them. Update the module
comment and `RegridSpace` docstring so they no longer call these methods private.

Do not introduce `AbstractChunkIndex`, `SpaceTraits`, a `PlanCore`, or a wrapper
object merely to make the list look uniform. Dispatch on the query operation is
the useful abstraction; each native index keeps its own representation.

Keep the semantic distinction between:

- `cellindices(space, chunk)`, which names the cells a chunk owns; and
- `chunkranges(space, chunk, spatialsize)`, which names rectangular array ranges
  that can be read from storage in one operation.

`chunktree` remains only as a compatibility bridge until Phase 4. At that point,
remove it if the repository and known downstream consumers use the query seam;
do not translate the new index protocol back into a legacy cap tree.

Exit criterion: one documented qualified interface explains every method
DiscreteGlobalGrids extends, declarations live with the space contract, and no
new interface type or forwarding wrapper has been added.

## Phase 2: build and own one authoritative dependency graph

Build destination-major graph rows through the Phase 1/1C chunk-query seam:

```julia
radius = support_radius(method, src_space)
index = chunkindex(src_space)
for (d, dstcap) in pairs(chunkextents(dst_space))
    row = Int[]
    candidatechunks!(row, index, dstcap; radius)
    rows[d] = row
end
dependencies = ChunkDependencyGraph(rows, nchunks(src_space))
```

Use one source index and independent per-destination cap queries rather than a
mixed cap/box dual-tree traversal. This aligns with the destination-row CSR,
parallelizes deterministically over destination chunks, and handles nonzero
support by dilating only the query cap. A generic Raster/R-tree dual traversal
would require a new angular-support rule for box/box node comparisons.

Before deleting the latitude-sorted cap join, establish the indexed builder's
no-false-negative property and compare whole plan-build time and peak memory —
not graph time alone; the recorded graph build is 0.12 s of a ~1 s plan build
whose remainder is the destination space (`scripts/copdem_production.jl:701-702`)
— on:

- the production Copernicus DEM to IGeo7 pair, whose recorded cap-join baseline
  is 0.122 seconds for 66,178 by 26,475 chunks and 326,392 edges, full graph,
  zero support radius, 4 threads
  (`regrid-notes/2026-08-21-copdem-lazy-source.md:110-113`). No note records
  graph bytes; that is a new measurement, not a comparison;
- a raster with very small spatial chunks;
- complete, rooted, and sparse DGG subsets;
- polar and antimeridian-heavy subsets; and
- highly nonuniform chunk coverage.

The first provisional production measurement on 2026-08-22 is a reason to keep
this gate. With four Julia threads, the existing latitude join produced 326,386
edges in a 3,352,520-byte graph with 0.0594 s median construction and 10,654,960
minimum allocated bytes. A direct Phase 1 index-row prototype produced 326,064
edges in 3,349,944 bytes but took 0.1500 s median and allocated at least
13,973,280 bytes. The 322 removed edges appear to be pairs retained only by cap
over-coverage, but that is not yet a correctness proof. The prototype was backed
out; Phase 2 must add actual-cell no-miss instrumentation and recover the graph
construction performance before cutover.

For the generic R-tree path, require adjacency identity with the flat
brute-force cap predicate. For raster and DGG hierarchies, treat that flat
relation as an upper bound: hierarchical traversal may safely remove edges
caused only by unused cap over-coverage. Verify small cases against actual cell
geometry and all cases through covering laws, numerical equivalence, and
graph-miss instrumentation.

Once the indexed builder passes those gates, remove `_caplatitude`, the sorted
latitude arrays, the global-radius band, and `_fillrow!`. There should be one
candidate query implementation, used by whole-graph construction.
`connectedchunks` (`discovery.jl:88-96`, not exported and not `public`) is a
third relation — it answers `sourcesof`'s question from different code with no
plan and no narrow phase — and goes in Phase 4; its oracle uses in
`lib/GlobalRegridding/test/test_lazy.jl:132-137, 148, 192, 215, 233, 341-342`
are rewritten against `sourcesof(dependencies(plan), d)`.

The resulting graph belongs to the logical plan because it depends on space
geometry and the method, not on source values or non-spatial slices. Reusing a
plan must reuse the graph. Production currently builds one global scheduling
graph while constructing a fresh one-destination lazy regrid for each column;
blindly adding graph construction to `ChunkedPlan` would therefore build
thousands of redundant small graphs. The whole per-column cost of the current
pattern is known: each of the 66,178 per-column regrids
(`copdem_production.jl:781-789, 981`) folds all 26,475 source caps into one
bounding cap (`src/regridding.jl:210-212`) and allocates a 26,475-byte
`emptymemo` (`lazy.jl:268`); `srccaps` is shared by reference.

Define explicit graph reuse before changing the constructor. Either production
uses one global plan/lazy array, or smaller plans accept a validated immutable
dependency object or row view. Validation must cover source and destination
chunk numbering, chunk counts, support radius, the narrow-phase tag, and the
space identity or geometry fingerprint needed to prevent accidental reuse with
different spaces.

The form chosen for Phase 3 is the row view. The global scheduling plan owns
the graph; each per-column plan is constructed with
`restrict(dependencies(globalplan), d)` — `ndst == 1`, the same `nsrc`, radius,
source fingerprint and narrow-phase tag, sharing the parent's CSR slice (one
row copy per column, mean about five entries, maximum 360). A per-column
rebuild would also call `refine` with the local destination number `1`
instead of the global `d`, so it is a correctness trap, not only waste. One
global lazy array — a `DestTiling` of 66,178 tiles (`lazy.jl:66-70`), one
shared `PerChunk` across workers where `budget` (`lazy.jl:76`) is per-worker
today — remains the later option and must be measured before adoption.

Plan construction may build chunk extents and adjacency, but must still:

- read no source data;
- perform no provider `HEAD` or other network metadata requests;
- build no weight blocks; and
- produce deterministic adjacency independent of thread count.

The graph remains a conservative relation: an edge means that the executor may
read the source chunk for that destination chunk. Before Phase 3, the
production driver demonstrated that a geometrically sound narrow phase applied
to the scheduling graph alone is narrower than the executor's independent
broad-phase reads (`regrid-notes/2026-08-21-dag-driver.md:81-119`). The remedy
is ownership, below, not disabling refinement.

### `refine` is set once

A narrow phase is an argument to plan construction and to nothing else.
`plan_regrid(...; lazy = true, refine = f)` stores the graph it builds on the
`ChunkedPlan`; `chunk_dependency_graph(dst_space, src_space; radius, prefilter)`
loses its `refine` keyword (`chunkgraph.jl:361-368`), and
`chunk_dependency_graph(plan; refine)` (`chunkgraph.jl:370-372`) becomes
`dependencies(plan)`, an accessor that builds nothing. Today the only
evaluation site is `chunkgraph.jl:461`; the struct records `(nsrc, ndst,
radius)` and no narrow phase (`chunkgraph.jl:52-64`).

Laws:

- Scheduler, executor, prefetcher and validator read one graph object per plan.
- The graph stamps chunk counts, radius, a space fingerprint, and a narrow-phase
  tag — the tag, not the callable, so it stays serializable and free of function
  type parameters (`chunkgraph.jl:39-41`, `Base.zero` at `:209-211`). Two
  objects claiming one relation must agree on all four, or plan construction
  fails.
- A `refine` is a pure, deterministic, thread-safe function of `(d, s)` — it
  runs inside the spawned row blocks (`chunkgraph.jl:396-408`). It may consult
  geometry, never data; emptiness stays in the executor (`lazy.jl:592-596`).
- `refine(d, s)` may return `false` only when no cell of destination chunk `d`
  lies within `support_radius(method, src_space)` of any cell of source chunk
  `s`. The lower bound is cell geometry dilated by that radius; the cap join is
  only the upper bound it narrows from.

Recording the narrow phase on the graph alone would make a violation
detectable but not impossible — two graphs with different stamps could still
reach different consumers, which is the production bug. A keyword that flows
into graph construction and nowhere else is what exists today. The plan is the
owner; the stamp is the witness; the keyword is the only door.

Exit criterion: no path supplies a narrow phase after a plan exists, a plan
exposes exactly one dependency graph, and any disagreement about the relation
is a construction-time error.

### Graph representation

Keep the bidirectional CSR representation and chunk-number queries needed by the
executor and production policy:

- `sourcesof(graph, destination_chunk)`;
- `consumersof(graph, source_chunk)`;
- `sourcedegree` and `consumerdegree`; and
- source/destination chunk counts.

Whether the complete `Graphs.AbstractGraph` compatibility surface remains in
core is a later packaging decision. It must not block executor integration. If
it is moved, move only the ecosystem adapter; retain the core relation and CSR
builder.

### Scaling

Graph vertices represent spatial chunks only. Temporal and other non-spatial
chunks remain in the executor's `othergroups` and do not multiply graph size.
A dataset with small spatial chunks and large temporal chunks therefore has many
spatial vertices but reuses the same graph over every temporal group.

Track graph construction time and bytes as first-class planning statistics.
CSR storage should remain `O(nsource + ndestination + nedges)`, keeping the
existing `Int32` edge indices and range check (`chunkgraph.jl:380-383`). Any extent metadata retained for wave
costing must be immutable and shared with the graph rather than replicated per
worker, destination column, or non-spatial slice.

Exit criterion: indexed construction has no geometric false negatives, the old
cap join is removed after its correctness and performance gates, and production
plus reused `ChunkedPlan`s share rather than replicate the authoritative
dependency relation.

## Phase 3: make lazy execution consume the graph

Replace `_connectedsource!`'s repeated native-index query with graph adjacency.

Before Phase 1 that query was a `1 × nchunks` linear scan in practice and cost
0.5 ms per production column at 360 tiles against roughly 10 core-seconds of
work (`regrid-notes/2026-08-21-polar-profile.md:124-136`). Phase 1 replaced it
with each space's native hierarchy, but it still repeats the query. The payoff
of this phase is one relation owned by the plan and consumed by scheduler and
executor alike, not merely query throughput.

For a destination tile aligned with one space chunk:

```julia
sources = sourcesof(plan.dependencies, destination_chunk)
```

For a derived executor tile spanning several destination chunks, use the
existing `DestTiling.capsof[t]` mapping and take the sorted union of the
corresponding graph rows:

```julia
sources = union(sourcesof(graph, d) for d in tiling.capsof[t])
```

This has the same conservative granularity as the current tile-to-space-chunk
mapping. Exact graph rows for executor-created subtiles can be considered later
if coarse destination chunks cause material over-reading. Add a benchmark and
reported over-read ratio for explicit and budget-derived destination tilings.

After graph integration:

- remove `LazyRegridArray.srctree`;
- remove per-read geometric dependency discovery;
- retain `srccaps` and `dstcaps` until Phase 4 moves wave costing
  (`_wavesize`, `_blockcosts!`, `lazy.jl:469-485, 502-518`) onto the
  dependency object; they are not held solely for discovery;
- retain `knownempty` filtering, because emptiness depends on data and
  non-spatial groups rather than geometry;
- retain deterministic ascending source-chunk application order; and
- ensure graph adjacency is reused by repeated reads and slices.

Once `_connectedsource!` reads graph rows, a dropped edge is no longer a late
demand — it is a source contribution that is never built, absent from
numerator and denominator alike, and the uncredited counter cannot see it. A
narrow phase therefore becomes a correctness gate rather than a residency
knob. Before enabling one, prove on small spaces that for every destination
chunk `d`, `sourcesof(dependencies(plan), d)` contains every source chunk
owning a structurally nonzero column of the eager whole-domain block
(`api.jl:210-217`) over `d`'s cells, and sandwich the refined graph between
that set and the brute-force `prefilter = false` cap reference
(`chunkgraph.jl:305-306`, `test_chunkgraph.jl:28`) over the `ToyLonLatSpace`
fixtures at several radii plus polar and antimeridian chunkings. The lon/lat
box narrow phase already holds this certificate against an independently
built exact adjacency — it kept all 186,069 edges
(`dag-driver.md:99-104`, `copdem_production.jl:650-651`) — and it must become
a runnable test rather than a note. Eager/lazy bit-identity over the same plan
remains the generic net: a `refine` that drops a true pair changes lazy
values.

Production obtains its relation from one plan object: the global scheduling
plan owns the graph, and each per-column regrid is constructed with a row view
of it, never a rebuild. Keep the run-time `uncredited` counter
(`scripts/copdem_policy.jl:223-231, 349-352`) and redocument it: after this
phase it is provably zero given that cache and executor hold the same graph
object, so it stops checking `refine` and becomes the structural assertion
that they do — across resume-retirement (`copdem_production.jl:1261-1271`),
the prefetcher, and the per-column plans. A long production run may still
serve an unexpected demand transiently, count and report it, and fail its
closing invariant (`copdem_production.jl:1047-1050`) instead of discarding
hours of otherwise diagnosable work.

Driver consequences: `copdem_policy.jl` is unchanged — its two-function seam
(`:9-14`) binds to `d -> sourcesof(dependencies(plan), d)` and
`s -> consumerdegree(dependencies(plan), s)` at `copdem_production.jl:1242-1243`
instead of `plan.graph`; `refinegraph` (`:81`) selects the closure passed to
`plan_regrid`; the comment block at `copdem_production.jl:640-668` is
rewritten, since its conclusion holds only while the executor discovers
independently; `boxesoverlap` and its helpers (`:669-687`) survive and carry
the certificate above as a test; `copdem_dag_validate.jl:31, 61` builds the
graph three times per validation run and must build the plan once in
`main_validate`; `chunk_dag_poc.jl:248-249` and
`copdem_download_count.jl:181-185` are studies, not schedulers, and move to
the plan constructor, gaining the lower bound next to their existing
"refinement only removes edges" check.

Exit criterion: lazy execution performs no chunk-tree search after the plan is
built, and eager/lazy values remain bit-identical under existing invariance
tests.

## Phase 4: retire duplicate discovery state and APIs

After indexed graph construction and graph-backed execution are established:

- `CapQuery`, dual-tree `_descend!`, `DGGChunkTree`, and the production
  chunk-lattice use of `RasterFlatTree` were already removed in Phase 1; the
  old raster `chunktree` wrapper remains temporarily for compatibility;
- remove `connectedchunks`, `connectedchunks!`, and `connectedchunkpairs`
  (`discovery.jl:88-154`); the multi-cap union survives as a private
  `_tilesources!(out, graph, capsof[t])`;
- remove `chunk_dependency_graph(plan; refine, prefilter)` in favour of
  `dependencies(plan)`, and `refine` from the space form; note
  `test/scripts/copdem_policy.jl:340-341` calls `_chunkgraph` positionally and
  breaks on any signature change;
- remove the latitude cap join after the Phase 2 performance gate;
- remove `LazyRegridArray.srcindex` and any complete cap vectors retained solely
  for per-read discovery; and
- decide whether the public `chunktree` surface has remaining consumers. Do not
  preserve it merely as a wrapper over the private `chunkindex` query seam.

Phase 1C makes `chunkindex`/`candidatechunks!` a documented qualified extension
seam rather than a private implementation detail. If `chunktree` has no
remaining downstream consumer, remove its export and declaration here, remove
`RasterFlatTree`'s final compatibility role, and make `chunkextent` plus the
query seam the only chunk-geometry interface. Phase 1B has already moved
arbitrary **cell** subsets to the packed fallback, so this phase must not retain
`RasterFlatTree` for a role it no longer serves.

Keep leaf extent access, and make it cheap:

```julia
chunkextent(space, chunk)
```

It is not cheap today: the generic method is `chunkextents(space)[chunk]`
(`discovery.jl:86`), which rebuilds every chunk cap on a `RasterGrid`
(`rastergrid.jl:912-919`), and `connectedchunks` pays that on every call
(`discovery.jl:94-96`). `DGGChunkTree`'s `Trees.ncells`/`getcell`
(`src/regridding.jl:226-229`) go with it; no consumer asks a chunk tree for
cell geometry.

`chunkextents(space)` may remain as a convenience materialization for graph
construction, diagnostics, or algorithms that genuinely need every cap. Wave
cost estimation currently consumes source and destination cap vectors; move the
metadata it needs into the immutable dependency object or make it request only
the extents it actually uses before removing those vectors.

`HierarchicalGridCursor` and the surviving ConservativeRegridding cell-tree
integrations remain. Phase 1A converges `RasterCellTree` and `RasterGridView`
onto one CR curvilinear-grid adapter, and Phase 1B replaces only the
generic/scattered fallback; Phase 4 does not reopen either geometry decision.

Exit criterion: one candidate-query implementation defines graph edges, and the
executor retains neither a second discovery algorithm nor redundant index
state; the documented space interface contains no compatibility tree that is
immediately flattened or translated back into the query protocol.

## Phase 5: make `WeightBlock` the block-builder result

Introduce one method-level construction seam:

```julia
weightblock(method, dst_space, dst_inds, src_space, src_inds) -> WeightBlock
```

The generic fallback:

1. creates `WeightCOO`;
2. calls the existing incremental weight emitter; and
3. finalizes the COO into `WeightBlock`.

Conservative construction overrides `weightblock` and adopts the assembled CSC
matrix directly. It computes denominators once and returns the final block.

Both eager whole-domain construction and chunked construction call this seam.
Remove the separate `wholeblock(::Conservative)` optimization after the unified
path preserves its allocation and timing behavior.

### Make the final block own its reference weights

Normalize the denominator/reference representation in the same phase:

1. Change `WeightCOO` from an always-allocated `denom` plus `hasdenom::Bool` to
   `denom::Union{Nothing,Vector{Float64}}`.
2. `markdenominated!` or the first `adddenom!` allocates the zero-filled vector;
   point methods leave it as `nothing`.
3. Make `WeightBlock` own the final data-independent reference vector used by
   accumulation. For a denominated block, `reference` aliases `denom`; for a
   non-denominated block, construct the row sums once when the block is made.
4. Remove `CachedBlock.ref`. Lazy application reads `entry.block.reference`,
   and repeated eager application reuses the same stored reference rather than
   recomputing it for every call.
5. Preserve `blockreference!` and `addreference!` as executor operations, but
   implement them from the block's stored reference rather than redispatching
   over sparse/dense weights and denominators.

This is one invariant, not another cache: every final block carries exactly the
reference weights its values are normalized against. `denom === nothing` still
records whether those weights were explicit denominators or row sums, so
`hasdenom` remains derivable without a second boolean.

`Spilled` need not serialize a second vector. Its existing matrix/denominator
format reconstructs `WeightBlock.reference` on read; bump the private spill
version only if the serialized bytes actually change. `_blockbytes` must count
aliased `denom`/`reference` storage once.

The motivation is memory, not tidiness. The round trip was measured at 0.17 %
of time but a 1335 MiB `maxrss` delta, 8.4× the result, from the
ConservativeRegridding matrix, the COO and the final matrix being live
together (`docs/design/regrid-notes/perf-P1.md:354-368`; that directory is
gitignored, so Phase 0 re-measures and commits the number). It is on the
chunked path, the one with no fast route today.

Facts the seam must respect:

- step 3 already exists as `WeightBlock(coo::WeightCOO, ndst, nsrc)`
  (`plans.jl:26-29`); the seam is a dispatch wrapper;
- the only required method for a regridder is `build_weights!`
  (`methods.jl:130`), so a `weightblock` with a generic fallback is compatible;
- type-dispatch overrides are lost by wrapper methods.
  `T6CountingMethod(Conservative())` (`test_integration.jl:133`) and the
  `P1CountingMethod` harness already miss `wholeblock`; this phase would extend
  that to the chunked path, so instrumented benchmarks would measure the slow
  path. Select the override by trait or require wrappers to forward it;
- carry the empty-side contract: `wholeblock(::Conservative)` `invoke`s the
  generic path when either side is empty (`conservative.jl:360-363`) so
  denominator presence matches. Pin the present asymmetry first — empty `dst_inds`
  returns before `markdenominated!` (`:330`), empty `src_inds` after (`:331`);
- `wholeblock` is private (`GlobalRegridding.jl:60-99`); its removal is
  internal;
- adjudicate the alternative: `BlockAreaOperator` already holds both index
  maps (`conservative.jl:280-286`) and could push into the `WeightCOO`
  directly (`perf-P1.md:364-366`), removing the CSC instead of adopting it.
  That is a smaller change than a new seam and reaches the same single
  materialization.

Required laws:

- eager and chunked blocks have identical sparse structure and values
  (`test_conservative.jl:385-397` already asserts this for the two existing
  paths; re-point it at `weightblock` rather than deleting it);
- differently partitioned source chunks sum to the whole-domain block
  (`test_conservative.jl:89-114`, `test_lazy.jl:691-717`);
- denominators are identical, including zero-coverage rows
  (`test_conservative.jl:231-242, 399-407`); denominators are per-block
  positive intersected area summed by the executor (`lazy.jl:429`,
  `executor.jl:432`, thresholded at `executor.jl:305-320`), so per-block
  computation is the only correct design;
- a denominated block's `reference` and `denom` are the same immutable-by-
  convention vector, while a point block has `denom === nothing` and a stored
  row-sum reference;
- repeated eager applications and lazy cache hits allocate no new reference
  vector;
- spill round trips reconstruct identical denominator/reference semantics and
  `_blockbytes` does not double-count aliased arrays;
- no CSC-to-COO-to-CSC round trip occurs for Conservative; and
- third-party incremental builders continue to work through the generic fallback.

Exit criterion: one block-building path serves eager and lazy plans without
regressing the measured eager Conservative optimization, and the chunked
path's peak RSS drops by the re-measured round-trip delta; one final
`WeightBlock` owns weights, optional denominator semantics, and its reusable
reference vector, while caches add no duplicate numerical representation.

## Phase 6: remove DGG-specific output application

Teach `DGGSpace` to describe its destination dimension through the generic
output interface:

```julia
GR.destinationdims(space::DGGSpace, sampling) =
    (Cells(CellLookup(space.grid)),)
```

Then remove:

- `_DirectToDGG`;
- `_ChunkedToDGG`;
- the specialized eager and lazy `GR.regrid` methods; and
- `_ascube`.

Remove the adjacent DGG target-conversion layer as well. `regridgrid` is an open
generic with four one-line conversions and a `RegridTarget` union whose only
consumer immediately constructs `DGGSpace`. Replace it with direct `DGGSpace`
forwarding constructors or direct `_asspace` methods for `AbstractGrid`,
`CellLookup`, `CellVector`, and `MultiOrderCellSet`. Retain the separate bare
`AbstractHierarchicalGridSystem` error and source-dependent destination-level
resolution; those encode real behavior rather than conversion syntax.

At the user API boundary, keep the complete documented keyword signature on
`plan_regrid`, the one owner of defaults and validation. Let the convenience
forms forward:

```julia
regrid(data; kwargs...) = regrid(data, plan_regrid(data; kwargs...))
regrid!(dest, data; kwargs...) = regrid!(dest, data, plan_regrid(data; kwargs...))
```

This removes two duplicate copies of the keyword defaults and type constraints.
Unknown-keyword and invalid-value errors still come from `plan_regrid`; test
their messages so forwarding does not weaken the API boundary. The docstrings
continue to spell out the supported keywords even though the forwarding methods
do not repeat them in code.

Generic `wrapoutput` and `wraplazy` must preserve the existing `Cells` lookup
and all non-spatial dimensions. Avoid constructing a `ShapedRegridArray` for a
destination that already has one cell axis.

This is not deletion-only. Measured with the method above: eager types, dims
and values are identical; lazy dims and values are identical, but `wraplazy`
wraps every lazy result in a `ShapedRegridArray` (`lazy.jl:843`), so the lazy
`parent` changes from `LazyRegridArray` to `ShapedRegridArray`. Add a one-axis
short-circuit in `wraplazy` and a `parent`-type assertion in
`test/systems/crosssystem/regrid.jl`; no existing test checks it
(`regrid_acceptance.jl:96`, `regrid.jl:362-365`). The eager DGG method is pure
duplication — it `invoke`s the generic path, which already labels axis 1, then
relabels (`src/regridding.jl:275-291`). Correct the message at
`api.jl:184-186` ("a lazy regrid returns an unlabelled disk array"), already
false for DGG destinations. `destinationdims(plan::ChunkedPlan)` ignores
`plan.sampling` (`api.jl:105-106`), consistent with the method above dropping it.

Exit criterion: all DGG target spellings return the same types, dimensions,
lookups, and values through the generic path; there is no DGG-specific plan
application or intermediate target-conversion generic, and `plan_regrid` is the
single executable owner of regridding keyword defaults and validation.

## Phase 7: consolidate conservative destination-geometry caches

Replace the wrapper stack used for destination tiles with one private prepared
object:

```text
PreparedDestinationTile
|- destination index set and local-position map
|- restricted destination tree, built once
`- optional prebuilt destination polygons
```

The prepared tile is created once per executor tile and shared by all source
blocks built for it.

Refactor Conservative block construction so `BlockAreaOperator` obtains
destination polygons directly from the prepared tile. This should permit
removal of:

- `TileCells <: RegridSpace` and its interface-forwarding methods;
- `CachedCellTree`; and
- the destination `CellMemo` when prebuilt polygons are present.

Retain a small per-task memo for a side whose polygons are not prebuilt, notably
the source. Large destination tiles that exceed the polygon-cache threshold may
also use that fallback.

Do not merge this polygon cache with `MemoRasterTree` or `CapCachedTree` merely
because all are called caches. Those types memoize spatial extents and leaf caps,
not clipping polygons. ConservativeRegridding's cached dual DFS
(`src/utils/CachedDualDepthFirstSearch.jl` — present at `66ed54c`, absent from
the stale local tree, so Phase 0 first) makes `MemoRasterTree`'s internal-node
extent memo redundant (`regrid-notes/2026-08-21-cr-cached-dfs.md:209-218`);
its leaf-entry memo may still be useful. Remove and benchmark the node-extent
half first. The DGG leaf-cap cache has already been re-measured and is null:
widening `_MEMO_EXTENT_SLOTS` flips sign across columns, and the leaf-cap memo
ported into `MemoBlockCursor` had a −6.4 % profile ceiling with timings inside
noise; both reverted at `e2f90d1`
(`regrid-notes/2026-08-21-polar-profile.md:281-315`). What remains to measure
is `_CACHED_BUCKET_SIZE = 49` now that side-2 extents are stack-cached
(`cr-cached-dfs.md:215-218`). `MemoRasterTree`'s own justification
(`raster_tree_memo.jl:3-13`) is an argument, not a measurement;
`CapCachedTree`'s is in-tree (`src/cap_cached_tree.jl:78-108`). Keep
`PreparedDestinationTile` only if it reduces wrappers and code without
regressing the Conservative production benchmark.

Required laws:

- a destination tile's restricted tree is built once (tested:
  `test_conservative.jl:312-343`, `test_lazy.jl:447-468`);
- cached and uncached blocks are bit-identical (tested:
  `test_conservative.jl:265-271, 337-342`);
- concurrent block builds never mutate shared prepared geometry (holds by
  construction, `ReentrantLock` at `conservative.jl:176-186`; untested);
- cache memory remains bounded by the existing tile threshold
  (`_TILE_CELL_CACHE_MAX`, `conservative.jl:119`, appears in no test); and
- the production Conservative benchmark does not regress (no such benchmark
  exists in-repo: `benchmark/` holds toys, geomorphometry, hex_border,
  maxneighbors; `scripts/bench/` holds `halo_subset_scaling.jl`).

The last two artifacts are Phase 0 deliverables; until they exist this phase
has no gate.

Exit criterion: destination geometry follows one lookup/cache path rather than
several nested wrappers.

## Phase 8: simplify private lazy residency and chunk metadata

Do not introduce a generic LRU cache abstraction. `PerChunk`
(`plans.jl:212-229`) and `SourceHold` (`lazy.jl:198-211`) differ in lifetime,
locking, limits, duplicate-build behavior, oversize handling, and protected
entry semantics. Their shared eviction scan is too small to justify a public or
parameter-heavy cache framework.

There are nevertheless two local representations worth simplifying.

### Give `SourceHold` one entry dictionary

`SourceHold` currently keeps source arrays in `held` and their recency stamps in
a parallel `used` dictionary. Replace those with one private entry:

```julia
mutable struct HeldSource{A}
    value::A
    bytes::Int
    used::Int
end
```

The hold owns `Dict{Tuple{Int,Int},HeldSource}`. Touching and eviction update or
inspect the entry directly; byte removal uses the stored size rather than
recomputing it. This removes the invariant that `held` and `used` contain the
same keys without trying to make source arrays and weight blocks share a value
type. A tiny private oldest-entry helper may be shared with `PerChunk` only if
the resulting call sites are shorter and preserve their different locking and
`keep` rules.

### Normalize DiskArrays chunk metadata once

Lazy construction currently asks `DiskArrays.eachchunk(source)` separately in
`_passthroughchunks` and `_sourceotherchunks`; deriving a source `RasterGrid`
may already have inspected the same grid for its spatial ranges. Normalize the
available `GridChunks` once per lazy-array construction and pass that immutable
description to:

- output pass-through chunk construction;
- non-spatial source read groups; and
- spatial source-space construction where the source is inferred in the same
  call and doing so does not couple `RasterGrid` to a `LazyRegridArray`.

Do not store a second chunk layout on the plan. Chunk metadata depends on the
applied data array, while the plan is reusable across compatible arrays.

Rename `isdiskbacked` to match its actual predicate: it tests
`DiskArrays.haschunks(data) isa DiskArrays.Chunked`, not whether storage is on
disk. A name such as `_haschunkedstorage` or `_ischunkedsource` keeps the current
default-lazy behavior explicit. Preserve the distinct `_isdisksource` check
where the executor genuinely needs to know whether a read must go through the
DiskArrays interface.

Required laws:

- source load, hit, skipped, dropped, current-byte and peak-byte statistics are
  identical for the same access sequence;
- the hold never retains an oversized source array and remains within its byte
  limit apart from the documented scratch buffer;
- output chunking and non-spatial read grouping are identical for regular,
  irregular, absent, and malformed/non-grid `eachchunk` results;
- lazy construction still reads no source values; and
- repeated `eachchunk` calls are absent from one lazy-array construction.

`Spilled` still wraps `PerChunk` and has a separate open fingerprint redesign;
that work remains outside this plan.

Exit criterion: source residency has one dictionary and one recency/byte record
per held array, DiskArrays chunk metadata is interpreted once per application,
and no generic cache hierarchy has been added.

## Verification matrix

Every phase should run the narrow tests first, then the full package and DGG
integration suites.

| Concern | Primary coverage |
|---|---|
| Raster DiskArrays ownership, transformed-grid quadtree coverage, and candidate identity | `lib/GlobalRegridding/test/test_rastergrid.jl`, `lib/GlobalRegridding/test/test_lazy.jl`, brute-force geometry property tests |
| Raster native-CRS transforms, Proj task isolation/lifetime, and unit-sphere naming | `lib/GlobalRegridding/test/test_rastergrid.jl`, Proj extension tests, concurrent-task projected-grid property tests |
| Public cap-cap/cap-box extent protocol and outward XYZ bounds | GeometryOps unit-spherical/Extents tests, ConservativeRegridding mixed dual-tree tests |
| Packed chunk and fallback-cell R-tree identity | GeometryOps FlexibleRTree tests, `lib/GlobalRegridding/test/test_lazy.jl`, `lib/GlobalRegridding/test/test_conservative.jl` |
| Complete qualified `RegridSpace` extension contract | interface conformance tests plus DGG method-coverage assertions |
| DGG chunk frontier, descendant coverage, and compact numbering | system `node_extent` tests, `test/systems/crosssystem/regrid.jl` |
| Graph construction and CSR laws | `lib/GlobalRegridding/test/test_chunkgraph.jl` |
| Executor graph use and graph-miss policy | `lib/GlobalRegridding/test/test_lazy.jl`, `test/scripts/copdem_policy.jl` |
| Weight-block identity, partition invariance, and one stored reference | `lib/GlobalRegridding/test/test_conservative.jl`, `lib/GlobalRegridding/test/test_interpolation.jl`, spill round trips |
| Eager/lazy equivalence and destination dimensions | `lib/GlobalRegridding/test/test_integration.jl`, `test/systems/crosssystem/regrid.jl` |
| Budget, residency, spilling, and empty chunks | `lib/GlobalRegridding/test/test_lazy.jl`, `test/systems/crosssystem/regrid_acceptance.jl` |
| Single-dictionary source hold and one-time DiskArrays metadata normalization | `lib/GlobalRegridding/test/test_lazy.jl`, counting `eachchunk` fixture |
| Production scheduling/cache compatibility | `scripts/copdem_dag_validate.jl`, production subset runs |
| Synthetic-source absoluteness | `test/scripts/copdem_source_mode.jl` |
| Every demand was a graph edge, against a real network | `scripts/copdem_prefetch_coldtest.jl:247-248` |

Gate hygiene:

- The suite's one broken test is
  `lib/GlobalRegridding/test/test_conservative.jl:147`, the GeometryOps
  non-convex destination clip; the same arm is
  `test/systems/crosssystem/regridding_conservation.jl:86-99`, three
  `@test_broken` per non-convex system. Phase 5's denominator law is asserted
  over a path with that known non-conserving arm.
- Eleven production columns are hybrid-source
  (`regrid-notes/2026-08-21-pole-artifact.md:196-208`: 122975, 122977, 122979,
  122981, 123147, 123172, 123176, 123199, 123202, 123203, 123204) and are not
  valid byte-identity references for any production-subset gate.
- Any peak-memory gate states which `shape` it runs under; `outer → inner` is
  the dominant RSS lever (`regrid-notes/2026-08-21-polar-profile.md:387-391`),
  and the unexplained epoch-C collapse (`polar-profile.md:373-386`) bounds
  what "no material regression" can mean on the shared box.
- Every number in this plan carries a note, script, or test reference; a
  number with no committed artifact is labelled as a new measurement.

Standing acceptance laws:

1. No dependency false negatives.
2. Source chunks are applied in deterministic order.
3. Eager, lazy, and differently chunked plans produce equivalent values.
4. Reusing a plan reuses its dependency graph and built weight blocks.
5. Plan construction reads no source values.
6. Memory use is bounded by declared budgets except for explicitly documented
   immutable planning structures.
7. Graph and executor agree exactly on the set of permitted reads, by
   construction: one graph object per plan, its narrow phase set once.
8. Provider/network metadata work remains driver policy, not regridding plan
   construction.
9. Every final weight block owns one reusable reference vector; cache wrappers
   do not copy it.
10. Every external `RegridSpace` specialization is named by the documented
    qualified extension contract.

## Commit strategy

Keep the work bisectable and avoid combining performance changes with API
cleanup:

1. Baseline: re-instantiate, repin both sources to reachable revisions,
   record tree SHAs in the commit message, create the missing baseline
   artifacts.
2. Tested DiskArrays-backed raster quadtree chunk index.
3. Tested `cap_xyz_extent` and generic `FlexibleRTrees` chunk index.
4. DGG original-tree chunk frontier and descendant-coverage tests.
5. GeometryOps transform reuse, precise native/unit-sphere naming, and the
   backwards-compatible callable adapter.
6. Proj extension with task-local context/clone ownership and projected-raster
   correctness tests; no direct PROJ calls.
7. Convergence of `RasterCellTree`/`RasterGridView` and cap machinery after its
   correctness and performance gate.
8. GeometryOps public spherical extent operations and ConservativeRegridding
   mixed cap/box dual-search conformance.
9. Packed fallback cell-tree adapter, followed by removal of `CellCapTree` and
   the scattered-cell role of `RasterFlatTree` after its performance gate.
10. Consolidated, documented qualified `RegridSpace` extension surface; no new
    interface wrapper.
11. Indexed graph-builder correctness and performance comparison.
12. Removal of the latitude join after that gate.
13. Plan-owned dependency graph: `refine` moves to `plan_regrid`, the graph
   gains its stamp, `restrict` row views; execution unchanged.
14. Lazy executor switched to graph rows and duplicate discovery removed,
    including the legacy `chunktree` surface if it has no external consumers.
15. Unified `weightblock` construction and one block-owned denominator/reference
    representation.
16. DGG output/target-conversion removal and keyword forwarding through
    `plan_regrid`.
17. Targeted stale-cache removal, then prepared destination geometry if earned.
18. Single-entry-dictionary `SourceHold`, one-time DiskArrays chunk-metadata
    normalization, and accurate chunked-storage naming.

Record allocation, graph-size, planning-time, and execution-time changes with
each performance-sensitive commit. Structural simplification is complete only
when the replacement has both fewer responsibilities and evidence that it
preserves the relevant scaling behavior.

## Definition of done

- Raster chunk queries use DiskArrays ranges over an existing CR quadtree;
  generic chunk queries use `FlexibleRTrees`; DGG queries expose the chunk-level
  frontier of `treeify(space.grid)` and its existing subtree-covering extents.
- `RasterGrid` accepts arbitrary native raster CRSs through explicitly named
  native-to-unit-sphere transforms; geographic conversion is supplied by
  GeometryOps and projected conversion by task-local Proj.jl clones.
- There is one structured raster lattice adapter, and generic projected-grid
  extent construction does not rely on geographic-only finite sampling.
- GeometryOps publicly owns cap-cap/cap-box intersection, Cartesian cap bounds,
  containment, and merging; ConservativeRegridding consumes the common extent
  predicate for mixed tree nodes.
- One packed `FlexibleRTrees` adapter handles arbitrary cell subsets; no custom
  general-purpose Raster, DGG, or fallback cell tree duplicates existing
  infrastructure.
- The documented qualified `RegridSpace` interface names every cross-package
  specialization, and the legacy `chunktree` contract is absent unless a real
  remaining consumer justifies it.
- A logical `ChunkedPlan` owns or validates a reference to one reusable
  dependency graph rather than rebuilding it per destination column.
- Lazy reads perform no geometric dependency discovery.
- Production scheduling and lazy execution consume the same relation.
- Generic cap-box broad-phase queries and exact cap filtering match the cap
  reference; raster-quadtree and DGG-frontier queries satisfy the stronger
  geometric no-false-negative law for zero and nonzero support radii.
- Conservative eager and chunked construction share one `WeightBlock` builder;
  the final block owns one reusable reference vector and cache wrappers do not
  copy it.
- DGG output uses the generic destination-dimension path, DGG targets resolve
  directly to `DGGSpace`, and `plan_regrid` alone owns keyword defaults.
- Destination polygon caching has one prepared-tile abstraction.
- Lazy source residency uses one entry dictionary, and each applied source's
  DiskArrays chunk metadata is normalized once.
- Point-method redesign remains untouched and independently actionable.
- The complete regridding test matrix is green.
- Representative production and small-spatial-chunk benchmarks show no material
  regression in time or peak memory.
- Both git sources are pinned to reachable revisions, the Manifest re-resolves
  to their trees, and the tree SHAs are recorded in the pin commit.
