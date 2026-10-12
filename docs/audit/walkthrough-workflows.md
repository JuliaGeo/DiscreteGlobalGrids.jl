# Newcomer walkthrough: workflow use cases

This report walks through UC-025 through UC-057 in
[`use-cases.md`](use-cases.md). I began each case at its listed public
documentation entry. I read implementation source only after the public path did
not answer a question. "Docstring" below means a public docstring rendered by an
`@docs` block on the named page.

## Execution record

I used the Julia daemon test environment named `workflow_walkthrough`. I ran five
small local scripts. No script downloaded data or rendered a plot.

| Run | Scope | Observed result |
|---|---|---|
| N | HEALPix level 2 neighborhood, needs, adjacency, boundary, growth, and error checks | Passed. `ring` lengths were 8 and 16; the radius-2 disk had 24 unique cells; an isolated subset had zero neighbors; threaded failure returned `NeighborCallbackError`. |
| R | HEALPix level transfer, plan reuse, regional target, and lazy plan | Passed except the intentional source-inference probe. A `DimArray` over `Cells` raised `ArgumentError` without `from`; explicit `from` passed. |
| C | Twelve forced chunks, custom chunk inspection, halo 2, and aliased output | Non-aliased output matched the eager result. Halo 2 still supplied one-ring callbacks. In-place aliasing produced different values. |
| S | Local Zarr round trip with two variables and a second dimension | Passed. Values and dimensions round-tripped lazily; the axis was a `ChunkedCellLookup`; `:auto` chose ranges; an existing path failed; an absent ID selection raised `BoundsError`. |
| P | Weighted contiguous partition, missing backend, and chunked regrid partition | Passed. Every task appeared exactly once; a missing Metis extension raised `PartitionBackendUnavailable`; all 48 regrid destination tasks were assigned. |

The scripts are temporary audit probes, not package tests. Static conclusions are
marked separately from executed conclusions.

## E. Compute with neighborhoods and boundaries

### UC-025: Find first-order neighbors or an exact higher-order ring

- **Discovery path:** `api/neighbors.md:7-34` → rendered `neighbors` and `ring`
  docstrings at `src/interface/grid.jl:316-464`.
- **Does `neighbors(..., k)` include all distances through `k`?** Yes. It
  concatenates rings 1 through `k` (`src/interface/grid.jl:332-338`).
- **Is the center included?** No (`docs/src/api/neighbors.md:12-16`).
- **How does `ring` differ, and what is order zero?** `ring` returns exactly one
  topological shell. `ring(..., 0)` returns the center; `neighbors(..., 0)` is
  empty (`src/interface/grid.jl:407-408,425-435`).
- **Is distance topological or metric?** Topological adjacency steps, selected by
  `Vertex()` or `Edge()` (`docs/src/api/neighbors.md:37-47`).
- **How are duplicates handled?** The documented ring concatenation makes each
  distance shell disjoint. Run N confirmed 24 unique cells in the two-ring disk.
- **Result:** **Discoverable and executed.** No source beyond rendered public
  docstrings was required.

### UC-026: Convolve a global dataset over second-order neighbors

- **Discovery path:** `tutorials/stencils.jl:106-133` →
  `api/neighbors.md:75-117` → rendered sweep docstrings.
- **Disk, exact ring, or diffusion?** The tutorial now distinguishes all three.
  A disk contains rings 1 and 2, an exact ring contains only distance 2, and two
  one-ring passes are diffusion with path-dependent multiplicity
  (`docs/src/tutorials/stencils.jl:129-133`).
- **Is the center included?** Direct disk/ring queries exclude it except
  `ring(..., 0)`; a sweep callback receives the center separately
  (`docs/src/api/neighbors.md:12-16,85-92`).
- **How should weights depend on distance?** This is the caller's kernel decision.
  The API defines membership and order, not a weighting law.
- **What changes at pentagons and seams?** Rows can have different lengths; the
  tutorial warns against fixed width at `docs/src/tutorials/stencils.jl:230-249`.
  System topology handles seams.
- **What natural API is missing?** A radius or exact-ring selector on all sweep
  forms. The docs state the limit at `docs/src/api/neighbors.md:77-83`.
- **Result:** **API GAP, statically confirmed and executed.** Run C showed
  `halo=2` still produces the one-ring result. Link: API-005, DOC-014, SW-001,
  SW-002. A larger halo is not counted as success.

### UC-027: Apply a first-order stencil to a field

- **Discovery path:** `tutorials/stencils.jl:106-133,366-375` →
  `api/neighbors.md:75-129`.
- **Callback arguments and center:** `Values()` calls
  `f(cell, center_value, neighbor_values)`; the center is the second argument.
  `Neighbors()` passes indexed handles, and `NeighborSlices()` passes whole
  non-spatial slices (`docs/src/api/neighbors.md:85-117`).
- **Callback return type:** One result per cell; a concrete tuple becomes one
  output per component (`src/engine/neighborhood.jl:505-506`).
- **Can the input be a plain vector?** Yes, when the cell collection and vector
  are passed separately (`docs/src/tutorials/stencils.jl:366-375`).
- **Are dimensions preserved?** `Values()` preserves the input dimensions;
  handle and slice modes return the cell dimension only
  (`docs/src/api/neighbors.md:85-92`).
- **Isolated cell:** Its neighbor sequence is empty. Run N returned center value
  7.0 and neighbor count 0 for a one-cell subset.
- **Result:** **Discoverable and executed.** The isolated-cell outcome follows
  clipping but is not called out with an example.

### UC-028: Read values, centroids, and identifiers in one stencil

- **Discovery path:** `api/neighbor-fields.md:7-64` → its complete worked example
  at `docs/src/api/neighbor-fields.md:66-121`.
- **Field-major or neighbor-major?** Field-major: `rings[j]` is field `j` for all
  neighbors. Use `zip(rings...)` for neighbor-major records
  (`docs/src/api/neighbor-fields.md:36-49`).
- **Which indices address which arrays?** `Index(Local())` names the swept
  collection; `Index(Global())` names the complete level. A `Value(data)` vector
  must follow the swept collection (`docs/src/api/neighbor-fields.md:50-58`).
- **How do fields stay aligned?** Slot `i` in every ring names the same neighbor
  (`docs/src/api/neighbor-fields.md:43-48`). Run N confirmed equal ring lengths
  for values, centroids, cells, and local indices.
- **Can a quantity be computed lazily?** `cellfield` computes missing values from
  a function and can use known values (`docs/src/api/neighbor-fields.md:122-173`).
- **Can a partial field supply only needed cells?** Yes. A labeled `known` subset
  resolves by cell identity; a dense vector must cover the collection
  (`src/engine/cellfield.jl:32-60`, rendered on the field page).
- **Result:** **Discoverable and executed.** No private source was required.

### UC-029: Use neighbor order for a directional kernel

- **Discovery path:** `api/neighbors.md:12-34,142-155` → rendered `neighbors`
  docstring at `src/interface/grid.jl:332-420`.
- **Is order counterclockwise?** Yes, seen from outside the sphere.
- **Where does the first slot begin?** It is deterministic within a system;
  rings 2 through `k` start on ring 1's spoke. The starting direction is not
  portable across systems (`src/interface/grid.jl:349-359`).
- **Does clipping preserve slot numbers?** It preserves surviving order, but
  removed absolute slots are not recoverable (`src/interface/grid.jl:383-403`).
- **Is a slot a compass direction?** Only if that system documents such a phase.
  The rotation direction is portable; compass labels are not.
- **Are IGeo7 direction codes portable?** No. `directioncode` belongs to IGeo7's
  system-specific contract; generic kernels must use the ordering contract.
- **Result:** **Discoverable, static.** The public docstring answers the generic
  questions; no execution was needed beyond Run N's order-preserving equality.

### UC-030: Run a stencil on a region with holes

- **Discovery path:** `tutorials/stencils.jl:278-288` →
  `api/neighbors.md:12-34` → rendered `neighbors` coverage contract.
- **Are neighbors clipped?** Yes, to subset membership.
- **Does a hole alter distance?** It removes returned cells. Distance remains the
  complete system's distance and does not route around the hole
  (`src/interface/grid.jl:383-403`).
- **How are outside values supplied?** Use `grow` for an in-memory margin or the
  chunk runner's halo for stored data (`docs/src/api/boundaries.md:79-94` and
  `docs/src/api/chunk-sweep.md:28-39`).
- **What if the center is outside the region?** Direct `neighbors` raises
  `ArgumentError` (`src/interface/grid.jl:405`). Sweeps visit region members only.
- **What boundary condition results?** A direct subset sweep omits outside values;
  the kernel therefore chooses the meaning of that clipped neighborhood. The
  docs show the numerical difference but do not name a standard boundary
  condition (`docs/src/tutorials/stencils.jl:280-288`).
- **Result:** **Discoverable and executed.** Run N confirmed radius-2 subset
  results equal complete-level results filtered by membership.

### UC-031: Cache adjacency for repeated graph calculations

- **Discovery path:** `api/neighbors.md:50-73` → rendered `adjacency` and
  `AdjacencyTable` docstrings.
- **Row and column index space:** Row `i` is local cell `i`. With `halo=0`, entries
  are local region positions; with `halo=1`, entries index `[region; halo]`;
  `halo=:mark` keeps complete slots and writes 0 for outside cells
  (`docs/src/api/neighbors.md:50-66`).
- **What does zero mark?** An outside neighbor in the marked complete row.
- **Why are wider halos rejected?** An adjacency row represents one-ring edges;
  wider reach is a different graph product. The public signature accepts only
  `0`, `1`, or `:mark`; invalid values raise.
- **Can the table be reused with another ordering?** No. Its local positions are
  tied to the region order used to build it.
- **Result:** **Discoverable and executed.** Run N matched every `halo=0` row to a
  direct query and validated the region-plus-halo index set.

### UC-032: Find the border, interior, and outside halo of a region

- **Discovery path:** `api/boundaries.md:7-77`.
- **Which side does each verb return?** `border` is inside and touches outside;
  `interior` is inside with all neighbors inside; `halo` is outside and touches
  the region (`docs/src/api/boundaries.md:12-22`).
- **Are holes included?** Yes. Cells across a hole are in the outside halo.
- **Lazy, sorted, and length?** The iterators are lazy and yield ascending indices
  by default. Some have no cheap exact `length`; collect to count
  (`docs/src/api/boundaries.md:54-71`).
- **Subtree versus arbitrary subset:** Both use the same result contract. The
  internals page explains that a rooted subtree can use a specialized traversal.
- **Result:** **Discoverable and executed.** Run N confirmed border and interior
  were disjoint and partitioned the region, and halo cells were external.

### UC-033: Grow a region for a multi-step computation

- **Discovery path:** `api/boundaries.md:79-94` → rendered `grow` docstring at
  `src/engine/region_algebra.jl:10-39`.
- **Does growth include the original region?** Yes, plus `n` same-level layers.
- **Which connectivity applies?** `Vertex()` by default; `Edge()` is accepted.
- **Does repeated growth match wider growth?** The implementation defines growth
  as repeated halo walks. Run N confirmed `grow(grow(region,1),1) == grow(region,2)`.
- **How are original output positions mapped?** Growth returns a newly ordered
  `CellVector`; the public workflow does not show the cell-ID lookup needed to map
  original outputs into it. Source reading confirmed there is no retained
  ownership map (`src/engine/region_algebra.jl:30-54`).
- **Result:** **DOC GAP, partly executed.** See UW-003.

### UC-034: Choose traversal order and handle callback failures

- **Discovery path:** `api/neighbors.md:116-129` → rendered sweep and
  `NeighborCallbackError` docstrings.
- **Which buffers may callbacks retain?** Public docs do not state ownership or
  lifetime for neighbor sequences, needs rings, or `ChunkCube` buffers. Source
  currently allocates value gathers and chunk buffers, but that is not a public
  promise.
- **Can callbacks mutate input or shared state?** Threaded callbacks must be
  order-independent (`src/engine/neighborhood.jl:509-513`). Shared mutation needs
  caller synchronization; input mutation safety is not specified.
- **Which iteration order is guaranteed?** `StorageOrder()` or an explicit
  permutation controls visits. Results always return in collection index order.
- **How are exceptions attributed?** Threaded failure wraps the cause with cell
  and local index; serial execution throws the original error.
- **How does threading affect determinism?** Pure, order-independent kernels have
  the same results. Side effects and order dependence are the caller's risk.
- **Result:** **DOC GAP and executed.** Run N received a correctly attributed
  `NeighborCallbackError`. Buffer ownership remains UW-002.

### UC-035: Compute a cost-distance field through neighboring cells

- **Discovery path:** `tutorials/stencils.jl:290-364`.
- **Connectivity:** The example uses default vertex connectivity. A caller can
  request edge connectivity in `neighbors`.
- **Edge distance and cost:** The example averages endpoint costs per graph step;
  it explicitly says a per-kilometer model must add center distance
  (`docs/src/tutorials/stencils.jl:292-310`).
- **Are missing cells barriers?** The example models barriers with `Inf`. It does
  not state how `missing` should be handled.
- **Which graph supplies shortest paths?** `neighbors(lookup, i)` uses the
  selected lookup, so absent cells remove edges. Unreachable cells remain `Inf`
  because the initialized distance is never changed.
- **Package API or example code?** Example code. `costdistance` is defined inside
  the tutorial and is not exported.
- **Result:** **DOC GAP, static.** Missing-value and disconnected-region behavior
  must be inferred from the example. See UW-004.

## F. Transfer data between spatial grids

### UC-036: Regrid a raster onto a DGGS

- **Discovery path:** `tutorials/regridding.jl:1-68` →
  `api/regridding.md:7-38` → `api/regridding-methods.md`.
- **Accepted raster, dimensions, and coordinates:** The tutorial uses a Rasters.jl
  raster with longitude, latitude, and time. GlobalRegridding's public README
  defines `RasterGrid` on longitude/latitude degrees and requires an explicit
  transform for projected coordinates (`lib/GlobalRegridding/README.md:17-18,60-64`).
- **Input footprints:** `Conservative()` treats values as cell-area averages and
  intersects raster-cell footprints (`docs/src/api/regridding-methods.md:24-30`).
- **Destination coverage:** A complete grid target covers the complete grid;
  source absence becomes missing output under the selected policy. Use a partial
  target to request only a region (`docs/src/api/regridding.md:29-32`).
- **Result wrapper and cell axis:** A dimensional input keeps its wrapper and
  non-spatial dimensions while replacing spatial axes with `Cells`
  (`docs/src/api/regridding.md:34-37`; tutorial example at lines 65-68).
- **Result:** **Discoverable and partly executed.** Run R confirmed constant
  preservation and target length. The raster-specific published tutorial was not
  executed because it downloads CPC data and renders plots.

### UC-037: Move a field between DGGS systems or levels

- **Discovery path:** `tutorials/between_grids.jl:39-86,154-164` →
  `api/regridding.md:13-38`.
- **When is source inference available?** Raster dimensions can infer a raster
  space. A plain vector requires `from`. The page says a dimensional source can
  describe its axes but does not state that a DGGS `Cells` axis still requires
  explicit `from` (`docs/src/api/regridding.md:13-35`).
- **When is `from` required?** Run R showed that omitting it from a DGGS
  `DimArray` raises `ArgumentError`; explicit `from=source_grid` succeeds. Source
  confirms the Cells dimension is rejected as a raster lattice at
  `lib/GlobalRegridding/src/api.jl:290-310`, despite the DGG adapter exposing
  `dimsource` at `src/regridding.jl:305-307`.
- **Can a regional source feed a complete target?** Yes. The source space can be
  a regional grid; uncovered target cells follow the missing policy.
- **Finer versus coarser:** Conservative transfer repeats values inside source
  footprints when refining and combines overlaps when coarsening; point methods
  sample/interpolate (`docs/src/api/regridding-methods.md:24-65`).
- **Does an equal level number mean equal resolution?** No. The between-grids
  tutorial compares physical cell size and chooses levels independently.
- **Result:** **DOC GAP with runtime evidence.** See UW-001. Other behavior is
  discoverable; Run R confirmed explicit and regional transfers.

### UC-038: Choose area conservation, interpolation, or containing-cell transfer

- **Discovery path:** `api/regridding-methods.md:7-67` →
  `tutorials/between_grids.jl:69-145`.
- **Data meanings:** `Conservative()` is for cell averages; `BarycentricPoint()`
  is for point samples; `NearestCell()` and `DirectNearest()` copy the containing
  source cell, suitable for categories. `Extensive()` handles totals through the
  missing-policy layer (`docs/src/api/regridding-methods.md:9-21,24-65`).
- **Which preserves integrals?** Conservative area weighting, subject to coverage
  and missing rules. Point methods do not.
- **What does nearest mean?** The source cell containing the destination center,
  not nearest centroid (`docs/src/api/regridding-methods.md:44-46`).
- **When is an operator useful?** `NearestCell()` supports plan reuse;
  `DirectNearest()` avoids building the operator for little reuse
  (`docs/src/api/regridding-methods.md:48-49`).
- **Which methods retain source values?** Containing-cell transfer copies one
  selected value exactly. Barycentric interpolation reproduces an exact source
  sample site and stays within contributing values because weights are
  nonnegative and sum to one (`docs/src/api/regridding-methods.md:36-46`).
- **Result:** **Discoverable and executed for containing-cell transfer.** Run R
  preserved a constant under `NearestCell()`.

### UC-039: Control missing values and partial source coverage

- **Discovery path:** `api/regridding-methods.md:69-114` → rendered `Weighted`
  and `Extensive` docstrings.
- **Threshold meaning:** `Weighted(t)` requires valid weight to reach fraction
  `t` of total source-covered weight, then normalizes by valid weight
  (`docs/src/api/regridding-methods.md:69-85`).
- **Uncovered versus missing source:** No source coverage and covered-but-invalid
  data both blank a weighted result, but they arise from different denominator
  and validity checks. The public page explains invalid weight but does not give
  a side-by-side uncovered example.
- **Output sentinel:** The table at `docs/src/api/regridding-methods.md:87-100`
  covers Raster `missingval`, `missing`, `NaN`, DimArray, and bare arrays.
- **Can `dest` hold it?** `regrid!` requires a destination element type that can
  hold its declared sentinel (`docs/src/api/regridding-methods.md:102-114`).
- **Totals versus averages:** `Weighted` normalizes; `Extensive` returns
  unnormalized weighted sums (`lib/GlobalRegridding/README.md:87-91`).
- **Result:** **Discoverable, static.** The published policy contract answers all
  questions; missing-data variants were not separately executed.

### UC-040: Reuse a regridding plan across time or variables

- **Discovery path:** `tutorials/regridding.jl:108-122` and
  `tutorials/between_grids.jl:154-164` → `api/regridding.md:40-59`.
- **What must remain fixed?** Source and destination spaces, cell order, method,
  and missing policy. Geometry or order changes require a new plan
  (`docs/src/api/regridding.md:40-44`).
- **Can the data mask change?** Values and invalid entries can change because
  validity is applied when the plan runs. The source sentinel fixed on the plan
  and its missing policy must remain appropriate. The page warns readers to
  review missing normalization for fields with different coverage
  (`docs/src/api/regridding.md:76-78`).
- **Destination allocation:** `regrid(data, plan)` allocates;
  `regrid!(dest, data, plan)` uses caller storage with checked shape and sentinel
  type (`docs/src/api/regridding.md:55-58`).
- **Same missing policy?** Yes, a plan fixes it. Build another plan to change it.
- **Concurrent application:** The public task page does not promise concurrent
  calls on one plan. Private storage uses locks, but this is not a user contract.
- **Result:** **Discoverable and executed except concurrency.** Run R reused one
  plan on two fields and matched independent linear behavior. Concurrency remains
  **untested** and should not be promised from implementation details.

### UC-041: Regrid onto a regional or mixed-level destination

- **Discovery path:** `tutorials/multiorder.jl` → `api/regridding.md:29-37` →
  rendered `DGGSpace` and target adapter docstrings.
- **Accepted target types:** Complete and partial grids, `CellVector`,
  `CellLookup`, and `MultiOrderCellSet` (`docs/src/api/regridding.md:29-32`).
- **Output order:** Destination local order; a DGG destination receives a
  `CellLookup` over those cells (`src/regridding.jl:311-323`, public docstring).
- **Mixed-level footprints:** The task page says mixed-level sets are accepted but
  does not explain that the adapter constructs `CellVector(set)` and then a
  `PartialGrid` (`src/regridding.jl:289-294`). This expands the represented area
  to one leaf level rather than returning one value per stored mixed-level member.
- **Empty destination:** Public docs do not state whether an empty target is valid.
  The adapter has empty-range handling, but no walkthrough or explicit contract.
- **Which collections are axes?** `CellLookup` is the data-axis form; grids,
  vectors, and mixed sets are target descriptions.
- **Result:** **DOC GAP, partly executed.** Run R confirmed a 15-cell regional
  target follows the corresponding full-target order. Mixed-level and empty
  behavior remain static/untested. See UW-005.

### UC-042: Regrid back to a raster or another external space

- **Discovery path:** `tutorials/regridding.jl:150-184` →
  `lib/GlobalRegridding/README.md:17-18,41-76`.
- **Destination forms:** An existing Rasters.jl object can serve as a template;
  GlobalRegridding also provides `RasterGrid(data_or_dims)`.
- **Resolution, extent, projection, dimension order:** Destination X/Y coordinates
  set resolution and extent. `RasterGrid` assumes lon/lat degrees; projected
  coordinates require transformation. Spatial dimensions come first
  (`lib/GlobalRegridding/README.md:41-64`).
- **Regional or chunked DGGS source:** A regional DGG source is accepted with an
  explicit `from`; chunked data can use the lazy plan path.
- **Which package owns the call?** GlobalRegridding owns the verbs and external
  spaces; DGG re-exports them and supplies `DGGSpace`
  (`docs/src/api/regridding.md:7-10`).
- **Result:** **Discoverable, not executed.** The published tutorial uses a remote
  raster, so the audit did not run it. No invented target constructor was needed.

### UC-043: Regrid large or tiled inputs within a memory budget

- **Discovery path:** `api/regridding.md:61-78` → rendered `PerChunk`, `Spilled`,
  `DGGSpace` docstrings → `api/partitioning.md`.
- **Weights versus source data:** `budget` is a target for both resident weights
  and loaded source data, not total process memory. `PerChunk` bounds cached
  weights; source reads use the plan's data share (`docs/src/api/regridding.md:61-71`;
  `lib/GlobalRegridding/src/plans.jl:114-133`).
- **Where do temporary weights go?** Memory under `PerChunk`; files in the caller's
  directory under `Spilled(dir)` (`lib/GlobalRegridding/src/plans.jl:225-260,435-456`).
- **Which chunk controls belong where?** `DGGSpace(...; chunklevel, chunkcells)`
  describes source/destination spatial chunks; `chunks` controls lazy destination
  output tiling (`docs/src/api/regridding.md:63-70`).
- **Peak memory:** The docs correctly call `budget` a target, not a strict process
  bound. They do not provide a complete estimator for runtime, geometry, and
  wrapper overhead.
- **Cleanup and restart:** The `Spilled` docstring says the caller owns directory
  deletion and spill files are tied to one storage instance, so they cannot
  restart another plan (`lib/GlobalRegridding/src/plans.jl:435-456`).
- **Result:** **Discoverable and partly executed.** Run R forced 12 source and 48
  destination chunks under a small lazy plan and matched eager output. Disk spill
  was not executed.

### UC-044: Interpret interpolation near boundaries, poles, and degeneracies

- **Discovery path:** `api/regridding-methods.md:116-134` and its point-method
  example at lines 136-171.
- **Outside a sample hull:** Unmapped.
- **Missing corners and regional rims:** A destination is unmapped when a required
  source stencil member is absent. Degenerate or folded polygons are also
  unmapped (`docs/src/api/regridding-methods.md:116-125`).
- **Pole cells:** Copernicus DEM uses the configured `poles` policy.
  `NearestCell()` is the default fallback; `nothing` leaves them unmapped
  (`docs/src/api/regridding-methods.md:125-134`).
- **Which fallback preserves a value?** `NearestCell()` selects one source value
  with weight one.
- **Can interpolation overshoot?** Barycentric weights are nonnegative and sum to
  one, so the result stays within contributing values
  (`docs/src/api/regridding-methods.md:36-42`).
- **Result:** **Discoverable, static.** The requested specialized fixtures were
  not rerun; the page states the tested contract directly.

## G. Read, write, and inspect stored cubes

### UC-045: Write and reopen a cell-indexed cube

- **Discovery path:** `tutorials/store_io.jl:1-53` → `api/store-io.md:7-54` →
  rendered `dggwrite` and `dggread` docstrings at `src/io/api.jl:26-74`.
- **Activation:** `using Zarr` loads the extension.
- **Inputs and dimensions:** `DimArray` or `DimStack` with a `Cells` dimension
  backed by a cell lookup. Other dimensions and coordinates are stored with each
  layer (`src/io/api.jl:51-67`).
- **Variables, metadata, missing values, chunks, and axis:** Layer metadata become
  array attributes; stack `metadata["attrs"]` becomes group attributes; chunks
  count cells and other dimensions; the chosen encoding stores or derives the
  cell axis. The full missing/fill details are in the rendered writer docstring.
- **Does reading load immediately?** Data stay lazy unless `lazy=false`; axis
  metadata and validation work occur at open (`src/io/api.jl:31-43`).
- **Existing path:** The writer does not overwrite an existing layer
  (`ext/DiscreteGlobalGridsZarrExt/write.jl:79-96`). Run S confirmed a second
  write to the same directory fails.
- **Result:** **Discoverable and executed.** Run S preserved two variables, a
  second dimension, values, and lazy storage.

### UC-046: Select a small region from a large or remote store

- **Discovery path:** `tutorials/store_io.jl:55-96,152-165` → rendered `dggread`
  docstring at `src/io/api.jl:26-48` → stored-axis section of `api/store-io.md`.
- **URL schemes:** Local paths and URLs are accepted. Public `gs://` uses HTTPS;
  `s3://` additionally needs `using AWSS3` (`src/io/api.jl:35-37`). HTTPS is shown
  in the tutorial. No remote access was executed.
- **What is lazy?** Data arrays remain Zarr-backed. Ranges and implicit axes use
  arithmetic; a dense axis may scan IDs depending on manifest and validation
  (`docs/src/api/store-io.md:26-47`). Selection reads the required data chunks.
- **Which chunks are fetched?** `ChunkManifest`, `chunkof`, and `chunkbounds`
  expose touched cell-axis chunks (`docs/src/tutorials/store_io.jl:81-96`).
- **Does the selected lookup stay chunked?** No. The selected result has a
  compressed in-memory `CellLookup`; the source remains chunked
  (`docs/src/tutorials/store_io.jl:71-79`).
- **Absent cells:** Exact absent-ID lookup fails rather than inventing a value.
  Run S observed `BoundsError`. The store tutorial does not state this beside its
  selector list.
- **Result:** **Discoverable and locally executed.** Remote behavior remains
  **untested** by design.

### UC-047: Choose a cell-ID encoding and preserve interoperability

- **Discovery path:** `tutorials/store_io.jl:124-150` →
  `api/store-io.md:114-150` → rendered encoding docstrings.
- **Readable and writable encodings:** Dense, ranges, and implicit are public.
  Dense is always the fallback; ranges needs a sorted unique one-level axis;
  implicit needs a complete level (`src/io/encodings.jl:313-430`).
- **What does `:auto` choose?** Eligible ranges, otherwise dense
  (`src/io/api.jl:59-60`). Run S observed ranges for a complete sorted axis.
- **Integer adjacency versus rank:** `:step` joins consecutive raw integers;
  `:rank` joins consecutive valid cells. Rank is zero-based canonical position
  (`docs/src/tutorials/store_io.jl:140-150`; `src/io/encodings.jl:53-99`).
- **Sorted and unique?** Required for ranges; malformed or overlapping intervals
  raise `DGGSFormatError`.
- **External readers:** Structural readers can safely count `:step` intervals;
  `:rank` requires a grid-aware rank/select reader
  (`docs/src/tutorials/store_io.jl:140-150`).
- **Result:** **Discoverable and partly executed.** Dense and automatic ranges
  round-tripped in Run S. Malformed fixtures were not rerun.

### UC-048: Recognize an existing store and diagnose format errors

- **Discovery path:** `api/store-io.md:15-24,75-112,153-164` → rendered
  `StoreDescription`, convention, `Detection`, and `DGGSFormatError` docstrings.
- **Recognized conventions:** Zarr DGGS, xdggs, legacy HEALPix, and DKRZ appear in
  the convention list at `docs/src/api/store-io.md:84-99`.
- **Identifying metadata:** `StoreDescription` carries grid, level, encoding,
  variable names, orientation, and ellipsoid information
  (`docs/src/api/store-io.md:15-20,75-82`).
- **Conflicts and incomplete declarations:** Several detections may fire only if
  their decoded descriptions agree. Invalid formats carry a named check in
  `DGGSFormatError` (`docs/src/api/store-io.md:84-104,153-164`).
- **Inspect without data:** `describe_store` works from `StoreSnapshot` metadata;
  reading source values is not part of detection.
- **Result:** **Discoverable, static.** Run S exercised successful detection; the
  malformed convention fixtures were not rerun.

### UC-049: Update or incrementally populate a store

- **Discovery path:** `api/store-io.md:49-54` →
  `api/subzone-layout.md:18-30` → rendered `subzonestore` and `dggwrite!`
  docstrings at `src/io/api.jl:76-108`.
- **Create versus update:** `dggwrite` creates a cell-axis store and refuses
  overwrite. `subzonestore` creates or reopens the ancestor layout;
  `dggwrite!` fills whole ancestor columns.
- **Variable, region, or cube:** The normal writer writes a `DimArray` or
  `DimStack`. Incremental `dggwrite!` accepts one ancestor's full values or a cube
  made of complete columns (`src/io/api.jl:95-106`).
- **Must the axis exist?** For incremental writes, the subzone store and its fixed
  layout exist first.
- **Incompatible dimensions:** Partial columns and mismatched axes are rejected;
  one-dimensional existing layers are not overwritten.
- **Append, transactions, concurrency:** Cell-axis append and general transactions
  are not supported. Disjoint subzone columns can be written by independent tasks
  because each column is a chunk; same-column coordination is the caller's
  responsibility (`src/io/api.jl:83-87`; rendered `SubzoneStore` docstring at
  `ext/DiscreteGlobalGridsZarrExt/subzones.jl:33-42`).
- **Result:** **Discoverable, static.** Run S confirmed overwrite refusal. The
  incremental subzone path was not executed.

### UC-050: Store fine cells by ancestor and subzone

- **Discovery path:** `api/subzone-layout.md:1-25` → rendered `SubzoneLayout` and
  coordinate helper docstrings.
- **Eligible systems and levels:** The system must have sorted subtrees, and
  `0 <= ancestor_level <= level` (`src/io/subzones.jl:78-140`).
- **Why two dimensions?** Rows are positions within a subtree; columns are
  ancestor cells. One column is one independently readable/writable chunk
  (`docs/src/api/subzone-layout.md:3-12`).
- **Unequal descendants:** Short pentagon-rooted columns have trailing fill;
  readers remove padding from the cell axis
  (`docs/src/api/subzone-layout.md:14-22`).
- **Cell-to-coordinate mapping:** `subzoneindex` and `columnrow` return one-based
  `(column,row)`; `columnindices` maps a column to its complete-grid range
  (`src/io/subzones.jl:240-290`).
- **Regional reads:** `dggread(...; ancestors=...)` selects columns, so unwritten
  ocean columns need no stored chunk and region reads follow subtree chunking
  (`ext/DiscreteGlobalGridsZarrExt/read.jl:70-77`).
- **Result:** **Discoverable, static.** Exceptional-column fixtures were not rerun.

### UC-051: Register a custom store convention, encoding, or grid reference

- **Discovery path:** `api/store-io.md:15-24,84-150` → rendered registration and
  qualified hook docstrings.
- **Required hooks and metadata:** A convention supplies `detect`, `decode`, and
  optional `encode!`; an encoding supplies `cellaxis`, validation,
  `write_eligible`, and its name; a grid reference maps a store name to system
  construction (`docs/src/api/store-io.md:101-150`).
- **Detection conflicts:** All fired conventions must decode to the same
  description; disagreement is an error (`docs/src/api/store-io.md:84-88`).
- **Does registration affect writes?** Yes for conventions with `encode!` and
  named encodings with writer hooks. Read-only registration does not make an
  implementation writable.
- **Identity validation:** The encoding/grid arithmetic validates IDs, level,
  sortedness, uniqueness, and declared counts; invalid stores raise a named
  `DGGSFormatError` check.
- **Which names stay qualified?** `detect`, `decode`, `encode!`, encoding hooks,
  and ID arithmetic are deliberately qualified (`docs/src/api/store-io.md:101-150`).
- **Result:** **Discoverable, static.** No registry mutation was performed because
  it would add little evidence beyond the existing isolated fixture tests.

## H. Compute out of core and distribute chunks

### UC-052: Apply a neighborhood kernel to a stored global dataset

- **Discovery path:** `tutorials/out_of_core.jl:18-65,78-103` →
  `api/chunk-sweep.md:7-113`.
- **Halo width:** One for the supported one-ring sweep. Wider halos provide input
  context but do not enlarge callbacks (`docs/src/api/chunk-sweep.md:77-85`).
- **Does the plan follow storage chunks?** Yes, including irregular stored chunk
  layouts (`docs/src/api/chunk-sweep.md:28-39`).
- **Output shape and chunking:** `dest` must hold one result per source cell and
  support range writes on the cell dimension. Stored output chunks should align
  with disjoint ownership for concurrent writers (`src/chunks.jl:686-707` and
  `docs/src/api/partitioning.md:69-82`).
- **Are halo values read from source?** Yes. Each chunk cube copies owned and halo
  values from the input, with at most one read per touched storage chunk
  (`src/chunks.jl:384-405,441-490`, rendered docstrings plus internals page).
- **Is in-place aliasing safe?** The public docs do not say. Run C used the same
  data as source and destination and got a different result from the eager
  reference because later chunks read earlier writes.
- **Result:** **DOC GAP with runtime failure.** Non-aliased Run C matched eager
  output across all chunk boundaries. Aliasing is UW-006.

### UC-053: Write a custom chunk callback

- **Discovery path:** `api/chunk-sweep.md:52-75` → rendered `ChunkCube`,
  `ownedindices`, `localindices`, `axisindices`, and `chunkhalo` docstrings.
- **Index spaces:** `ownedindices` are caller-axis positions;
  `localindices` are owned positions inside the chunk cube;
  `axisindices[k]` maps every chunk-buffer cell back to the caller axis
  (`docs/src/api/chunk-sweep.md:54-66`).
- **Halo placement:** The buffer is sorted as lower halo, contiguous owned run,
  upper halo. `chunkhalo(cc)` is the complement of `localindices(cc)`
  (`src/chunks.jl:312-367`).
- **Can callbacks retain buffers?** The public docs do not state lifetime or
  ownership. Current code allocates a new in-memory block for each callback, but
  callers should not rely on that private behavior. Link UW-002.
- **Extra dimensions and empty chunks:** Extra dimensions are retained in each
  chunk cube. Plans contain nonempty owned chunk ranges; public docs do not call
  out empty input behavior.
- **Result:** **DOC GAP and executed.** Run C verified all owned cells were written
  once and the three index spaces agreed.

### UC-054: Run wider or repeated stencils across chunk boundaries

- **Discovery path:** `api/chunk-sweep.md:77-85` →
  `api/neighbors.md:75-83` → `tutorials/out_of_core.jl:78-103`.
- **Natural wider request:** Missing radius/exact-ring selector on sweep methods.
  Link API-005.
- **Required input halo:** At least the requested graph reach for a future direct
  sweep. The current supported sweep needs width one.
- **Diffusion versus radius two:** Two passes include intermediate path
  multiplicity and reintroduce centers; they are not one uniform radius-2 disk
  operation (`docs/src/tutorials/stencils.jl:129-133`).
- **Intermediate halo values:** A second pass needs first-pass values computed on
  the surrounding region, so a one-time input halo alone cannot replace a pass or
  exchange of intermediate values. No built-in multi-pass orchestration exists.
- **Result:** **API GAP, statically confirmed and executed.** Run C proved
  `halo=2` still returns the one-ring result. Link API-005, DOC-014, SW-001,
  SW-002.

### UC-055: Assign chunks to threads or distributed workers

- **Discovery path:** `api/partitioning.md:7-87` →
  `tutorials/out_of_core.jl:78-103`.
- **Does partitioning launch workers?** No. It creates logical assignments; an
  executor maps them to tasks or processes (`docs/src/api/partitioning.md:7-15,58-60`).
- **Worker setup:** Each worker opens its own data handles and needs packages used
  by its kernel. Optional partition libraries are required only on the
  coordinator (`docs/src/api/partitioning.md:143-168`).
- **Disjoint writes:** Plan pieces own disjoint cell-axis ranges. Stored output
  must also align those with physical write chunks (`docs/src/api/partitioning.md:69-82`).
- **Serialization:** `ChunkPartition` contains ordinary numeric arrays and can be
  serialized without the plan, store, or backend library
  (`docs/src/api/partitioning.md:11-15`).
- **Nested threading and duplicate reads:** The example sets `threaded=false`
  inside spawned chunk tasks. Each chunk currently reads its halo separately;
  grouping does not itself create a shared cache (`docs/src/api/partitioning.md:74-87`).
- **Result:** **Discoverable and partly executed.** Run P verified exact logical
  task coverage. Process-level execution was not repeated.

### UC-056: Balance work and reduce shared source reads

- **Discovery path:** `api/partitioning.md:17-60,153-252`.
- **Objectives:** `WeightedContiguous` balances contiguous traversal work;
  Metis reduces shared-source graph edges; KaHyPar minimizes weighted source
  connectivity; Scotch maps the shared-source graph to capacity-weighted workers
  (`docs/src/api/partitioning.md:153-188,207-246`).
- **Optional packages:** None, Metis, KaHyPar_jll, and Scotch respectively
  (`docs/src/api/partitioning.md:155-168`).
- **Weights, capacities, imbalance:** Work weights estimate compute; source
  weights estimate read/transfer cost; capacities are relative compute rates,
  not memory; imbalance relaxes native load bounds
  (`docs/src/api/partitioning.md:17-42,184-205`).
- **Seed:** Each native backend's algorithm object has a seed. Same inputs and
  build are repeatable; labels may change across library versions
  (`docs/src/api/partitioning.md:223-252`).
- **Edge limiting:** `maxedges` bounds Metis/Scotch projected graph expansion and
  can reject too-large expansions; KaHyPar keeps hyperedges directly.
- **Result:** **Discoverable and partly executed.** Run P confirmed weighted
  contiguous coverage and the missing-backend error. Native heuristics were not
  run.

### UC-057: Partition a chunked regridding plan

- **Discovery path:** `api/partitioning.md:89-151` →
  `api/regridding.md:61-78` → GlobalRegridding dependency docstrings.
- **Which chunks are tasks and dependencies?** Destination chunks are work rows;
  source chunks are dependency resources (`docs/src/api/partitioning.md:89-111`).
- **How are assignments executed?** Partitioning returns IDs and row positions;
  the caller slices or schedules the plan. The executor opens sources and writes
  destinations (`docs/src/api/partitioning.md:138-151`).
- **Plan reuse and spill behavior:** Partitioning reads the existing dependency
  relation and does not rebuild weights. It does not change the plan's storage
  object, cache, or spill directory (`docs/src/api/partitioning.md:108-111`).
- **Who coordinates writes?** The caller's executor. It owns process placement,
  retries, data exchange, completion, and output writes
  (`docs/src/api/partitioning.md:143-151`).
- **Result:** **Discoverable and executed at planning level.** Run P partitioned a
  48-destination-chunk lazy plan, covered every task once, and exposed each
  partition's source sets. Applying separately partitioned regrid pieces is not a
  public execution API and was not presented as one.

## New walkthrough findings

### UW-001: State when a DGGS dimensional source still needs `from`

- **Use case:** UC-037
- **Kind / severity:** DOC GAP / High
- **Evidence:** `docs/src/api/regridding.md:29-35` says a dimensional source can
  describe its spatial axes but does not distinguish raster axes from `Cells`.
  Run R omitted `from` from a `DimArray` over `Cells` and received
  `ArgumentError`; explicit `from` passed. The rejection is at
  `lib/GlobalRegridding/src/api.jl:290-310`.
- **Correction:** Say that raster X/Y axes can infer `RasterGrid`, while a DGGS
  value array currently needs an explicit source grid or space. Add one DGGS
  source example with `from`.

### UW-002: Specify callback buffer ownership and lifetime

- **Use cases:** UC-034, UC-053
- **Kind / severity:** DOC GAP / Medium
- **Evidence:** Sweep docs define callback contents and chunk docs define index
  spaces, but neither says whether a callback may retain or mutate neighbor
  sequences, needs rings, or `chunkcube(cc)` data. Current allocations at
  `src/engine/neighborhood.jl:463-480` and `src/chunks.jl:458-468` are private
  implementation details.
- **Correction:** Give each callback argument an explicit borrowed or owned
  contract. State whether retention and mutation are supported across callbacks
  and threads.

### UW-003: Show how owned cells map into a grown region

- **Use case:** UC-033
- **Kind / severity:** DOC GAP / Medium
- **Evidence:** `docs/src/api/boundaries.md:79-94` explains membership but not
  value remapping. `grow` returns a newly sorted `CellVector` and no ownership map
  (`src/engine/region_algebra.jl:30-54`).
- **Correction:** Add a small example that builds a grown input field, maps the
  original cells by identity, computes only owned outputs, and restores original
  order.

### UW-004: State cost-distance missing and unreachable behavior

- **Use case:** UC-035
- **Kind / severity:** DOC GAP / Medium
- **Evidence:** `docs/src/tutorials/stencils.jl:290-364` models a barrier with
  `Inf` but does not define `missing` input or disconnected output. The example
  initializes all distances to `Inf`, so unreachable cells remain `Inf`.
- **Correction:** State that the tutorial function is example code, define valid
  costs, say how to preconvert missing cells to barriers, and state that
  unreachable cells remain `Inf`.

### UW-005: Explain how mixed-level targets become a regridding space

- **Use case:** UC-041
- **Kind / severity:** DOC GAP / High
- **Evidence:** `docs/src/api/regridding.md:29-32` lists `MultiOrderCellSet` as a
  target without describing output granularity. The adapter converts it through
  `CellVector(set)` and `PartialGrid` at `src/regridding.jl:289-294`, producing a
  one-level leaf space rather than one value per stored mixed-level member.
- **Correction:** State the expansion level, output order, and value meaning. Add
  a mixed-level example that compares stored members with returned leaf cells.
  State empty-target behavior after validating it.

### UW-006: Forbid source/destination aliasing in chunked sweeps

- **Use case:** UC-052
- **Kind / severity:** DOC GAP / High
- **Evidence:** `src/chunks.jl:686-707` describes destination shape but not
  aliasing. Run C used `parent(A)` as `dest`; the result differed from the eager
  one-ring sweep because later chunk loads observed earlier writes. A separate
  destination matched exactly.
- **Correction:** State that `dest` must not alias the swept data or any requested
  value field unless the implementation adds snapshot semantics. Add a small
  rejection or non-alias test.

## Coverage and status

The `Source read` column records whether answering the case required implementation
source beyond the public manual and rendered docstrings. A `Yes` means that public
documentation first left a contract question unanswered.

| Use case | Source read | Result |
|---|---:|---|
| UC-025 | No | Discoverable |
| UC-026 | No | API GAP (`API-005`) |
| UC-027 | No | Discoverable |
| UC-028 | No | Discoverable |
| UC-029 | No | Discoverable |
| UC-030 | No | Discoverable |
| UC-031 | No | Discoverable |
| UC-032 | No | Discoverable |
| UC-033 | Yes | DOC GAP (`UW-003`) |
| UC-034 | Yes | DOC GAP (`UW-002`) |
| UC-035 | Yes | DOC GAP (`UW-004`) |
| UC-036 | No | Discoverable |
| UC-037 | Yes | DOC GAP (`UW-001`) |
| UC-038 | No | Discoverable |
| UC-039 | No | Discoverable |
| UC-040 | Yes | Discoverable; concurrent use remains untested |
| UC-041 | Yes | DOC GAP (`UW-005`) |
| UC-042 | No | Discoverable; runtime untested |
| UC-043 | Yes | Discoverable; spill behavior untested |
| UC-044 | No | Discoverable; specialized cases untested |
| UC-045 | Yes | Discoverable and executed |
| UC-046 | No | Discoverable; remote behavior untested |
| UC-047 | No | Discoverable |
| UC-048 | No | Discoverable |
| UC-049 | No | Discoverable |
| UC-050 | No | Discoverable |
| UC-051 | No | Discoverable |
| UC-052 | Yes | DOC GAP and runtime failure (`UW-006`) |
| UC-053 | Yes | DOC GAP (`UW-002`) |
| UC-054 | No | API GAP (`API-005`) |
| UC-055 | No | Discoverable; process execution untested |
| UC-056 | No | Discoverable; native heuristics untested |
| UC-057 | No | Discoverable at planning level |

| Range | Cases | Question-level walkthrough | Runtime evidence |
|---|---:|---|---|
| Neighborhoods and boundaries | UC-025–035 | Complete | Runs N and C; UC-035 static |
| Regridding | UC-036–044 | Complete | Run R; raster download, spill, and specialized edge fixtures untested |
| Store I/O | UC-045–051 | Complete | Run S; remote and subzone cases static |
| Chunk execution and partitioning | UC-052–057 | Complete | Runs C and P; native backends and process execution untested |

Every assigned use case has a concrete public-doc path, an answer for each
inventory question, a source-reading marker where needed, and a result label.
Confirmed unsupported wider sweeps remain API-005. This report adds no new API
gap: the other residuals are documentation contracts or explicitly untested
behaviors.
