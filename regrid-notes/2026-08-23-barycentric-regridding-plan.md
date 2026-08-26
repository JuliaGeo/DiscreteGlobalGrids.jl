# Generic barycentric point regridding implementation plan

- Date: 2026-08-23
- Status: in progress, authoritative for this work
- Scope: `lib/GlobalRegridding` and the qualified `DGGSpace` extensions in this
  repository

This plan consolidates the geometry research and the subsequent API,
Copernicus DEM, chunk-discovery, and performance discussion. The supporting
notes remain useful evidence:

- `regrid-notes/generic-barycentric-patch-regridding.md`
- `regrid-notes/esmf-regridding-methods.md`

Where either note proposes `sampleelement`, `patchsites`, or one universal
`stencilreach`, this plan supersedes it. The existing regridding simplification
plan remains authoritative for conservative regridding and shared spatial
infrastructure; this plan owns only the point-method redesign. Its fused tile
weights are an explicit point-method exception to the older plan's assumption
that an exact chunk dependency graph can always be materialized before weights.

## Outcome

Add a generic `BarycentricPoint` method for fields sampled at source-cell
sample sites: a point sample at each destination point, interpolated between
the source sample sites around it. It must:

1. give true tensor-product Q1 bilinear interpolation on a raster;
2. give P1 barycentric interpolation on the triangular dual of a hexagonal
   source grid;
3. support explicitly selected coordinates on other convex dual polygons;
4. construct every destination stencil once, independently of source chunking;
5. load exactly the source chunks used by those stencils, without a one-chunk
   halo; and
6. make source rims, partial grids, degeneracies, and poles explicit policy
   cases.

The intended data flow is:

```text
destination tile
      |
      | one pass over destination sample points
      v
weightsat!(row, sampler, point)
      |
      | source local indices + weights
      v
partition by chunkat(source_space, local_index)
      |
      +--> source chunk 4  --> WeightBlock
      +--> source chunk 9  --> WeightBlock
      `--> source chunk 10 --> WeightBlock
      |
      v
one cached TileWeights with an exact, sorted dependency list
```

The unavoidable work is proportional to the number of destination cells and
the number of nonzeros ultimately emitted. It must not also be multiplied by
the number of candidate source chunks.

## Decisions this work must preserve

1. **Conservative regridding stays separate.** `Conservative()` continues to
   discover polygon overlaps and build area weights through its present path.
   It is not expressed as a point stencil and is not folded into the APIs below.
2. **The caller chooses the sampling semantics.** `BarycentricPoint` means a
   point sample at each destination sample site. `Conservative` means
   an area operation. No metadata-driven choice between them is part of this
   work.
3. **Copernicus DEM is intended to use point interpolation eventually.** Its
   published pixels are pixel-is-point. The current production conservative
   path remains unchanged until the point method passes the Copernicus-specific
   correctness and performance gates below.
4. **Other DEMs that encode area means remain conservative.** This plan does not
   invent a `valuesampling`, `samplepoint`, or point-to-area conversion API.
5. **The stable method-facing seam is a stencil, not a polygon.** A source space
   may expose or internally construct dual cells, but the regridding executor
   consumes source local indices and scalar weights.
6. **Chunking cannot affect mathematics.** The full stencil is determined before
   it is partitioned into source chunks. No source-chunk boundary may clamp,
   truncate, or otherwise change a stencil.
7. **No neighbor-chunk buffering.** Source dependencies are the exact sorted
   union of chunks owning nonzero stencil entries.
8. **The first implementation is non-extrapolating.** A valid stencil uses
   nonnegative P1, Q1, or convex generalized barycentric weights. The default
   at a geometric rim or unsupported degeneracy is unmapped/missing. Nearest
   fallback may be explicit; silent clamping is not barycentric interpolation.
9. **Patch recovery is a later, distinct method.** `BarycentricPoint` will not
   grow patch rings or fit polynomials.

## Vocabulary

Three index spaces appear below. The words are `DiscreteGlobalGrids`' own
(`localindex`/`globalindex`, `Index(Local())`/`Index(Global())`).

- **Local index** — a cell's place in the collection a space wraps,
  `1:ncells(src_space)`. It is what `cellat(space, p)` answers, what
  `cellcentroid(space, i)` and `chunkat(space, i)` take, and the space every
  stencil is written in. A `DGGSpace` may wrap a `PartialGrid`, so this is not
  the complete level's numbering.
- **Chunk-local index** — the local index of one chunk read as a partial grid.
  A chunk is a smaller collection, not another kind of index: chunk-local `j`
  satisfies `ownedindices(space, k)[j] == i`, and when the chunk is a
  contiguous range, as every `DGGSpace` chunk is, that is
  `i - first(ownedindices(space, k)) + 1`. `WeightCOO` rows and columns are
  chunk-local, rows within the destination tile and columns within the source
  chunk.
- **Global index** — the cell's place in the complete level (`globalindex`).
  It equals the local index only on a complete grid. This plan uses it in one
  place: CopDEM topology is constructed on the complete level and translated to
  local indices before any stencil is emitted.

Four further words are used precisely:

- **Tile vs chunk** — a destination *chunk* is the space's unit,
  `ownedindices(space, k)`; a *tile* is the lazy array's write block, a run of
  destination chunks (`DestTiling` in `lib/GlobalRegridding/src/lazy.jl`). This
  plan builds weights per tile and reads sources per chunk.
- **Sample site** — where a source value is taken to sit:
  `cellcentroid(space, i)` unless a space says otherwise. Sample site is the
  word used throughout; "centroid" appears only where the centroid specifically
  is meant.
- **Rim vs `border`** — `border(cv)` is a set of cells and `halo(cv)` the cells
  outside them. A *rim* is an area: the part of the border cells lying outside
  the dual complex, where no dual cell exists. "Rim policy" is this plan's term
  for what happens on a rim.
- **Partial grid** — a collection holding fewer cells than its level. A required
  sample site outside the collection has no local index at all, which is a rim
  case rather than a reason to substitute another cell.

## Geometry and mathematical contract

Source cell values are treated as values at source sample sites. Normally the
site is the cell centroid. Interpolation happens on the dual complex whose
nodes are those sites:

| primal source grid | usual primal-vertex valence | dual cell | basis |
|---|---:|---|---|
| quadrilateral raster | 4 | dual quadrilateral | tensor Q1 |
| hexagons | 3 | dual triangle | mean value |
| triangles | about 6 | dual polygon | named generalized coordinates |
| nonconforming CopDEM transition | 3 or 4 | triangle or quadrilateral | mean value or Q1 |

For a valid destination point `p`, the method emits

```math
\hat f(p)=\sum_i w_i(p)f_i, \qquad \sum_i w_i(p)=1.
```

In a local two-dimensional chart it also reproduces affine coordinates. Q1 on
a rectangular raster additionally reproduces `a + bx + cy + dxy`. The dual
cell's kind, rather than merely its vertex count, chooses the basis:

- `Bilinear`: inverse isoparametric coordinates followed by tensor Q1 weights;
- `MeanValue`: mean-value coordinates on a supported convex polygon, which on a
  triangle are that triangle's barycentric coordinates.

A four-node dual cell is not automatically `Bilinear`. A raster or source-space
specialization must assert that kind. Likewise, an arbitrary n-gon is
not silently triangulated unless a future method names and documents that
triangulation policy.

On the sphere, containment remains unit-spherical. Raster Q1 uses the raster's
declared native chart. Unstructured dual-cell weights use a stable local
tangent chart and reject folded, antipodal, or otherwise nonlocal cells. Raw
global longitude/latitude is not the generic coordinate system.

## Core API direction

### Weight rows and samplers

The hot API must be allocation-free in steady state:

```julia
weightsat!(row, sampler, p) -> status
```

`row` is a reusable `WeightRow`: it owns parallel vectors of source local
indices and `Float64` weights, and is cleared on entry. Each worker owns one.
`status` distinguishes at least mapped from unmapped; diagnostic builds may
distinguish outside coverage, rim, and degenerate geometry.

Preparation is once per plan or source space, not once per destination:

```julia
sampler(method, space) -> Sampler
```

A `Sampler` is the source space made ready to be sampled at points. It may
cache raster axis locators, topology tables, or other immutable query state,
and a `cellfield` of source sample sites (`src/engine/cellfield.jl`, whose
`known=` accepts sites already computed) is the natural thing for it to hold.
Sampler state must be safe for concurrent reads; task-local scratch belongs to
the worker's `WeightRow`. The exact names may remain internal initially,
but the separation between preparation and the per-point hot loop is required.

The qualified extension fallback is conceptually:

```julia
weightsat!(row, method::BarycentricPoint, src_space, p)
```

Source spaces with a faster fused algorithm specialize it. The returned
indices are the source space's local indices, `1:ncells(src_space)`; source
chunking is not an input.

### Dual cells

Do not add `sampleelement`. If a reusable geometry layer proves useful, call
the point lookup `dualcellat`, consistent with `cellat`:

```julia
dualcellat(sampler, p) -> DualCell
```

A `DualCell` is the polygon whose vertices are source sample sites and which
contains `p`, the dual cell of a primal vertex: a small structured object
containing ordered source local indices, sample-site geometry, and an explicit
basis kind. It is not required to be a `GeoInterface.Polygon`, and Q1/P1
identity must not be erased by converting everything to a generic polygon. A
natural polar construction could have many nodes, so the representation must
not promise a fixed size.

This is initially an internal factoring seam. The stable extension point is
`weightsat!`; only export `dualcellat` after raster, a conforming
DGGS, and CopDEM demonstrate that one contract is both sufficient and fast.

### Method types

Point execution dispatches on `outputsampling(method) isa Points`, the sampling
trait already carried by every method, so it never tests concrete method names.
`BarycentricPoint` reports `Points()`. `NearestCell` reports `Points()` too and
may migrate to the same execution path after the new path is proven. Area
methods report `Intervals(Center())`; `Conservative` and custom existing
methods remain on `buildweights!` unless they explicitly opt in. No new trait
is introduced.

Keep `BilinearPoint` during migration. It may share the raster Q1 kernel, but do
not remove or silently alias it until its current rim-clamping behavior and
downstream uses have an explicit compatibility decision. The new generic API
is named `BarycentricPoint` because “bilinear” is only correct for its raster
Q1 specialization.

## Raster behavior

The raster specialization is the reference fast path:

1. prepare both monotonic chart axes and periods once;
2. convert a destination point through `chartcoords` once;
3. locate its two bracketing sample coordinates on each axis;
4. emit the Cartesian product's at most four nonzero weights; and
5. wrap periodic axes exactly.

Unlike the current `_locate`, the new method does not clamp an arbitrarily
distant point to an edge sample. The dual Q1 domain ends at the outermost source
sample sites on a nonperiodic axis. The half-cell strip between those sites and
the source cell boundary is a geometric rim and follows the explicit rim
policy.

Repeated indices on a one-cell axis are combined before emission. Zero weights
are dropped. Identity at every source sample site and exact Q1 reproduction are
hard acceptance laws.

## Conforming DGGS behavior

For cell-centred DGGS data, each primal vertex produces a dual cell
from the source cells incident on that vertex. A host source cell contributes
one candidate dual cell per incident primal vertex; the dual cell containing
`p` supplies the stencil.

The current `RegridSpace` contract does not expose primal vertex incidence, and
an ordered `Vertex()` one-ring alone does not identify every incident-cell set
at arbitrary valence. Therefore this phase must prototype the minimum topology
hook behind `DGGSpace` on at least:

- one quadrilateral or raster-like grid;
- one hexagonal grid, including a pentagonal defect if present; and
- one higher-valence dual polygon.

Prefer a high-level, query-oriented implementation over manufacturing a
numbering of primal vertices. Possible internal factorizations are `dualcellat`
or an iterator of dual cells incident on a host cell. Do not export a
lower-level vertex-incidence API until those systems show that it removes work
rather than merely moving it.

For a regular hexagonal source, the required result is a three-source P1
triangle and O(1) or bounded local topology work per destination. Coordinate
matching of complete cell polygons is acceptable only as a correctness oracle,
not as the production hot path for million-cell tiles.

## Copernicus DEM specialization

CopDEM must specialize the point query. Its row-dependent longitude counts and
hanging vertices are real topology, not storage artifacts.

### Ordinary rows and tile seams

Within one latitude band, adjacent sample rows share the same longitude
lattice. The dual cells are four-site Q1 rectangles. A 1-degree COG seam
inside that band changes chunk ownership only; it does not change the stencil.

### Band transitions

At a change in `ncols`, construct the strip between adjacent sample rows from
the union of their exact longitude breakpoint sets:

- a breakpoint present on only one row yields a three-site P1 triangle;
- a breakpoint common to both rows yields a four-site Q1
  quadrilateral/trapezoid.

Use CopDEM's integer row coordinates and `_facing` arithmetic. Do not infer
these dual cells by intersecting floating cell rings: a coarse cell can have a
hanging vertex in the interior of its edge, so walking only its published
corners misses valid dual cells. Off the poles, lookup and stencil size must
stay O(1), with three or four source entries.

### Partial grids and rims

Construct topology in the complete CopDEM level's index space — the global
index, the one place this plan uses it — then translate each node to its local
index through the actual `DGGSpace`/partial-grid membership. If any required
node has no local index, it is a required sample site outside the collection:
treat it as a rim and apply the selected explicit rim policy.
Do not replace the missing node with a cell from a different chunk or silently
renormalize geometry.

Missing *field values* are a separate executor concern. For nonnegative
barycentric weights, the existing `Weighted(threshold)` policy can continue to
measure valid weight and renormalize it. `Weighted(1)` is the strict complete
stencil choice. Geometry-rim fallback happens before weights exist.

### Poles are a decision gate

The natural dual cell at a pole contains the entire polar source row:
129,600 GLO-30 sites or 43,200 GLO-90 sites. Evaluating or storing such a
coordinate set is necessarily linear in that row length and violates the
ordinary constant-stencil performance model.

There is also a semantic issue: `cell_centroid` is the midpoint of the clamped
polar cell box, while the original DEM post coordinate and that midpoint differ
in the pole rows. Before enabling CopDEM point regridding at the poles, choose
and test one explicit policy:

1. natural all-row polar coordinates, accepting their cost;
2. an explicitly named nearest/local polar fallback; or
3. another documented reconstruction that defines both its virtual geometry
   and its sample-coordinate meaning.

Do not make this choice implicitly in the generic method. Until it is made,
CopDEM polar queries are unmapped under `BarycentricPoint` and the current
production path remains available.

## Exact source-chunk planning

### Current problem

`buildblock` currently builds one `(destination tile, source chunk)` pair at a
time. Each pair calls `buildweights!` over all destination indices. A point
method consequently repeats point location for every candidate source chunk.
With roughly a million destination cells this is unacceptable even when the
per-point stencil itself is cheap.

The current cap index plus `supportradius` is a sound broad phase when its
radius is conservative, but it may return false-positive chunks. It cannot by
itself give exact point-stencil dependencies.

### Cost model

For `N` destination cells and maximum stencil size `k`, the target is `O(Nk)`
arithmetic and `O(Nk)` final sparse storage. Ordinarily `k = 4` for rasters and
CopDEM and `k = 3` for a hexagonal source. Thus a million-cell destination tile
means roughly three or four million emitted entries, not a million searches per
candidate source chunk. The exact dependency list is obtained during this same
pass; it is not a second scan.

This still deserves careful engineering, but it is bounded and parallelizable.
The hot loop must avoid cell-polygon construction, spatial-tree searches where
native arithmetic exists, dictionaries per destination, and heap allocation per
stencil. Threading should use bounded destination ranges with worker-local
buffers and deterministic merging. Measure the serial kernel first, then retain
threading only where it improves the real large-tile benchmark.

### Tile weights as the build unit

Add point-specific tile weights whose atomic build unit is one destination
tile, not one chunk pair:

```julia
TileWeights(
    sourcechunks::Vector{Int},      # sorted, exact
    blocks::Vector{WeightBlock},    # same order
)
```

Construction performs one destination pass. For every nonzero stencil entry it:

1. obtains the owning source chunk with `chunkat(src_space, i)` from the
   source local index `i`;
2. obtains the source's chunk-local index, its local index in that chunk read
   as a partial grid;
3. appends `(dst_tilelocal, src_chunklocal, weight)` to that chunk's builder;
4. combines duplicate source cells deterministically; and
5. finalizes nonempty blocks in ascending source-chunk order.

A contiguous `ownedindices` range makes the chunk-local index arithmetic, as
the vocabulary section states. Other spaces may use one cached `indexmap` per
contributing chunk. Add a more specialized chunk-local hook only if profiles
show the existing maps dominate.

The implementation should process destination rows in bounded batches and use
reusable buffers. It need not retain a second whole-tile copy of all stencils. The
final sparse weights are inherently O(nonzeros); temporary memory must also be
O(nonzeros in the tile) or smaller, not O(destination cells × candidate chunks).

### Lazy execution and caching

For point methods, lazy execution obtains a tile's weights before source reads.
Their `sourcechunks` replaces `_connectedsource!`'s cap candidates for that
tile, so only exact dependencies are read. Applying blocks remains
deterministic in ascending source-chunk order.

The exact dependency manifest and the weights are deliberately built together.
For a generic point method, discovering exact dependencies requires evaluating
the same stencils that define the weights; discarding those weights would force
a second million-cell pass. `ChunkedPlan` and lazy-array construction still read
no source values and need not build all destination tiles up front. The first
request for a destination tile builds that tile's geometry-only weights before
loading any source chunk.

A tile's weights are the cache/locking unit: concurrent requests for one tile
must build its stencils once. Their complete weight footprint counts against
the existing weight budget. If spill storage remains pair-oriented, store one
small tile manifest plus its blocks; do not re-run the destination pass merely
to reconstruct an evicted dependency list.

`DirectPlan` may use the same stencil loop with one whole-domain source index.
The conservative chunk-pair cache and executor remain unchanged.

If a whole-run scheduler requires a graph before any destination is requested, it
has two honest choices: accept a no-false-negative support superset but delay
source reads/prefetch until the exact tile manifest exists, or explicitly build
the relevant tiles' weights eagerly. It must not obtain an “exact” graph
by buffering adjacent chunks or repeat stencil construction later.

### Optional support-query acceleration

Do not make a support query a prerequisite for correctness. Exact dependencies
already fall out of the fused stencil pass, and every requested destination
still needs a value.

If later measurements show that ahead-of-time graph construction, scheduling,
or prefetch needs a cheaper geometric summary, add an optional opaque query.
It is not named until P3-P5 profiles justify it; if it comes, it sits beside
the existing `supportradius` as `chunksupport`.

The support object is not necessarily a polygon. A raster may use index ranges;
CopDEM may use band and rational-breakpoint arithmetic. It answers how one
source chunk relates to a queried region — `Outside`, `Inside`, or `Unknown` —
and, exactly, whether a single destination point uses that chunk. `Outside`
must be a proof of no contribution, `Inside` a proof for the whole queried
region, and `Unknown` must fall back to exact stencil evaluation. False
negatives are forbidden. This optional API applies to point methods only;
conservative regridding keeps overlap discovery.

## Missing data, boundaries, and degeneracies

The first `BarycentricPoint` implementation supports only interpolation inside
a valid convex dual cell. It must diagnose or map to the selected
rim policy for:

- source coverage outside the dual complex;
- a partial grid missing a required node;
- repeated sample sites or zero-area triangles;
- failed or nonunique inverse Q1 maps;
- unsupported nonconvex/self-intersecting polygons;
- tangent-chart folding or nonlocal spherical cells; and
- the unresolved CopDEM polar construction.

Default: emit no row, so the destination is missing under the existing executor
semantics. Initial optional fallback: nearest source sample, explicitly named in
the method configuration and reported in diagnostics. Defer one-sided signed
extrapolation and ghost-node schemes.

All accepted v1 barycentric weights are nonnegative. This makes the current
valid-weight accounting meaningful. Signed methods must not reuse that coverage
logic without a separate design.

## Patch-like regridding: later plan

ESMF patch is quadratic recovery over topological stars followed by blending
with the containing dual cell's low-order basis. It is not “more barycentric
neighbors.” A faithful implementation needs:

- stars of incident dual cells around every corner of the containing dual cell;
- scaled local tangent coordinates;
- rank-revealing QR/SVD for a degree-two fit;
- an explicit quadratic-to-linear-to-missing/nearest fallback; and
- missing-data semantics that remain valid for signed weights.

A generic moving-least-squares method is also possible, but must be named as
MLS/patch-like rather than ESMF patch. Neither belongs in the first
`BarycentricPoint` implementation. The ESMF report is the mathematical basis
for a separate follow-up plan after the low-order dual cell and topology APIs
are proven.

## Implementation tasks and gates

```text
P0 baseline
  -> P1 point contracts
  -> P2 raster Q1 correctness
  -> P3 fused tile-weight execution
  -> P4 conforming DGGS dual cells
  -> P5 CopDEM nonconforming specialization
  -> P6 production opt-in and documentation

Patch/MLS and the optional support query start only after P3-P5 measurements.
```

## Progress

- P0 `df1d8b8` — the baseline: `BilinearPoint`'s values pinned, the repeat of
  point location per candidate chunk counted, and one large tile priced.
  `regrid-notes/2026-08-25-p0-point-baseline.md`.
- P1 `e57bd82` — the point seam: `WeightRow`, `sampler`, `weightsat!`, dual
  cells and the coordinate kernels.
  `regrid-notes/2026-08-25-p1-point-contracts.md`.
- P2 `874bf74` — the raster specialization: Q1 stencils on prepared chart axes,
  with every intentional difference from `BilinearPoint` documented.
  `regrid-notes/2026-08-25-p2-raster-q1.md`.
- P3 — realised by the S2 and S3 commits `8c05a2f` and `2800ba3`, and recorded
  closed in `d227de2`: the fused build unit and its exact reads: one
  `TileWeights` per destination tile, cached and locked by tile number, reading
  the source chunks its stencils name and no others.
  `regrid-notes/2026-08-25-s2-tile-weights.md`,
  `regrid-notes/2026-08-25-s3-exact-reads.md`, and the gate's measurement in
  `regrid-notes/2026-08-25-p3-fused-tile-execution.md`.
- P4 `2634061` — dual cells behind `DGGSpace`: the host's ring-ordered vertex
  fans, an azimuthal-equidistant chart centred on the queried point,
  `MeanValue` on every cell, the `hasdualcells` trait and the `sampler`
  refusal of a space with neither chart nor dual cells.
  `regrid-notes/2026-08-26-p4-dggs-dual-cells.md`.
- P5 `2744b85` — Copernicus DEM row stencils: Q1 rectangles within a band, P1
  triangles and Q1 trapezoids at band edges from exact integer arithmetic,
  rims on holdings, poles `WeightsDegenerate` with the polar policy unchosen.
  `regrid-notes/2026-08-26-p5-copdem-point-stencils.md`.
- P6 next. The polar decision gate is open.

### P0 — baseline and instrumentation

- Record current `BilinearPoint` values at raster interiors, periodic seams,
  nonperiodic rims, and outside coverage.
- Add a counting test that demonstrates current chunk-pair point lookup scales
  as destination cells times candidate chunks.
- Record time, allocations, peak weight memory, and source reads for a
  representative large destination tile, without changing production output.

**Gate:** the baseline makes both the repeated-work problem and intended
behavior changes observable.

### P1 — point contracts and geometry kernels

- Add the `BarycentricPoint` method reporting `Points()`, the reusable
  `WeightRow`, and the `sampler` preparation seam.
- Implement and test P1 triangle, inverse Q1 quadrilateral, and convex
  mean-value coordinate kernels independently of spaces and chunks.
- Add the internal `DualCell` representation only if the three
  kernels benefit from it; do not export it yet.
- Preserve `Conservative` and custom `buildweights!` dispatch unchanged.

**Gate:** kernels reproduce constants and their promised coordinate fields,
drop zero weights, reject degeneracies, and allocate nothing once the
`WeightRow` is warm.

### P2 — raster specialization

- Prepare chart axes once and implement at-most-four-entry Q1 stencils.
- Cover ascending/descending axes, periodic seams, one-cell axes, source sample
  identity, and explicit nonperiodic rims.
- Compare `BarycentricPoint` with `BilinearPoint` in the interior and document
  every intentional edge difference.

**Gate:** raster Q1 laws pass in direct and chunked correctness tests, before
the performance cutover.

### P3 — fused tile-weight execution

- Add the tile-level build/cache unit and exact dependency manifest.
- Route `BarycentricPoint` through it in lazy execution; retain the old path for
  conservative and unopted custom methods.
- Partition once by source local index and apply exact source blocks in
  ascending order.
- Preserve eager/lazy and source-rechunking equivalence.
- Add call-count and read-count assertions: one stencil query per requested
  destination, and no source chunk outside the stencil union is read.

**Gate:** work is O(destination cells + nonzeros), not O(destination cells ×
candidate chunks), and budget/spill tests pass.

### P4 — conforming DGGS support

- Prototype the smallest fast topology seam behind `DGGSpace`.
- Implement a hexagonal source fast path producing dual triangles.
- Exercise a variable-valence defect and a higher-valence polygon with its
  explicitly selected basis.
- Test continuity across primal edges, face seams, and chunk boundaries.

**Gate:** topology lookup is bounded local work, uses no per-query polygon
matching in production, and is invariant to chunking.

### P5 — CopDEM support

- Implement ordinary Q1 and transition P1/Q1 stencils using exact band
  arithmetic and `_facing`.
- Test all band ratios and both hemispheres on a small structurally equivalent
  CopDEM twin, including antimeridian and 1-degree storage seams.
- Test partial grids as explicit rims.
- Benchmark a real-sized row/tile hot loop.
- Stop at the polar decision gate and record the chosen policy before enabling
  polar output.

**Gate:** every nonpolar query is O(1), emits at most four entries, reads exactly
its owning chunks, and agrees across tile seams and band transitions.

### P6 — production opt-in

- Add an explicit CopDEM production configuration selecting
  `BarycentricPoint`; do not change the default silently.
- Compare the existing conservative result and the point result on known DEM
  samples, smooth synthetic surfaces, band transitions, partial-grid rims, and
  polar policy regions.
- Record end-to-end time, source reads, peak residency, and output differences.
- Document that point output is not conservative and that area-mean DEMs should
  continue to use `Conservative`.

**Gate:** the user can choose either semantic path explicitly, with no regression
to the current conservative production run.

## Standing acceptance laws

1. Every mapped row sums to one within a stated floating tolerance.
2. Raster Q1, triangle P1, and the selected polygon coordinates pass their
   reproduction laws.
3. A destination sample exactly at a source sample site reproduces that source.
4. No geometric or discovery false negative is permitted.
5. Direct, eager, lazy, and differently chunked execution produce equivalent
   operators and values.
6. Source chunks read for a point tile equal the sorted union owning its
   nonzero stencil entries.
7. Weight construction is independent of source field values and reads no
   source data.
8. Steady-state point lookup allocates nothing per destination.
9. Large-tile planning computes each point stencil once and keeps temporary
   memory bounded by the final nonzeros plus documented worker scratch.
10. Conservative tests and production behavior remain unchanged.

## Explicitly deferred

- automatic inference of point versus area-mean source semantics;
- a public universal dual-mesh object;
- vector/tensor basis transport;
- signed extrapolation, ghost nodes, or implicit rim clamping;
- ESMF patch and generic MLS;
- second-order conservative reconstruction; and
- the optional support query, `chunksupport`, unless P3-P5 profiles justify it.

## Current code seams

At the time of this plan:

- method, `outputsampling`, and `WeightCOO` contracts are in
  `lib/GlobalRegridding/src/methods.jl`;
- the old chart-only implementation is in
  `lib/GlobalRegridding/src/interpolation.jl`;
- chunk-pair rebuilding occurs in `lib/GlobalRegridding/src/plans.jl`;
- cap-based candidate discovery and lazy source selection are in
  `lib/GlobalRegridding/src/discovery.jl` and
  `lib/GlobalRegridding/src/lazy.jl`;
- the qualified `DGGSpace` bridge is in `src/regridding.jl`; and
- CopDEM's exact row transition topology is in
  `src/systems/CopernicusDEM/system.jl` around `_facing` and `one_ring`.

These locations are evidence, not a requirement to preserve the current file
layout. The architectural boundaries and acceptance laws above are the
requirements.
