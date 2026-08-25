# P1 — the point contract and its geometry kernels

- Date: 2026-08-25
- Plan: `regrid-notes/2026-08-23-barycentric-regridding-plan.md`, task P1
- Scope: `lib/GlobalRegridding/src/barycentric.jl`, the method in `methods.jl`,
  the `DGGSpace` sample-site hook in `src/regridding.jl`

## 1. The contract

`BarycentricPoint()` is a point method: `outputsampling` reports `Points()`, so
the shared builders put it on the point path without naming its type.
`Conservative`, `NearestCell`, `BilinearPoint` and every custom
`buildweights!` are untouched.

The hot call is one destination point at a time:

```julia
weightsat!(row::WeightRow, sampler, p) -> WeightStatus
```

`WeightRow` owns parallel `Vector{Int}` source local indices and
`Vector{Float64}` weights. It is cleared on entry, refilled, and left empty for
every status but `WeightsMapped`, so one row serves a whole sweep and belongs to
one task. `WeightStatus` is an `@enum` — one concrete return type, compared with
`===`, so the per-point path stays inferable and allocation-free — with
`WeightsMapped` against the three diagnostic values `WeightsOutside`,
`WeightsRim` and `WeightsDegenerate`. Execution reads `ismapped` and nothing
else. Nothing on this path is told about chunks: a row carries the source
space's local indices, `1:ncells(space)`, and the caller partitions them.

Preparation happens once per source space, never once per destination:

```julia
sampler(::BarycentricPoint, space) -> Sampler   # space, sample sites, prepared state
samplesites(space::RegridSpace) -> AbstractVector      # local index -> point
samplerstate(space::RegridSpace) -> state              # default `nothing`
```

`samplesites(space)[i] == cellcentroid(space, i)` for every space. The fallback
is a lazy vector that reads `cellcentroid` on demand and holds nothing; a
`DGGSpace` answers with the collection's centroid field, the same pure vector a
sweep asked for `Centroid()` reads. A `Sampler` is immutable and read
concurrently; per-point scratch is the calling task's row.

The default `weightsat!` locates a dual cell and weights its nodes:

```julia
dualcellat(sampler, p) -> DualCell    # fallback: a cell with no nodes
chartat(sampler, p) -> Union{NTuple{2,Float64},Nothing}   # fallback: the cell chart
```

A `DualCell` is the polygon of source sample sites containing the point:
ordered source local indices, their coordinates in one two-dimensional chart,
and an explicit `BasisKind` (`Linear`, `Bilinear`, `MeanValue`). The kind is the
cell's own statement, not a reading of its node count, and the cell promises no
fixed node count. It is read and never written, so a prepared cell is shared by
concurrent queries and costs a query nothing; a space that must build one per
point owns `weightsat!` whole and uses its task's own buffers. The generic
`dualcellat` constructs nothing and answers no nodes, so a source space that
implements neither maps nothing rather than inventing a stencil.

The three kernels take bare geometry — no space, no chunk, no sampler:

```julia
linearweights!(row, indices, nodes, p::NTuple{2,Float64}) -> WeightStatus
bilinearweights!(row, indices, nodes, p::NTuple{2,Float64}) -> WeightStatus
meanvalueweights!(row, indices, nodes, p::NTuple{2,Float64}) -> WeightStatus
```

`nodes` holds node coordinates in one chart in cyclic order, `p` is a point in
that chart, `indices` names the nodes' source cells in the same order. Every
kernel: clears the row on entry; emits nonnegative weights that sum to one;
drops zero weights, so a point on an edge emits two entries and a point at a
node one; reproduces every affine field, and `a + bx + cy + dxy` for Q1 on an
axis-aligned rectangle; answers `WeightsOutside` for a point outside the cell,
which is never clamped in or extrapolated from; and answers
`WeightsDegenerate` for repeated nodes, a collinear triangle, a folded,
zero-area or reflex quadrilateral, a non-convex polygon, and an inverse map that
does not converge. `Bilinear` inverts the isoparametric map by Newton iteration
from the centre of the parameter square, to a step below `1e-12` within 20
iterations and a residual below `1e-8` of the cell's size; coordinates within
`1e-12` of the square's edge are snapped to it, which is what makes a node
exact and an edge point two entries. Rounding margins are `1e-12` throughout,
relative to the cell's own size; they are roundoff handling, not a rim policy.

`buildweights!(coo, ::BarycentricPoint, …)` is the chunk-pair route: one
`weightsat!` per destination in the pair, keeping the entries the source chunk
owns and writing them as `WeightCOO` documents. It reads no source data and no
chunk boundary changes a stencil, so a chunked operator equals the eager one —
but every `(destination tile, source chunk)` pair re-runs every point query,
which is exactly the repeated work the fused tile-weight build removes.

`supportradius` is not defined for `BarycentricPoint`, so it is `0.0` and chunk
discovery prunes on cap overlap alone.

## 2. What it verifies

`lib/GlobalRegridding/test/test_barycentric.jl`, 219 assertions:

- each kernel, on hand-written node coordinates whose source indices are
  deliberately not their node positions: interior points, edge points and node
  points; constants, affine fields and the coordinate field reproduced; row
  sums; the entry count at nodes and edges; every degeneracy class rejected;
  points outside on every side unmapped and leaving no row;
- Q1's own field `a + bx + cy + dxy` on an axis-aligned rectangle, and on a
  general convex quadrilateral the forward map of the coordinates read back off
  the weights returning the query point;
- mean-value coordinates equal to barycentric ones on a triangle, and equal
  weights at the centre of a regular hexagon;
- the row cleared on entry and left empty when unmapped, and zero bytes
  allocated by each kernel and by `weightsat!` once warm;
- the sample sites equal to the space's centroids, a sampler that prepares once,
  and the generic source that maps nothing;
- a toy source space answering one hand-built four-node `Bilinear` dual cell
  over its own sample sites in two chunks: every stencil crosses the source-chunk
  seam, the eager whole-domain operator equals the one assembled from all chunk
  pairs, the operator reproduces the bilinear field, eager and chunked values
  agree, and a destination outside the dual cell takes no weights at all.

`test/systems/crosssystem/regrid.jl` pins the bridge: a `DGGSpace`'s sample
sites are its centroid field, and its entries are `cellcentroid`.

## 3. What P2 needs from it

- A raster answers `dualcellat` and `chartat`, or replaces `weightsat!`. The
  four-node cell it builds must say `Bilinear`; nothing infers that from the
  node count.
- `supportradius(::BarycentricPoint, ::RasterGrid)`. Without one the method's
  radius is `0.0`, and a chunked plan finds a source chunk only where its cap
  already meets the destination's. A Q1 stencil reaches out to the bracketing
  sample sites, up to a cell beyond a destination's own extent, so the raster
  bound is what keeps chunk discovery free of false negatives.
- The rim policy for the half-cell strip between the outermost sample sites and
  the source's boundary. `WeightsRim` exists for it and no kernel emits it: a
  rim is a statement about a source collection, not about a cell's geometry, so
  it can only come from a space's own `dualcellat` or `weightsat!`.
- Whatever a raster sampler must precompute — axis locators, periods — belongs
  in `samplerstate`, which the sampler already carries and reads concurrently.
