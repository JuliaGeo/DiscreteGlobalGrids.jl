# The space contract. A space is a cell collection addressed by dense positions
# `1:ncells(space)` and grouped into chunks `1:nchunks(space)`. Both trees are
# `GeometryOps.SpatialTreeInterface` trees whose node extents are
# `GO.UnitSpherical.SphericalCap`s; that is the whole predicate vocabulary the
# dual descent needs, and it is what makes the two sides of a regrid comparable
# without a shared coordinate system.

"""
    RegridSpace

A cell collection that can take part in a regrid, on either side.

Implementations provide [`celltree`](@ref), [`chunktree`](@ref),
[`ncells`](@ref), [`getcell`](@ref), [`cellindices`](@ref),
[`nchunks`](@ref), and [`manifold`](@ref). [`cellat`](@ref) and
[`cellcentroid`](@ref) are optional fast paths that some methods require;
[`hascellchart`](@ref) is the trait that gates the interpolating ones.

A space is geometry and structure only. It is not asked for data, never sees
the array being regridded, and is expected to be cheap to construct — cell
polygons are synthesized on demand, not stored.
"""
abstract type RegridSpace end

"""
    celltree(space::RegridSpace)

A spatial tree over `space`'s cells whose leaf indices are cell positions in
`1:ncells(space)`.

**Required.**

The result implements `GeometryOps.SpatialTreeInterface`, with
`GO.UnitSpherical.SphericalCap` node extents at every level. A tree that
derives extents rather than storing them should also define
`STI.node_extent_is_expensive` on its node type, so the dual descent caches
them.

Construction should be cheap relative to `ncells`; callers may build one per
chunk pair.
"""
function celltree end

"""
    chunktree(space::RegridSpace)

A spatial tree over `space`'s chunks whose leaf indices are chunk numbers in
`1:nchunks(space)`.

**Required.**

Same conventions as [`celltree`](@ref): `SpatialTreeInterface`, `SphericalCap`
node extents. This is the tree the lazy path descends against the destination's
to find which source chunks a destination chunk needs, so a chunk's extent must
cover every cell [`cellindices`](@ref) assigns to it.
"""
function chunktree end

"""
    ncells(space::RegridSpace) -> Int

The number of cells in `space`. Positions run over `1:ncells(space)`.

**Required.**

This is `ConservativeRegridding.Trees.ncells` — the same binding, extended here
— so a space is a `Trees` source without a wrapper.

Must be O(1) and must not change over the lifetime of the space.
"""
function ncells end

"""
    getcell(space::RegridSpace, i::Int) -> GI.Polygon

The boundary of the cell at position `i`, as a GeoInterface polygon of one
explicitly closed `LinearRing` with unit-sphere `(x, y, z)` coordinates.

**Required.**

This is `ConservativeRegridding.Trees.getcell`, extended here.

The ring is counter-clockwise seen from outside the sphere, so that area,
containment, and clipping all read the same winding. Edges are great-circle
arcs between consecutive vertices; a space whose cell edges are not geodesics
(a parallel of latitude, say) densifies them itself.

Cells are synthesized on demand. `i` outside `1:ncells(space)` throws a
`BoundsError`.
"""
function getcell end

"""
    nchunks(space::RegridSpace) -> Int

The number of chunks in `space`. Chunk numbers run over `1:nchunks(space)` and
are the leaf indices of [`chunktree`](@ref).

**Required.**

A space with no natural chunking answers `1`, with
`cellindices(space, 1) == 1:ncells(space)`.
"""
function nchunks end

"""
    cellindices(space::RegridSpace, chunk::Int) -> AbstractVector{Int}

The cell positions belonging to `chunk`, ascending.

**Required.**

Chunks partition `1:ncells(space)`: every position belongs to exactly one
chunk. Spaces whose cell order makes a chunk contiguous — a chunked 1-D cell
list, a DGGS subtree in canonical order — return an `AbstractUnitRange`, and
callers that can exploit contiguity should test for one rather than assume it.

Weight builders receive these vectors as `dst_inds` and `src_inds` and address
their `WeightCOO` entries by position *within* them, not by the cell positions
themselves; see [`build_weights!`](@ref).
"""
function cellindices end

"""
    manifold(space::RegridSpace) -> GeometryOpsCore.Manifold

The manifold `space`'s geometry lives on, which for every space here is a
sphere.

**Required.**

This is `GeometryOpsCore.manifold`, extended here. It is what areas, clipping,
and point-in-cell tests are computed on; the two sides of a regrid must agree,
and a mismatch is an error rather than a silent reprojection.
"""
function manifold end

"""
    cellat(space::RegridSpace, p::GO.UnitSphericalPoint) -> Union{Int,Nothing}

The position of the cell of `space` containing `p`, or `nothing` when `p` lies
outside the space's coverage.

**Optional**, and a fast path rather than a convenience: methods that locate
points ([`NearestCell`](@ref), the interpolating family) require it of their
source space and error without it. `nothing` is a real answer — a space need
not cover the sphere.

A shared-boundary point is assigned deterministically to one incident cell;
which one is the space's business, but it must never be a non-incident cell and
must never be `nothing` inside coverage.
"""
function cellat end

"""
    cellcentroid(space::RegridSpace, i::Int) -> GO.UnitSphericalPoint

The centroid of the cell at position `i`, strictly interior to the cell.

**Optional.** Methods that sample the source at a destination point
([`NearestCell`](@ref), [`BilinearPoint`](@ref)) require it of their
destination space.

Cheaper than deriving it from [`getcell`](@ref) is the point of having it:
spaces with a closed-form cell centre should not build a polygon to answer.
"""
function cellcentroid end

"""
    hascellchart(space::RegridSpace) -> Bool

Whether `space` carries a structured chart over its cells — a local coordinate
system in which neighbouring cell centres form a regular lattice that an
interpolation stencil can be written against.

Defaults to `false`. A raster space is `true`; an arbitrary cell collection is
`false` until it can name such a chart.

This gates the interpolating methods: [`BilinearPoint`](@ref) and its relatives
error when their **source** space answers `false`, because their stencil is
defined on the chart and not on the cell polygons.
"""
hascellchart(::RegridSpace) = false
