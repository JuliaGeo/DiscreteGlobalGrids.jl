# Regridding spaces use dense cell indices and spherical-cap spatial trees.

"""
    RegridSpace

A source or destination cell collection.

The qualified extension contract is grouped below by responsibility. Every
space provides cell geometry, chunk ownership, and a manifold. Restricted cell
trees, native chunk indexes, rectangular storage reads, point lookup, charts,
and labelled output have generic fallbacks or are required only by the methods
that use them.

Spaces contain geometry and structure, not field data. Construction should be
cheap, with cell polygons generated on demand.
"""
abstract type RegridSpace end

# --------------------------------------------------------------------------
# Cell geometry
# --------------------------------------------------------------------------

"""
    celltree(space::RegridSpace)

Return a `SpatialTreeInterface` tree over cell indices `1:ncells(space)`.
Every node extent must be a `SphericalCap`. Define
`STI.node_extent_is_expensive` when extents are computed on demand.
"""
function celltree end

"""
    subtree(space::RegridSpace, inds) -> tree

Return a spatial tree over `inds`, with leaves addressed by the space's
local index. The fallback packs GeometryOps' Cartesian cell extents in an
R-tree.
Spaces with a cheaper restricted tree should specialize this function.
"""
function subtree end

"""
    ncells(space::RegridSpace) -> Int

Return the stable cell count in `O(1)`. Cell indices are `1:ncells(space)`.
"""
function ncells end

"""
    getcell(space::RegridSpace, i::Int) -> GI.Polygon

Return cell `i` as a GeoInterface polygon with one explicitly closed ring of
unit-sphere `(x, y, z)` coordinates. The ring is counter-clockwise from outside
the sphere and its segments are great-circle arcs. Densify non-geodesic edges.
Throw `BoundsError` for an invalid index.
"""
function getcell end

"""
    expensivecellgeometry(space::RegridSpace) -> Bool

Whether one [`getcell`](@ref) costs enough that an area method should keep the
polygons of a destination tile rather than synthesize each one again for every
source leaf that overlaps it.

  - `true` by default: a space that derives cell boundaries from an index pays
    for every call.
  - A space whose cells are a few coordinates read off a lattice should answer
    `false`; keeping them there costs more memory traffic than it saves.
  - A property of the space, not of the tree an index set produces. The
    spatial-tree trait `STI.node_extent_is_expensive` describes extents, and the
    packed R-tree fallback answers `false` for it whatever space it wraps.
"""
expensivecellgeometry(::RegridSpace) = true

# --------------------------------------------------------------------------
# Cell adjacency
# --------------------------------------------------------------------------

"""
    cellneighbors(space::RegridSpace, i::Int) -> AbstractVector{Int}

The local indices of the cells adjacent to cell `i`, in no particular order,
without `i` itself.

Adjacent means sharing an edge or a vertex. A space with lattice or grid
topology answers from that, where edge neighbours are enough; the generic
fallback answers from geometry, every cell sharing a vertex with `i`'s,
diagonals included.

Only cells of the collection are returned, so a partial space has fewer
neighbours at its rim. The result may alias the space's own storage and is valid
until the next call.

Gradient-recovering methods use this, and their least-squares fit depends on
neither the order nor the count, so 4, 6 or 8 are all fine.
"""
function cellneighbors end

"""
    celldiameter(space::RegridSpace) -> Float64

An upper bound, in radians, on the angular diameter of every cell of `space` —
the largest angular distance between two points of one cell.

A bound, not a measurement: overestimating costs discovery work and nothing
else. A method whose stencil reaches one cell past its source — a cell and its
[`cellneighbors`](@ref) — declares this as its [`supportradius`](@ref).

The generic fallback answers twice the widest leaf cap of [`celltree`](@ref),
`O(ncells(space))`. A structured space should answer from its resolution.
"""
function celldiameter end

# --------------------------------------------------------------------------
# Chunk ownership and spatial discovery
# --------------------------------------------------------------------------

"""
    nchunks(space::RegridSpace) -> Int

Return the number of chunks. A space without natural chunking returns one chunk
containing `1:ncells(space)`.
"""
function nchunks end

"""
    ownedindices(space::RegridSpace, chunk::Int) -> AbstractVector{Int}

The space's local indices of the cells `chunk` owns — the cells it produces
results for — ascending. Chunks must partition `1:ncells(space)`. Return an
`AbstractUnitRange` when those indices are contiguous. Weight builders address
entries by chunk-local index within this result.
"""
function ownedindices end

# `cellindices` is the old name of `ownedindices` and forwards to it, so a call
# of the old name answers the same with a deprecation warning. In this package
# "cell index" is the local index a space numbers its cells by, and the name
# collided with the typed cell id it means elsewhere. Only callers are carried:
# a space that defines the old name supplies no chunk ownership, and the
# generic dispatches on the new one.

"""
    cellindices(space::RegridSpace, chunk::Int) -> AbstractVector{Int}

Deprecated. Use [`ownedindices`](@ref), which this forwards to, so existing
calls keep their old behaviour exactly.
"""
function cellindices end

@deprecate cellindices(space::RegridSpace, chunk::Int) ownedindices(space, chunk) false

"""
    chunkextent(space::RegridSpace, chunk::Integer) -> SphericalCap

Return the spherical cap covering every cell owned by `chunk`. The fallback
indexes [`chunkextents`](@ref); a space specializes it when one extent is
cheaper to obtain than the complete vector, as [`RasterGrid`](@ref) does.

For the cap of a chunk a *relation* was built over, prefer
[`destinationextent`](@ref) or `sourceextent` on the relation: it already holds
the caps it was built from, so reading them back off the space recomputes work
and risks answering with a cap the relation never saw.
"""
function chunkextent end

"""
    chunkextents(space::RegridSpace) -> Vector{SphericalCap}

Return the chunk extents in chunk-number order. Required of every space: there
is no fallback.

These are the caps as *values*, not a query. They stamp a relation's identity
([`spacestamp`](@ref)), they are the destination caps a relation is built by
querying with, and the generic [`chunkindex`](@ref) packs them. A chunk query
goes to [`candidatechunks!`](@ref) on the space's own index, so a native space
may report caps here that its index does not itself test.
"""
function chunkextents end

"""
    chunkindex(space::RegridSpace) -> index

Build the source-chunk query object consumed by [`candidatechunks!`](@ref).
The fallback packs [`chunkextents`](@ref) in a GeometryOps `FlexibleRTree`.
Structured spaces may return any native hierarchy; indexes need not share a
type or expose one common node-extent representation.
"""
function chunkindex end

"""
    candidatechunks!(out::Vector{Int}, index, dstcap; radius = 0.0) -> out

Replace `out` with ascending unique chunk numbers from `index` that may lie
within `radius` radians of `dstcap`. False positives are allowed; false
negatives are not. A native index specializes this query operation directly.
"""
function candidatechunks! end

"""
    chunkat(space::RegridSpace, i::Integer) -> Int
    chunkat(space::RegridSpace, p::GO.UnitSphericalPoint) -> Union{Int,Nothing}

Return the chunk containing cell index `i` or point `p`. The fallback scans
all chunks; structured spaces should provide an `O(1)` or `O(log nchunks)`
method. The point form returns `nothing` outside the space's coverage.
"""
function chunkat end

function chunkat(space::RegridSpace, i::Integer)
    i = Int(i)
    1 <= i <= ncells(space) || throw(BoundsError(space, i))
    for c in 1:nchunks(space)
        i in ownedindices(space, c) && return c
    end
    throw(ArgumentError(
        "cell index $i of $(typeof(space)) belongs to no chunk; chunks must " *
        "partition 1:ncells(space)"))
end

function chunkat(space::RegridSpace, p::US.UnitSphericalPoint)
    i = cellat(space, p)
    i === nothing && return nothing
    return chunkat(space, i)
end

# --------------------------------------------------------------------------
# Array storage
# --------------------------------------------------------------------------

"""
    chunkranges(space::RegridSpace, chunk, spatialsize::NTuple{NS,Int})
        -> NTuple{NS,UnitRange{Int}}

Return the rectangular array ranges that storage can read for `chunk` in one
operation, in spatial-dimension order. Flattening that block must enumerate
[`ownedindices`](@ref) in the same order, but the two contracts are distinct:
`ownedindices` describes cell ownership and need not be a storage rectangle.
Non-rectangular spaces must specialize this function.
"""
function chunkranges end

# --------------------------------------------------------------------------
# Manifold, point lookup, and cell charts
# --------------------------------------------------------------------------

"""
    manifold(space::RegridSpace) -> GeometryOpsCore.Manifold

Return the geometry manifold. Source and destination manifolds must match;
regridding does not reproject coordinates.
"""
function manifold end

"""
    cellat(space::RegridSpace, p::GO.UnitSphericalPoint) -> Union{Int,Nothing}

Return the index containing `p`, or `nothing` outside coverage. Point-based
methods require this optional interface. Assign boundary points consistently to
an incident cell.
"""
function cellat end

"""
    cellcentroid(space::RegridSpace, i::Int) -> GO.UnitSphericalPoint

Return an interior centroid for cell `i`. Point-sampling methods require this
optional interface on the destination. Prefer a direct calculation over
building the polygon.
"""
function cellcentroid end

"""
    hascellchart(space::RegridSpace) -> Bool

Return whether cells have a structured chart suitable for interpolation.
Defaults to `false`; chart-based methods require `true` from the source space.
"""
hascellchart(::RegridSpace) = false

"""
    chartaxes(space::RegridSpace) -> (xs, ys)

Return strictly monotonic cell-centre coordinates for each separable lattice
axis. Required when [`hascellchart`](@ref) is `true`.
"""
function chartaxes end

"""
    chartcoords(space::RegridSpace, p) -> Union{Tuple{Real,Real},Nothing}

Convert `p` to native chart coordinates, or return `nothing` outside the chart.
Coordinates must use the same branch as [`chartaxes`](@ref), except on periodic
axes.
"""
function chartcoords end

"""
    chartlocalindex(space::RegridSpace, ix::Int, iy::Int) -> Int

Return the space's local index for the cell at lattice index `(ix, iy)`.
Required when [`hascellchart`](@ref) is `true`.
"""
function chartlocalindex end

"""
    chartperiod(space::RegridSpace) -> (px, py)

Return each axis period in native coordinates, or `nothing` for no wrap.
Defaults to `(nothing, nothing)`.
"""
function chartperiod end

"""
    chartspacing(space::RegridSpace) -> (Δx, Δy)

Return upper bounds, in radians, on adjacent-centre distance along each axis.
Required when [`hascellchart`](@ref) is `true`.
"""
function chartspacing end

# `chartposition` is the old name of `chartlocalindex` and forwards to it, so a
# call of the old name answers the same with a deprecation warning. Only
# callers are carried: a space that defines the old name supplies no chart
# hook, and `_chart_required` names the new one.

"""
    chartposition(space::RegridSpace, ix::Int, iy::Int) -> Int

Deprecated. Use [`chartlocalindex`](@ref), which this forwards to, so existing
calls keep their old behaviour exactly.
"""
function chartposition end

@deprecate chartposition(space::RegridSpace, ix::Int, iy::Int) chartlocalindex(space, ix, iy) false

# --------------------------------------------------------------------------
# Output labelling and target resolution
# --------------------------------------------------------------------------

"""
    destinationdims(space::RegridSpace, sampling) -> Tuple or nothing

Return the `DimensionalData` dimensions that label a result over this space,
fastest dimension first, with lookups carrying `sampling`. The default is
`nothing`: the result keeps one flat `Cell` axis over `1:ncells(space)`.
"""
destinationdims(::RegridSpace, ::DD.Lookups.Sampling) = nothing

"""
    dimsource(lookup) -> `from` target or nothing

Return the source a lookup already names, or `nothing`. A lookup that carries
its own cells is not a raster axis, so a package that supplies one extends this
and a regrid given no `from` names it instead of asking for `xdim`.
"""
dimsource(::Any) = nothing

"""
    _asspace(space, name) -> RegridSpace
    _asspace(space, name, src_space) -> RegridSpace

Resolve a `to` or `from` argument into a [`RegridSpace`](@ref). Packages that
supply spaces extend the two-argument form for their own target spellings, and
the three-argument form when the destination depends on the resolved source
space. `name` names the keyword in error messages.
"""
function _asspace end
