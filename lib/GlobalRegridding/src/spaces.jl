# Regridding spaces use dense cell positions and spherical-cap spatial trees.

"""
    RegridSpace

A source or destination cell collection.

Implementations provide [`celltree`](@ref), [`chunktree`](@ref),
[`ncells`](@ref), [`getcell`](@ref), [`cellindices`](@ref),
[`nchunks`](@ref), and [`manifold`](@ref). [`cellat`](@ref) and
[`cellcentroid`](@ref) are optional fast paths that some methods require;
[`hascellchart`](@ref) is the trait that gates the interpolating ones.

Spaces contain geometry and structure, not field data. Construction should be
cheap, with cell polygons generated on demand.
"""
abstract type RegridSpace end

"""
    celltree(space::RegridSpace)

Return a `SpatialTreeInterface` tree over cell positions `1:ncells(space)`.
Every node extent must be a `SphericalCap`. Define
`STI.node_extent_is_expensive` when extents are computed on demand.
"""
function celltree end

"""
    chunktree(space::RegridSpace)

Return a spatial tree over chunk numbers `1:nchunks(space)`. Each chunk extent
must cover every cell returned by [`cellindices`](@ref) for that chunk.
"""
function chunktree end

"""
    ncells(space::RegridSpace) -> Int

Return the stable cell count in `O(1)`. Cell positions are `1:ncells(space)`.
"""
function ncells end

"""
    getcell(space::RegridSpace, i::Int) -> GI.Polygon

Return cell `i` as a GeoInterface polygon with one explicitly closed ring of
unit-sphere `(x, y, z)` coordinates. The ring is counter-clockwise from outside
the sphere and its segments are great-circle arcs. Densify non-geodesic edges.
Throw `BoundsError` for an invalid position.
"""
function getcell end

"""
    nchunks(space::RegridSpace) -> Int

Return the number of chunks. A space without natural chunking returns one chunk
containing `1:ncells(space)`.
"""
function nchunks end

"""
    cellindices(space::RegridSpace, chunk::Int) -> AbstractVector{Int}

Return a chunk's ascending cell positions. Chunks must partition
`1:ncells(space)`. Return an `AbstractUnitRange` when positions are contiguous.
Weight builders address entries by local position within this result.
"""
function cellindices end

"""
    chunkat(space::RegridSpace, i::Integer) -> Int
    chunkat(space::RegridSpace, p::GO.UnitSphericalPoint) -> Union{Int,Nothing}

Return the chunk containing cell position `i` or point `p`. The fallback scans
all chunks; structured spaces should provide an `O(1)` or `O(log nchunks)`
method. The point form returns `nothing` outside the space's coverage.
"""
function chunkat end

function chunkat(space::RegridSpace, i::Integer)
    p = Int(i)
    1 <= p <= ncells(space) || throw(BoundsError(space, p))
    for c in 1:nchunks(space)
        p in cellindices(space, c) && return c
    end
    throw(ArgumentError(
        "cell position $p of $(typeof(space)) belongs to no chunk; chunks must " *
        "partition 1:ncells(space)"))
end

function chunkat(space::RegridSpace, p::US.UnitSphericalPoint)
    i = cellat(space, p)
    i === nothing && return nothing
    return chunkat(space, i)
end

"""
    manifold(space::RegridSpace) -> GeometryOpsCore.Manifold

Return the geometry manifold. Source and destination manifolds must match;
regridding does not reproject coordinates.
"""
function manifold end

"""
    cellat(space::RegridSpace, p::GO.UnitSphericalPoint) -> Union{Int,Nothing}

Return the position containing `p`, or `nothing` outside coverage. Point-based
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
