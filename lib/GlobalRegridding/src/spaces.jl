# Regridding spaces use dense cell indices and spherical-cap spatial trees.

"""
    RegridSpace

Base type for source and destination cell geometry. Every space provides cell
geometry, chunk ownership, and a manifold. Optional interfaces add restricted
trees, native chunk indexes, rectangular reads, point lookup, charts, and
labelled output.

A space holds geometry and structure; data arrays hold field values. Constructors
should remain cheap and generate cell polygons on demand.
"""
abstract type RegridSpace end

# Cell geometry

"""
    celltree(space::RegridSpace)

Return a `SpatialTreeInterface` tree over cell indices `1:ncells(space)`.
Every node extent must be a `SphericalCap`. Define
`STI.node_extent_is_expensive` when extents are computed on demand.
"""
function celltree end

"""
    subtree(space::RegridSpace, inds) -> tree

Return a spatial tree over `inds` whose leaves use the space's local indices.
The fallback packs GeometryOps Cartesian cell extents in an R-tree. Spaces with
a native restricted tree should specialize this function.
"""
function subtree end

"""
    ncells(space::RegridSpace) -> Int

Return the stable cell count in `O(1)`. Cell indices are `1:ncells(space)`.
"""
function ncells end

"""
    getcell(space::RegridSpace, i::Int) -> GI.Polygon

Return cell `i` as a GeoInterface polygon. The polygon contract requires:

  - one explicitly closed ring of unit-sphere `(x, y, z)` coordinates;
  - counter-clockwise order when viewed from outside the sphere;
  - great-circle segments, with non-geodesic edges densified.

Invalid indices throw `BoundsError`.
"""
function getcell end

"""
    expensivecellgeometry(space::RegridSpace) -> Bool

Return whether repeated [`getcell`](@ref) calls justify caching destination
polygons across overlapping source leaves. The default is `true` for spaces
that derive boundaries from cell identifiers. Lattice spaces should return
`false` when coordinate reads are cheaper than cache traffic.

This trait describes cell geometry. `STI.node_extent_is_expensive` separately
describes spatial-tree extents.
"""
expensivecellgeometry(::RegridSpace) = true

# Chunk ownership and spatial discovery

"""
    nchunks(space::RegridSpace) -> Int

Return the number of chunks. A space without natural chunking returns one chunk
containing `1:ncells(space)`.
"""
function nchunks end

"""
    ownedindices(space::RegridSpace, chunk::Int) -> AbstractVector{Int}

Return the ascending local indices owned by `chunk`. Chunks must partition
`1:ncells(space)`. Contiguous ownership should use `AbstractUnitRange`. Weight
builders address entries by position within this result.
"""
function ownedindices end

# Deprecation forwards calls only; extensions must implement `ownedindices`.

"""
    cellindices(space::RegridSpace, chunk::Int) -> AbstractVector{Int}

Deprecated. Use [`ownedindices`](@ref), which this forwards to, so existing
calls keep their old behaviour exactly.
"""
function cellindices end

@deprecate cellindices(space::RegridSpace, chunk::Int) ownedindices(space, chunk) false

"""
    chunkextent(space::RegridSpace, chunk::Integer) -> SphericalCap

Return a spherical cap covering every cell owned by `chunk`. The fallback
indexes [`chunkextents`](@ref). Spaces such as [`RasterGrid`](@ref) specialize
this method when one cap is cheaper than the complete vector.

A dependency relation retains the caps used at construction. Read those through
[`destinationextent`](@ref) or [`sourceextent`](@ref) to preserve its identity
and avoid recomputation.
"""
function chunkextent end

"""
    chunkextents(space::RegridSpace) -> Vector{SphericalCap}

Return chunk caps in chunk-number order. Every space must implement this method.

The returned values serve three roles: they define [`spacestamp`](@ref), provide
destination queries during relation construction, and populate the generic
[`chunkindex`](@ref). Native indexes may use a different internal extent
representation in [`candidatechunks!`](@ref).
"""
function chunkextents end

"""
    chunkindex(space::RegridSpace) -> index

Build the source-chunk query object consumed by [`candidatechunks!`](@ref). The
fallback packs [`chunkextents`](@ref) in a GeometryOps `FlexibleRTree`.
Structured spaces may return a native hierarchy and extent representation.
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

Return the chunk containing cell index `i` or point `p`. The fallback scans all
chunks; structured spaces should provide an `O(1)` or `O(log nchunks)` method.
The point form returns `nothing` outside coverage.
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

# Array storage

"""
    chunkranges(space::RegridSpace, chunk, spatialsize::NTuple{NS,Int})
        -> NTuple{NS,UnitRange{Int}}

Return rectangular spatial ranges for reading `chunk` in one storage operation.
Flattening the block must match the order of [`ownedindices`](@ref). Spaces with
non-rectangular ownership must specialize this function.
"""
function chunkranges end

# Manifold, point lookup, and cell charts

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

# Deprecation forwards calls only; extensions must implement `chartlocalindex`.

"""
    chartposition(space::RegridSpace, ix::Int, iy::Int) -> Int

Deprecated. Use [`chartlocalindex`](@ref), which this forwards to, so existing
calls keep their old behaviour exactly.
"""
function chartposition end

@deprecate chartposition(space::RegridSpace, ix::Int, iy::Int) chartlocalindex(space, ix, iy) false

# Output labelling and target resolution

"""
    destinationdims(space::RegridSpace, sampling) -> Tuple or nothing

Return the `DimensionalData` dimensions that label a result over this space,
fastest dimension first, with lookups carrying `sampling`. The default is
`nothing`: the result keeps one flat `Cell` axis over `1:ncells(space)`.
"""
destinationdims(::RegridSpace, ::DD.Lookups.Sampling) = nothing

"""
    dimsource(lookup) -> `from` target or nothing

Return the source target a lookup names, or `nothing`.

Packages extend this for lookups over explicit cells. Source inference resolves
the target through [`sourcespacefor`](@ref); lookups returning `nothing` remain
eligible for raster inference. An explicit `from` takes precedence over both
paths.
"""
dimsource(::Any) = nothing

"""
    sourceview(lookup, data, method) -> array or nothing
    sourceview(data, method) -> array

Return the source array presented to `method`.

  - The axis-level form returns a specialized view or `nothing`.
  - The array-level form returns the first specialized view or `data` unchanged.

Compressed axes use this hook to map stored values onto the cells named by
[`dimsource`](@ref). The returned view must name its space through its own
`dimsource` and order values like that space. Implementations may reject methods
that fail [`refinementinvariant`](@ref). [`checksource`](@ref) handles conflicts
with an explicit `from`.
"""
sourceview(::Any, ::Any, ::Any) = nothing

sourceview(data, method) = data

function sourceview(data::DD.AbstractDimArray, method)
    for d in DD.dims(data)
        view = sourceview(DD.lookup(d), data, method)
        view === nothing || return view
    end
    return data
end

"""
    checksource(from, data, space) -> nothing

Validate that an explicit `from` describes the layout of `data`.

Plans call this once after resolving `space`. Sources with compressed or
method-specific layouts should throw an `ArgumentError` here when the target
conflicts with the stored values.
"""
checksource(::Any, ::Any, ::RegridSpace) = nothing

"""
    _asspace(space, name) -> RegridSpace
    _asspace(space, name, src_space) -> RegridSpace

Resolve a `to` or `from` target into a [`RegridSpace`](@ref). Packages extend
the two-argument form for target spellings and the three-argument form for
destinations that depend on the resolved source. `name` identifies the keyword
in errors.
"""
function _asspace end

"""
    sourcespacefor(target, method) -> RegridSpace

Resolve a source target into the [`RegridSpace`](@ref) that `method` reads.

The default delegates to [`_asspace`](@ref)`(target, "from")`. Targets with
multiple presentations specialize this method and choose through
[`sourcesampling`](@ref). The returned space must match the cells and ordering
of [`sourceview`](@ref) for the same method.
"""
sourcespacefor(target, method) = _asspace(target, "from")
