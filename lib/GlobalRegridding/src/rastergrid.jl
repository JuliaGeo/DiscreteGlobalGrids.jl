# `RasterGrid`, the dimensional-raster `RegridSpace`. Owned by task T2.
#
# The whole space is two cell-edge vectors plus a chart mapping native raster
# coordinates onto the unit sphere. Polygons, centroids and caps are synthesized
# from those on demand, so construction costs nothing that scales with the cell
# count and no array element is ever read.

# ===========================================================================
# The chart
# ===========================================================================

"""
    LonLatToSphere()

The chart of a plain geographic raster: longitude and latitude in **degrees**
to a `UnitSphericalPoint`.

A projected raster supplies its own callable instead; see [`RasterGrid`](@ref).
"""
struct LonLatToSphere end

function (::LonLatToSphere)(lon::Real, lat::Real)
    coslat = cosd(Float64(lat))
    return USPoint(coslat * cosd(Float64(lon)), coslat * sind(Float64(lon)),
        sind(Float64(lat)))
end

"""
    SphereToLonLat()

The inverse of [`LonLatToSphere`](@ref): a unit-sphere point to longitude in
`(-180, 180]` and latitude in `[-90, 90]`, both in degrees.
"""
struct SphereToLonLat end

(::SphereToLonLat)(p) = (atand(p[2], p[1]), asind(clamp(p[3], -1.0, 1.0)))

"""
    chartinverse(transform) -> callable or nothing

The inverse of a [`RasterGrid`](@ref) chart, mapping a unit-sphere point back to
the pair of native coordinates. `nothing` — the default — means the chart cannot
be inverted, and [`cellat`](@ref) on a space carrying it throws.
"""
chartinverse(::Any) = nothing
chartinverse(::LonLatToSphere) = SphereToLonLat()

"""
    chartlimits(transform) -> (xperiod, ybounds)

The limits a chart puts on its native coordinates: the period of the first
coordinate, or `nothing` when it does not wrap, and the closed bounds of the
second, or `nothing` when it is unbounded.

[`RasterGrid`](@ref) reduces a [`cellat`](@ref) query by the period before
searching, and clamps a midpointed edge vector into the bounds. Defaults to
`(nothing, nothing)`, so a projected chart is unconstrained unless it says
otherwise.
"""
chartlimits(::Any) = (nothing, nothing)
chartlimits(::LonLatToSphere) = (360.0, (-90.0, 90.0))

"""
    chartarcs(transform, dx, dy) -> (Δx, Δy)

An upper bound, in radians of arc, on the sphere distance a chart step of `dx`
in the first native coordinate or `dy` in the second can span.

The one place a chart's native units are converted to angle. `chartspacing` of a
[`RasterGrid`](@ref) is this applied to its largest edge steps, and its only
consumer is [`support_radius`](@ref), so an over-estimate costs chunk-discovery
work and an under-estimate silently truncates stencils.
"""
chartarcs(t::Any, dx, dy) = throw(ArgumentError(
    "a RasterGrid on a $(typeof(t)) chart cannot bound its cell spacing in " *
    "radians of arc; define GlobalRegridding.chartarcs(::$(typeof(t)), dx, dy) " *
    "to use BilinearPoint with it."))

# A step of Δλ along a parallel spans Δλ·cos φ ≤ Δλ radians, and a step of Δφ
# along a meridian spans exactly Δφ.
chartarcs(::LonLatToSphere, dx, dy) = (deg2rad(Float64(dx)), deg2rad(Float64(dy)))

"""
    LonLatEdgeTables(xedges, yedges)

`cos` and `sin` of every cell edge of a [`LonLatToSphere`](@ref) raster, plus the
widest longitude step.

The corner at edge coordinates `(xedges[i], yedges[j])` is
`(cy[j]·cx[i], cy[j]·sx[i], sy[j])` — the same three products the chart itself
performs on the same arguments, so a table lookup returns the chart's own
`Float64`s. `maxdx` bounds how far a cell's east–west edge bows poleward of its
parallel, which is what lets a wide box be bounded by something short of the
whole sphere.
"""
struct LonLatEdgeTables
    cx::Vector{Float64}
    sx::Vector{Float64}
    cy::Vector{Float64}
    sy::Vector{Float64}
    maxdx::Float64
end

LonLatEdgeTables(xe::Vector{Float64}, ye::Vector{Float64}) =
    LonLatEdgeTables(cosd.(xe), sind.(xe), cosd.(ye), sind.(ye), _maxstep(xe))

@inline _tablecorner(t::LonLatEdgeTables, i::Int, j::Int) =
    @inbounds USPoint(t.cy[j] * t.cx[i], t.cy[j] * t.sx[i], t.sy[j])

"""
    chartedgetables(transform, xedges, yedges) -> tables or nothing

Per-edge trigonometric tables for a chart whose corners factor through the two
coordinates separately, or `nothing` for a chart that has to be evaluated.

Every corner, cap and polygon a [`RasterGrid`](@ref) synthesizes stands on a pair
of cell edges, so a chart that can be tabulated once per edge — `O(nx + ny)` —
never has to be evaluated per cell. `nothing` keeps the general path, in which
`transform` is called for every point.
"""
chartedgetables(::Any, xedges, yedges) = nothing
chartedgetables(::LonLatToSphere, xedges::Vector{Float64}, yedges::Vector{Float64}) =
    LonLatEdgeTables(xedges, yedges)

# ===========================================================================
# The space
# ===========================================================================

"""
    RasterGrid(A::DimensionalData.AbstractDimArray; kwargs...)
    RasterGrid(dims::Tuple; chunks = nothing, kwargs...)

A [`RegridSpace`](@ref) over the cells of a dimensional raster.

Only geometry is taken from `A`: the two spatial dimensions' lookups become
cell-**edge** vectors, the parent array's chunking becomes the chunk structure,
and no element is ever read. Construction is `O(nx + ny + nchunks)` — cell
polygons, centroids and spherical caps are synthesized from the edges when they
are asked for, never stored.

The spatial dimensions are found by `DimensionalData`'s `XDim`/`YDim` traits and
then by name (`x`, `lon`, `long`, `longitude`; `y`, `lat`, `latitude`); pass
`xdim`/`ydim` to name them explicitly. Every other dimension is ignored.

Edges come from interval bounds where the lookup has them, and from midpointing
the cell centres where it does not (the outer two edges extrapolate the first
and last half-widths). The edge vectors are kept in **lookup order**, so a
reverse-ordered lookup gives a descending edge vector, and every consumer here
derives the sense of the axis from that rather than assuming it.

# Cell order

Cell positions enumerate the two spatial dimensions in the **array's own
dimension order, fastest first**. With the usual `(X, Y)` order, position `i` of
cell `(ix, iy)` is `ix + (iy - 1) * nx`; with `(Y, X)` it is
`iy + (ix - 1) * ny`. Flattening a spatial slice of the raster with `vec` is
therefore exactly position order, with no permutation. [`cellsubscript`](@ref)
and [`cellposition`](@ref) are the mapping. Chunk numbers follow the same
convention over the chunk lattice.

# Keywords

  - `xdim`, `ydim`: the spatial dimensions, when the traits and names cannot
    find them.
  - `transform`: the chart, a callable `(x, y) -> UnitSphericalPoint`. Defaults
    to [`LonLatToSphere`](@ref); a projected raster passes its own, and the
    orientation of the cell rings adapts to it.
  - `inverse`: the chart's inverse, `p -> (x, y)`. Defaults to
    [`chartinverse`](@ref) of `transform`; `nothing` disables [`cellat`](@ref).
  - `xperiod`, `ybounds`: the native limits, defaulting to
    [`chartlimits`](@ref) of `transform`.
  - `chunks`: for the dimension-tuple form, `(xchunks, ychunks)` as vectors of
    index ranges; `nothing` is one whole-domain chunk.

# Geometry

Cells are geodesic quadrilaterals through their four edge-vector corners, so
neighbouring cells share an edge exactly and the polygons of a global raster
tile the sphere. Rings are counter-clockwise seen from outside whatever the
lookup order and whatever the handedness of the chart.
"""
struct RasterGrid{XD,YD,F,G,T} <: RegridSpace
    xdim::XD
    ydim::YD
    "cell edges along the X dimension, in lookup order, strictly monotone"
    xedges::Vector{Float64}
    "cell edges along the Y dimension, in lookup order, strictly monotone"
    yedges::Vector{Float64}
    "the X dimension's chunk index ranges, ascending and partitioning `1:nx`"
    xchunks::Vector{UnitRange{Int}}
    "the Y dimension's chunk index ranges, ascending and partitioning `1:ny`"
    ychunks::Vector{UnitRange{Int}}
    "native `(x, y)` to the unit sphere"
    transform::F
    "the unit sphere back to native `(x, y)`, or `nothing`"
    inverse::G
    "period of the native X coordinate, or `nothing`"
    xperiod::Union{Nothing,Float64}
    "whether the X dimension varies fastest in cell positions and chunk numbers"
    xfast::Bool
    "whether `(xlo, ylo), (xhi, ylo), (xhi, yhi), (xlo, yhi)` is already counter-clockwise from outside"
    ccw::Bool
    "the chart's per-edge tables ([`chartedgetables`](@ref)), or `nothing`"
    tables::T
end

const _XNAMES = (:x, :lon, :long, :longitude)
const _YNAMES = (:y, :lat, :latitude)

# Corner samples for a cell cap; a lon–lat cell box is farthest from any centre
# at a corner, so four points are the exact bound there.
const _CELL_CAP_SAMPLES = 0
# Extra samples per edge for a multi-cell box, headroom for a chart that is not
# a lon–lat one.
const _BOX_CAP_SAMPLES = 3
# Cells per leaf of the recursive cell tree — small enough that a leaf is
# genuinely local, large enough that the descent overhead is amortized.
const _CELL_TREE_LEAF = 16

const _WHOLE_SPHERE = SphericalCap(USPoint(0.0, 0.0, 1.0), nextfloat(Float64(pi)))

function RasterGrid(A::DD.AbstractDimArray;
    xdim = nothing, ydim = nothing, chunks = nothing, kwargs...)
    ds = DD.dims(A)
    xd = _resolvedim(ds, xdim, DD.XDim, _XNAMES, "X")
    yd = _resolvedim(ds, ydim, DD.YDim, _YNAMES, "Y")
    xnum, ynum = _dimnums(ds, xd, yd)
    ch = chunks === nothing ? _spatialchunks(parent(A), xnum, ynum) : chunks
    return _rastergrid(xd, yd, ch, xnum < ynum; kwargs...)
end

function RasterGrid(ds::Tuple{Vararg{DD.Dimension}};
    xdim = nothing, ydim = nothing, chunks = nothing, kwargs...)
    xd = _resolvedim(ds, xdim, DD.XDim, _XNAMES, "X")
    yd = _resolvedim(ds, ydim, DD.YDim, _YNAMES, "Y")
    xnum, ynum = _dimnums(ds, xd, yd)
    return _rastergrid(xd, yd, chunks, xnum < ynum; kwargs...)
end

RasterGrid(ds::DD.Dimension...; kwargs...) = RasterGrid(ds; kwargs...)

function _rastergrid(xd, yd, chunks, xfast::Bool;
    transform = LonLatToSphere(),
    inverse = chartinverse(transform),
    xperiod = chartlimits(transform)[1],
    ybounds = chartlimits(transform)[2])
    xe = _edges(xd, nothing)
    ye = _edges(yd, ybounds)
    nx, ny = length(xe) - 1, length(ye) - 1
    xc, yc = chunks === nothing ? ([1:nx], [1:ny]) :
             (_chunkranges(chunks[1], nx, "X"), _chunkranges(chunks[2], ny, "Y"))
    return RasterGrid(xd, yd, xe, ye, xc, yc, transform, inverse,
        xperiod === nothing ? nothing : Float64(xperiod), xfast,
        _chartorientation(transform, xe, ye), chartedgetables(transform, xe, ye))
end

function _dimnums(ds, xd, yd)
    xnum = DD.dimnum(ds, DD.basetypeof(xd))
    ynum = DD.dimnum(ds, DD.basetypeof(yd))
    xnum == ynum && throw(ArgumentError(
        "RasterGrid: the X and Y dimensions resolved to the same dimension"))
    return (xnum, ynum)
end

function Base.show(io::IO, space::RasterGrid)
    print(io, "RasterGrid(", DD.name(space.xdim), "×", DD.name(space.ydim), " ",
        _nx(space), "×", _ny(space), ", chunks=", length(space.xchunks), "×",
        length(space.ychunks), ")")
end

_nx(space::RasterGrid) = length(space.xedges) - 1
_ny(space::RasterGrid) = length(space.yedges) - 1

"""
    rastersize(space::RasterGrid) -> (n1, n2)

The shape of a spatial slice of the raster, in the array's own dimension order.
`vec` of such a slice is indexed by cell position.
"""
rastersize(space::RasterGrid) =
    space.xfast ? (_nx(space), _ny(space)) : (_ny(space), _nx(space))

"""
    DimensionalData.dims(space::RasterGrid)

The space's two spatial dimensions, in the array's own dimension order.
"""
DD.dims(space::RasterGrid) =
    space.xfast ? (space.xdim, space.ydim) : (space.ydim, space.xdim)

# --- Dimension and edge extraction -----------------------------------------

function _resolvedim(ds, given, D::Type, names, what)
    if given !== nothing
        d = DD.dims(ds, given)
        d === nothing && throw(ArgumentError(
            "RasterGrid: no dimension matching $given among $(map(DD.name, ds))"))
        return d
    end
    d = DD.dims(ds, D)
    d === nothing || return d
    for dim in ds
        Symbol(lowercase(String(DD.name(dim)))) in names && return dim
    end
    throw(ArgumentError("""
    RasterGrid could not find the $what dimension among $(map(DD.name, ds)).
    Expected a `DimensionalData.$(what)Dim` or a dimension named one of $names;
    pass `$(lowercase(what))dim = ...` to name it explicitly.
    """))
end

# Cell edges in lookup order. Interval lookups already carry them; point
# lookups are midpointed, with the outer two edges extrapolating the first and
# last half-widths.
function _edges(dim::DD.Dimension, ybounds)
    lk = DD.lookup(dim)
    lk isa DD.Lookups.AbstractSampled || throw(ArgumentError(
        "RasterGrid needs a `Sampled` lookup on $(DD.name(dim)); got $(typeof(lk))"))
    n = length(lk)
    n >= 1 || throw(ArgumentError("RasterGrid: dimension $(DD.name(dim)) is empty"))
    e = DD.sampling(lk) isa DD.Lookups.Intervals ?
        _interval_edges(lk, dim) : _point_edges(lk, dim)
    ybounds === nothing || map!(v -> clamp(v, ybounds[1], ybounds[2]), e, e)
    _checkmonotone(e, dim)
    return e
end

function _interval_edges(lk, dim)
    bnds = DD.intervalbounds(lk)
    n = length(bnds)
    rev = DD.order(lk) isa DD.Lookups.ReverseOrdered
    e = Vector{Float64}(undef, n + 1)
    lo(k) = Float64(min(bnds[k][1], bnds[k][2]))
    hi(k) = Float64(max(bnds[k][1], bnds[k][2]))
    e[1] = rev ? hi(1) : lo(1)
    for k in 1:n
        e[k+1] = rev ? lo(k) : hi(k)
    end
    # Non-abutting intervals are a gappy axis, which no single edge vector can
    # describe; catching it here beats silently widening every second cell.
    for k in 1:(n-1)
        width = hi(k) - lo(k)
        shared = rev ? lo(k) - hi(k + 1) : lo(k + 1) - hi(k)
        abs(shared) <= 1e-8 * max(width, 1.0) || throw(ArgumentError(
            "RasterGrid: the intervals of $(DD.name(dim)) do not abut at cell $k"))
    end
    return e
end

function _point_edges(lk, dim)
    v = Float64.(collect(DD.val(lk)))
    n = length(v)
    e = Vector{Float64}(undef, n + 1)
    if n == 1
        span = DD.span(lk)
        span isa DD.Lookups.Regular || throw(ArgumentError(
            "RasterGrid cannot infer a cell width for the single point of $(DD.name(dim));" *
            " give the lookup `Intervals` sampling or a `Regular` span"))
        half = Float64(DD.val(span)) / 2
        return [v[1] - half, v[1] + half]
    end
    for k in 2:n
        e[k] = (v[k-1] + v[k]) / 2
    end
    e[1] = 2 * v[1] - e[2]
    e[n+1] = 2 * v[n] - e[n]
    return e
end

function _checkmonotone(e::Vector{Float64}, dim)
    up = e[2] > e[1]
    ok = all(k -> up ? e[k+1] > e[k] : e[k+1] < e[k], 1:(length(e)-1))
    ok || throw(ArgumentError(
        "RasterGrid: the cell edges of $(DD.name(dim)) are not strictly monotone"))
    return nothing
end

# --- Chunk extraction ------------------------------------------------------

function _spatialchunks(data, xnum::Int, ynum::Int)
    sz = size(data)
    if DiskArrays.haschunks(data) isa DiskArrays.Chunked
        cs = DiskArrays.eachchunk(data)
        if cs isa DiskArrays.GridChunks
            return (collect(cs.chunks[xnum]), collect(cs.chunks[ynum]))
        end
    end
    return ([1:sz[xnum]], [1:sz[ynum]])
end

function _chunkranges(ranges, n::Int, what)
    out = UnitRange{Int}[UnitRange{Int}(first(r), last(r)) for r in ranges]
    expected = 1
    for r in out
        first(r) == expected || throw(ArgumentError(
            "RasterGrid: the $what chunks do not partition 1:$n in ascending order"))
        expected = last(r) + 1
    end
    expected == n + 1 || throw(ArgumentError(
        "RasterGrid: the $what chunks cover 1:$(expected - 1), not 1:$n"))
    return out
end

# ===========================================================================
# The contract
# ===========================================================================

ncells(space::RasterGrid) = _nx(space) * _ny(space)

manifold(::RasterGrid) = GOCore.Spherical(; radius = 1.0)

# The chart an interpolating method writes its stencil against is the lattice of
# cell centres, and placing a destination centroid on it needs the inverse
# transform. A raster built without one still has cells, just no chart.
hascellchart(space::RasterGrid) = space.inverse !== nothing

"""
    cellsubscript(space::RasterGrid, i::Int) -> (ix, iy)
    cellposition(space::RasterGrid, ix::Integer, iy::Integer) -> Int

The lattice coordinates of a cell position, and their inverse — the space's
chart. `ix` indexes the X dimension and `iy` the Y dimension whichever way round
they sit in the array; the flattening follows the array's own dimension order,
fastest first.
"""
function cellsubscript(space::RasterGrid, i::Int)
    1 <= i <= ncells(space) || throw(BoundsError(space, i))
    if space.xfast
        iy, ix = fldmod1(i, _nx(space))
    else
        ix, iy = fldmod1(i, _ny(space))
    end
    return (ix, iy)
end

function cellposition(space::RasterGrid, ix::Integer, iy::Integer)
    1 <= ix <= _nx(space) && 1 <= iy <= _ny(space) ||
        throw(BoundsError(space, (ix, iy)))
    return space.xfast ? Int(ix) + (Int(iy) - 1) * _nx(space) :
           Int(iy) + (Int(ix) - 1) * _ny(space)
end

"""
    cellbox(space::RasterGrid, ix::Integer, iy::Integer) -> (xlo, xhi, ylo, yhi)

The cell's native bounding box, always ascending — the edge vectors' own order
is discarded here, which is what keeps ring winding independent of it.
"""
function cellbox(space::RasterGrid, ix::Integer, iy::Integer)
    xlo, xhi = minmax(space.xedges[ix], space.xedges[ix+1])
    ylo, yhi = minmax(space.yedges[iy], space.yedges[iy+1])
    return (xlo, xhi, ylo, yhi)
end

function getcell(space::RasterGrid, i::Int)
    ix, iy = cellsubscript(space, i)
    c = _cellcorners(space, ix, iy)
    ring = _cellring(space.ccw ? c : (c[4], c[3], c[2], c[1]))
    return GI.Polygon([GI.LinearRing(ring)])
end

"""
    _cellcorners(space::RasterGrid, ix, iy) -> NTuple{4,UnitSphericalPoint}

The cell's four corners in ascending-native order —
`(xlo, ylo), (xhi, ylo), (xhi, yhi), (xlo, yhi)` — from the chart's edge tables
where it has them and from the chart itself otherwise.

The order discards the edge vectors' own sense, which is what keeps ring winding
independent of lookup order; `space.ccw` says whether it needs reversing.
"""
@inline function _cellcorners(space::RasterGrid, ix::Integer, iy::Integer)
    t = space.tables
    if t === nothing
        xlo, xhi, ylo, yhi = cellbox(space, ix, iy)
        f = space.transform
        return (f(xlo, ylo), f(xhi, ylo), f(xhi, yhi), f(xlo, yhi))
    end
    ilo, ihi = _edgeorder(space.xedges, ix)
    jlo, jhi = _edgeorder(space.yedges, iy)
    return (_tablecorner(t, ilo, jlo), _tablecorner(t, ihi, jlo),
        _tablecorner(t, ihi, jhi), _tablecorner(t, ilo, jhi))
end

# The edge-vector positions of a cell's lower and upper native bound: the index
# form of the `minmax` in `cellbox`, so the two agree on a descending axis.
@inline function _edgeorder(e::Vector{Float64}, k::Integer)
    i = Int(k)
    return @inbounds e[i] <= e[i+1] ? (i, i + 1) : (i + 1, i)
end

# A cell against a pole has two coincident corners; dropping the repeat keeps
# every edge of the ring non-degenerate.
#
# The distinct corners are counted before anything is allocated, so the ring is
# built at its final length: a cell polygon is synthesized once per candidate
# pair in a plan build, and growing a vector by `push!` there is most of what the
# build allocates.
function _cellring(corners::NTuple{4,USPoint})
    n = _ringlength(corners)
    ring = Vector{USPoint}(undef, n + 1)
    @inbounds begin
        ring[1] = corners[1]
        k = 1
        for j in 2:4
            k == n && break
            corners[j] == ring[k] && continue
            ring[k+=1] = corners[j]
        end
        ring[n+1] = corners[1]
    end
    return ring
end

# How many of the four corners survive dropping consecutive repeats and a tail
# equal to the head — the length of the open ring.
@inline function _ringlength(c::NTuple{4,USPoint})
    n = 1
    last = c[1]
    for j in 2:4
        c[j] == last && continue
        n += 1
        last = c[j]
    end
    return n > 1 && last == c[1] ? n - 1 : n
end

function cellcentroid(space::RasterGrid, i::Int)
    ix, iy = cellsubscript(space, i)
    return space.transform((space.xedges[ix] + space.xedges[ix+1]) / 2,
        (space.yedges[iy] + space.yedges[iy+1]) / 2)
end

"""
    cellat(space::RasterGrid, p) -> Union{Int,Nothing}

The cell containing `p`, by inverting the chart and binary-searching both edge
vectors — `O(log nx + log ny)`, and `nothing` outside the raster's coverage. A
periodic native X coordinate is reduced into the raster's span first, so a
global raster answers for every longitude however its edges are numbered.

A point on a shared edge goes to the cell on the larger-native side; a point on
the raster's own outer edge goes inward. Throws when the space carries no
inverse chart.
"""
function cellat(space::RasterGrid, p)
    space.inverse === nothing && throw(ArgumentError(
        "cellat needs an inverse chart; this RasterGrid was built with `inverse = nothing`"))
    x, y = space.inverse(p)
    ix = _edgeindex(space.xedges, Float64(x), space.xperiod)
    ix === nothing && return nothing
    iy = _edgeindex(space.yedges, Float64(y), nothing)
    iy === nothing && return nothing
    return cellposition(space, ix, iy)
end

function _edgeindex(edges::Vector{Float64}, v::Float64, period)
    n = length(edges) - 1
    lo, hi = minmax(edges[1], edges[end])
    if period !== nothing
        v = lo + mod(v - lo, period)
    end
    (lo <= v <= hi) || return nothing
    k = edges[1] < edges[end] ? searchsortedlast(edges, v) :
        searchsortedfirst(edges, v; rev = true) - 1
    return clamp(k, 1, n)
end

nchunks(space::RasterGrid) = length(space.xchunks) * length(space.ychunks)

"""
    chunksubscript(space::RasterGrid, chunk::Int) -> (cx, cy)

The chunk's lattice coordinates along the X and Y dimensions. Chunk numbers use
the same fastest-first convention as cell positions.
"""
function chunksubscript(space::RasterGrid, chunk::Int)
    1 <= chunk <= nchunks(space) || throw(BoundsError(space, chunk))
    if space.xfast
        cy, cx = fldmod1(chunk, length(space.xchunks))
    else
        cx, cy = fldmod1(chunk, length(space.ychunks))
    end
    return (cx, cy)
end

"""
    chunkbox(space::RasterGrid, chunk::Int) -> (xrange, yrange)

The chunk's cell index ranges along the X and Y dimensions.
"""
function chunkbox(space::RasterGrid, chunk::Int)
    cx, cy = chunksubscript(space, chunk)
    return (space.xchunks[cx], space.ychunks[cy])
end

function cellindices(space::RasterGrid, chunk::Int)
    xr, yr = chunkbox(space, chunk)
    nx, ny = _nx(space), _ny(space)
    if space.xfast
        # A chunk spanning the whole fastest dimension is contiguous in
        # position space, and callers are allowed to exploit that.
        length(xr) == nx && return ((first(yr)-1)*nx+1):(last(yr)*nx)
        out = Vector{Int}(undef, length(xr) * length(yr))
        k = 0
        for iy in yr, ix in xr
            out[k+=1] = ix + (iy - 1) * nx
        end
        return out
    else
        length(yr) == ny && return ((first(xr)-1)*ny+1):(last(xr)*ny)
        out = Vector{Int}(undef, length(xr) * length(yr))
        k = 0
        for ix in xr, iy in yr
            out[k+=1] = iy + (ix - 1) * ny
        end
        return out
    end
end

# ===========================================================================
# Caps
# ===========================================================================

# The `j`-th of `4m` points walking the native box boundary counter-clockwise,
# starting at `(xlo, ylo)`; `j = 0, m, 2m, 3m` are the four corners.
@inline function _boxpoint(t, xlo, xhi, ylo, yhi, m::Int, j::Int)
    e, k = divrem(j, m)
    u = k / m
    return e == 0 ? t(xlo + u * (xhi - xlo), ylo) :
           e == 1 ? t(xhi, ylo + u * (yhi - ylo)) :
           e == 2 ? t(xhi - u * (xhi - xlo), yhi) :
           t(xlo, yhi - u * (yhi - ylo))
end

"""
    _boxcap(space, xlo, xhi, ylo, yhi, nsamp) -> SphericalCap

A cap covering every cell whose native box lies inside `[xlo, xhi] × [ylo, yhi]`.

The cap is centred on the mean of `4(nsamp + 1)` boundary samples and sized to
reach the farthest of them. Covering follows from convexity: under a lon–lat
chart the point of the box farthest from any centre is a corner, so the cap
contains every cell corner in the box, and a cap of radius at most `π/2` is
convex and therefore also contains the geodesic arcs between them — that is, the
cell polygons themselves, which bow outside the box.

Past `π/2` that argument lapses, as it does when the samples average to nothing
and leave no centre; both cases report the whole sphere here, and [`_rectcap`](@ref)
is where a box too wide for this construction gets a bound worth having.

Kept free of that fallback deliberately: this is the innermost function of a plan
build, and a call to anything more than arithmetic in either bail branch costs the
600 000 nodes that never take it more than the 22 that do ever save.
"""
function _boxcap(space::RasterGrid, xlo, xhi, ylo, yhi, nsamp::Int)
    t = space.transform
    m = nsamp + 1
    n = 4m
    sx = sy = sz = 0.0
    for j in 0:(n-1)
        p = _boxpoint(t, xlo, xhi, ylo, yhi, m, j)
        sx += p[1]
        sy += p[2]
        sz += p[3]
    end
    nrm = sqrt(sx^2 + sy^2 + sz^2)
    nrm <= eps(Float64) && return _WHOLE_SPHERE
    centre = USPoint(sx / nrm, sy / nrm, sz / nrm)
    r = 0.0
    for j in 0:(n-1)
        r = max(r, US.spherical_distance(centre,
            _boxpoint(t, xlo, xhi, ylo, yhi, m, j)))
    end
    # Nudged outward past dot-product rounding noise, so containment stays
    # closed; a cap past π/2 loses convexity and with it the covering argument.
    r = nextfloat(r * 1.0001 + 1e-12)
    r > Float64(pi) / 2 && return _WHOLE_SPHERE
    return SphericalCap(centre, r)
end

# The cell cap of `_boxcap(space, cellbox(space, ix, iy)..., 0)`, taken from the
# chart's edge tables when it has them: `_CELL_CAP_SAMPLES == 0` means all four
# sample points are corners, and a corner is a pair of edge coordinates.
function _rastercellcap(space::RasterGrid, ix::Integer, iy::Integer)
    space.tables === nothing &&
        return _boxcap(space, cellbox(space, ix, iy)..., _CELL_CAP_SAMPLES)
    return _cornercap(_cellcorners(space, ix, iy))
end

# `_boxcap`'s construction over four corners already in hand — the same sums in
# the same order over the same points, and the same bail, so the cap is the
# chart-evaluating path's cap bit for bit.
function _cornercap(c::NTuple{4,USPoint})
    sx = sy = sz = 0.0
    for p in c
        sx += p[1]
        sy += p[2]
        sz += p[3]
    end
    nrm = sqrt(sx^2 + sy^2 + sz^2)
    nrm <= eps(Float64) && return _WHOLE_SPHERE
    centre = USPoint(sx / nrm, sy / nrm, sz / nrm)
    r = 0.0
    for p in c
        r = max(r, US.spherical_distance(centre, p))
    end
    r = nextfloat(r * 1.0001 + 1e-12)
    r > Float64(pi) / 2 && return _WHOLE_SPHERE
    return SphericalCap(centre, r)
end

"""
    _widecap(space::RasterGrid, xlo, xhi, ylo, yhi) -> SphericalCap

A cap covering the box's cells without the convexity argument [`_boxcap`](@ref)
rests on — the answer for a box too wide to be bounded by a cap through its own
boundary samples.

A chart that cannot bound its own boxes reports the whole sphere, which is what
a general `transform` gets. Under [`LonLatToSphere`](@ref) the box's cells lie in
a latitude band: an east–west cell edge is a great-circle arc, so it bows to
`atan(tan φ / cos(Δλ/2))` for a cell `Δλ` wide, while the meridians bounding the
box do not bow at all. Two constructions contain that band — a cap about a pole
([`_polarcap`](@ref)), which holds for a box of any longitude span including a
full row, and, for a box narrower than 180° of longitude, a cap on the box's
mid-meridian, whose farthest point is a corner because the distance from a centre
on that meridian falls off monotonically in both coordinates. The tighter is
taken. The whole sphere comes back only when neither beats it, which a band
reaching from pole to pole genuinely does not.
"""
_widecap(space::RasterGrid, xlo, xhi, ylo, yhi) =
    _widecap(space.tables, xlo, xhi, ylo, yhi)

_widecap(::Nothing, xlo, xhi, ylo, yhi) = _WHOLE_SPHERE

function _widecap(t::LonLatEdgeTables, xlo, xhi, ylo, yhi)
    polar = _polarcap(t, ylo, yhi)
    abs(xhi - xlo) / 2 < 90.0 || return polar
    a, b = _bowedband(Float64(ylo), Float64(yhi), t.maxdx)
    chart = LonLatToSphere()
    centre = chart((xlo + xhi) / 2, (a + b) / 2)
    r = 0.0
    for x in (xlo, xhi), y in (a, b)
        r = max(r, US.spherical_distance(centre, chart(x, y)))
    end
    r = nextfloat(r * 1.0001 + 1e-12)
    (r >= Float64(pi) || polar.radius <= r) && return polar
    return SphericalCap(centre, r)
end

"""
    _polarcap(tables, ylo, yhi) -> SphericalCap

The tighter of the two caps that bound the box's latitude band by a pole, or the
whole sphere when neither is tighter.

`{lat ≥ a}` **is** a cap of radius `90° − a` about the north pole and `{lat ≤ b}`
one of radius `90° + b` about the south, so this covers a box of any longitude
span at all — including a full row, which is exactly where a cap through the
box's own boundary has to reach around the sphere. `a` and `b` come from
[`_bowedband`](@ref), so the poleward bow of the cells' east–west edges is inside
them. Cheap enough to offer against every node extent: two `tan`s and a `cos`.
"""
function _polarcap(t::LonLatEdgeTables, ylo, yhi)
    a, b = _bowedband(Float64(ylo), Float64(yhi), t.maxdx)
    centre = USPoint(0.0, 0.0, 1.0)
    r = deg2rad(90.0 - a)
    rs = deg2rad(90.0 + b)
    if rs < r
        centre = USPoint(0.0, 0.0, -1.0)
        r = rs
    end
    r = nextfloat(r * 1.0001 + 1e-12)
    r >= Float64(pi) && return _WHOLE_SPHERE
    return SphericalCap(centre, r)
end

# The latitude band the cells of a box `[ylo, yhi]` tall actually occupy, given
# the widest cell in the raster. A great-circle arc between two points of a
# parallel bows toward the nearer pole, so the band grows outward at whichever
# end faces away from the equator, and not at all at an end on it.
function _bowedband(ylo::Float64, yhi::Float64, dxmax::Float64)
    h = min(abs(dxmax), 360.0) / 2
    c = h >= 90.0 ? 0.0 : cosd(h)
    hi = yhi > 0 ? (c <= 0.0 ? 90.0 : min(90.0, atand(tand(yhi) / c))) : yhi
    lo = ylo < 0 ? (c <= 0.0 ? -90.0 : max(-90.0, atand(tand(ylo) / c))) : ylo
    return (lo, hi)
end

"""
    _rectcap(space, ix0, ix1, iy0, iy1) -> SphericalCap

The extent of an index rectangle — a tree node, a chunk.

Both constructions cover, so the tighter one is taken: sampling the box boundary
([`_boxcap`](@ref)) wins on a compact box, and the chart's own bound
([`_widecap`](@ref)) wins on one wide enough that a cap through its boundary has
to reach around the sphere — a full-longitude band above all, which is how a
global raster is normally chunked and which `_boxcap` alone answers with the
whole sphere.
"""
function _rectcap(space::RasterGrid, ix0::Int, ix1::Int, iy0::Int, iy1::Int)
    xlo, xhi = minmax(space.xedges[ix0], space.xedges[ix1+1])
    ylo, yhi = minmax(space.yedges[iy0], space.yedges[iy1+1])
    sampled = _boxcap(space, xlo, xhi, ylo, yhi, _BOX_CAP_SAMPLES)
    # A cap about a pole reaches at least to the near end of the box's latitude
    # span and the bow only widens that, while the mid-meridian cap is built on
    # the bowed box and so never beats a sampled cap that survived at all — a
    # surviving cap tighter than `90° − |φ|` therefore settles it for nothing.
    # Every deep node of a large raster takes this exit, and there are hundreds
    # of thousands of them against a couple of dozen that do not.
    sampled.radius < Float64(pi) &&
        sampled.radius <= deg2rad(min(90.0 - ylo, 90.0 + yhi)) && return sampled
    wide = _widecap(space, xlo, xhi, ylo, yhi)
    return wide.radius < sampled.radius ? wide : sampled
end

# The winding of `(xlo, ylo), (xhi, ylo), (xhi, yhi), (xlo, yhi)` under the
# chart, from the Newell normal of the widest probe cell: positive against the
# outward radial means counter-clockwise seen from outside. A lon–lat chart is
# always positive; a chart that flips handedness is caught here once, at
# construction, rather than per cell.
function _chartorientation(transform, xedges::Vector{Float64}, yedges::Vector{Float64})
    nx, ny = length(xedges) - 1, length(yedges) - 1
    best = 0.0
    for iy in _probes(ny), ix in _probes(nx)
        xlo, xhi = minmax(xedges[ix], xedges[ix+1])
        ylo, yhi = minmax(yedges[iy], yedges[iy+1])
        s = _quadsign(transform, xlo, xhi, ylo, yhi)
        abs(s) > abs(best) && (best = s)
    end
    best == 0.0 && throw(ArgumentError(
        "RasterGrid: the chart collapses every probe cell to a degenerate ring, " *
        "so the winding of a cell cannot be determined"))
    return best > 0.0
end

_probes(n::Int) = n <= 3 ? (1:n) : (1, cld(n, 2), n)

function _quadsign(t, xlo, xhi, ylo, yhi)
    c = (t(xlo, ylo), t(xhi, ylo), t(xhi, yhi), t(xlo, yhi))
    nx = ny = nz = 0.0
    for k in 1:4
        a, b = c[k], c[k == 4 ? 1 : k+1]
        nx += a[2] * b[3] - a[3] * b[2]
        ny += a[3] * b[1] - a[1] * b[3]
        nz += a[1] * b[2] - a[2] * b[1]
    end
    cx = c[1][1] + c[2][1] + c[3][1] + c[4][1]
    cy = c[1][2] + c[2][2] + c[3][2] + c[4][2]
    cz = c[1][3] + c[2][3] + c[3][3] + c[4][3]
    return nx * cx + ny * cy + nz * cz
end

# ===========================================================================
# Trees
# ===========================================================================

"""
    RasterCellTree(space, ix0, ix1, iy0, iy1)

A `SpatialTreeInterface` tree over the cells of an index rectangle of `space`,
halving the longer side at each level down to a leaf of a few cells.

The tree stores nothing but the rectangle: node extents are derived from the
native box on demand (and marked expensive, so the dual descent caches them),
and leaf entries are cell positions with their own caps. That is what lets a
raster of any size be treed in `O(1)`.
"""
struct RasterCellTree{S<:RasterGrid}
    space::S
    ix0::Int
    ix1::Int
    iy0::Int
    iy1::Int
end

Base.show(io::IO, t::RasterCellTree) = print(io, "RasterCellTree(",
    t.ix0, ":", t.ix1, ", ", t.iy0, ":", t.iy1, ")")

_treecells(t::RasterCellTree) = (t.ix1 - t.ix0 + 1) * (t.iy1 - t.iy0 + 1)

STI.isspatialtree(::Type{<:RasterCellTree}) = true
STI.node_extent_is_expensive(::Type{<:RasterCellTree}) = true
STI.isleaf(t::RasterCellTree) = _treecells(t) <= _CELL_TREE_LEAF
STI.nchild(t::RasterCellTree) = STI.isleaf(t) ? 0 : 2
STI.node_extent(t::RasterCellTree) = _rectcap(t.space, t.ix0, t.ix1, t.iy0, t.iy1)

function STI.getchild(t::RasterCellTree)
    if t.ix1 - t.ix0 >= t.iy1 - t.iy0
        mid = (t.ix0 + t.ix1) ÷ 2
        return (RasterCellTree(t.space, t.ix0, mid, t.iy0, t.iy1),
            RasterCellTree(t.space, mid + 1, t.ix1, t.iy0, t.iy1))
    else
        mid = (t.iy0 + t.iy1) ÷ 2
        return (RasterCellTree(t.space, t.ix0, t.ix1, t.iy0, mid),
            RasterCellTree(t.space, t.ix0, t.ix1, mid + 1, t.iy1))
    end
end

STI.getchild(t::RasterCellTree, i::Int) = STI.getchild(t)[i]

STI.child_indices_extents(t::RasterCellTree) =
    ((cellposition(t.space, ix, iy), _rastercellcap(t.space, ix, iy))
     for iy in t.iy0:t.iy1, ix in t.ix0:t.ix1)

GOCore.best_manifold(t::RasterCellTree) = manifold(t.space)
Trees.ncells(t::RasterCellTree) = ncells(t.space)
Trees.getcell(t::RasterCellTree, i::Int) = getcell(t.space, i)
Trees.getcell(t::RasterCellTree) =
    (getcell(t.space, cellposition(t.space, ix, iy))
     for iy in t.iy0:t.iy1, ix in t.ix0:t.ix1)

"""
    RasterFlatTree(space, indices, caps)

A one-node `SpatialTreeInterface` tree with stored extents, for an index set
with no rectangular shape to exploit — an arbitrary cell subset, or the chunk
lattice.

`indices` are cell positions for a cell tree and chunk numbers for a chunk tree;
`ConservativeRegridding`'s cell accessors read cells either way, since a chunk
tree is never handed to it.
"""
struct RasterFlatTree{S<:RasterGrid}
    space::S
    indices::Vector{Int}
    caps::Vector{Cap}
    extent::Cap
end

function RasterFlatTree(space::RasterGrid, indices, caps)
    ix = collect(Int, indices)
    cs = collect(Cap, caps)
    extent = isempty(cs) ? _WHOLE_SPHERE : foldl(US._merge, cs)
    return RasterFlatTree{typeof(space)}(space, ix, cs, extent)
end

Base.show(io::IO, t::RasterFlatTree) =
    print(io, "RasterFlatTree(", length(t.indices), " entries)")

STI.isspatialtree(::Type{<:RasterFlatTree}) = true
STI.node_extent_is_expensive(::Type{<:RasterFlatTree}) = false
STI.isleaf(::RasterFlatTree) = true
STI.nchild(::RasterFlatTree) = 0
STI.getchild(::RasterFlatTree) = ()
STI.node_extent(t::RasterFlatTree) = t.extent
STI.child_indices_extents(t::RasterFlatTree) = zip(t.indices, t.caps)

GOCore.best_manifold(t::RasterFlatTree) = manifold(t.space)
Trees.ncells(t::RasterFlatTree) = ncells(t.space)
Trees.getcell(t::RasterFlatTree, i::Int) = getcell(t.space, i)
Trees.getcell(t::RasterFlatTree) = (getcell(t.space, i) for i in t.indices)

celltree(space::RasterGrid) = RasterCellTree(space, 1, _nx(space), 1, _ny(space))

"""
    celltree(space::RasterGrid, chunk::Int)
    celltree(space::RasterGrid, indices::AbstractVector{<:Integer})

Cell trees restricted to part of the raster: a chunk keeps its index rectangle
and so keeps the `O(1)` recursive tree; an arbitrary position list gets a flat
tree with one stored cap per cell.

Leaf indices are cell positions in both cases, so a builder localizing them
against its `src_inds`/`dst_inds` sees the same numbering as the whole-space
tree.
"""
function celltree(space::RasterGrid, chunk::Int)
    xr, yr = chunkbox(space, chunk)
    return RasterCellTree(space, first(xr), last(xr), first(yr), last(yr))
end

function celltree(space::RasterGrid, indices::AbstractVector{<:Integer})
    caps = [_rastercellcap(space, cellsubscript(space, Int(i))...) for i in indices]
    return RasterFlatTree(space, indices, caps)
end

function chunktree(space::RasterGrid)
    n = nchunks(space)
    caps = Vector{Cap}(undef, n)
    for c in 1:n
        xr, yr = chunkbox(space, c)
        caps[c] = _rectcap(space, first(xr), last(xr), first(yr), last(yr))
    end
    return RasterFlatTree(space, 1:n, caps)
end

"""
    subtree(space::RasterGrid, inds)

The raster's own tree restricted to `inds`.

An index set that is a rectangle of the lattice — every chunk, and the whole
space — keeps the `O(1)` [`RasterCellTree`](@ref); anything else falls back to a
flat tree with one analytic cap per listed cell. Leaf indices are cell positions
either way, as the generic [`subtree`](@ref) requires.
"""
function subtree(space::RasterGrid, inds)
    rect = _indexrect(space, inds)
    rect === nothing && return celltree(space, inds)
    return RasterCellTree(space, rect...)
end

# The number of cells along the dimension that varies fastest in cell positions.
_nfast(space::RasterGrid) = space.xfast ? _nx(space) : _ny(space)

"""
    _indexrect(space::RasterGrid, inds) -> (ix0, ix1, iy0, iy1) or nothing

The index rectangle `inds` enumerates, or `nothing` when it is not one.

`inds` is a rectangle when it lists exactly the cells of a lattice box in the
space's own position order, which is what [`cellindices`](@ref) produces for a
chunk and what a whole-space `1:ncells` is.
"""
function _indexrect(space::RasterGrid, inds)
    n = ncells(space)
    (isempty(inds) || length(inds) > n) && return nothing
    nfast = _nfast(space)
    lo, hi = Int(first(inds)), Int(last(inds))
    (1 <= lo <= n && 1 <= hi <= n) || return nothing
    b0, a0 = fldmod1(lo, nfast)
    b1, a1 = fldmod1(hi, nfast)
    (a1 >= a0 && b1 >= b0) || return nothing
    length(inds) == (a1 - a0 + 1) * (b1 - b0 + 1) || return nothing
    k = 0
    for b in b0:b1, a in a0:a1
        Int(inds[k+=1]) == a + (b - 1) * nfast || return nothing
    end
    return space.xfast ? (a0, a1, b0, b1) : (b0, b1, a0, a1)
end

# ===========================================================================
# The cell chart
# ===========================================================================
#
# `chartaxes` and its companions are declared in `interpolation.jl`, which the
# module includes after this file; these methods therefore precede the generic
# declarations they extend. Nothing is called at load time, so the order is a
# reading inconvenience only.

"""
    chartaxes(space::RasterGrid)

The cell-centre coordinates along each spatial dimension, in **lookup order** —
descending for a reverse-ordered lookup, which the stencil locator handles.
Cell `(ix, iy)` is centred at `(xs[ix], ys[iy])`, the same subscripting
[`cellposition`](@ref) inverts.
"""
chartaxes(space::RasterGrid) = (_centres(space.xedges), _centres(space.yedges))

_centres(e::Vector{Float64}) = [(e[k] + e[k+1]) / 2 for k in 1:(length(e)-1)]

"""
    chartcoords(space::RasterGrid, p)

`p` in the raster's native coordinates, or `nothing` when the chart cannot be
inverted. A periodic native X is answered on the branch the edge vector is
written on, so a raster spanning `0:360` places a point at 350 rather than -10.
"""
function chartcoords(space::RasterGrid, p)
    space.inverse === nothing && return nothing
    x, y = space.inverse(p)
    return (_onbranch(space.xedges, Float64(x), space.xperiod), Float64(y))
end

function _onbranch(edges::Vector{Float64}, v::Float64, period)
    period === nothing && return v
    lo = min(edges[1], edges[end])
    return lo + mod(v - lo, period::Float64)
end

chartposition(space::RasterGrid, ix::Int, iy::Int) = cellposition(space, ix, iy)

"""
    chartperiod(space::RasterGrid)

`(360.0, nothing)` for a raster whose X edges span the chart's whole period, and
`(nothing, nothing)` otherwise.

The chart being periodic is not enough: a regional raster on a longitude chart
does not close, and reporting a period would wrap its stencils onto the far edge
of the domain.
"""
function chartperiod(space::RasterGrid)
    p = space.xperiod
    p === nothing && return (nothing, nothing)
    span = abs(space.xedges[end] - space.xedges[1])
    return (isapprox(span, p; rtol = 1e-9) ? p : nothing, nothing)
end

"""
    chartspacing(space::RasterGrid)

The largest edge step along each dimension, converted to radians of arc by
[`chartarcs`](@ref). Adjacent cell centres are `(e[k+2] - e[k]) / 2` apart in
native units, which the largest edge step bounds.
"""
chartspacing(space::RasterGrid) =
    chartarcs(space.transform, _maxstep(space.xedges), _maxstep(space.yedges))

_maxstep(e::Vector{Float64}) =
    maximum(abs(e[k+1] - e[k]) for k in 1:(length(e)-1))
