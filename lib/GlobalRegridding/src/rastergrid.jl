# Dimensional-raster implementation of `RegridSpace`.

# Chart interface

"""
    LonLatToSphere()

Convert longitude and latitude in degrees to a `UnitSphericalPoint`.
"""
struct LonLatToSphere end

function (::LonLatToSphere)(lon::Real, lat::Real)
    coslat = cosd(Float64(lat))
    return USPoint(coslat * cosd(Float64(lon)), coslat * sind(Float64(lon)),
        sind(Float64(lat)))
end

"""
    SphereToLonLat()

Convert a unit-sphere point to longitude and latitude in degrees.
"""
struct SphereToLonLat end

(::SphereToLonLat)(p) = (atand(p[2], p[1]), asind(clamp(p[3], -1.0, 1.0)))

"""
    chartinverse(transform) -> callable or nothing

Return a chart's inverse, or `nothing` when unavailable.
"""
chartinverse(::Any) = nothing
chartinverse(::LonLatToSphere) = SphereToLonLat()

"""
    chartlimits(transform) -> (xperiod, ybounds)

Return the first-coordinate period and second-coordinate bounds. Either may be
`nothing`. The default is unconstrained and non-periodic.
"""
chartlimits(::Any) = (nothing, nothing)
chartlimits(::LonLatToSphere) = (360.0, (-90.0, 90.0))

"""
    chartarcs(transform, dx, dy) -> (Δx, Δy)

Return angular upper bounds, in radians, for native chart steps `dx` and `dy`.
These bounds are used for interpolation support discovery.
"""
chartarcs(t::Any, dx, dy) = throw(ArgumentError(
    "a RasterGrid on a $(typeof(t)) chart cannot bound its cell spacing in " *
    "radians of arc; define GlobalRegridding.chartarcs(::$(typeof(t)), dx, dy) " *
    "to use BilinearPoint with it."))

# Longitude distance is bounded by Δλ; latitude distance equals Δφ.
chartarcs(::LonLatToSphere, dx, dy) = (deg2rad(Float64(dx)), deg2rad(Float64(dy)))

"""
    LonLatEdgeTables(xedges, yedges)

Cache sine and cosine values for longitude/latitude edges, plus the largest
longitude step. This avoids repeated chart evaluation for cell corners.
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

Return optional per-edge lookup tables for a chart. `nothing` uses direct chart
evaluation.
"""
chartedgetables(::Any, xedges, yedges) = nothing
chartedgetables(::LonLatToSphere, xedges::Vector{Float64}, yedges::Vector{Float64}) =
    LonLatEdgeTables(xedges, yedges)

# Raster space

"""
    RasterGrid(A::DimensionalData.AbstractDimArray; kwargs...)
    RasterGrid(dims::Tuple; chunks = nothing, kwargs...)

Build a [`RegridSpace`](@ref) from a dimensional raster or dimension tuple.
Construction reads no array values. Spatial lookups provide cell edges and
array chunking provides spatial chunks. X/Y traits are preferred, then common
dimension names; `xdim` and `ydim` override detection. Interval lookups use
their bounds, while point lookups are midpointed. Edge order matches lookup order.

# Cell order

Cells follow the array's dimension order, fastest dimension first. Thus `vec`
of a spatial slice matches cell position order. Chunk numbers use the same rule.

# Keywords

  - `xdim`, `ydim`: explicit spatial dimensions.
  - `transform`: native `(x, y)` to the unit sphere.
  - `inverse`: unit sphere to native coordinates; `nothing` disables `cellat`.
  - `xperiod`, `ybounds`: native chart limits.
  - `chunks`: `(xchunks, ychunks)` index ranges for the dimension-tuple form.

# Geometry

Cells are geodesic quadrilaterals with counter-clockwise exterior rings.
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

# Four corners are sufficient for a lon/lat cell cap.
const _CELL_CAP_SAMPLES = 0
# Extra boundary samples for multi-cell boxes.
const _BOX_CAP_SAMPLES = 3

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

Return spatial slice shape in array dimension order. `vec` follows cell positions.
"""
rastersize(space::RasterGrid) =
    space.xfast ? (_nx(space), _ny(space)) : (_ny(space), _nx(space))

"""
    DimensionalData.dims(space::RasterGrid)

Return spatial dimensions in array order — the order the space was
constructed with, whichever of X and Y comes first.
"""
DD.dims(space::RasterGrid) =
    space.xfast ? (space.xdim, space.ydim) : (space.ydim, space.xdim)

"""
    destinationdims(space::RasterGrid, sampling)

Return [`DimensionalData.dims`](@ref) relabelled to carry `sampling`, so
regrid output echoes the space's construction order. A dimension whose
lookup already samples that way is returned untouched;
otherwise its lookup is rebuilt, with `Intervals` bounds taken from the cell
edges and `Points` values from the cell centres.
"""
destinationdims(space::RasterGrid, sampling::DD.Lookups.Sampling) =
    space.xfast ?
    (_withsampling(space.xdim, space.xedges, sampling),
        _withsampling(space.ydim, space.yedges, sampling)) :
    (_withsampling(space.ydim, space.yedges, sampling),
        _withsampling(space.xdim, space.xedges, sampling))

function _withsampling(dim, edges::Vector{Float64}, sampling::DD.Lookups.Points)
    lk = DD.lookup(dim)
    DD.sampling(lk) isa DD.Lookups.Points && return dim
    return DD.rebuild(dim, DD.Lookups.Sampled(_centres(edges); order = DD.order(lk),
        sampling, metadata = DD.metadata(lk)))
end

function _withsampling(dim, edges::Vector{Float64}, sampling::DD.Lookups.Intervals)
    lk = DD.lookup(dim)
    DD.sampling(lk) isa DD.Lookups.Intervals && return dim
    n = length(edges) - 1
    bounds = Matrix{Float64}(undef, 2, n)
    for k in 1:n
        bounds[1, k] = edges[k]
        bounds[2, k] = edges[k+1]
    end
    return DD.rebuild(dim, DD.Lookups.Sampled(DD.val(lk); order = DD.order(lk),
        span = DD.Lookups.Explicit(bounds), sampling, metadata = DD.metadata(lk)))
end

# Dimension and edge extraction

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

# Build cell edges in lookup order from intervals or point midpoints.
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
    # A single edge vector cannot represent gaps between intervals.
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

# Chunk extraction

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

# `RegridSpace` interface

ncells(space::RasterGrid) = _nx(space) * _ny(space)

manifold(::RasterGrid) = GOCore.Spherical(; radius = 1.0)

# Interpolation requires an inverse transform into the cell-centre chart.
hascellchart(space::RasterGrid) = space.inverse !== nothing

"""
    cellsubscript(space::RasterGrid, i::Int) -> (ix, iy)
    cellposition(space::RasterGrid, ix::Integer, iy::Integer) -> Int

Convert between cell positions and `(ix, iy)` lattice coordinates. Position
order follows the array's fastest dimension.
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

Return the cell's ascending native-coordinate bounds.
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

Return corners in ascending native order, using edge tables when available.
`space.ccw` records whether polygon construction must reverse them.
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

# Return edge indices in ascending native-coordinate order.
@inline function _edgeorder(e::Vector{Float64}, k::Integer)
    i = Int(k)
    return @inbounds e[i] <= e[i+1] ? (i, i + 1) : (i + 1, i)
end

# Drop coincident polar corners and allocate the closed ring at its final size.
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

# Count distinct consecutive corners before closure.
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

Return the cell containing `p` in `O(log nx + log ny)`, or `nothing` outside
coverage. Periodic X coordinates wrap into the raster span. Shared edges select
the larger-native cell; outer edges select the interior cell. Requires an inverse chart.
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

Return a chunk's `(cx, cy)` lattice coordinates.
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

Return a chunk's X and Y cell-index ranges.
"""
function chunkbox(space::RasterGrid, chunk::Int)
    cx, cy = chunksubscript(space, chunk)
    return (space.xchunks[cx], space.ychunks[cy])
end

"""
    chunkposition(space::RasterGrid, cx::Integer, cy::Integer) -> Int

Return the chunk number at `(cx, cy)`.
"""
function chunkposition(space::RasterGrid, cx::Integer, cy::Integer)
    1 <= cx <= length(space.xchunks) && 1 <= cy <= length(space.ychunks) ||
        throw(BoundsError(space, (cx, cy)))
    return space.xfast ? Int(cx) + (Int(cy) - 1) * length(space.xchunks) :
           Int(cy) + (Int(cx) - 1) * length(space.ychunks)
end

# Locate an axis index by binary search without allocating chunk starts.
function _chunkofindex(ranges::Vector{UnitRange{Int}}, k::Int)
    lo, hi = 1, length(ranges)
    while lo < hi
        mid = (lo + hi + 1) >> 1
        first(ranges[mid]) <= k ? (lo = mid) : (hi = mid - 1)
    end
    return lo
end

# Locate a cell's chunk with one binary search per axis.
function chunkat(space::RasterGrid, i::Integer)
    ix, iy = cellsubscript(space, Int(i))
    return chunkposition(space, _chunkofindex(space.xchunks, ix),
        _chunkofindex(space.ychunks, iy))
end

function cellindices(space::RasterGrid, chunk::Int)
    xr, yr = chunkbox(space, chunk)
    nx, ny = _nx(space), _ny(space)
    if space.xfast
        # Full-width chunks are contiguous in position order.
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

# Spherical caps

# Return boundary sample `j`; multiples of `m` are corners.
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

Return a cap around `4(nsamp + 1)` boundary samples. For lon/lat boxes, a cap no
wider than `π/2` contains the corners and their geodesic edges. Return the whole
sphere when no stable convex cap exists; [`_rectcap`](@ref) handles wide boxes.
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
    # Caps beyond π/2 lose the convexity bound.
    r = _padcap(r)
    r > Float64(pi) / 2 && return _WHOLE_SPHERE
    return SphericalCap(centre, r)
end

# Use tabulated corners for a cell cap when available.
function _rastercellcap(space::RasterGrid, ix::Integer, iy::Integer)
    space.tables === nothing &&
        return _boxcap(space, cellbox(space, ix, iy)..., _CELL_CAP_SAMPLES)
    return _cornercap(_cellcorners(space, ix, iy))
end

# Build the same cap as `_boxcap(..., 0)` from precomputed corners.
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
    r = _padcap(r)
    r > Float64(pi) / 2 && return _WHOLE_SPHERE
    return SphericalCap(centre, r)
end

"""
    _widecap(space::RasterGrid, xlo, xhi, ylo, yhi) -> SphericalCap

Return a safe cap for boxes too wide for [`_boxcap`](@ref). Lon/lat charts use
the tighter of a polar latitude-band cap and, below 180° width, a mid-meridian
cap. Other charts return the whole sphere.
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
    r = _padcap(r)
    (r >= Float64(pi) || polar.radius <= r) && return polar
    return SphericalCap(centre, r)
end

"""
    _polarcap(tables, ylo, yhi) -> SphericalCap

Return the tighter north- or south-pole cap containing the bowed latitude band.
This remains valid for any longitude span.
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
    r = _padcap(r)
    r >= Float64(pi) && return _WHOLE_SPHERE
    return SphericalCap(centre, r)
end

# Expand a latitude band to include poleward-bowing great-circle edges.
function _bowedband(ylo::Float64, yhi::Float64, dxmax::Float64)
    h = min(abs(dxmax), 360.0) / 2
    c = h >= 90.0 ? 0.0 : cosd(h)
    hi = yhi > 0 ? (c <= 0.0 ? 90.0 : min(90.0, atand(tand(yhi) / c))) : yhi
    lo = ylo < 0 ? (c <= 0.0 ? -90.0 : max(-90.0, atand(tand(ylo) / c))) : ylo
    return (lo, hi)
end

# Read an index rectangle's corners from edge tables.
@inline _boxcorners(t::LonLatEdgeTables, ix0::Int, ix1::Int, iy0::Int, iy1::Int) =
    (_tablecorner(t, ix0, iy0), _tablecorner(t, ix1 + 1, iy0),
        _tablecorner(t, ix1 + 1, iy1 + 1), _tablecorner(t, ix0, iy1 + 1))

"""
    _sampledcap(tables, space, ix0, ix1, iy0, iy1, xlo, xhi, ylo, yhi) -> SphericalCap

Return the sampled cap used by [`_rectcap`](@ref). Lon/lat boxes narrower than
180° use four tabulated corners; wider or untabulated boxes sample the boundary.
"""
@inline _sampledcap(::Nothing, space::RasterGrid, ix0::Int, ix1::Int, iy0::Int,
    iy1::Int, xlo, xhi, ylo, yhi) = _boxcap(space, xlo, xhi, ylo, yhi, _BOX_CAP_SAMPLES)

@inline function _sampledcap(t::LonLatEdgeTables, space::RasterGrid, ix0::Int, ix1::Int,
    iy0::Int, iy1::Int, xlo, xhi, ylo, yhi)
    xhi - xlo < 180.0 || return _boxcap(space, xlo, xhi, ylo, yhi, _BOX_CAP_SAMPLES)
    return _cornercap(_boxcorners(t, ix0, ix1, iy0, iy1))
end

"""
    _rectcap(space, ix0, ix1, iy0, iy1) -> SphericalCap

Return the tighter sampled or wide-box cap for a tree-node rectangle.
"""
function _rectcap(space::RasterGrid, ix0::Int, ix1::Int, iy0::Int, iy1::Int)
    xlo, xhi = minmax(space.xedges[ix0], space.xedges[ix1+1])
    ylo, yhi = minmax(space.yedges[iy0], space.yedges[iy1+1])
    sampled = _sampledcap(space.tables, space, ix0, ix1, iy0, iy1, xlo, xhi, ylo, yhi)
    return _tighterof(space, xlo, xhi, ylo, yhi, sampled)
end

"""
    _chunkcap(space, ix0, ix1, iy0, iy1) -> SphericalCap

Return a chunk rectangle's cap. Extra boundary samples improve the centre and
reduce false-positive chunk connections. This is computed once per chunk.
"""
function _chunkcap(space::RasterGrid, ix0::Int, ix1::Int, iy0::Int, iy1::Int)
    xlo, xhi = minmax(space.xedges[ix0], space.xedges[ix1+1])
    ylo, yhi = minmax(space.yedges[iy0], space.yedges[iy1+1])
    sampled = _boxcap(space, xlo, xhi, ylo, yhi, _BOX_CAP_SAMPLES)
    return _tighterof(space, xlo, xhi, ylo, yhi, sampled)
end

# Skip the wide-box calculation when the sampled cap already beats its lower bound.
@inline function _tighterof(space::RasterGrid, xlo, xhi, ylo, yhi, sampled::Cap)
    sampled.radius < Float64(pi) &&
        sampled.radius <= deg2rad(min(90.0 - ylo, 90.0 + yhi)) && return sampled
    wide = _widecap(space, xlo, xhi, ylo, yhi)
    return wide.radius < sampled.radius ? wide : sampled
end

# Determine chart handedness once from probe-cell Newell normals.
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

# Spatial trees

"""
    RasterCellTree(space, ix0, ix1, iy0, iy1)

A recursive spatial tree over a raster index rectangle. It splits the longer
side and derives extents on demand, so construction is `O(1)`.
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

# ConservativeRegridding fetches matrix sizes and polygons through these
# bindings during `intersection_areas`.
GOCore.best_manifold(t::RasterCellTree) = manifold(t.space)
Trees.ncells(t::RasterCellTree) = ncells(t.space)
# `Trees.ncells` answers for the whole space, so the frontier's default estimate
# would be wrong here; the node's index rectangle is exact.
Trees.split_weight(t::RasterCellTree) = _treecells(t)
Trees.getcell(t::RasterCellTree, i::Int) = getcell(t.space, i)
Trees.getcell(t::RasterCellTree) =
    (getcell(t.space, cellposition(t.space, ix, iy))
     for iy in t.iy0:t.iy1, ix in t.ix0:t.ix1)

"""
    RasterGridView(space)

A zero-copy `ConservativeRegridding` curvilinear-grid view of a `RasterGrid`.
Its first index is whichever spatial dimension is fastest in the source array,
so the existing `TopDownQuadtreeCursor` reports the same linear cell positions
as the regridding space. The wrapped space retains the actual task-local chart
transform; raster storage chunking remains solely in `space.xchunks` and
`space.ychunks`, populated from `DiskArrays.eachchunk`.
"""
struct RasterGridView{M<:GOCore.Manifold,S<:RasterGrid} <: Trees.AbstractCurvilinearGrid{M}
    manifold::M
    space::S
end

RasterGridView(space::RasterGrid) = RasterGridView(manifold(space), space)

GOCore.manifold(grid::RasterGridView) = grid.manifold

Trees.ncells(grid::RasterGridView, dim::Int) = if grid.space.xfast
    dim == 1 ? _nx(grid.space) : dim == 2 ? _ny(grid.space) : throw(BoundsError(grid, dim))
else
    dim == 1 ? _ny(grid.space) : dim == 2 ? _nx(grid.space) : throw(BoundsError(grid, dim))
end

@inline function _viewsubscript(grid::RasterGridView, i::Int, j::Int)
    return grid.space.xfast ? (i, j) : (j, i)
end

function Trees.getcell(grid::RasterGridView, i::Int, j::Int)
    ix, iy = _viewsubscript(grid, i, j)
    return getcell(grid.space, cellposition(grid.space, ix, iy))
end

function Trees.getvertex(grid::RasterGridView, i::Int, j::Int)
    ix, iy = _viewsubscript(grid, i, j)
    space = grid.space
    if space.tables === nothing
        return space.transform(space.xedges[ix], space.yedges[iy])
    end
    return _tablecorner(space.tables, ix, iy)
end

function Trees.cell_range_extent(grid::RasterGridView, irange::UnitRange{Int},
        jrange::UnitRange{Int})
    xr, yr = grid.space.xfast ? (irange, jrange) : (jrange, irange)
    return _rectcap(grid.space, first(xr), last(xr), first(yr), last(yr))
end

# The range cap is computed from the chart rather than read directly from axis
# coordinates. This lets ConservativeRegridding's own traversal policy cache
# child extents when it is used in a dual search.
Trees.extent_is_expensive(::Type{<:RasterGridView}) = true

_rasterchunkcursor(space::RasterGrid) = Trees.TopDownQuadtreeCursor(RasterGridView(space))

"""
    RasterFlatTree(space, indices, caps)

A one-node spatial tree with stored extents for arbitrary cells or the chunk
lattice. `indices` are cell positions or chunk numbers, respectively.
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
# Same whole-space `Trees.ncells` caveat; the stored entries are the node.
Trees.split_weight(t::RasterFlatTree) = length(t.indices)
Trees.getcell(t::RasterFlatTree, i::Int) = getcell(t.space, i)
Trees.getcell(t::RasterFlatTree) = (getcell(t.space, i) for i in t.indices)

celltree(space::RasterGrid) = RasterCellTree(space, 1, _nx(space), 1, _ny(space))

"""
    celltree(space::RasterGrid, chunk::Int)
    celltree(space::RasterGrid, indices::AbstractVector{<:Integer})

Return a recursive tree for a chunk rectangle or a flat tree for arbitrary cell
positions. Leaves always use global cell positions.
"""
function celltree(space::RasterGrid, chunk::Int)
    xr, yr = chunkbox(space, chunk)
    return RasterCellTree(space, first(xr), last(xr), first(yr), last(yr))
end

function celltree(space::RasterGrid, indices::AbstractVector{<:Integer})
    caps = [_rastercellcap(space, cellsubscript(space, Int(i))...) for i in indices]
    return RasterFlatTree(space, indices, caps)
end

function chunkextents(space::RasterGrid)
    n = nchunks(space)
    caps = Vector{Cap}(undef, n)
    for c in 1:n
        xr, yr = chunkbox(space, c)
        caps[c] = _chunkcap(space, first(xr), last(xr), first(yr), last(yr))
    end
    return caps
end

# Compatibility for callers that still request the old public tree. Production
# chunk discovery uses `_rasterchunkcursor` and never constructs this flat node.
chunktree(space::RasterGrid) = RasterFlatTree(space, 1:nchunks(space), chunkextents(space))

"""
    subtree(space::RasterGrid, inds)

Return a memoized recursive tree when `inds` form a lattice rectangle,
otherwise a flat tree with one cap per cell.
"""
function subtree(space::RasterGrid, inds)
    rect = _indexrect(space, inds)
    rect === nothing && return celltree(space, inds)
    return MemoRasterTree(RasterCellTree(space, rect...))
end

# The number of cells along the dimension that varies fastest in cell positions.
_nfast(space::RasterGrid) = space.xfast ? _nx(space) : _ny(space)

"""
    _indexrect(space::RasterGrid, inds) -> (ix0, ix1, iy0, iy1) or nothing

Return the lattice rectangle enumerated by `inds`, or `nothing`. A contiguous
range needs no scan: matching the corner subscripts' cell count already proves
the rectangle, so only scattered index sets are walked.
"""
function _indexrect(space::RasterGrid, inds)
    corners = _cornersubscripts(space, inds)
    corners === nothing && return nothing
    a0, a1, b0, b1 = corners
    nfast = _nfast(space)
    k = 0
    for b in b0:b1, a in a0:a1
        Int(inds[k+=1]) == a + (b - 1) * nfast || return nothing
    end
    return _lattice(space, corners)
end

function _indexrect(space::RasterGrid, inds::AbstractUnitRange{<:Integer})
    corners = _cornersubscripts(space, inds)
    corners === nothing && return nothing
    return _lattice(space, corners)
end

# Return the candidate rectangle's fast/slow bounds, or `nothing`.
function _cornersubscripts(space::RasterGrid, inds)
    n = ncells(space)
    (isempty(inds) || length(inds) > n) && return nothing
    nfast = _nfast(space)
    lo, hi = Int(first(inds)), Int(last(inds))
    (1 <= lo <= n && 1 <= hi <= n) || return nothing
    b0, a0 = fldmod1(lo, nfast)
    b1, a1 = fldmod1(hi, nfast)
    (a1 >= a0 && b1 >= b0) || return nothing
    length(inds) == (a1 - a0 + 1) * (b1 - b0 + 1) || return nothing
    return (a0, a1, b0, b1)
end

_lattice(space::RasterGrid, (a0, a1, b0, b1)) =
    space.xfast ? (a0, a1, b0, b1) : (b0, b1, a0, a1)

# Cell chart

"""
    chartaxes(space::RasterGrid)

Return cell-centre coordinates in lookup order. Reverse-ordered dimensions
therefore return descending coordinates.
"""
chartaxes(space::RasterGrid) = (_centres(space.xedges), _centres(space.yedges))

_centres(e::Vector{Float64}) = [(e[k] + e[k+1]) / 2 for k in 1:(length(e)-1)]

"""
    chartcoords(space::RasterGrid, p)

Return `p` in native coordinates, or `nothing` without an inverse. Periodic X
returns the representative nearest the edge span.
"""
function chartcoords(space::RasterGrid, p)
    space.inverse === nothing && return nothing
    x, y = space.inverse(p)
    return (_onbranch(space.xedges, Float64(x), space.xperiod), Float64(y))
end

# nearest periodic representative: folding to [lo, lo+p) throws points just west of lo a period east
function _onbranch(edges::Vector{Float64}, v::Float64, period)
    period === nothing && return v
    p = period::Float64
    mid = (edges[1] + edges[end]) / 2
    return v - p * round((v - mid) / p)
end

chartposition(space::RasterGrid, ix::Int, iy::Int) = cellposition(space, ix, iy)

"""
    chartperiod(space::RasterGrid)

Return the X period only when the raster spans the chart's full period. Regional
rasters do not wrap.
"""
function chartperiod(space::RasterGrid)
    p = space.xperiod
    p === nothing && return (nothing, nothing)
    span = abs(space.xedges[end] - space.xedges[1])
    return (isapprox(span, p; rtol = 1e-9) ? p : nothing, nothing)
end

"""
    chartspacing(space::RasterGrid)

Return the largest edge step on each axis in radians of arc.
"""
chartspacing(space::RasterGrid) =
    chartarcs(space.transform, _maxstep(space.xedges), _maxstep(space.yedges))

_maxstep(e::Vector{Float64}) =
    maximum(abs(e[k+1] - e[k]) for k in 1:(length(e)-1))
