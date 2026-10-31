# A dimensional raster as a GridInterface grid: a one-face `Structured` grid
# whose positions are the array's storage order and whose face lattice is the
# native axes ascending, ordered so the lattice is counter-clockwise.

# Native-coordinate transform interface

# Constructor-only sentinel, needed to distinguish an omitted inverse from an
# explicit `unit_sphere_to_native = nothing`.
struct _UnsetRasterKeyword end
const _UNSET_RASTER_KEYWORD = _UnsetRasterKeyword()

# Compatibility adapter for the original `transform = (x, y) -> ...` contract.
# Raster geometry calls only the one-coordinate form after construction.
struct _TwoArgumentNativeToUnitSphere{F}
    f::F
end

@inline (t::_TwoArgumentNativeToUnitSphere)(xy) = t.f(xy[1], xy[2])

function _one_coordinate_native_to_unit_sphere(f)
    applicable(f, (0.0, 0.0)) && return f
    applicable(f, 0.0, 0.0) && return _TwoArgumentNativeToUnitSphere(f)
    throw(ArgumentError(
        "native_to_unit_sphere must accept one `(x, y)` coordinate; the legacy " *
        "two-argument form `(x, y) -> ...` is also supported"))
end

function _one_coordinate_unit_sphere_to_native(f)
    f === nothing && return nothing
    applicable(f, USPoint(1.0, 0.0, 0.0)) && return f
    throw(ArgumentError(
        "unit_sphere_to_native must accept one unit-sphere point coordinate"))
end

"""
    _prepare_raster_transform_pair(native_to_unit_sphere, unit_sphere_to_native)

Prepare the two callable transformations stored by `RasterGrid`.
`unit_sphere_to_native` may be `_UNSET_RASTER_KEYWORD`, in which case the
forward transformation's default inverse is used. Extensions may specialize
this hook to return lightweight views which obtain task-local transformation
resources when called.
"""
function _prepare_raster_transform_pair(forward, backward)
    prepared_forward = _one_coordinate_native_to_unit_sphere(forward)
    raw_backward = backward === _UNSET_RASTER_KEYWORD ?
        _default_unit_sphere_to_native(forward) : backward
    prepared_backward = _one_coordinate_unit_sphere_to_native(raw_backward)
    return (prepared_forward, prepared_backward)
end

"""
    _task_prepared_raster_transform(transform)

Return a callable prepared for repeated use by the current task. Ordinary Julia
transformations need no preparation. Extensions may hoist task-owned resources
out of multi-point loops without storing them on a shared `RasterGrid`.
"""
_task_prepared_raster_transform(transform) = transform

"""
    _default_unit_sphere_to_native(native_to_unit_sphere) -> callable or nothing

Return the default inverse of a native-to-unit-sphere transformation, or
`nothing` when no inverse is known.
"""
_default_unit_sphere_to_native(::Any) = nothing
_default_unit_sphere_to_native(t::US.UnitSphereFromGeographic) = inv(t)

"""
    _crs_native_to_unit_sphere(A) -> transform or _UNSET_RASTER_KEYWORD

Return the forward chart implied by the CRS metadata of a dimensional array or
dimension tuple, or `_UNSET_RASTER_KEYWORD` when it carries none or names
longitude/latitude in Greenwich degrees. Plain dimensional data carries no CRS. The
Rasters extension reads `Projected` and `Mapped` lookups and calls
[`_projected_crs_native_to_unit_sphere`](@ref) for any CRS it cannot classify
as Greenwich degrees on its own.
"""
_crs_native_to_unit_sphere(::Any) = _UNSET_RASTER_KEYWORD

"""
    _projected_crs_native_to_unit_sphere(crs) -> transform or _UNSET_RASTER_KEYWORD

Return the forward chart for a CRS value requiring coordinate normalization.
Proj reads geographic prime meridians and angular units or transforms projected
coordinates to geographic longitude/latitude, so the generic method throws: without Proj the CRS cannot
be interpreted, and silently treating it as longitude/latitude would be wrong.
"""
_projected_crs_native_to_unit_sphere(crs) = throw(ArgumentError(
    "the raster's CRS $(repr(crs)) is not EPSG:4326 or another form that " *
    "GlobalRegridding recognizes as geographic longitude/latitude without " *
    "Proj. Load Proj (`import Proj`) so the CRS can be classified and a " *
    "transformation derived from it, or pass `native_to_unit_sphere` " *
    "explicitly. If Proj is already loaded, the GlobalRegriddingRastersProjExt " *
    "extension did not load; check " *
    "`Base.get_extension(GlobalRegridding, :GlobalRegriddingRastersProjExt)` " *
    "and rebuild stale precompilation with `Pkg.precompile()`."))

function _with_crs_transform(A, kwargs)
    (haskey(kwargs, :native_to_unit_sphere) || haskey(kwargs, :transform)) &&
        return kwargs
    t = _crs_native_to_unit_sphere(A)
    t === _UNSET_RASTER_KEYWORD && return kwargs
    return (; kwargs..., native_to_unit_sphere = t)
end

"""
    _native_coordinate_limits(native_to_unit_sphere) -> (xperiod, ybounds)

Return the first native coordinate's period and the second native coordinate's
bounds. Either may be `nothing`. Known longitude/latitude charts supply defaults.
"""
_native_coordinate_limits(::Any) = (nothing, nothing)
_native_coordinate_limits(::US.UnitSphereFromGeographic) = (360.0, (-90.0, 90.0))

"""
    _spherical_step_bounds_radians(native_to_unit_sphere, dx, dy) -> (Δx, Δy)

Return angular upper bounds, in radians, for native coordinate steps `dx` and
`dy`. These bounds are used for interpolation support discovery.
"""
_spherical_step_bounds_radians(t::Any, dx, dy) = throw(ArgumentError(
    "a RasterGrid using $(typeof(t)) cannot bound its native cell spacing in " *
    "radians of unit-sphere arc; define " *
    "GlobalRegridding._spherical_step_bounds_radians(::$(typeof(t)), dx, dy) " *
    "to interpolate on its chart."))

# Longitude distance is bounded by Δλ; latitude distance equals Δφ.
_spherical_step_bounds_radians(::US.UnitSphereFromGeographic, dx, dy) =
    (deg2rad(Float64(dx)), deg2rad(Float64(dy)))
_spherical_step_bounds_radians(t::_TwoArgumentNativeToUnitSphere, dx, dy) =
    _spherical_step_bounds_radians(t.f, dx, dy)

"""
    GeographicEdgeTables(xedges, yedges)

Cache sine and cosine values for longitude/latitude edges, plus the largest
longitude step. This avoids repeated chart evaluation for cell corners.
"""
struct GeographicEdgeTables
    cx::Vector{Float64}
    sx::Vector{Float64}
    cy::Vector{Float64}
    sy::Vector{Float64}
    maxdx::Float64
end

GeographicEdgeTables(xe::Vector{Float64}, ye::Vector{Float64}) =
    GeographicEdgeTables(cosd.(xe), sind.(xe), cosd.(ye), sind.(ye), _maxstep(xe))

@inline _tablecorner(t::GeographicEdgeTables, i::Int, j::Int) =
    @inbounds USPoint(t.cy[j] * t.cx[i], t.cy[j] * t.sx[i], t.sy[j])

"""
    _geographic_edge_tables(native_to_unit_sphere, xedges, yedges) -> tables or nothing

Return optional per-edge lookup tables for a chart. `nothing` uses direct chart
evaluation.
"""
_geographic_edge_tables(::Any, xedges, yedges) = nothing
_geographic_edge_tables(::US.UnitSphereFromGeographic,
    xedges::Vector{Float64}, yedges::Vector{Float64}) =
    GeographicEdgeTables(xedges, yedges)

# Raster space

"""
    RasterGrid(A::DimensionalData.AbstractDimArray; kwargs...)
    RasterGrid(dims::Tuple; chunks = nothing, kwargs...)

Build a GridInterface grid from a dimensional raster or dimension tuple: a
one-face `GridInterface.Structured` grid whose positions are cell indices in
the array's storage order. Construction reads no array values. Spatial
lookups provide cell edges and array chunking provides spatial chunks. X/Y
traits are preferred, then common dimension names; `xdim` and `ydim` override
detection. Interval lookups use their bounds, while point lookups are
midpointed. Edge order matches lookup order.

# Cell order

Cells follow the array's dimension order, fastest dimension first. Thus `vec`
of a spatial slice matches cell index order. Chunk numbers use the same rule.

# Face lattice

The face lattice is the native axes ascending, whatever the lookup order and
storage layout: lattice axis 1 is X and axis 2 is Y when the ascending-native
quad is counter-clockwise from outside, and Y then X when the chart reverses
orientation, so every lattice cell is counter-clockwise as the `Structured`
family requires. An X axis spanning `xperiod` wraps (the first and last
columns are neighbours); a lattice side whose vertices coincide (a pole) is
collapsed: it has no neighbour across it, its cells are triangles, and every
cell of the polar row meets at the pole.

# Keywords

  - `xdim`, `ydim`: explicit spatial dimensions.
  - `native_to_unit_sphere`: transform one native `(x, y)` coordinate to the
    unit sphere, returning a `UnitSphericalPoint{Float64}`. With the Proj
    extension loaded, a `Proj.Transformation` is interpreted as native
    coordinates to geographic longitude/latitude and composed with the
    unit-sphere conversion; construct such templates with `always_xy = true`.
    The selected pipeline is cloned for each task.
  - `unit_sphere_to_native`: transform one unit-sphere point back to native
    coordinates; `nothing` disables `cellat`.
  - `transform`, `inverse`: compatibility spellings for the preceding two
    keywords. A legacy two-argument forward closure is adapted at construction.
  - `xperiod`, `ybounds`: native chart limits.
  - `chunks`: `(xchunks, ychunks)` index ranges for the dimension-tuple form.

# Chart resolution

The forward chart is chosen in this order:

 1. An explicit `native_to_unit_sphere` (or `transform`) keyword.
 2. CRS metadata on the X/Y lookups. With Rasters loaded, `Projected` lookups
    supply their `crs` and `Mapped` lookups their `mappedcrs`. Greenwich degrees
    keep the built-in chart. Other geographic coordinates need Proj to convert
    angular units and shift the prime meridian, preserving the datum and latitude
    convention. A projected CRS needs Proj:
    with both loaded, the chart is
    `Proj.Transformation(crs, "EPSG:4326"; always_xy = true)` and `cellat`
    works through its inverse; with Rasters alone, construction throws an
    `ArgumentError` naming the CRS.
 3. Geographic longitude/latitude in degrees.

# Geometry

Cells are geodesic quadrilaterals with counter-clockwise exterior rings.
"""
struct RasterGrid{XD,YD,F,G,T} <: AbstractStructuredGrid
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
    "one native `(x, y)` coordinate to the unit sphere"
    native_to_unit_sphere::F
    "one unit-sphere point back to native `(x, y)`, or `nothing`"
    unit_sphere_to_native::G
    "period of the native X coordinate, or `nothing`"
    xperiod::Union{Nothing,Float64}
    "whether the X dimension varies fastest in cell indices and chunk numbers"
    xfast::Bool
    "whether `(xlo, ylo), (xhi, ylo), (xhi, yhi), (xlo, yhi)` is counter-clockwise from outside, so lattice axis 1 is X; otherwise it is Y"
    ccw::Bool
    "geographic-degrees per-edge tables, or `nothing`"
    tables::T
    "whether `xedges` ascend"
    xasc::Bool
    "whether `yedges` ascend"
    yasc::Bool
    "`xedges` ascending"
    xaxis::Vector{Float64}
    "`yedges` ascending"
    yaxis::Vector{Float64}
    "whether the X edges span `xperiod`, so the lattice wraps across the X axis"
    xwraps::Bool
    "per lattice side 1:4, whether all its vertices are one point (a pole)"
    collapsed::NTuple{4,Bool}
end

const _XNAMES = (:x, :lon, :long, :longitude)
const _YNAMES = (:y, :lat, :latitude)

function RasterGrid(A::DD.AbstractDimArray;
    xdim = nothing, ydim = nothing, chunks = nothing, kwargs...)
    ds = DD.dims(A)
    xd = _resolvedim(ds, xdim, DD.XDim, _XNAMES, "X")
    yd = _resolvedim(ds, ydim, DD.YDim, _YNAMES, "Y")
    xnum, ynum = _dimnums(ds, xd, yd)
    ch = chunks === nothing ? _spatialchunks(parent(A), xnum, ynum) : chunks
    return _rastergrid(xd, yd, ch, xnum < ynum; _with_crs_transform(A, kwargs)...)
end

function RasterGrid(ds::Tuple{Vararg{DD.Dimension}};
    xdim = nothing, ydim = nothing, chunks = nothing, kwargs...)
    xd = _resolvedim(ds, xdim, DD.XDim, _XNAMES, "X")
    yd = _resolvedim(ds, ydim, DD.YDim, _YNAMES, "Y")
    xnum, ynum = _dimnums(ds, xd, yd)
    return _rastergrid(xd, yd, chunks, xnum < ynum; _with_crs_transform(ds, kwargs)...)
end

RasterGrid(ds::DD.Dimension...; kwargs...) = RasterGrid(ds; kwargs...)

function _rastergrid(xd, yd, chunks, xfast::Bool;
    native_to_unit_sphere = _UNSET_RASTER_KEYWORD,
    unit_sphere_to_native = _UNSET_RASTER_KEYWORD,
    transform = _UNSET_RASTER_KEYWORD,
    inverse = _UNSET_RASTER_KEYWORD,
    xperiod = _UNSET_RASTER_KEYWORD,
    ybounds = _UNSET_RASTER_KEYWORD)

    native_to_unit_sphere !== _UNSET_RASTER_KEYWORD &&
        transform !== _UNSET_RASTER_KEYWORD && throw(ArgumentError(
            "pass only `native_to_unit_sphere` or its compatibility alias `transform`"))
    unit_sphere_to_native !== _UNSET_RASTER_KEYWORD &&
        inverse !== _UNSET_RASTER_KEYWORD && throw(ArgumentError(
            "pass only `unit_sphere_to_native` or its compatibility alias `inverse`"))

    raw_forward = native_to_unit_sphere !== _UNSET_RASTER_KEYWORD ?
        native_to_unit_sphere : transform !== _UNSET_RASTER_KEYWORD ?
        transform : US.UnitSphereFromGeographic()
    raw_backward = unit_sphere_to_native !== _UNSET_RASTER_KEYWORD ?
        unit_sphere_to_native : inverse !== _UNSET_RASTER_KEYWORD ?
        inverse : _UNSET_RASTER_KEYWORD
    forward, backward = _prepare_raster_transform_pair(raw_forward, raw_backward)
    default_xperiod, default_ybounds = _native_coordinate_limits(raw_forward)
    xp = xperiod === _UNSET_RASTER_KEYWORD ? default_xperiod : xperiod
    yb = ybounds === _UNSET_RASTER_KEYWORD ? default_ybounds : ybounds

    xe = _edges(xd, nothing)
    ye = _edges(yd, yb)
    forward((xe[1], ye[1])) isa USPoint || throw(ArgumentError(
        "native_to_unit_sphere must return a UnitSphericalPoint{Float64}"))
    nx, ny = length(xe) - 1, length(ye) - 1
    xc, yc = chunks === nothing ? ([1:nx], [1:ny]) :
             (_chunkranges(chunks[1], nx, "X"), _chunkranges(chunks[2], ny, "Y"))
    xperiod = xp === nothing ? nothing : Float64(xp)
    xasc, yasc = xe[end] > xe[1], ye[end] > ye[1]
    ccw = _native_chart_orientation(forward, xe, ye)
    xwraps = xperiod !== nothing && isapprox(abs(xe[end] - xe[1]), xperiod; rtol = 1e-9)
    return RasterGrid(xd, yd, xe, ye, xc, yc, forward, backward, xperiod, xfast, ccw,
        _geographic_edge_tables(raw_forward, xe, ye),
        xasc, yasc, xasc ? xe : reverse(xe), yasc ? ye : reverse(ye), xwraps,
        _collapsedsides(forward, xe, ye, ccw))
end

# Whether each native end of the edge grid (xlo, xhi, ylo, yhi) is one point
# on the sphere, as lattice sides 1:4: with X first the Y ends are sides 1 and
# 3 and the X ends sides 4 and 2; with Y first the roles swap.
function _collapsedsides(native_to_unit_sphere, xe::Vector{Float64}, ye::Vector{Float64}, ccw::Bool)
    f = _task_prepared_raster_transform(native_to_unit_sphere)
    xlo, xhi = extrema(xe)
    ylo, yhi = extrema(ye)
    onepoint(pts) = (p0 = first(pts); all(p -> GridInterface._near(p, p0), pts))
    xlo_end = onepoint(f((xlo, y)) for y in ye)
    xhi_end = onepoint(f((xhi, y)) for y in ye)
    ylo_end = onepoint(f((x, ylo)) for x in xe)
    yhi_end = onepoint(f((x, yhi)) for x in xe)
    return ccw ? (ylo_end, xhi_end, yhi_end, xlo_end) : (xlo_end, yhi_end, xhi_end, ylo_end)
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

Return spatial slice shape in array dimension order. `vec` follows cell indices.
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
regrid output echoes the grid's construction order. A dimension whose
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

# A raster's own chunking, read from the array it is built over. This is a
# space's geometry, fixed when the space is constructed and independent of any
# data a plan is later applied to, so it stays separate from the chunking a lazy
# array reads off the source it regrids.
function _spatialchunks(data, xnum::Int, ynum::Int)
    sz = size(data)
    if declareschunks(data)
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

# Cells, in GridInterface's spelling

ncells(space::RasterGrid) = _nx(space) * _ny(space)

GOCore.best_manifold(::RasterGrid) = GOCore.Spherical(; radius = 1.0)

# A raster's cells are single rings the engine clips with its convex clipper
# (see `intersection_area.jl`): not `GeneralCells`, whose general clipper is a
# different (and, on the sphere, less exact) operation.
GridInterface.cellshape(::RasterGrid) = GridInterface.SimpleCells()

# A raster locates a point through its inverse chart (`cellat`) and names its
# neighbours by lattice arithmetic, so an `AnalyticLocator` can walk it.
GridInterface.location_is_analytic(space::RasterGrid) = space.unit_sphere_to_native !== nothing

# Seams: the lattice wraps across a periodic axis and stops at a collapsed
# side. The lattice is counter-clockwise, so `Structured`'s derivations (ordered
# neighbours, fans, vertex keys) read these directly.
@inline function _periodicaxes(space::RasterGrid)
    return space.ccw ? (space.xwraps, false) : (false, space.xwraps)
end

function facecross(space::RasterGrid, f, i, j, s)
    space.collapsed[s] && return nothing
    ni, nj = _latticesize(space)
    px, py = _periodicaxes(space)
    s == 1 && return py ? (1, i, nj, 3) : nothing
    s == 2 && return px ? (1, 1, j, 4) : nothing
    s == 3 && return py ? (1, i, 1, 1) : nothing
    return px ? (1, ni, j, 2) : nothing
end

# Every cell along a collapsed side meets at its point, counter-clockwise
# about it: westward along a south pole, eastward along a north pole, and the
# rotations of those for a collapsed column.
function facefan(space::RasterGrid, f, a, b)
    ni, nj = _latticesize(space)
    c = space.collapsed
    c[1] && b == 1 && return Int[faceposition(space, 1, i, 1) for i in ni:-1:1]
    c[3] && b == nj + 1 && return Int[faceposition(space, 1, i, nj) for i in 1:ni]
    c[2] && a == ni + 1 && return Int[faceposition(space, 1, ni, j) for j in nj:-1:1]
    c[4] && a == 1 && return Int[faceposition(space, 1, 1, j) for j in 1:nj]
    return GridInterface.Structured._fanwalk(space, f, a, b)
end

GridInterface.maxneighbors(space::RasterGrid, ::GridInterface.Vertex) =
    any(space.collapsed) ? nothing : 8

# A lattice closed along one axis and collapsed at both ends of the other is
# the sphere.
GridInterface.isglobal(space::RasterGrid) =
    (space.xwraps && (space.ccw ? (space.collapsed[1] && space.collapsed[3]) :
                                  (space.collapsed[2] && space.collapsed[4])))

# The lattice corners through a chart hoisted once per cell, as `Structured`'s
# default would read them one `facevertex` at a time.
function GridInterface.cell_corners(space::RasterGrid, k::Integer)
    i, j = _latticeindex(space, cellsubscript(space, Int(k))...)
    return GridInterface.Structured.Corners(_latticecorners(space, i, j))
end
GridInterface.cell_boundary(space::RasterGrid, k::Integer) = GridInterface.cell_corners(space, k)

# The clipper measures the ring `cell_geometry` returns: the four corners with
# polar duplicates dropped, so the cap over the corners bounds it; the relative
# margin absorbs rounding in the cap arithmetic.
function GridInterface.cell_extent(space::RasterGrid, k::Integer)
    i, j = _latticeindex(space, cellsubscript(space, Int(k))...)
    c = _latticecorners(space, i, j)
    cap = SphericalCap(c[1], 0.0)
    for k in 2:4
        cap = Extents.union(cap, SphericalCap(c[k], 0.0))
    end
    cap.radius > Float64(pi) / 2 && return _WHOLE_SPHERE
    return Extents.grow(cap, 5e-5)
end

"""
    cellsubscript(space::RasterGrid, i::Int) -> (ix, iy)
    localindex(space::RasterGrid, ix::Integer, iy::Integer) -> Int

Convert between the raster's cell indices and `(ix, iy)` storage subscripts,
the positions along the X and Y lookups. Index order follows the array's
fastest dimension.
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

function localindex(space::RasterGrid, ix::Integer, iy::Integer)
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

"""
    GridInterface.cell_geometry(space::RasterGrid, i) -> GI.Polygon

Cell `i` as a GeoInterface polygon with one explicitly closed ring of
unit-sphere `(x, y, z)` coordinates, counter-clockwise from outside the sphere,
its segments great-circle arcs.
"""
function GridInterface.cell_geometry(space::RasterGrid, k::Integer)
    i, j = _latticeindex(space, cellsubscript(space, Int(k))...)
    return GI.Polygon([GI.LinearRing(_cellring(_latticecorners(space, i, j)))])
end

# A raster cell is four chart evaluations off two coordinate axes, which costs
# about as much as reading a stored polygon back and allocates no more.
GridInterface.cell_geometry_is_expensive(::RasterGrid) = false

"""
    _latticecorners(space::RasterGrid, i, j[, native_to_unit_sphere]) -> NTuple{4,UnitSphericalPoint}

The lattice corners of cell `(i, j)`, counter-clockwise: `(i, j), (i+1, j),
(i+1, j+1), (i, j+1)`, through edge tables when the chart has them.
"""
@inline _latticecorners(space::RasterGrid, i::Integer, j::Integer) =
    _latticecorners(space, i, j, _task_prepared_raster_transform(space.native_to_unit_sphere))

@inline function _latticecorners(space::RasterGrid, i::Integer, j::Integer, f)
    a, b = Int(i), Int(j)
    return (_rastervertex(space, f, a, b), _rastervertex(space, f, a + 1, b),
        _rastervertex(space, f, a + 1, b + 1), _rastervertex(space, f, a, b + 1))
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

function GridInterface.cell_centroid(space::RasterGrid, i::Integer)
    ix, iy = cellsubscript(space, Int(i))
    return space.native_to_unit_sphere(((space.xedges[ix] + space.xedges[ix+1]) / 2,
        (space.yedges[iy] + space.yedges[iy+1]) / 2))
end

"""
    GridInterface.cellat(space::RasterGrid, p) -> Union{Int,Nothing}

The cell containing `p`, or `nothing` outside the raster. `Structured`'s axis
floor on the face axes, repaired against the geodesic edges: a point just
across a parallel can lie under the bulge of the cell on the other side, and a
point just outside the chart can lie under the bulge of a boundary cell, which
the floor refuses and the nearest boundary cell (or one of its neighbours)
claims. Without an inverse chart the point is located by the lattice walk.
"""
function GridInterface.cellat(space::RasterGrid, p)
    space.unit_sphere_to_native === nothing &&
        return invoke(GridInterface.cellat, Tuple{AbstractStructuredGrid,Any}, space, p)
    m = GOCore.best_manifold(space)
    p = GridInterface.check_point(m, p)
    u, v = facecoords(space, 1, p)
    ij = facelocalindex(space, 1, u, v)
    ij === nothing || return GridInterface.Structured._repair(m, space, 1, ij[1], ij[2], p)
    xs, ys = faceaxes(space, 1)
    k = faceposition(space, 1, _nearestcell(xs, Float64(u)), _nearestcell(ys, Float64(v)))
    GridInterface.cell_contains(space, k, p) && return k
    for q in GridInterface.neighbors(space, k)
        GridInterface.cell_contains(space, q, p) && return q
    end
    return nothing
end

# The cell of an ascending axis nearest coordinate `v`: the floor inside the
# axis, the end cell beyond it. `facecoords` already put `v` on the branch
# nearest the axis, so a coordinate an ulp west of a wrapping axis's origin
# stays at that end.
function _nearestcell(xs::Vector{Float64}, v::Float64)
    v = clamp(v, first(xs), last(xs))
    return min(searchsortedlast(xs, v), length(xs) - 1)
end

# Storage chunks

_nrasterchunks(space::RasterGrid) = length(space.xchunks) * length(space.ychunks)

"""
    chunksubscript(space::RasterGrid, chunk::Int) -> (cx, cy)

Return a chunk's `(cx, cy)` lattice coordinates.
"""
function chunksubscript(space::RasterGrid, chunk::Int)
    1 <= chunk <= _nrasterchunks(space) || throw(BoundsError(space, chunk))
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
    chunknumber(space::RasterGrid, cx::Integer, cy::Integer) -> Int

Return the number of the chunk at `(cx, cy)`.
"""
function chunknumber(space::RasterGrid, cx::Integer, cy::Integer)
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
function _rasterchunkat(space::RasterGrid, i::Integer)
    ix, iy = cellsubscript(space, Int(i))
    return chunknumber(space, _chunkofindex(space.xchunks, ix),
        _chunkofindex(space.ychunks, iy))
end

"""
    RasterBlock(space::RasterGrid, xr, yr)

The cell indices of the storage block `xr × yr`, lazily and strictly
ascending (fast dimension first, as the array stores them), whatever the face
lattice's orientation. What a storage chunk owns.
"""
struct RasterBlock{S<:RasterGrid} <: AbstractVector{Int}
    space::S
    xr::UnitRange{Int}
    yr::UnitRange{Int}
end
Base.size(t::RasterBlock) = (length(t.xr) * length(t.yr),)
Base.IndexStyle(::Type{<:RasterBlock}) = IndexLinear()
@inline function Base.getindex(t::RasterBlock, k::Int)
    @boundscheck checkbounds(t, k)
    fr, sr = _fastslow(t.space, t.xr, t.yr)
    a, b = divrem(k - 1, length(fr))
    return (first(fr) + b) + (first(sr) + a - 1) * _nfast(t.space)
end
Base.first(t::RasterBlock) = localindex(t.space, first(t.xr), first(t.yr))
Base.last(t::RasterBlock) = localindex(t.space, last(t.xr), last(t.yr))
GridInterface._strictly_ascending(::RasterBlock) = true
function Base.in(p::Integer, t::RasterBlock)
    1 <= p <= ncells(t.space) || return false
    ix, iy = cellsubscript(t.space, Int(p))
    return ix in t.xr && iy in t.yr
end

const RasterBlockSubset{S} = GridInterface.Subset{S,RasterBlock{S}}

_rasterowned(space::RasterGrid, chunk::Int) = RasterBlock(space, chunkbox(space, chunk)...)

# A storage block's extent and tree are the lattice block's, O(1).
GridInterface.grid_extent(s::RasterBlockSubset) =
    (t = s.cells; blockextent(t.space, 1, _latticeranges(t.space, t.xr, t.yr)...))
GridInterface.subtree(space::RasterGrid, t::RasterBlock; leaf_capacity = GridInterface.Balanced()) =
    (t.space === space ? _rastercursor(space, t.xr, t.yr) :
     invoke(GridInterface.subtree, Tuple{RasterGrid,AbstractVector{<:Integer}}, space, t; leaf_capacity))
GridInterface.subtree(space::RasterGrid, s::RasterBlockSubset; leaf_capacity = GridInterface.Balanced()) =
    GridInterface.subtree(space, s.cells; leaf_capacity)

# `chunkbox` is `(X, Y)`; array dimensions follow `xfast` order.
chunkranges(space::RasterGrid, chunk::Integer, ::NTuple{2,Int}) =
    space.xfast ? chunkbox(space, Int(chunk)) : reverse(chunkbox(space, Int(chunk)))

# Spherical caps

"""
    _geographic_range_extent(tables, space, ix0, ix1, iy0, iy1)

Return an analytic cap for a rectangle in the built-in geographic-degrees
chart. Every raster cell edge is a declared shortest geodesic. The maximum
single-cell longitude step therefore bounds the poleward bow of every such
edge. A polar latitude-band cap is valid for all longitude widths; rectangles
narrower than 180 degrees may use the tighter mid-meridian cap.

This specialization is deliberately selected only by `GeographicEdgeTables`,
which are installed only for `GeometryOps.UnitSphereFromGeographic`. Arbitrary
curvilinear and projected charts use ConservativeRegridding's perimeter walk.
"""
function _geographic_range_extent(t::GeographicEdgeTables, space::RasterGrid,
        ix0::Int, ix1::Int, iy0::Int, iy1::Int)
    xlo, xhi = minmax(space.xedges[ix0], space.xedges[ix1+1])
    ylo, yhi = minmax(space.yedges[iy0], space.yedges[iy1+1])
    polar = _geographic_polar_cap(t, ylo, yhi)
    abs(xhi - xlo) < 180.0 || return polar
    meridian = _geographic_meridian_cap(t, xlo, xhi, ylo, yhi)
    return meridian.radius < polar.radius ? meridian : polar
end

"""
    _geographic_polar_cap(tables, ylo, yhi) -> SphericalCap

Return the tighter north- or south-pole cap containing the bowed geographic
latitude band. This remains valid for any longitude span, including full
longitude and antimeridian-crossing ranges.
"""
function _geographic_polar_cap(t::GeographicEdgeTables, ylo, yhi)
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

"""
    _geographic_meridian_cap(tables, xlo, xhi, ylo, yhi) -> SphericalCap

Return a cap around a geographic box narrower than 180 degrees longitude. The
box is first expanded to the proven poleward latitude limits of its geodesic
cell edges. Distance from the mid-meridian centre is maximized at a box corner,
so the resulting cap covers every point of every leaf boundary directly; it
does not depend on convexity of a wider-than-hemisphere cap.
"""
function _geographic_meridian_cap(t::GeographicEdgeTables, xlo, xhi, ylo, yhi)
    a, b = _bowedband(Float64(ylo), Float64(yhi), t.maxdx)
    geographic_to_sphere = US.UnitSphereFromGeographic()
    centre = geographic_to_sphere(((xlo + xhi) / 2, (a + b) / 2))
    r = 0.0
    for x in (xlo, xhi), y in (a, b)
        r = max(r, US.spherical_distance(centre, geographic_to_sphere((x, y))))
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

# Determine chart handedness once from probe-cell Newell normals.
function _native_chart_orientation(native_to_unit_sphere,
    xedges::Vector{Float64}, yedges::Vector{Float64})
    native_to_unit_sphere = _task_prepared_raster_transform(native_to_unit_sphere)
    nx, ny = length(xedges) - 1, length(yedges) - 1
    best = 0.0
    for iy in _probes(ny), ix in _probes(nx)
        xlo, xhi = minmax(xedges[ix], xedges[ix+1])
        ylo, yhi = minmax(yedges[iy], yedges[iy+1])
        s = _quadsign(native_to_unit_sphere, xlo, xhi, ylo, yhi)
        abs(s) > abs(best) && (best = s)
    end
    best == 0.0 && throw(ArgumentError(
        "RasterGrid: the chart collapses every probe cell to a degenerate ring, " *
        "so the winding of a cell cannot be determined"))
    return best > 0.0
end

_probes(n::Int) = n <= 3 ? (1:n) : (1, cld(n, 2), n)

function _quadsign(t, xlo, xhi, ylo, yhi)
    c = (t((xlo, ylo)), t((xhi, ylo)), t((xhi, yhi)), t((xlo, yhi)))
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

# The face: the native axes ascending, counter-clockwise

# The raster is a one-face `Structured` grid. Positions are the cell indices
# in storage order; the lattice is the native axes ascending, X first when
# that quad is counter-clockwise from outside and Y first otherwise, so every
# lattice cell winds as the family requires. The maps below convert between
# the two, and the quadtree cursor reports positions through them.

nfaces(::RasterGrid) = 1
facesize(space::RasterGrid, f) = _latticesize(space)
@inline faceposition(space::RasterGrid, f, i, j) = localindex(space, _storagecell(space, Int(i), Int(j))...)
@inline facecell(space::RasterGrid, k) = (1, _latticeindex(space, cellsubscript(space, Int(k))...)...)
@inline facevertex(space::RasterGrid, f, a, b) =
    _rastervertex(space, _task_prepared_raster_transform(space.native_to_unit_sphere),
        Int(a), Int(b))
# Positions ascend with lattice `i` fastest only when lattice axis 1 is the
# storage-fast axis and neither lookup descends.
GridInterface.Structured.facecolumnmajor(space::RasterGrid, f) =
    space.ccw == space.xfast && space.xasc && space.yasc

_latticesize(space::RasterGrid) = _latticeorder(space, _nx(space), _ny(space))

# (X, Y) pairs to lattice-axis pairs and back: the swap is its own inverse.
@inline _latticeorder(space::RasterGrid, x, y) = space.ccw ? (x, y) : (y, x)

# Lattice cell `(i, j)` as storage subscripts `(ix, iy)`, and back.
@inline function _storagecell(space::RasterGrid, i::Int, j::Int)
    a, b = _latticeorder(space, i, j)
    return (space.xasc ? a : _nx(space) + 1 - a, space.yasc ? b : _ny(space) + 1 - b)
end
@inline function _latticeindex(space::RasterGrid, ix::Int, iy::Int)
    return _latticeorder(space, space.xasc ? ix : _nx(space) + 1 - ix,
        space.yasc ? iy : _ny(space) + 1 - iy)
end

# Lattice vertex `(a, b)` as storage edge indices `(ex, ey)`.
@inline function _storagevertex(space::RasterGrid, a::Int, b::Int)
    p, q = _latticeorder(space, a, b)
    return (space.xasc ? p : _nx(space) + 2 - p, space.yasc ? q : _ny(space) + 2 - q)
end

# Lattice cell ranges as storage ranges `(xr, yr)`, and back.
function _storageranges(space::RasterGrid, ir::UnitRange{Int}, jr::UnitRange{Int})
    ar, br = _latticeorder(space, ir, jr)
    nx, ny = _nx(space), _ny(space)
    return (space.xasc ? ar : (nx + 1 - last(ar)):(nx + 1 - first(ar)),
        space.yasc ? br : (ny + 1 - last(br)):(ny + 1 - first(br)))
end
function _latticeranges(space::RasterGrid, xr::UnitRange{Int}, yr::UnitRange{Int})
    nx, ny = _nx(space), _ny(space)
    return _latticeorder(space, space.xasc ? xr : (nx + 1 - last(xr)):(nx + 1 - first(xr)),
        space.yasc ? yr : (ny + 1 - last(yr)):(ny + 1 - first(yr)))
end

# Lattice vertex `(a, b)` through the given chart callable.
@inline function _rastervertex(space::RasterGrid, native_to_unit_sphere, a::Int, b::Int)
    ex, ey = _storagevertex(space, a, b)
    space.tables === nothing &&
        return native_to_unit_sphere((space.xedges[ex], space.yedges[ey]))
    return _tablecorner(space.tables, ex, ey)
end

# The number of cells along the dimension that varies fastest in cell indices.
_nfast(space::RasterGrid) = space.xfast ? _nx(space) : _ny(space)

# Storage-order pairs from X/Y pairs and back: the swap is its own inverse.
@inline _fastslow(space::RasterGrid, x, y) = space.xfast ? (x, y) : (y, x)

# A cursor node's positions as runs along the storage-fast dimension, exact
# for every lattice orientation.
function GridInterface.position_span(q::TopDownQuadtreeCursor{<:RasterGrid,false})
    space = q.grid
    fr, sr = _fastslow(space, _storageranges(space, q.irange, q.jrange)...)
    return GridInterface.LatticeBlock(first(fr) + (first(sr) - 1) * _nfast(space),
        length(fr), _nfast(space), length(sr))
end

# The structured grid's extent seam (`Structured.blockextent`): every cursor
# node, leaf cell, scattered-cell cap and storage chunk extent comes through
# here. Only the exact built-in geographic chart takes the analytic
# specialization; every other chart delegates to GridInterface's conservative
# perimeter walk over the face vertices, on a vertex source that has hoisted
# the chart once for this task (one TLS lookup per block, not per vertex).
function blockextent(space::RasterGrid, f, irange::UnitRange{Int},
        jrange::UnitRange{Int})
    tables = space.tables
    if tables === nothing
        return GridInterface.Structured._blockextent(GOCore.best_manifold(space),
            _RasterPerimeter(space), 1, irange, jrange)
    end
    xr, yr = _storageranges(space, irange, jrange)
    return _geographic_range_extent(tables, space, first(xr), last(xr), first(yr), last(yr))
end

# The vertex source the perimeter walk reads for a curvilinear chart: the
# raster's lattice through a task-prepared chart callable. Ephemeral, never
# stored, never a grid of its own.
struct _RasterPerimeter{S<:RasterGrid,F} <: AbstractStructuredGrid
    space::S
    native_to_unit_sphere::F
end

_RasterPerimeter(space::RasterGrid) =
    _RasterPerimeter(space, _task_prepared_raster_transform(space.native_to_unit_sphere))

GOCore.best_manifold(p::_RasterPerimeter) = GOCore.best_manifold(p.space)
nfaces(::_RasterPerimeter) = 1
facesize(p::_RasterPerimeter, f) = _latticesize(p.space)
@inline facevertex(p::_RasterPerimeter, f, i, j) =
    _rastervertex(p.space, p.native_to_unit_sphere, Int(i), Int(j))

# Spatial trees

const _RASTER_CURSOR_LEAFSIZE = (4, 4)

# The cursor over the storage block `xr × yr`, as the lattice block it is.
function _rastercursor(space::RasterGrid, xr::UnitRange{Int}, yr::UnitRange{Int})
    ranges = _latticeranges(space, xr, yr)
    return TopDownQuadtreeCursor(space, 1, ranges...; leafsize = _RASTER_CURSOR_LEAFSIZE)
end

_rastercursor(space::RasterGrid) = _rastercursor(space, 1:_nx(space), 1:_ny(space))

# The raster's one tree is its quadtree cursor at every intent and leaf size.
GridInterface.treeify(space::RasterGrid; leaf_capacity = GridInterface.Balanced()) =
    _rastercursor(space)
GridInterface.treeify(space::RasterGrid, ::GridInterface.LeafHint) = _rastercursor(space)
GridInterface.treeify(space::RasterGrid, ::Integer) = _rastercursor(space)

"""
    GridInterface.subtree(space::RasterGrid, inds::AbstractVector{<:Integer}; leaf_capacity = Balanced())

The quadtree cursor over the lattice rectangle `inds` enumerates, when it is
one, otherwise GridInterface's packed R-tree over those cells. Leaves emit the
raster's cell indices. The raster's cursor has the one leaf size its `treeify`
has (`_RASTER_CURSOR_LEAFSIZE`), so the keyword reaches only the R-tree.

Both branches beat GridInterface's generic restriction of the raster's tree:
a rectangle not aligned with the bisection would root above itself with
looser extents, and a scattered set would be pruned on rows only (a block's
[`GridInterface.position_span`](@ref) is a run of whole rows), visiting every
block of a member's row.
"""
function GridInterface.subtree(space::RasterGrid, inds::AbstractVector{<:Integer};
    leaf_capacity = GridInterface.Balanced())
    rect = _indexrect(space, inds)
    if rect !== nothing
        ix0, ix1, iy0, iy1 = rect
        return _rastercursor(space, ix0:ix1, iy0:iy1)
    end
    return GridInterface.fallback_subtree(GOCore.best_manifold(space), space, inds,
        GridInterface._capacity(space, leaf_capacity))
end

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

# Cell chart: the face axes

# The face axes are the cell edges ascending, in lattice order; the chart
# coordinates are the native ones on the branch nearest the raster's span
# (`_rastercoords`), so `Structured`'s axis floor locates a point's cell.
hasfaceaxes(space::RasterGrid, f) = space.unit_sphere_to_native !== nothing
faceaxes(space::RasterGrid, f) = _latticeorder(space, space.xaxis, space.yaxis)
function facecoords(space::RasterGrid, f, p)
    xy = _rastercoords(space, p)
    xy === nothing && return nothing
    return _latticeorder(space, xy[1], xy[2])
end
# The chart's period, whatever the raster's span; the chart wraps only where
# the axis spans it (`_chartperiod`).
faceperiod(space::RasterGrid, f) = _latticeorder(space, space.xperiod, nothing)

_centres(e::Vector{Float64}) = [(e[k] + e[k+1]) / 2 for k in 1:(length(e)-1)]

# `p` in native coordinates, or `nothing` without an inverse. Periodic X
# returns the representative nearest the edge span.
function _rastercoords(space::RasterGrid, p)
    space.unit_sphere_to_native === nothing && return nothing
    x, y = space.unit_sphere_to_native(p)
    return (_onbranch(space.xedges, Float64(x), space.xperiod), Float64(y))
end

# nearest periodic representative: folding to [lo, lo+p) throws points just west of lo a period east
function _onbranch(edges::Vector{Float64}, v::Float64, period)
    period === nothing && return v
    p = period::Float64
    mid = (edges[1] + edges[end]) / 2
    return v - p * round((v - mid) / p)
end

# The largest edge step on each axis in radians of arc.
_rasterspacing(space::RasterGrid) =
    _spherical_step_bounds_radians(space.native_to_unit_sphere,
        _maxstep(space.xedges), _maxstep(space.yedges))

_maxstep(e::Vector{Float64}) =
    maximum(abs(e[k+1] - e[k]) for k in 1:(length(e)-1))

_facespacing(space::RasterGrid) = _latticeorder(space, _rasterspacing(space)...)

# Storage chunks as GridInterface's chunks

"""
    RasterChunks(space::RasterGrid)

The raster's storage chunks as GridInterface's `chunks(space)`: chunk `c` is the
`Subset` over the chunk's storage block (a [`RasterBlock`](@ref), held
lazily), with the chunk number as its root. Building one, its extent (the
analytic range cap, through `blockextent`) and its tree (the chunk's own
cursor) are all O(1). Lazy: nothing per chunk or per cell is stored.
"""
struct RasterChunks{S<:RasterGrid} <: AbstractVector{GridInterface.Subset{S}}
    space::S
end

Base.size(c::RasterChunks) = (_nrasterchunks(c.space),)
Base.IndexStyle(::Type{<:RasterChunks}) = IndexLinear()
function Base.getindex(c::RasterChunks, k::Int)
    @boundscheck checkbounds(c, k)
    return GridInterface.Subset(c.space, _rasterowned(c.space, k); root = k)
end

GridInterface.chunks(space::RasterGrid) = RasterChunks(space)
GridInterface.chunkat(c::RasterChunks, i::Integer) = _rasterchunkat(c.space, Int(i))

# A block of the raster's own lattice (a cursor node's `block_indices`) is the
# raster's cursor over it, with the raster's leaf size.
function GridInterface.subtree(space::RasterGrid, t::TileIndices;
    leaf_capacity = GridInterface.Balanced())
    t.grid === space || return invoke(GridInterface.subtree,
        Tuple{RasterGrid,AbstractVector{<:Integer}}, space, t; leaf_capacity)
    return TopDownQuadtreeCursor(space, 1, t.ir, t.jr; leafsize = _RASTER_CURSOR_LEAFSIZE)
end
