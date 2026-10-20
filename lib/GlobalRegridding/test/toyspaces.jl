# Analytic test grids and methods.
#
# A regrid side is a GridInterface grid. `ToyLonLatSpace` is a one-face
# `Structured` grid whose vertices, chart, sample sites and point location are
# the toy's own closed-form lon/lat formulas; the wrappers below are
# `AbstractWrappedGrid`s that count or relabel. The short names after the
# imports are test-local spellings of the GridInterface verbs and
# GlobalRegridding's chunk helpers that the suites read a grid through.

using GlobalRegridding
import GlobalRegridding as GR
import GlobalRegridding: AbstractRegriddingMethod, WeightCOO, ncells,
    buildweights!, addweight!, adddenom!

import GeometryOps as GO
import GeometryOpsCore as GOCore
import GeoInterface as GI
import DimensionalData as DD
import GridInterface as GrI
import GridInterface.Structured as GrS

const US = GO.UnitSpherical
const USPoint = GO.UnitSpherical.UnitSphericalPoint{Float64}
# Generic constructor and concrete annotation aliases.
const SphericalCap = GO.UnitSpherical.SphericalCap
const Cap = GO.UnitSpherical.SphericalCap{Float64}

# The verbs a test reads a regrid side through, in the names the suites use.
# Grid-shaped questions are GridInterface's; chunk questions are
# GlobalRegridding's private helpers over `chunks(grid)`.
getcell(space, i::Integer) = GrI.cell_geometry(space, Int(i))
cellcentroid(space, i::Integer) = GrI.cell_centroid(space, Int(i))
cellcorners(space, i::Integer) = GrI.cell_corners(space, Int(i))
cellat(space, p) = GrI.cellat(space, p)
manifold(space) = GOCore.best_manifold(space)
cellneighbors(space, i::Integer) = GrI.neighbors(space, Int(i))
cellcap(space, i::Integer) = convert(Cap, GrI.cell_extent(space, Int(i)))
hasanalyticlocation(space) = GrI.location_is_analytic(space)
hascellchart(space) = GR._hascellchart(space)
hasdualcells(space) = GR._hasdualcells(space)
nchunks(space) = GR._nchunks(space)
ownedindices(space, chunk::Integer) = GR._ownedindices(space, Int(chunk))
chunkextents(space) = GR._chunkextents(space)
chunkextent(space, chunk::Integer) = GR._chunkextent(space, Int(chunk))
chunkat(space, i) = GR._chunkat(space, i)
celltree(space) = GR._celltree(space)
celltree(space, chunk::Integer) = GR._subtree(space, ownedindices(space, chunk))
celltree(space, inds::AbstractVector) = GR._subtree(space, inds)
subtree(space, inds) = GR._subtree(space, inds)
chartaxes(space) = GR._chartaxes(space)
chartcoords(space, p) = GR._chartcoords(space, p)
chartlocalindex(space, ix, iy) = GR._chartlocalindex(space, Int(ix), Int(iy))
chartperiod(space) = GR._chartperiod(space)
chartspacing(space) = GR._chartspacing(space)

# Unit-sphere helpers

"""
    toy_point(lon, lat) -> UnitSphericalPoint

Convert longitude and latitude in degrees to the unit sphere.
"""
function toy_point(lon::Real, lat::Real)
    coslat = cosd(Float64(lat))
    return USPoint(coslat * cosd(Float64(lon)), coslat * sind(Float64(lon)),
        sind(Float64(lat)))
end

"""
    toy_lonlat(p) -> (lon, lat)

Convert a unit-sphere point to longitude and latitude in degrees.
"""
toy_lonlat(p) = (atand(p[2], p[1]), asind(clamp(p[3], -1.0, 1.0)))

const TOY_FULL_SPHERE = SphericalCap(USPoint(0.0, 0.0, 1.0), nextfloat(Float64(pi)))

"""
    toy_cap(points) -> SphericalCap

Return a cap around `points`, or the full sphere when the radius exceeds `π/2`.
"""
function toy_cap(points)
    sx = sy = sz = 0.0
    n = 0
    for p in points
        sx += p[1]; sy += p[2]; sz += p[3]
        n += 1
    end
    n == 0 && return TOY_FULL_SPHERE
    norm = sqrt(sx^2 + sy^2 + sz^2)
    norm <= eps(Float64) && return TOY_FULL_SPHERE
    centre = USPoint(sx / norm, sy / norm, sz / norm)
    radius = 0.0
    for p in points
        radius = max(radius, US.spherical_distance(centre, USPoint(p[1], p[2], p[3])))
    end
    radius > Float64(pi) / 2 && return TOY_FULL_SPHERE
    # Nudged outward past dot-product rounding noise, so containment stays closed.
    return SphericalCap(centre, nextfloat(radius * 1.0001 + 1e-12))
end

"""
    toy_bandcap(lat0, lat1, dlon) -> SphericalCap

Return a polar cap covering the bowed latitude band `[lat0, lat1]` for cells of
width `dlon`. This remains valid for full-longitude bands.
"""
function toy_bandcap(lat0::Real, lat1::Real, dlon::Real)
    h = min(abs(Float64(dlon)), 360.0) / 2
    c = h >= 90 ? 0.0 : cosd(h)
    a = lat0 < 0 ? (c <= 0 ? -90.0 : max(-90.0, atand(tand(Float64(lat0)) / c))) :
        Float64(lat0)
    b = lat1 > 0 ? (c <= 0 ? 90.0 : min(90.0, atand(tand(Float64(lat1)) / c))) :
        Float64(lat1)
    centre = USPoint(0.0, 0.0, 1.0)
    radius = deg2rad(90.0 - a)
    if deg2rad(90.0 + b) < radius
        centre = USPoint(0.0, 0.0, -1.0)
        radius = deg2rad(90.0 + b)
    end
    radius = nextfloat(radius * 1.0001 + 1e-12)
    radius >= Float64(pi) && return TOY_FULL_SPHERE
    return SphericalCap(centre, radius)
end

# Lon/lat test space

"""
    ToyLonLatSpace(nlon, nlat; lon = (-180, 180), lat = (-90, 90), chunks = (nlon, nlat))

A regular `nlon × nlat` geodesic lon/lat grid on the unit sphere, as a one-face
`GridInterface.Structured` grid. Longitude varies fastest. `chunks` gives
`(longitude, latitude)` tile sizes, and `chunks(grid)` is the tile partition.

The vertices, cell polygons, centroids, chart and point location are the toy's
own closed-form formulas; the seams (a periodic longitude, collapsed poles) are
a `RegularGrid`'s over the same lattice, and everything else — rings, extents,
the quadtree, tiles — is derived by `GridInterface.Structured`.
"""
struct ToyLonLatSpace{R} <: GrS.AbstractStructuredGrid
    nlon::Int
    nlat::Int
    lon0::Float64
    lon1::Float64
    lat0::Float64
    lat1::Float64
    chunklon::Int
    chunklat::Int
    reg::R
end

function ToyLonLatSpace(nlon::Integer, nlat::Integer;
    lon::Tuple{Real,Real} = (-180.0, 180.0),
    lat::Tuple{Real,Real} = (-90.0, 90.0),
    chunks::Tuple{Integer,Integer} = (nlon, nlat))
    nlon >= 1 && nlat >= 1 || throw(ArgumentError("a toy space needs at least one cell"))
    lon[1] < lon[2] || throw(ArgumentError("longitude span must be increasing"))
    lat[1] < lat[2] || throw(ArgumentError("latitude span must be increasing"))
    -90 <= lat[1] && lat[2] <= 90 || throw(ArgumentError("latitude span must lie in [-90, 90]"))
    chunks[1] >= 1 && chunks[2] >= 1 || throw(ArgumentError("chunk sizes must be positive"))
    lons = range(Float64(lon[1]), Float64(lon[2]); length = Int(nlon) + 1)
    lats = range(Float64(lat[1]), Float64(lat[2]); length = Int(nlat) + 1)
    reg = GrS.RegularGrid(GOCore.Spherical(; radius = 1.0), lons, lats;
        periodic = (lon[2] - lon[1] >= 360, false),
        collapsed = (lat[1] <= -90, lat[2] >= 90))
    return ToyLonLatSpace(Int(nlon), Int(nlat),
        Float64(lon[1]), Float64(lon[2]), Float64(lat[1]), Float64(lat[2]),
        min(Int(chunks[1]), Int(nlon)), min(Int(chunks[2]), Int(nlat)), reg)
end

Base.show(io::IO, space::ToyLonLatSpace) =
    print(io, "ToyLonLatSpace(", space.nlon, "×", space.nlat,
        ", chunks=", nchunklon(space), "×", nchunklat(space), ")")

dlon(space::ToyLonLatSpace) = (space.lon1 - space.lon0) / space.nlon
dlat(space::ToyLonLatSpace) = (space.lat1 - space.lat0) / space.nlat

nchunklon(space::ToyLonLatSpace) = cld(space.nlon, space.chunklon)
nchunklat(space::ToyLonLatSpace) = cld(space.nlat, space.chunklat)

"""
    cellsubscript(space::ToyLonLatSpace, i::Int) -> (ix, iy)
    localindex(space::ToyLonLatSpace, ix::Int, iy::Int) -> Int

Convert between the space's local cell indices and lattice coordinates.
"""
function cellsubscript(space::ToyLonLatSpace, i::Int)
    1 <= i <= ncells(space) || throw(BoundsError(space, i))
    iy, ix = fldmod1(i, space.nlon)
    return (ix, iy)
end

function localindex(space::ToyLonLatSpace, ix::Integer, iy::Integer)
    1 <= ix <= space.nlon && 1 <= iy <= space.nlat ||
        throw(BoundsError(space, (ix, iy)))
    return Int(ix) + (Int(iy) - 1) * space.nlon
end

"""
    cellbounds(space::ToyLonLatSpace, i::Int) -> (lon_lo, lon_hi, lat_lo, lat_hi)

Return the cell's lon/lat bounds in degrees.
"""
function cellbounds(space::ToyLonLatSpace, i::Int)
    ix, iy = cellsubscript(space, i)
    return (space.lon0 + (ix - 1) * dlon(space), space.lon0 + ix * dlon(space),
        space.lat0 + (iy - 1) * dlat(space), space.lat0 + iy * dlat(space))
end

"""
    graticule_area(space::ToyLonLatSpace, i::Int) -> Float64

Return the parallel-bounded area `Δλ (sin φ_hi − sin φ_lo)` in steradians.
This differs from the geodesic polygon area used by [`getcell`](@ref).
"""
function graticule_area(space::ToyLonLatSpace, i::Int)
    _, _, lat_lo, lat_hi = cellbounds(space, i)
    return deg2rad(dlon(space)) * (sind(lat_hi) - sind(lat_lo))
end

# The Structured face: one face, the toy's lattice, its vertices by formula.

GOCore.best_manifold(::ToyLonLatSpace) = GOCore.Spherical(; radius = 1.0)
GrS.nfaces(::ToyLonLatSpace) = 1
GrS.facesize(g::ToyLonLatSpace, f) = (g.nlon, g.nlat)
GrS.facevertex(g::ToyLonLatSpace, f, a, b) =
    toy_point(g.lon0 + (a - 1) * dlon(g), g.lat0 + (b - 1) * dlat(g))
GrS.faceposition(g::ToyLonLatSpace, f, i, j) = localindex(g, i, j)
GrS.facecell(g::ToyLonLatSpace, k) = (1, cellsubscript(g, Int(k))...)
GrS.facecross(g::ToyLonLatSpace, f, i, j, s) = GrS.facecross(g.reg, f, i, j, s)
GrS.facefan(g::ToyLonLatSpace, f, a, b) = GrS.facefan(g.reg, f, a, b)
GrI.maxneighbors(g::ToyLonLatSpace, c::GrI.Vertex) = GrI.maxneighbors(g.reg, c)

# The cell chart: face axes at the lattice's edges, in native degrees, with
# longitude read on the axes' branch and a period only where the span is global.
GrS.hasfaceaxes(::ToyLonLatSpace, f) = true
GrS.faceaxes(g::ToyLonLatSpace, f) =
    ([g.lon0 + (a - 1) * dlon(g) for a in 1:(g.nlon + 1)],
        [g.lat0 + (b - 1) * dlat(g) for b in 1:(g.nlat + 1)])
function GrS.facecoords(g::ToyLonLatSpace, f, p)
    lon, lat = toy_lonlat(p)
    return (g.lon0 + mod(lon - g.lon0, 360.0), lat)
end
GrS.faceperiod(g::ToyLonLatSpace, f) = (360.0, nothing)

# The chart's spacing: Δλ bounds adjacent distance along a parallel.
GR._facespacing(g::ToyLonLatSpace) = (deg2rad(dlon(g)), deg2rad(dlat(g)))

# The toy's cells are synthesized polygons an area method keeps per tile.
GrI.cell_geometry_is_expensive(::ToyLonLatSpace) = true

# A full-longitude block's perimeter cap is the whole sphere; the bowed
# latitude band bounds it instead, so row chunks keep bounded extents.
function GrS.blockextent(g::ToyLonLatSpace, f, ir::UnitRange{Int}, jr::UnitRange{Int})
    cap = invoke(GrS.blockextent, Tuple{Any,Any,UnitRange{Int},UnitRange{Int}}, g, f, ir, jr)
    cap.radius < Float64(pi) && return cap
    return toy_bandcap(g.lat0 + (first(jr) - 1) * dlat(g), g.lat0 + last(jr) * dlat(g),
        dlon(g))
end

# The tiles are the toy's `(longitude, latitude)` chunk sizes.
GrI.chunks(g::ToyLonLatSpace) = GrI.chunks(g, (g.chunklon, g.chunklat))

"""
    cellcorners(space::ToyLonLatSpace, i::Int) -> NTuple{4,UnitSphericalPoint}

Return counter-clockwise graticule corners. Polar cells repeat one corner.
"""
function GrI.cell_corners(space::ToyLonLatSpace, i::Integer)
    lon_lo, lon_hi, lat_lo, lat_hi = cellbounds(space, Int(i))
    return (toy_point(lon_lo, lat_lo), toy_point(lon_hi, lat_lo),
        toy_point(lon_hi, lat_hi), toy_point(lon_lo, lat_hi))
end

GrI.cell_geometry(space::ToyLonLatSpace, i::Integer) =
    GI.Polygon([GI.LinearRing(_closed_ring(cellcorners(space, Int(i))))])

# Drop repeated polar corners.
function _closed_ring(corners)
    ring = USPoint[]
    for p in corners
        (isempty(ring) || ring[end] != p) && push!(ring, p)
    end
    length(ring) > 1 && ring[end] == ring[1] && pop!(ring)
    push!(ring, ring[1])
    return ring
end

function GrI.cell_centroid(space::ToyLonLatSpace, i::Integer)
    lon_lo, lon_hi, lat_lo, lat_hi = cellbounds(space, Int(i))
    return toy_point((lon_lo + lon_hi) / 2, (lat_lo + lat_hi) / 2)
end

"""
    cellat(space::ToyLonLatSpace, p) -> Union{Int,Nothing}

Return the cell containing `p`, wrapping global longitude. Return `nothing`
outside partial coverage. Shared boundaries select the east/north cell.
"""
function GrI.cellat(space::ToyLonLatSpace, p)
    lon, lat = toy_lonlat(p)
    wrapped = _wrap_lon(lon, space.lon0, space.lon1)
    wrapped === nothing && return nothing
    space.lat0 <= lat <= space.lat1 || return nothing
    ix = clamp(floor(Int, (wrapped - space.lon0) / dlon(space)) + 1, 1, space.nlon)
    iy = clamp(floor(Int, (lat - space.lat0) / dlat(space)) + 1, 1, space.nlat)
    return localindex(space, ix, iy)
end

function _wrap_lon(lon::Float64, lo::Float64, hi::Float64)
    hi - lo >= 360 && return lo + mod(lon - lo, 360.0)
    wrapped = lo + mod(lon - lo, 360.0)
    return wrapped <= hi ? wrapped : nothing
end

function chunksubscript(space::ToyLonLatSpace, chunk::Int)
    1 <= chunk <= nchunks(space) || throw(BoundsError(space, chunk))
    cy, cx = fldmod1(chunk, nchunklon(space))
    return (cx, cy)
end

# Tree-build counting wrapper

"""
    CountingSpace(space)

Wrap a grid and count restricted-tree builds at the `GridInterface.subtree` /
`treeify` seam GlobalRegridding's `_subtree` reads, delegating every other
grid verb unchanged; `cs.chunkreads[]` counts `chunks(cs)` resolutions. Read
`cs.builds[]`; block builders may run on
`Threads.@spawn` tasks, hence the atomic counter.
"""
struct CountingSpace{S} <: GrI.AbstractWrappedGrid
    space::S
    builds::Threads.Atomic{Int}
    chunkreads::Threads.Atomic{Int}
end

CountingSpace(space) = CountingSpace(space, Threads.Atomic{Int}(0), Threads.Atomic{Int}(0))
Base.parent(cs::CountingSpace) = cs.space
GrI.chunks(cs::CountingSpace) = (Threads.atomic_add!(cs.chunkreads, 1); GrI.chunks(cs.space))

GrI.subtree(cs::CountingSpace, inds::AbstractVector{<:Integer}; kw...) =
    (Threads.atomic_add!(cs.builds, 1); GrI.subtree(cs.space, inds; kw...))
GrI.subtree(cs::CountingSpace, s::GrI.Subset; kw...) =
    (Threads.atomic_add!(cs.builds, 1); GrI.subtree(cs.space, s; kw...))
GrI.treeify(cs::CountingSpace; leaf_capacity = GrI.Balanced()) =
    (Threads.atomic_add!(cs.builds, 1); GrI.treeify(cs.space; leaf_capacity))
GrI.treeify(cs::CountingSpace, hint::GrI.LeafHint) =
    (Threads.atomic_add!(cs.builds, 1); GrI.treeify(cs.space, hint))
GrI.treeify(cs::CountingSpace, n::Integer) =
    (Threads.atomic_add!(cs.builds, 1); GrI.treeify(cs.space, n))

# One-axis destination wrapper

"""
    ToyCellAxisSpace(space)

Wrap a grid and label results over it with one axis — `Dim{:toycell}` over the
cells' local indices — where the wrapped grid labels with its own lattice.
This is the destination shape a cell collection has: the cells are the whole of
it, so a result needs no reshape to sit on that axis.
"""
struct ToyCellAxisSpace{S} <: GrI.AbstractWrappedGrid
    space::S
end

Base.parent(cs::ToyCellAxisSpace) = cs.space
Base.show(io::IO, cs::ToyCellAxisSpace) =
    print(io, "ToyCellAxisSpace(", cs.space, ")")

GR.destinationdims(cs::ToyCellAxisSpace, ::DD.Lookups.Sampling) =
    (DD.Dim{:toycell}(1:ncells(cs.space)),)

# Geometry-free test method

"""
    ToyDiagonalMethod(; scale = 1.0, withdenom = true)

Build diagonal weights of `scale` for shared cell indices. `withdenom = false`
omits denominators. This isolates executor behavior from geometry.
"""
struct ToyDiagonalMethod <: AbstractRegriddingMethod
    scale::Float64
    withdenom::Bool
end

ToyDiagonalMethod(; scale::Real = 1.0, withdenom::Bool = true) =
    ToyDiagonalMethod(Float64(scale), withdenom)

"""
    countbuild!(method)

Increment `method.builds` under a lock for concurrent test builds.
"""
const TOY_COUNT_LOCK = ReentrantLock()

countbuild!(method) = @lock TOY_COUNT_LOCK method.builds += 1

function buildweights!(coo::WeightCOO, method::ToyDiagonalMethod,
    dst, dst_inds, src, src_inds)
    chunklocal_of = Dict{Int,Int}(p => k for (k, p) in enumerate(src_inds))
    for (j, p) in enumerate(dst_inds)
        k = get(chunklocal_of, p, 0)
        k == 0 && continue
        addweight!(coo, j, k, method.scale)
        method.withdenom && adddenom!(coo, j, method.scale)
    end
    return coo
end

"""
    WaveFailMethod(bad, delay)

Fail the build of the source chunk whose first cell index is `bad`, and make
every other build take `delay` seconds before recording itself in `finished`.

This exists to pin one thing: a wave that loses a task must still wait for the
rest. Put the failing chunk first and the survivors sleep, so a `_fillwave!` that
raises on the first `fetch` and walks away leaves `finished` short.
"""
struct WaveFailMethod <: AbstractRegriddingMethod
    bad::Int
    delay::Float64
    finished::Threads.Atomic{Int}
end

WaveFailMethod(bad::Integer, delay::Real) =
    WaveFailMethod(Int(bad), Float64(delay), Threads.Atomic{Int}(0))

function buildweights!(coo::WeightCOO, method::WaveFailMethod,
    dst, dst_inds, src, src_inds)
    Int(first(src_inds)) == method.bad &&
        error("WaveFailMethod: the chunk at index $(method.bad) fails by design")
    sleep(method.delay)
    chunklocal_of = Dict{Int,Int}(p => k for (k, p) in enumerate(src_inds))
    for (j, p) in enumerate(dst_inds)
        k = get(chunklocal_of, p, 0)
        k == 0 && continue
        addweight!(coo, j, k, 1.0)
        adddenom!(coo, j, 1.0)
    end
    Threads.atomic_add!(method.finished, 1)
    return coo
end
