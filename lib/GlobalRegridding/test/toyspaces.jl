# Analytic test spaces and methods.

using GlobalRegridding
import GlobalRegridding as GR
import GlobalRegridding: RegridSpace, AbstractRegriddingMethod, WeightCOO,
    celltree, chunktree, ncells, getcell, nchunks, cellindices,
    cellcentroid, cellat, hascellchart, manifold,
    build_weights!, addweight!, adddenom!

import GeometryOps as GO
import GeometryOpsCore as GOCore
import GeoInterface as GI
import ConservativeRegridding: Trees

const US = GO.UnitSpherical
const USPoint = GO.UnitSphericalPoint{Float64}
# Generic constructor and concrete annotation aliases.
const SphericalCap = GO.UnitSpherical.SphericalCap
const Cap = GO.UnitSpherical.SphericalCap{Float64}

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

# Flat spatial tree

"""
    ToyCapTree(space, indices)

A one-node spatial tree with stored caps. Indices are cell positions for cell
trees and chunk numbers for chunk trees.
"""
struct ToyCapTree{S}
    space::S
    indices::Vector{Int}
    caps::Vector{Cap}
    extent::Cap
end

function ToyCapTree(space, indices, caps)
    ix = collect(Int, indices)
    cs = collect(Cap, caps)
    extent = isempty(cs) ? TOY_FULL_SPHERE :
             foldl(US._merge, cs)
    return ToyCapTree{typeof(space)}(space, ix, cs, extent)
end

Base.show(io::IO, tree::ToyCapTree) =
    print(io, "ToyCapTree(", length(tree.indices), " entries)")

GO.SpatialTreeInterface.isspatialtree(::Type{<:ToyCapTree}) = true
GO.SpatialTreeInterface.node_extent_is_expensive(::Type{<:ToyCapTree}) = false
GO.SpatialTreeInterface.isleaf(::ToyCapTree) = true
GO.SpatialTreeInterface.nchild(::ToyCapTree) = 0
GO.SpatialTreeInterface.getchild(::ToyCapTree) = ()
GO.SpatialTreeInterface.node_extent(tree::ToyCapTree) = tree.extent
GO.SpatialTreeInterface.child_indices_extents(tree::ToyCapTree) =
    zip(tree.indices, tree.caps)

# Cell-tree access uses the wrapped space's global positions.
GOCore.best_manifold(tree::ToyCapTree) = manifold(tree.space)
Trees.ncells(tree::ToyCapTree) = ncells(tree.space)
Trees.getcell(tree::ToyCapTree, i::Int) = getcell(tree.space, i)
Trees.getcell(tree::ToyCapTree) = (getcell(tree.space, i) for i in tree.indices)

# Lon/lat test space

"""
    ToyLonLatSpace(nlon, nlat; lon = (-180, 180), lat = (-90, 90), chunks = (nlon, nlat))

A regular `nlon × nlat` geodesic lon/lat grid on the unit sphere. Longitude
varies fastest. `chunks` gives `(longitude, latitude)` chunk sizes.
"""
struct ToyLonLatSpace <: RegridSpace
    nlon::Int
    nlat::Int
    lon0::Float64
    lon1::Float64
    lat0::Float64
    lat1::Float64
    chunklon::Int
    chunklat::Int
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
    return ToyLonLatSpace(Int(nlon), Int(nlat),
        Float64(lon[1]), Float64(lon[2]), Float64(lat[1]), Float64(lat[2]),
        min(Int(chunks[1]), Int(nlon)), min(Int(chunks[2]), Int(nlat)))
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
    cellposition(space::ToyLonLatSpace, ix::Int, iy::Int) -> Int

Convert between cell positions and lattice coordinates.
"""
function cellsubscript(space::ToyLonLatSpace, i::Int)
    1 <= i <= ncells(space) || throw(BoundsError(space, i))
    iy, ix = fldmod1(i, space.nlon)
    return (ix, iy)
end

function cellposition(space::ToyLonLatSpace, ix::Integer, iy::Integer)
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

# `RegridSpace` interface

ncells(space::ToyLonLatSpace) = space.nlon * space.nlat

manifold(::ToyLonLatSpace) = GOCore.Spherical(; radius = 1.0)

hascellchart(::ToyLonLatSpace) = true

"""
    cellcorners(space::ToyLonLatSpace, i::Int) -> NTuple{4,UnitSphericalPoint}

Return counter-clockwise graticule corners. Polar cells repeat one corner.
"""
function cellcorners(space::ToyLonLatSpace, i::Int)
    lon_lo, lon_hi, lat_lo, lat_hi = cellbounds(space, i)
    return (toy_point(lon_lo, lat_lo), toy_point(lon_hi, lat_lo),
        toy_point(lon_hi, lat_hi), toy_point(lon_lo, lat_hi))
end

getcell(space::ToyLonLatSpace, i::Int) =
    GI.Polygon([GI.LinearRing(_closed_ring(cellcorners(space, i)))])

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

function cellcentroid(space::ToyLonLatSpace, i::Int)
    lon_lo, lon_hi, lat_lo, lat_hi = cellbounds(space, i)
    return toy_point((lon_lo + lon_hi) / 2, (lat_lo + lat_hi) / 2)
end

"""
    cellat(space::ToyLonLatSpace, p) -> Union{Int,Nothing}

Return the cell containing `p`, wrapping global longitude. Return `nothing`
outside partial coverage. Shared boundaries select the east/north cell.
"""
function cellat(space::ToyLonLatSpace, p)
    lon, lat = toy_lonlat(p)
    wrapped = _wrap_lon(lon, space.lon0, space.lon1)
    wrapped === nothing && return nothing
    space.lat0 <= lat <= space.lat1 || return nothing
    ix = clamp(floor(Int, (wrapped - space.lon0) / dlon(space)) + 1, 1, space.nlon)
    iy = clamp(floor(Int, (lat - space.lat0) / dlat(space)) + 1, 1, space.nlat)
    return cellposition(space, ix, iy)
end

function _wrap_lon(lon::Float64, lo::Float64, hi::Float64)
    hi - lo >= 360 && return lo + mod(lon - lo, 360.0)
    wrapped = lo + mod(lon - lo, 360.0)
    return wrapped <= hi ? wrapped : nothing
end

nchunks(space::ToyLonLatSpace) = nchunklon(space) * nchunklat(space)

function chunksubscript(space::ToyLonLatSpace, chunk::Int)
    1 <= chunk <= nchunks(space) || throw(BoundsError(space, chunk))
    cy, cx = fldmod1(chunk, nchunklon(space))
    return (cx, cy)
end

function cellindices(space::ToyLonLatSpace, chunk::Int)
    cx, cy = chunksubscript(space, chunk)
    ix0 = (cx - 1) * space.chunklon + 1
    ix1 = min(space.nlon, cx * space.chunklon)
    iy0 = (cy - 1) * space.chunklat + 1
    iy1 = min(space.nlat, cy * space.chunklat)
    # Full-width chunks are contiguous in position order.
    nchunklon(space) == 1 &&
        return cellposition(space, 1, iy0):cellposition(space, space.nlon, iy1)
    out = Vector{Int}(undef, (ix1 - ix0 + 1) * (iy1 - iy0 + 1))
    k = 0
    for iy in iy0:iy1, ix in ix0:ix1
        out[k += 1] = cellposition(space, ix, iy)
    end
    return out
end

celltree(space::ToyLonLatSpace) = ToyCapTree(space, 1:ncells(space))

ToyCapTree(space::ToyLonLatSpace, indices) =
    ToyCapTree(space, collect(Int, indices),
        [toy_cap(cellcorners(space, i)) for i in indices])

# Fall back to a latitude-band cap when a corner cap becomes non-convex.
function chunktree(space::ToyLonLatSpace)
    n = nchunks(space)
    caps = Vector{Cap}(undef, n)
    points = USPoint[]
    for c in 1:n
        empty!(points)
        lat0, lat1 = 90.0, -90.0
        for i in cellindices(space, c)
            append!(points, cellcorners(space, i))
            _, _, a, b = cellbounds(space, i)
            lat0, lat1 = min(lat0, a), max(lat1, b)
        end
        cap = toy_cap(points)
        caps[c] = cap.radius >= Float64(pi) ? toy_bandcap(lat0, lat1, dlon(space)) : cap
    end
    return ToyCapTree(space, collect(1:n), caps)
end

# Tree-build counting wrapper

"""
    CountingSpace(space)

Wrap a `RegridSpace` and count restricted-tree builds at the `GR.subtree`
seam, delegating every other space verb unchanged. Read `cs.builds[]`; block
builders may run on `Threads.@spawn` tasks, hence the atomic counter.
"""
struct CountingSpace{S<:RegridSpace} <: RegridSpace
    space::S
    builds::Threads.Atomic{Int}
end

CountingSpace(space::RegridSpace) = CountingSpace(space, Threads.Atomic{Int}(0))

GR.subtree(cs::CountingSpace, inds) =
    (Threads.atomic_add!(cs.builds, 1); GR.subtree(cs.space, inds))

ncells(cs::CountingSpace) = ncells(cs.space)
getcell(cs::CountingSpace, i::Int) = getcell(cs.space, i)
manifold(cs::CountingSpace) = manifold(cs.space)
hascellchart(cs::CountingSpace) = hascellchart(cs.space)
cellcentroid(cs::CountingSpace, i::Int) = cellcentroid(cs.space, i)
cellat(cs::CountingSpace, p) = cellat(cs.space, p)
nchunks(cs::CountingSpace) = nchunks(cs.space)
cellindices(cs::CountingSpace, chunk::Int) = cellindices(cs.space, chunk)
celltree(cs::CountingSpace) = celltree(cs.space)
chunktree(cs::CountingSpace) = chunktree(cs.space)

# Geometry-free test method

"""
    ToyDiagonalMethod(; scale = 1.0, withdenom = true)

Build diagonal weights of `scale` for shared cell positions. `withdenom = false`
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

function build_weights!(coo::WeightCOO, method::ToyDiagonalMethod,
    ::RegridSpace, dst_inds, ::RegridSpace, src_inds)
    local_of = Dict{Int,Int}(p => k for (k, p) in enumerate(src_inds))
    for (j, p) in enumerate(dst_inds)
        k = get(local_of, p, 0)
        k == 0 && continue
        addweight!(coo, j, k, method.scale)
        method.withdenom && adddenom!(coo, j, method.scale)
    end
    return coo
end

"""
    WaveFailMethod(bad, delay)

Fail the build of the source chunk whose first cell position is `bad`, and make
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

function build_weights!(coo::WeightCOO, method::WaveFailMethod,
    ::RegridSpace, dst_inds, ::RegridSpace, src_inds)
    Int(first(src_inds)) == method.bad &&
        error("WaveFailMethod: the chunk at position $(method.bad) fails by design")
    sleep(method.delay)
    local_of = Dict{Int,Int}(p => k for (k, p) in enumerate(src_inds))
    for (j, p) in enumerate(dst_inds)
        k = get(local_of, p, 0)
        k == 0 && continue
        addweight!(coo, j, k, 1.0)
        adddenom!(coo, j, 1.0)
    end
    Threads.atomic_add!(method.finished, 1)
    return coo
end
