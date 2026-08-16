# Fully analytic test doubles: a coarse lon–lat space implementing the whole
# `RegridSpace` contract, and a method whose weights are a diagonal. Everything
# here is closed-form and small enough to check by hand — this is a stand-in for
# a real space, not a grid implementation.

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
# The UnionAll is what constructs (its two-argument outer constructor derives
# `radiuslike`); the parametrized alias is what annotates.
const SphericalCap = GO.UnitSpherical.SphericalCap
const Cap = GO.UnitSpherical.SphericalCap{Float64}

# ===========================================================================
# Unit-sphere helpers
# ===========================================================================

"""
    toy_point(lon, lat) -> UnitSphericalPoint

Longitude/latitude in **degrees** onto the unit sphere.
"""
function toy_point(lon::Real, lat::Real)
    coslat = cosd(Float64(lat))
    return USPoint(coslat * cosd(Float64(lon)), coslat * sind(Float64(lon)),
        sind(Float64(lat)))
end

"""
    toy_lonlat(p) -> (lon, lat)

A unit-sphere point back to longitude/latitude in **degrees**.
"""
toy_lonlat(p) = (atand(p[2], p[1]), asind(clamp(p[3], -1.0, 1.0)))

const TOY_FULL_SPHERE = SphericalCap(USPoint(0.0, 0.0, 1.0), nextfloat(Float64(pi)))

"""
    toy_cap(points) -> SphericalCap

A cap covering `points`, centred on their normalized mean. Returns the full
sphere past radius `π/2`, where the cap stops being convex and would no longer
be guaranteed to contain the great-circle arcs between the points.
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

A cap covering every cell of the latitude band `[lat0, lat1]`, whatever longitude
span the band covers.

`{lat ≥ a}` **is** a cap of radius `90° − a` about the north pole and `{lat ≤ b}`
one of radius `90° + b` about the south, so no convexity argument is needed and a
full-longitude band is bounded as tightly as a narrow one. `a` and `b` widen
`lat0`/`lat1` by the poleward bow of a `dlon`-wide cell's great-circle east–west
edge, `atan(tan φ / cos(Δλ/2))`, at whichever end faces away from the equator.
The full sphere comes back only for a band reaching both poles.

This is the bound [`toy_cap`](@ref) cannot give: a cap through the cells' corners
stops being convex past `π/2`, so a band of every longitude degenerates to the
whole sphere there and every chunk of a banded space appears to touch every other.
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

# ===========================================================================
# A flat spatial tree
# ===========================================================================

"""
    ToyCapTree(space, indices)

A one-node `SpatialTreeInterface` tree over `indices` of `space`, with stored
`SphericalCap` extents — a brute-force leaf, which is all a handful of cells
needs.

`indices` are cell positions for a cell tree and chunk numbers for a chunk
tree; `caps` are the corresponding extents. The subset constructor is the
chunk-restricted cell tree a weight builder works against.
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

# `ConservativeRegridding` addresses a cell tree as a cell source. A chunk tree
# is never handed to it, so these read the space's cells unconditionally.
GOCore.best_manifold(tree::ToyCapTree) = manifold(tree.space)
Trees.ncells(tree::ToyCapTree) = ncells(tree.space)
Trees.getcell(tree::ToyCapTree, i::Int) = getcell(tree.space, i)
Trees.getcell(tree::ToyCapTree) = (getcell(tree.space, i) for i in tree.indices)

# ===========================================================================
# The lon–lat space
# ===========================================================================

"""
    ToyLonLatSpace(nlon, nlat; lon = (-180, 180), lat = (-90, 90), chunks = (nlon, nlat))

A regular `nlon × nlat` lon–lat cell grid on the unit sphere.

Cells are geodesic quadrilaterals through the four graticule corners, in
**degrees** over the half-open longitude span `lon` and the closed latitude
span `lat`. Position `i` of cell `(ix, iy)` is `ix + (iy - 1) * nlon`:
longitude varies fastest, so a full-width chunk is a contiguous range.

`chunks` is the chunk **size** in cells, `(along lon, along lat)`; the default
is a single chunk over everything. Chunk `c = cx + (cy - 1) * nchunklon` holds
the cells of the corresponding cell block, and `cellindices` returns a
`UnitRange` exactly when chunks span the full longitude row.

Everything is closed-form: [`graticule_area`](@ref) is the exact area of the
cell's parallel-bounded footprint, and a global space's graticule areas sum to
`4π` steradians.
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

The lattice coordinates of a cell position, and their inverse. This is the
space's chart, which an interpolation stencil is written against.
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

The cell's graticule box, in degrees.
"""
function cellbounds(space::ToyLonLatSpace, i::Int)
    ix, iy = cellsubscript(space, i)
    return (space.lon0 + (ix - 1) * dlon(space), space.lon0 + ix * dlon(space),
        space.lat0 + (iy - 1) * dlat(space), space.lat0 + iy * dlat(space))
end

"""
    graticule_area(space::ToyLonLatSpace, i::Int) -> Float64

The exact area in steradians of the cell's parallel-bounded footprint,
`Δλ (sin φ_hi − sin φ_lo)`. A global space's graticule areas sum to `4π`.

This is **not** the area of [`getcell`](@ref)'s polygon, whose east–west edges
are great-circle arcs through the corners rather than parallels. Such an arc
bows toward the nearer pole, adding area along the cell's poleward edge and
removing it along the equatorward one, so the difference has no fixed sign and
vanishes only for an edge on the equator or a degenerate polar edge. The
polygons still tile the sphere exactly — neighbouring cells share an edge — so
their areas sum to `4π` as well.
"""
function graticule_area(space::ToyLonLatSpace, i::Int)
    _, _, lat_lo, lat_hi = cellbounds(space, i)
    return deg2rad(dlon(space)) * (sind(lat_hi) - sind(lat_lo))
end

# --- The contract ----------------------------------------------------------

ncells(space::ToyLonLatSpace) = space.nlon * space.nlat

manifold(::ToyLonLatSpace) = GOCore.Spherical(; radius = 1.0)

hascellchart(::ToyLonLatSpace) = true

"""
    cellcorners(space::ToyLonLatSpace, i::Int) -> NTuple{4,UnitSphericalPoint}

The cell's four graticule corners, counter-clockwise seen from outside: south-
west, south-east, north-east, north-west. A polar cell repeats a corner.
"""
function cellcorners(space::ToyLonLatSpace, i::Int)
    lon_lo, lon_hi, lat_lo, lat_hi = cellbounds(space, i)
    return (toy_point(lon_lo, lat_lo), toy_point(lon_hi, lat_lo),
        toy_point(lon_hi, lat_hi), toy_point(lon_lo, lat_hi))
end

getcell(space::ToyLonLatSpace, i::Int) =
    GI.Polygon([GI.LinearRing(_closed_ring(cellcorners(space, i)))])

# A polar cell's two upper (or lower) corners coincide; dropping the repeat
# keeps every edge non-degenerate.
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

The cell containing `p`, by closed-form inversion. Longitude is wrapped into
the space's span; a point outside the latitude span, or outside a partial
longitude span, is `nothing`.

A point on a shared boundary goes to the cell east and north of it, except on
the space's own outer boundary, which goes inward.
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
    # A full-width chunk is contiguous in position space; say so, since callers
    # are allowed to exploit it.
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

# A chunk's extent is bounded by its cells' corners rather than its own box
# corners, so it covers every cell `cellindices` assigns to it — and by its
# latitude band instead where that construction degenerates, which is any chunk
# more than a quadrant across, a full-longitude row above all.
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

# ===========================================================================
# A geometry-free method
# ===========================================================================

"""
    ToyDiagonalMethod(; scale = 1.0, withdenom = true)

A method whose weight matrix is a diagonal: source cell `p` contributes `scale`
to destination cell `p`, for every position `p` the two chunks share.

It pairs cells by position and reads no geometry at all, so it exercises the
executor — accumulation, missing handling, finalize, N-D slicing — without
depending on any weight-construction task. Over a space regridded onto itself
it is `scale` times the identity: with [`Weighted`](@ref) the result is the
source field unchanged, with [`Extensive`](@ref) it is the field times `scale`.

`withdenom = false` emits no denominators, which is how a block that finalizes
as a raw sum under either policy is produced.
"""
struct ToyDiagonalMethod <: AbstractRegriddingMethod
    scale::Float64
    withdenom::Bool
end

ToyDiagonalMethod(; scale::Real = 1.0, withdenom::Bool = true) =
    ToyDiagonalMethod(Float64(scale), withdenom)

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
