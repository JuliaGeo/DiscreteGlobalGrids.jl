"""
    HealpixLookups

Prototype of a HEALPix DGGS lookup for DimensionalData/Rasters-style vector data cubes,
following the EOPF/ESA (GRID4EARTH) conventions:

- a 1-D vector of **sorted, 0-based, nested-order** `Int64` cell ids at a fixed zoom
  `level` (`nside = 2^level`, `12 * 4^level` global cells);
- **partial coverage**: only covered cells are stored, nothing is materialized for
  uncovered cells;
- metadata attrs `grid_name="healpix"`, `level`, `indexing_scheme="nested"` (xdggs/EOPF).

Spatial queries run the package-level spherical tree descent
(`_query_positions`, `src/core/lookup_ops.jl`) over the nested hierarchy:
pixel `p` at level `l` owns the leaf id range `[p*4^Δ, (p+1)*4^Δ)`, so with
sorted ids "the sort order is the index" — the tree's binary searches prune
on range emptiness, HEALPix's O(1) exact subtree caps on geometry. All
predicates are spherical, so query geometries may cross the antimeridian or
enclose a pole; their ring edges are great-circle arcs.

Limitations (prototype): spherical HEALPix (not the WGS84 authalic-sphere variant).
"""
module HealpixLookups

import ..Helpers
using ..Helpers: strictly_increasing
# The package namespace, reached the way the kernel wirings reach it — this
# module is a grandchild of it, hence the third dot. It supplies the lookup
# supertype and `DGGSGlobeIds`.
import ...DiscreteGlobalGrids as DGG
# `zonal` and `stencil` started life here and are still exported here, but as
# the package-level generics (`src/core/lookup_ops.jl`) — the same bindings,
# imported back, so a `using` of both namespaces cannot make them ambiguous.
# The spatial query underneath them is package-level too, the spherical tree
# descent this module's selectors route through (`_query_positions`).
import ...DiscreteGlobalGrids: zonal, stencil, neighbor_indices
import DimensionalData as DD
import GeoInterface as GI
import Extents
import Healpix
using DimensionalData: Lookups

export HealpixLookup, Cells, Touching, nested_neighbors, zonal, stencil, cell_centers, cell_polygons

#=
## 8-neighbors in the nested scheme

Transcription of the HEALPix C++ `neighbors()` algorithm (healpix_base.cc).
Healpix.jl exposes the xyf <-> nested conversions but no neighbor function.
=#
const NB_XOFFSET = (-1, -1, 0, 1, 1, 1, 0, -1)
const NB_YOFFSET = (0, 1, 1, 1, 0, -1, -1, -1)
const NB_FACEARRAY = (
    (8, 9, 10, 11, -1, -1, -1, -1, 10, 11, 8, 9),   # S
    (5, 6, 7, 4, 8, 9, 10, 11, 9, 10, 11, 8),       # SE
    (-1, -1, -1, -1, 5, 6, 7, 4, -1, -1, -1, -1),   # E
    (4, 5, 6, 7, 11, 8, 9, 10, 11, 8, 9, 10),       # SW
    (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11),         # center
    (1, 2, 3, 0, 0, 1, 2, 3, 5, 6, 7, 4),           # NE
    (-1, -1, -1, -1, 7, 4, 5, 6, -1, -1, -1, -1),   # W
    (3, 0, 1, 2, 3, 0, 1, 2, 4, 5, 6, 7),           # NW
    (2, 3, 0, 1, -1, -1, -1, -1, 0, 1, 2, 3),       # N
)
const NB_SWAPARRAY = (
    (0, 0, 3), (0, 0, 6), (0, 0, 0), (0, 0, 5), (0, 0, 0),
    (5, 0, 0), (0, 0, 0), (6, 0, 0), (3, 0, 0),
)

"""
    nested_neighbors(res::Healpix.Resolution, pix::Integer) -> NTuple{8,Int}

8 neighbors of 0-based nested pixel `pix`, in order SW, W, NW, N, NE, E, SE, S.
Nonexistent neighbors (24 pixels per grid sit at a degree-3 base-tiling vertex and have
only 7) are -1.
"""
function nested_neighbors(res::Healpix.Resolution, pix::Integer)
    nside = res.nside
    x, y, f = Healpix.pix2xyfNest(res, Int(pix) + 1)   # x, y 0-based; pixel arg 1-based
    ntuple(8) do m
        xx = x + NB_XOFFSET[m]
        yy = y + NB_YOFFSET[m]
        nbnum = 5                                       # 1-based "center"
        if xx < 0
            xx += nside; nbnum -= 1
        elseif xx >= nside
            xx -= nside; nbnum += 1
        end
        if yy < 0
            yy += nside; nbnum -= 3
        elseif yy >= nside
            yy -= nside; nbnum += 3
        end
        nf = NB_FACEARRAY[nbnum][f + 1]
        nf < 0 && return -1
        bits = NB_SWAPARRAY[nbnum][(f >> 2) + 1]
        (bits & 1) != 0 && (xx = nside - xx - 1)
        (bits & 2) != 0 && (yy = nside - yy - 1)
        (bits & 4) != 0 && ((xx, yy) = (yy, xx))
        return Healpix.xyf2pixNest(res, xx, yy, nf) - 1
    end
end

#=
## The lookup type
=#

DD.@dim Cells "cells"

"""
    HealpixLookup(cell_ids::AbstractVector{<:Integer}; level::Int)
    HealpixLookup(globe_ids::DGGSGlobeIds)

EOPF/ESA-convention HEALPix lookup: sorted, 0-based, **nested-order** Int64 cell ids at
zoom `level` (nside = 2^level); only covered cells are stored (partial coverage).
Spherical HEALPix (not the WGS84 authalic-sphere variant).

Every construction path — this one, `DD.rebuild` on a slice, the raw positional form —
runs the same three checks, and `nside` is always `2^level` (see the struct).

The second form is the complementary case, complete coverage:
`HealpixLookup(DGGSGlobeIds(HEALPixDGGS(), 12))` holds all 8e8 level-12 cells in two
words, defaulting `level` to the ids' own.
"""
struct HealpixLookup{A<:AbstractVector{Int64},M} <: DGG.AbstractDGGSLookup{Int64}
    data::A
    level::Int
    nside::Int
    metadata::M
    # `nside` is derived, not supplied: it is `2^level` and nothing else, so it
    # is computed here and there is no constructor that takes it. (An `nside`
    # decoupled from `level` would send `resolution(l)` — hence every center,
    # polygon, neighbor and point lookup — to a different grid than the ids
    # were indexed in, silently.) The three checks below travel with it for the
    # same reason: they are the invariants every method on this type assumes,
    # and the raw positional constructor used to sit open beside the keyword
    # form that ran them.
    function HealpixLookup(data::AbstractVector{Int64}, level::Integer, metadata)
        lvl = Int(level)
        nside = Healpix.order2nside(lvl)   # DomainError outside 0:ORDER_MAX
        strictly_increasing(data) ||
            throw(ArgumentError("cell_ids must be sorted ascending and unique (EOPF nested-order convention)"))
        npix = Healpix.nside2npix(nside)
        (isempty(data) || (first(data) >= 0 && last(data) < npix)) ||
            throw(ArgumentError("cell_ids must be in [0, $npix) for level $lvl"))
        return new{typeof(data),typeof(metadata)}(data, lvl, nside, metadata)
    end
end
function HealpixLookup(cell_ids::AbstractVector{<:Integer}; level::Int, metadata=Dict{String,Any}())
    defaults = Dict{String,Any}(
        "grid_name" => "healpix",
        "level" => level,
        "indexing_scheme" => "nested",
    )
    return HealpixLookup(collect(Int64, cell_ids), level, merge(defaults, metadata))
end

# The globe fast path. The `collect` above would materialize every id (8e8 at
# level 12, 3.5e18 at level 29) and the inner constructor's `strictly_increasing`
# would walk them all; the latter is skipped by
# `Helpers.strictly_increasing(::DGGSGlobeIds)` in `core/globe_ids.jl`, the
# former by this method. The remaining `[0, npix)` endpoint check is O(1) and
# still runs, so a `level` disagreeing with the ids' own is still rejected.
#
# This is load-bearing rather than an optimization: `DD.rebuild` defaults
# `data=l.data` and routes through the keyword constructor, so without this
# method *any* rebuild — `format`, `set`, a metadata change — would silently
# materialize the whole globe.
function HealpixLookup(cell_ids::DGG.DGGSGlobeIds; level::Int=cell_ids.level,
        metadata=Dict{String,Any}())
    # Rejected here rather than left to the inner constructor's `[0, npix)`
    # range check, which a foreign system's ids could pass by accident — and
    # because falling through to the generic method is the materialization this
    # path exists to prevent.
    cell_ids.system isa DGG.HEALPixDGGS || throw(ArgumentError(
        "a HealpixLookup cannot be built from a $(DGG.system_name(cell_ids.system)) globe"))
    defaults = Dict{String,Any}(
        "grid_name" => "healpix",
        "level" => level,
        "indexing_scheme" => "nested",
    )
    return HealpixLookup(cell_ids, level, merge(defaults, metadata))
end

resolution(l::HealpixLookup) = Healpix.Resolution(l.nside)

DD.order(::HealpixLookup) = Lookups.ForwardOrdered()
DD.parent(l::HealpixLookup) = l.data
Base.size(l::HealpixLookup) = size(l.data)
Base.getindex(l::HealpixLookup, i::Int) = l.data[i]
Lookups.metadata(l::HealpixLookup) = l.metadata
DD.Dimensions.format(l::HealpixLookup, D::Type, values, axis::AbstractRange) = l
DD.format(dim::DD.Dimension{<:HealpixLookup}, axis::AbstractRange) = dim

function DD.rebuild(l::HealpixLookup; data=l.data, metadata=l.metadata, kw...)
    # Route slices and other rebuilds through the public constructor: it fills the
    # metadata defaults, and the ordering, uniqueness and level-specific range checks
    # then run in the inner constructor, which no path can go around.
    return HealpixLookup(data; level=l.level, metadata)
end
@inline Lookups.reducelookup(l::HealpixLookup) = Lookups.NoLookup(Base.OneTo(1))

# Without these, `show` of any DimArray holding this lookup throws a MethodError
# (DD 0.30 calls both arities; GeometryLookup only defines the 3-arg one).
Lookups.show_properties(io::IO, l::HealpixLookup) =
    print(io, " level: ", l.level, " ncells: ", length(l.data))
Lookups.show_properties(io::IO, mime, l::HealpixLookup) = Lookups.show_properties(io, l)

#=
## Coordinate helpers (lon/lat degrees, lon wrapped to (-180, 180])
=#
_wraplon(lon) = lon > 180 ? lon - 360 : lon
function _cell_center_lonlat(res, pix)               # pix 0-based nested
    theta, phi = Healpix.pix2angNest(res, pix + 1)
    (_wraplon(rad2deg(phi)), 90 - rad2deg(theta))
end
cell_center(l::HealpixLookup, cell_id::Integer) = _cell_center_lonlat(resolution(l), cell_id)
cell_centers(l::HealpixLookup) = [cell_center(l, c) for c in l.data]

# boundary of nested pixel as a closed lon/lat ring with 4*step points, lons unwrapped
# near the center lon. step > 1 densifies the (curved) edges. This is the planar
# *presentation* polygon (`cell_polygons` for plotting and planar consumers); the
# spatial queries run on the unit-sphere geometry of the kernel instead.
function _cell_polygon_lonlat(res, pix; step::Int=1)
    ringpix = Healpix.nest2ring(res, pix + 1)
    b = Healpix.boundariesRing(res, ringpix, step, Float64)   # (4*step) x 3 cartesian
    clon, _ = _cell_center_lonlat(res, pix)
    n = size(b, 1)
    points = Vector{Tuple{Float64,Float64}}(undef, n + 1)
    @inbounds for i in 1:n
        x, y, z = b[i, 1], b[i, 2], b[i, 3]
        lon = rad2deg(atan(y, x))
        lat = rad2deg(asin(clamp(z, -1, 1)))
        lon -= 360 * round((lon - clon) / 360)              # keep contiguous around center
        points[i] = (lon, lat)
    end
    points[end] = points[1]
    return GI.Polygon([GI.LinearRing(points)])
end
cell_polygon(l::HealpixLookup, cell_id::Integer) = _cell_polygon_lonlat(resolution(l), cell_id)
cell_polygons(l::HealpixLookup) = [cell_polygon(l, c) for c in l.data]

#=
## Selectors
=#
Lookups.selectindices(l::HealpixLookup, sel::Lookups.StandardIndices) = sel
function Lookups.selectindices(l::HealpixLookup, sel::Lookups.At{<:Integer})
    i = _cellid_index(l, Lookups.val(sel))
    iszero(i) &&
        throw(ArgumentError("cell id $(Lookups.val(sel)) not stored in lookup"))
    return i
end
# The lookup position holding a cell id, 0 when absent — this file's sentinel,
# kept because `_point_index` below reads it. The two branches behind it —
# binary search over stored ids, `cell_to_ordinal` over a globe — are one
# generic method pair in core (`core/lookups.jl`), chosen by the type of the
# id vector and so invisible at every call site here.
_cellid_index(l, cid) = something(DGG.cell_position(l.data, cid), 0)
function Lookups.selectindices(l::HealpixLookup, sel::Lookups.Contains)
    v = Lookups.val(sel)
    if GI.isgeometry(v) && GI.trait(v) isa GI.PointTrait
        return _point_index(l, GI.x(v), GI.y(v))
    elseif v isa Tuple{Real,Real}
        return _point_index(l, v[1], v[2])
    elseif GI.isgeometry(v)
        return DGG._query_positions(l, v, :center)
    end
    throw(ArgumentError("unsupported Contains selector value $v"))
end
function _point_index(l, lon, lat)
    pix = Healpix.ang2pixNest(resolution(l), deg2rad(90 - lat), deg2rad(mod(lon, 360))) - 1
    i = _cellid_index(l, pix)
    i == 0 && throw(ArgumentError("cell containing ($lon, $lat) (id $pix) not stored in lookup"))
    i
end
"""
    Touching(geom)

Selector for stored cells whose (spherical) polygon intersects `geom`. This is the
geometry analogue of `DD.Lookups.Touches`, which only admits tuples/`Extents.Extent`s
in its type parameter and therefore cannot carry a polygon. Use as
`A[Cells(Touching(geom))]`.

The query is the package-level spherical tree descent (`src/core/lookup_ops.jl`):
`geom` may cross the antimeridian or enclose a pole, and its ring edges are
great-circle arcs — see [`zonal`](@ref) for the full geometry conventions.
"""
struct Touching{T} <: Lookups.ArraySelector{T}
    val::T
end
Lookups.val(sel::Touching) = sel.val
Lookups.selectindices(l::HealpixLookup, sel::Touching) =
    DGG._query_positions(l, sel.val, :touches)

# DD-native Touches carries tuples or extents; treat an extent as the
# lon/lat-aligned box *region*.
function Lookups.selectindices(l::HealpixLookup, sel::Lookups.Touches)
    ext = Lookups.val(sel)
    ext isa Extents.Extent ||
        throw(ArgumentError("only Touches(::Extents.Extent) is supported; use Touching(geom) for geometries"))
    return DGG._query_positions(l, _extent_polygon(ext), :touches)
end

#=
An `Extents.Extent` is bounded by parallels above and below, and the spherical
query engine draws ring edges as great-circle arcs — which parallels are not —
so the box polygon densifies its two constant-latitude edges (a great-circle
chord across 0.25° of longitude stays within ~10⁻⁶ rad of the parallel,
far under any cell size this prototype reaches). The meridian sides ARE
great circles and stay two-point edges.
=#
function _extent_polygon(ext::Extents.Extent)
    (x1, x2), (y1, y2) = Float64.(ext.X), Float64.(ext.Y)
    x2 - x1 < 360 || throw(ArgumentError(
        "a full-longitude extent has no boundary meridians; query the two polar sides separately"))
    lons = range(x1, x2; length=max(2, ceil(Int, (x2 - x1) / 0.25) + 1))
    points = Vector{Tuple{Float64,Float64}}(undef, 2 * length(lons) + 1)
    k = 0
    for lon in lons
        points[k += 1] = (lon, y1)
    end
    for lon in reverse(lons)
        points[k += 1] = (lon, y2)
    end
    points[k += 1] = points[1]
    return GI.Polygon([GI.LinearRing(points)])
end

#=
## Zonal / stencil wiring

`zonal`, `stencil` and `neighbor_indices` are the package-level generics
(`src/core/lookup_ops.jl`), and so is the spatial query they and the
selectors above resolve HEALPix cells through — the spherical tree descent
over `treeify(DGGSPartialGrid(l))`. What this module contributes to that
query lives in the kernel wiring (`HealpixKernel.jl`): the O(1) exact
subtree caps the descent prunes with — which only a system whose parents
geographically contain their children can offer — plus the densified pixel
outlines (`subtree_polygon_unitsphere`) available to traversals that want
polygon-level subtree classification.
=#

end # module
