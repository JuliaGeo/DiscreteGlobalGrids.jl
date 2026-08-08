"""
    HealpixLookups

Prototype of a HEALPix DGGS lookup for DimensionalData/Rasters-style vector data cubes,
following the EOPF/ESA (GRID4EARTH) conventions:

- a 1-D vector of **sorted, 0-based, nested-order** `Int64` cell ids at a fixed zoom
  `level` (`nside = 2^level`, `12 * 4^level` global cells);
- **partial coverage**: only covered cells are stored, nothing is materialized for
  uncovered cells;
- metadata attrs `grid_name="healpix"`, `level`, `indexing_scheme="nested"` (xdggs/EOPF).

Spatial queries use the nested hierarchy directly (the same tree as
`ConservativeRegriddingHealpixExt.HealpixTreeNode`): pixel `p` at level `l` owns the leaf
id range `[p*4^Δ, (p+1)*4^Δ)`, so with sorted ids "the sort order is the index" — tree
descent prunes on `searchsorted` range emptiness plus spatial disjointness.

Limitations (prototype): spherical HEALPix (not the WGS84 authalic-sphere variant);
planar lon/lat predicates, so query geometries crossing the antimeridian or enclosing a
pole are unsupported; `stencil` handles 1-D (Cells-only) arrays.
"""
module HealpixLookups

import ..Helpers
using ..Helpers: strictly_increasing
# The package namespace, reached the way the kernel wirings reach it — this
# module is a grandchild of it, hence the third dot. It supplies the lookup
# supertype and `DGGSGlobeIds`.
import ...DiscreteGlobalGrids as DGG
import DimensionalData as DD
import GeometryOps as GO
import GeoInterface as GI
import Extents
import Healpix
using DimensionalData: Lookups
using Statistics

export HealpixLookup, Cells, Touching, nested_neighbors, zonal, stencil, cell_centers, cell_polygons

const MAX_QUERY_DENSIFY_DELTA = 8

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
# near the center lon. step > 1 densifies the (curved) edges; internal tree nodes MUST
# pass step = 2^Δ so the planar outline exactly bounds the descendant leaf polygons.
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
# kept because `_point_index` and `neighbor_indices` below both read it. The two
# branches behind it — binary search over stored ids, `cell_to_ordinal` over a
# globe — are one generic method pair in core (`core/lookups.jl`), chosen by the
# type of the id vector and so invisible at every call site here.
_cellid_index(l, cid) = something(DGG.cell_position(l.data, cid), 0)
function Lookups.selectindices(l::HealpixLookup, sel::Lookups.Contains)
    v = Lookups.val(sel)
    if GI.isgeometry(v) && GI.trait(v) isa GI.PointTrait
        return _point_index(l, GI.x(v), GI.y(v))
    elseif v isa Tuple{Real,Real}
        return _point_index(l, v[1], v[2])
    elseif GI.isgeometry(v)
        return _query_indices(l, v, :center)
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

Selector for stored cells whose polygon intersects `geom`. This is the geometry analogue
of `DD.Lookups.Touches`, which only admits tuples/`Extents.Extent`s in its type parameter
and therefore cannot carry a polygon. Use as `A[Cells(Touching(geom))]`.
"""
struct Touching{T} <: Lookups.ArraySelector{T}
    val::T
end
Lookups.val(sel::Touching) = sel.val
Lookups.selectindices(l::HealpixLookup, sel::Touching) = _query_indices(l, sel.val, :touches)

# DD-native Touches carries tuples or extents; treat them as a box polygon.
function Lookups.selectindices(l::HealpixLookup, sel::Lookups.Touches)
    ext = Lookups.val(sel)
    ext isa Extents.Extent ||
        throw(ArgumentError("only Touches(::Extents.Extent) is supported; use Touching(geom) for geometries"))
    (x1, x2), (y1, y2) = ext.X, ext.Y
    box = GI.Polygon([GI.LinearRing([(x1, y1), (x2, y1), (x2, y2), (x1, y2), (x1, y1)])])
    return _query_indices(l, box, :touches)
end

#=
## Hierarchical cover query
=#

"""
Tree-descent cover. `mode=:center`: stored cells whose center is contained in `geom`.
`mode=:touches`: stored cells whose (planar lon/lat) polygon intersects `geom`.
Prunes on (a) emptiness of the stored-id range under a node (binary search — the sort
order is the index) and (b) spatial disjointness of the node's *densified* polygon from
`geom`. The densification (`step = 2^Δ`) is load-bearing: with only 4 corners the planar
diamond under-covers the curved node region and the prune drops real results.
Limitation: planar lon/lat predicates — geometries crossing the antimeridian or
enclosing a pole are not handled.
"""
function _query_indices(l::HealpixLookup, geom, mode::Symbol)
    out = Int[]
    geomext = GI.extent(geom)
    for f in 0:11
        _query!(out, l, geom, geomext, 0, f, mode)
    end
    return out
end

_query_polygon_step(Δ::Integer) =
    Δ <= MAX_QUERY_DENSIFY_DELTA ? (1 << Int(Δ)) : nothing

function _query!(out, l, geom, geomext, level, pix, mode)
    Δ = l.level - level
    lo = Int64(pix) * 4^Δ
    hi = lo + 4^Δ
    i1 = searchsortedfirst(l.data, lo)
    (i1 > length(l.data) || l.data[i1] >= hi) && return        # nothing stored under this node

    step = _query_polygon_step(Δ)
    if step === nothing
        for c in 0:3
            _query!(out, l, geom, geomext, level + 1, 4 * pix + c, mode)
        end
        return
    end

    # densified node outline: exactly bounds the descendant leaf planar polygons
    nodepoly = _cell_polygon_lonlat(Healpix.Resolution(2^level), pix; step)
    Extents.disjoint(GI.extent(nodepoly), geomext) && return
    if level == l.level                                        # leaf: i1 is the stored index
        hit = mode === :center ?
            GO.contains(geom, _cell_center_lonlat(resolution(l), pix)) :
            GO.intersects(geom, nodepoly)
        hit && push!(out, i1)
        return
    end
    GO.disjoint(geom, nodepoly) && return
    if GO.covers(geom, nodepoly)   # whole subtree inside: valid for both modes
        i2 = searchsortedlast(l.data, hi - 1)
        append!(out, i1:i2)
        return
    end
    for c in 0:3
        _query!(out, l, geom, geomext, level + 1, 4 * pix + c, mode)
    end
end

#=
## Zonal statistics
=#

"""
    zonal(f, A::DD.AbstractDimArray; of, boundary=:center, skipmissing=true)

Zonal statistics over the `Cells` dimension of `A` (which must hold a `HealpixLookup`).
`of` is a geometry, feature(collection), or vector thereof. HEALPix cells are equal-area,
so e.g. `zonal(mean, ...)` is the true (unweighted) areal mean — no latitude weighting.
Returns one value per geometry; `missing` where no stored cell matches.
"""
function zonal(f, A::DD.AbstractDimArray; of, boundary::Symbol=:center, skipmissing::Bool=true)
    l = DD.lookup(A, Cells)
    geoms = _geometries(of)
    mode = boundary === :center ? :center : :touches
    map(geoms) do g
        idx = _query_indices(l, g, mode)
        isempty(idx) && return missing
        sub = A[Cells(idx)]
        vals = skipmissing ? Base.skipmissing(sub) : sub
        isempty(vals) ? missing : f(vals)
    end
end

_geometries(of) = GI.isgeometry(of) ? [of] :
    GI.trait(of) isa GI.AbstractFeatureCollectionTrait ? [GI.geometry(f) for f in GI.getfeature(of)] :
    GI.trait(of) isa GI.AbstractFeatureTrait ? [GI.geometry(of)] :
    of isa AbstractVector ? map(g -> GI.trait(g) isa GI.AbstractFeatureTrait ? GI.geometry(g) : g, of) :
    throw(ArgumentError("cannot extract geometries from $(typeof(of))"))

#=
## Stencil operations
=#

"""
    neighbor_indices(l::HealpixLookup) -> Vector{NTuple{8,Int}}

For each stored cell, positions (into the lookup) of its 8 HEALPix neighbors;
0 where the neighbor cell is not stored (coverage boundary) or does not exist.
Computed once; this is the static "halo table" (DLWP-HPX pattern from the research doc).
"""
function neighbor_indices(l::HealpixLookup)
    res = resolution(l)
    map(l.data) do cid
        nbs = nested_neighbors(res, cid)
        ntuple(m -> nbs[m] < 0 ? 0 : _cellid_index(l, nbs[m]), 8)
    end
end

"""
    stencil(f, A::DD.AbstractDimArray; nbidx=nothing)

Apply `f(center_value, neighbor_values::Vector)` over every stored cell of the (1-D)
`Cells` array `A`. Neighbors outside the stored coverage are simply absent from the
vector (partial-coverage semantics: reductions skip, like `nanmean`). Pass a precomputed
`neighbor_indices(lookup)` as `nbidx` to amortize the halo table across many stencils.
"""
function stencil(f, A::DD.AbstractDimArray; nbidx=nothing)
    l = DD.lookup(A, Cells)
    nbi = nbidx === nothing ? neighbor_indices(l) : nbidx
    data = parent(A)
    length(nbi) == length(data) ||
        throw(DimensionMismatch("neighbor index table and array must have the same length"))
    out = similar(data, Base.promote_op(f, eltype(data), Vector{eltype(data)}))
    for i in eachindex(data)
        indices = nbi[i]
        nneighbors = count(>(0), indices)
        values = Vector{eltype(data)}(undef, nneighbors)
        k = 0
        @inbounds for j in indices
            if j > 0
                k += 1
                values[k] = data[j]
            end
        end
        @inbounds out[i] = f(data[i], values)
    end
    DD.rebuild(A; data=out)
end

end # module
