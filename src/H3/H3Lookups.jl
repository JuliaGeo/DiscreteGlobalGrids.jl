module H3Lookups

# `H3Native` is a sibling submodule of `H3`, bound here with `import` (not
# `using`) so its generic exports (`cell_area`, `cell_boundary`,
# `MAX_RESOLUTION`, ...) do not collide with this module's own definitions.
import ..H3Native
import ..Helpers
using ..Helpers: strictly_increasing, to_uint64_id
# The package namespace, reached the way the kernel wirings reach it — this
# module is a grandchild of it, hence the third dot. It supplies the lookup
# supertype and `DGGSGlobeIds`.
import ...DiscreteGlobalGrids as DGG
using DimensionalData
using DimensionalData: Lookups
using GeoInterface
using GeometryOps

const DD = DimensionalData
const GI = GeoInterface

export H3Lookup,
    H3Cells,
    H3Native,
    Touching,
    cell_area,
    cell_boundary,
    cell_center,
    cell_centers,
    cell_polygon,
    cell_polygons,
    cell_to_children,
    cell_to_parent,
    lonlat_to_cell,
    res0_cells

DD.@dim H3Cells "H3 cell ids"

"""
    H3Lookup(cell_ids; resolution, metadata=Dict{String,Any}(), validate=false)
    H3Lookup(globe_ids::DGGSGlobeIds)

A `DimensionalData` lookup over H3 cell ids (`UInt64`s) at a single resolution,
sorted strictly ascending — the order every selector's binary search assumes.

The cheap structural invariants are enforced on *every* construction path (this
one, the raw positional form, `DD.rebuild` on a slice): sorted and unique ids,
`resolution` in range, and — since an H3 index encodes its own resolution in
bits 52:55 — the first and last id's encoded resolution against the declared
one. `validate=true` adds the expensive per-id pass: `is_valid_cell` plus the
resolution of *every* id. It is off by default because that pass is O(n) native
calls, which is real money on a continent-sized lookup.

The second form is the globe-complete dimension: `H3Lookup(DGGSGlobeIds(H3DGGS(),
9))` holds all 4.8e9 res-9 cells in two words, defaulting `resolution` to the
ids' own level. There `validate` is a documented no-op — a
[`DGGSGlobeIds`](@ref) enumerates `ordinal_to_cell`, so its ids are valid at
`resolution` by construction, and running the per-id pass would be an O(n) trap
that establishes nothing.
"""
struct H3Lookup{A<:AbstractVector{UInt64},M} <: DGG.AbstractDGGSLookup{UInt64}
    data::A
    resolution::Int
    metadata::M
    # The keyword form below used to be the only door these checks stood at,
    # with the default positional constructor open beside it — so
    # `H3Lookup(res5_ids, 6, md)` built a lookup whose every `Contains`
    # selector answers "not stored" for cells that are stored. The endpoint
    # resolution check is two bit-ops and closes exactly that: with the ids
    # sorted, a whole vector at the wrong resolution cannot hide behind them.
    # A mixed-resolution *interior* still needs `validate=true` — an id is an
    # opaque encoding and only the native layer can judge it.
    function H3Lookup(data::AbstractVector{UInt64}, resolution::Integer, metadata)
        strictly_increasing(data) ||
            throw(ArgumentError("H3 cell ids must be sorted ascending and unique"))
        res = Int(resolution)
        0 <= res <= H3Native.MAX_RESOLUTION ||
            throw(ArgumentError("H3 resolution must be in 0:$(H3Native.MAX_RESOLUTION)"))
        isempty(data) || (H3Native.get_resolution(first(data)) == res &&
                          H3Native.get_resolution(last(data)) == res) ||
            throw(ArgumentError("H3 cell ids do not all match resolution $res"))
        return new{typeof(data),typeof(metadata)}(data, res, metadata)
    end
end

const _to_id = to_uint64_id

function H3Lookup(cell_ids; resolution::Integer, metadata=Dict{String,Any}(), validate::Bool=false)
    ids = [_to_id(id) for id in cell_ids]
    res = Int(resolution)
    md = merge(
        Dict{String,Any}("grid_name" => "h3", "resolution" => res, "indexing_scheme" => "h3-u64"),
        metadata,
    )
    # The cheap structural checks run first, in the inner constructor.
    lookup = H3Lookup(ids, res, md)
    if validate
        all(id -> H3Native.is_valid_cell(id) && H3Native.get_resolution(id) == res, ids) ||
            throw(ArgumentError("H3 cell ids do not all match resolution $res"))
    end
    return lookup
end

# The globe fast path. Both O(n) steps of the constructor above have to be
# skipped, not merely made faster: the comprehension would materialize every id
# (~39 GB at res 9) and `validate`'s per-id pass would run 4.8e9 native calls to
# re-establish what the ordinal contract already guarantees. The inner
# constructor's own O(n) step is skipped by `strictly_increasing(::DGGSGlobeIds)`
# in `core/globe_ids.jl`; the endpoint resolution check still runs, and still
# rejects a level that disagrees with `resolution`, at O(1).
#
# This is load-bearing rather than an optimization: `DD.rebuild` defaults
# `data=l.data` and routes through the keyword constructor, so without this
# method *any* rebuild — `format`, `set`, a metadata change — would silently
# materialize the whole globe.
function H3Lookup(cell_ids::DGG.DGGSGlobeIds; resolution::Integer=cell_ids.level,
        metadata=Dict{String,Any}(), validate::Bool=false)
    # Rejected here rather than left to the inner constructor's endpoint
    # check, which reads a foreign id's bits and could conceivably agree by
    # accident — and because falling through to the generic method is the
    # materialization this whole path exists to prevent.
    cell_ids.system isa DGG.H3DGGS || throw(ArgumentError(
        "an H3Lookup cannot be built from a $(DGG.system_name(cell_ids.system)) globe"))
    res = Int(resolution)
    md = merge(
        Dict{String,Any}("grid_name" => "h3", "resolution" => res, "indexing_scheme" => "h3-u64"),
        metadata,
    )
    return H3Lookup(cell_ids, res, md)
end

H3Lookup(l::H3Lookup; data=l.data, resolution=l.resolution, metadata=l.metadata) =
    H3Lookup(data; resolution, metadata)

Lookups.metadata(l::H3Lookup) = l.metadata
DD.parent(l::H3Lookup) = l.data
DD.order(::H3Lookup) = Lookups.ForwardOrdered()
Base.size(l::H3Lookup) = size(l.data)
Base.length(l::H3Lookup) = length(l.data)
Base.getindex(l::H3Lookup, i::Int) = l.data[i]
Base.firstindex(l::H3Lookup) = firstindex(l.data)
Base.lastindex(l::H3Lookup) = lastindex(l.data)
Base.iterate(l::H3Lookup, state...) = iterate(l.data, state...)

function DD.rebuild(l::H3Lookup; data=l.data, metadata=l.metadata, kw...)
    return H3Lookup(data; resolution=l.resolution, metadata, kw...)
end

DD.format(dim::DD.Dimension{<:H3Lookup}, ::AbstractRange) = dim
DD.Dimensions.format(l::H3Lookup, ::Type, values, axis::AbstractRange) = l
Lookups.reducelookup(l::H3Lookup) = Lookups.NoLookup(Base.OneTo(1))
Lookups.show_properties(io::IO, l::H3Lookup) =
    print(io, " resolution: ", l.resolution, " ncells: ", length(l.data))
Lookups.show_properties(io::IO, mime, l::H3Lookup) = Lookups.show_properties(io, l)

"""
    Touching(geometry)

Selector for every stored cell whose (spherical) polygon intersects
`geometry` (`Contains` selects by cell *center* instead). Answered by the
package-level spherical tree descent (`src/core/lookup_ops.jl`): `geometry`
may cross the antimeridian or enclose a pole, and its ring edges are
great-circle arcs — see [`DiscreteGlobalGrids.zonal`](@ref) for the full
geometry conventions.
"""
struct Touching{G} <: Lookups.ArraySelector{G}
    val::G
end

Lookups.val(sel::Touching) = sel.val

# The lookup position holding a cell id, or `nothing`. One line here because
# the two branches — binary search over stored ids, `cell_to_ordinal` over a
# globe — are one generic method pair in core (`core/lookups.jl`), chosen by
# the type of the id vector and so invisible at every call site below.
_cell_position(l::H3Lookup, id::UInt64) = DGG.cell_position(l.data, id)

Lookups.selectindices(l::H3Lookup, sel::Lookups.StandardIndices) = sel

function Lookups.selectindices(l::H3Lookup, sel::Lookups.At)
    id = _to_id(Lookups.val(sel))
    i = _cell_position(l, id)
    isnothing(i) && throw(ArgumentError("H3 cell $(string(id; base=16)) is not present in the lookup"))
    return i
end

function _lonlat(value)
    if value isa Tuple && length(value) == 2
        return (Float64(value[1]), Float64(value[2]))
    elseif GI.isgeometry(value)
        trait = GI.trait(value)
        trait isa GI.PointTrait && return (Float64(GI.x(value)), Float64(GI.y(value)))
    end
    return nothing
end

function Lookups.selectindices(l::H3Lookup, sel::Lookups.Contains)
    value = Lookups.val(sel)
    point = _lonlat(value)
    if !isnothing(point)
        id = H3Native.lonlat_to_cell(point[1], point[2], l.resolution)
        i = _cell_position(l, id)
        return isnothing(i) ? Int[] : [i]
    end

    return DGG._query_positions(l, value, :center)
end

Lookups.selectindices(l::H3Lookup, sel::Touching) =
    DGG._query_positions(l, Lookups.val(sel), :touches)

lonlat_to_cell(lon::Real, lat::Real, resolution::Integer) =
    H3Native.lonlat_to_cell(lon, lat, resolution)

res0_cells() = H3Native.res0_cells()
cell_area(id) = H3Native.cell_area(_to_id(id))
cell_boundary(id; closed_ring::Bool=true) = H3Native.cell_boundary(_to_id(id); closed_ring)
cell_center(id) = H3Native.cell_center(_to_id(id))
cell_to_children(id, resolution::Union{Nothing,Integer}=nothing) =
    H3Native.cell_to_children(_to_id(id), resolution)
cell_to_parent(id, resolution::Integer) =
    H3Native.cell_to_parent(_to_id(id), resolution)

cell_centers(l::H3Lookup) = [cell_center(id) for id in l.data]

function cell_polygon(id)
    ring = cell_boundary(id; closed_ring=true)
    isempty(ring) && throw(ArgumentError("H3 cell boundary is empty"))
    return GI.Polygon([GI.LinearRing(ring)])
end

cell_polygons(l::H3Lookup) = [cell_polygon(id) for id in l.data]

end
