module A5Lookups

# `A5Native` is a sibling submodule of `A5`, bound here with `import` (not
# `using`) so its generic exports (`cell_area`, `cell_boundary`,
# `MAX_RESOLUTION`, ...) do not collide with this module's own definitions.
import ..A5Native
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

export A5Lookup,
    A5Cells,
    A5Native,
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

DD.@dim A5Cells "A5 cell ids"

"""
    A5Lookup(cell_ids; resolution, metadata=Dict{String,Any}(), validate=false)
    A5Lookup(globe_ids::DGGSGlobeIds)

A `DimensionalData` lookup over A5 cell ids (`UInt64`s) at a single resolution,
sorted strictly ascending — the order every selector's binary search assumes.

The cheap structural invariants are enforced on *every* construction path (this
one, the raw positional form, `DD.rebuild` on a slice): sorted and unique ids,
`resolution` in range, and the first and last id's encoded resolution against
the declared one — an A5 index carries its resolution in the position of its
low marker bit, which `A5Native.get_resolution` reads in a short shift loop.
`validate=true` adds that same resolution test over *every* id; it is off by
default because it is O(n) rather than O(1).

Neither pass proves an id is a *reachable* cell: A5's encoding has holes (the
res-30 gap, unsupported origins), and only the native deserializer can speak to
those. They surface at the first `A5Native` call that decodes the id.

The second form is the globe-complete dimension: `A5Lookup(DGGSGlobeIds(A5DGGS(),
5))` holds every res-5 cell in two words, defaulting `resolution` to the ids'
own level. There `validate` is a documented no-op — a [`DGGSGlobeIds`](@ref)
enumerates `ordinal_to_cell`, so its ids are valid at `resolution` by
construction, and running the per-id pass would be an O(n) trap that
establishes nothing.
"""
struct A5Lookup{A<:AbstractVector{UInt64},M} <: DGG.AbstractDGGSLookup{UInt64}
    data::A
    resolution::Int
    metadata::M
    # The keyword form below used to be the only door these checks stood at,
    # with the default positional constructor open beside it — so
    # `A5Lookup(res3_ids, 4, md)` built a lookup whose every `Contains`
    # selector answers "not stored" for cells that are stored. With the ids
    # sorted, checking the two endpoints means a whole vector at the wrong
    # resolution cannot hide behind them; a mixed-resolution *interior* still
    # needs `validate=true`.
    function A5Lookup(data::AbstractVector{UInt64}, resolution::Integer, metadata)
        strictly_increasing(data) ||
            throw(ArgumentError("A5 cell ids must be sorted ascending and unique"))
        res = Int(resolution)
        0 <= res <= A5Native.MAX_RESOLUTION ||
            throw(ArgumentError("A5 resolution must be in 0:$(A5Native.MAX_RESOLUTION)"))
        isempty(data) || (A5Native.get_resolution(first(data)) == res &&
                          A5Native.get_resolution(last(data)) == res) ||
            throw(ArgumentError("A5 cell ids do not all match resolution $res"))
        return new{typeof(data),typeof(metadata)}(data, res, metadata)
    end
end

const _to_id = to_uint64_id

function A5Lookup(cell_ids; resolution::Integer, metadata=Dict{String,Any}(), validate::Bool=false)
    ids = [_to_id(id) for id in cell_ids]
    res = Int(resolution)
    md = merge(
        Dict{String,Any}("grid_name" => "a5", "resolution" => res, "indexing_scheme" => "u64-hex"),
        metadata,
    )
    # The cheap structural checks run first, in the inner constructor.
    lookup = A5Lookup(ids, res, md)
    if validate
        all(id -> A5Native.get_resolution(id) == res, ids) ||
            throw(ArgumentError("A5 cell ids do not all match resolution $res"))
    end
    return lookup
end

# The globe fast path. Both O(n) steps of the constructor above have to be
# skipped, not merely made faster: the comprehension would materialize every id
# and `validate`'s per-id pass would re-establish what the ordinal contract
# already guarantees. The inner constructor's own O(n) step is skipped by
# `strictly_increasing(::DGGSGlobeIds)` in `core/globe_ids.jl`; the endpoint
# resolution check still runs, and still rejects a level that disagrees with
# `resolution`, at O(1).
#
# This is load-bearing rather than an optimization: `DD.rebuild` defaults
# `data=l.data` and routes through the keyword constructor, so without this
# method *any* rebuild — `format`, `set`, a metadata change — would silently
# materialize the whole globe.
function A5Lookup(cell_ids::DGG.DGGSGlobeIds; resolution::Integer=cell_ids.level,
        metadata=Dict{String,Any}(), validate::Bool=false)
    # Rejected here rather than left to the inner constructor's endpoint
    # check, which reads a foreign id's bits and could conceivably agree by
    # accident — and because falling through to the generic method is the
    # materialization this whole path exists to prevent.
    cell_ids.system isa DGG.A5DGGS || throw(ArgumentError(
        "an A5Lookup cannot be built from a $(DGG.system_name(cell_ids.system)) globe"))
    res = Int(resolution)
    md = merge(
        Dict{String,Any}("grid_name" => "a5", "resolution" => res, "indexing_scheme" => "u64-hex"),
        metadata,
    )
    return A5Lookup(cell_ids, res, md)
end

A5Lookup(l::A5Lookup; data=l.data, resolution=l.resolution, metadata=l.metadata) =
    A5Lookup(data; resolution, metadata)

Lookups.metadata(l::A5Lookup) = l.metadata
DD.parent(l::A5Lookup) = l.data
DD.order(::A5Lookup) = Lookups.ForwardOrdered()
Base.size(l::A5Lookup) = size(l.data)
Base.length(l::A5Lookup) = length(l.data)
Base.getindex(l::A5Lookup, i::Int) = l.data[i]
Base.firstindex(l::A5Lookup) = firstindex(l.data)
Base.lastindex(l::A5Lookup) = lastindex(l.data)
Base.iterate(l::A5Lookup, state...) = iterate(l.data, state...)

function DD.rebuild(l::A5Lookup; data=l.data, metadata=l.metadata, kw...)
    return A5Lookup(data; resolution=l.resolution, metadata, kw...)
end

DD.format(dim::DD.Dimension{<:A5Lookup}, ::AbstractRange) = dim
DD.Dimensions.format(l::A5Lookup, ::Type, values, axis::AbstractRange) = l
Lookups.reducelookup(l::A5Lookup) = Lookups.NoLookup(Base.OneTo(1))
Lookups.show_properties(io::IO, l::A5Lookup) =
    print(io, " resolution: ", l.resolution, " ncells: ", length(l.data))
Lookups.show_properties(io::IO, mime, l::A5Lookup) = Lookups.show_properties(io, l)

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
_cell_position(l::A5Lookup, id::UInt64) = DGG.cell_position(l.data, id)

function Lookups.selectindices(l::A5Lookup, sel::Lookups.At)
    id = _to_id(Lookups.val(sel))
    i = _cell_position(l, id)
    isnothing(i) && throw(ArgumentError("A5 cell $(string(id; base=16)) is not present in the lookup"))
    return i
end

function _lonlat(value)
    if value isa Tuple && length(value) == 2
        return (Float64(value[1]), Float64(value[2]))
    elseif GI.isgeometry(value) && GI.geomtrait(value) isa GI.PointTrait
        return (Float64(GI.x(value)), Float64(GI.y(value)))
    end
    return nothing
end

function Lookups.selectindices(l::A5Lookup, sel::Lookups.Contains)
    value = Lookups.val(sel)
    point = _lonlat(value)
    if !isnothing(point)
        id = A5Native.lonlat_to_cell(point[1], point[2], l.resolution)
        i = _cell_position(l, id)
        return isnothing(i) ? Int[] : [i]
    end

    return DGG._query_positions(l, value, :center)
end

Lookups.selectindices(l::A5Lookup, sel::Touching) =
    DGG._query_positions(l, Lookups.val(sel), :touches)

lonlat_to_cell(lon::Real, lat::Real, resolution::Integer) =
    A5Native.lonlat_to_cell(lon, lat, resolution)

res0_cells() = A5Native.res0_cells()
cell_area(resolution::Integer) = A5Native.cell_area(resolution)
cell_boundary(id; segments=:auto) = A5Native.cell_boundary(_to_id(id); segments)
cell_center(id) = A5Native.cell_to_lonlat(_to_id(id))
cell_to_children(id, resolution::Union{Nothing,Integer}=nothing) =
    A5Native.cell_to_children(_to_id(id), resolution)
cell_to_parent(id, resolution::Union{Nothing,Integer}=nothing) =
    A5Native.cell_to_parent(_to_id(id), resolution)

cell_centers(l::A5Lookup) = [cell_center(id) for id in l.data]

function cell_polygon(id; segments=:auto)
    ring = cell_boundary(id; segments)
    isempty(ring) && throw(ArgumentError("A5 cell boundary is empty"))
    return GI.Polygon([GI.LinearRing(ring)])
end

cell_polygons(l::A5Lookup; segments=:auto) =
    [cell_polygon(id; segments) for id in l.data]

end
