# IGeo7Lookups.jl — DimensionalData integration for the clean IGEO7 core
# (design.md Section 9 task 8; spec/interface-contract.md "Lookup/Trees").
#
# Ported from the untainted IGeo7Lookups.jl of the prototype tree — repo wiring
# code, an allowed source under CLEANROOM.md — with the clean native module
# swapped in: every `IGeo7Native.*` call is now `IGeo7.*`. Struct name, exported
# names, method signatures and selector behavior are unchanged.
#
# Differences from the sibling, all forced by design.md Section 3 (the Z7
# UInt64 *is* the cell id):
#   * the `z7_*` "passthroughs" are direct calls, not a conversion layer;
#     `cell_to_z7` / `z7_to_cell` are validation-only identities;
#   * metadata advertises `indexing_scheme = "z7-u64"` with only `"z7-string"`
#     external, and the clean projection name.
#
# CLEANROOM: no sealed source was opened; the layer sees only the public
# native API listed in `IGeo7.jl`'s export block.

module IGeo7Lookups

using DimensionalData
using DimensionalData: Lookups
using GeoInterface
using GeometryOps
import SmallCollections
using SmallCollections: SmallVector

# The clean native core (this module's parent) and the package-wide helpers.
import ..IGeo7
import ..Helpers
# The package namespace, reached the way the kernel wirings reach it — this
# module is a grandchild of it, hence the third dot. It supplies the lookup
# supertype and `DGGSGlobeIds`.
import ...DiscreteGlobalGrids as DGG

const _to_id = Helpers.to_uint64_id
const _strictly_increasing = Helpers.strictly_increasing

const DD = DimensionalData
const GI = GeoInterface

export IGeo7Lookup,
    IGeo7,
    Touching,
    cell_area,
    cell_boundary,
    cell_center,
    cell_centers,
    cell_polygon,
    cell_polygons,
    cell_to_index,
    cell_to_children,
    cell_to_parent,
    cell_to_z7,
    index_to_cell,
    lonlat_to_cell,
    lonlat_to_index,
    lonlat_to_z7,
    num_cells,
    res0_cells,
    z7_base_cell,
    z7_child,
    z7_children,
    z7_digit,
    z7_from_hex,
    z7_from_string,
    z7_is_descendant,
    z7_is_pentagon,
    z7_parent,
    z7_resolution,
    z7_to_cell,
    z7_to_hex,
    z7_to_string

"""
    IGeo7Lookup(cell_ids; resolution, metadata=Dict{String,Any}(), validate=false)
    IGeo7Lookup(globe_ids::DGGSGlobeIds)

A `DimensionalData` lookup over IGEO7 cell ids (Z7 `UInt64`s, design Section
3) at a single resolution, sorted strictly ascending — which is also the
canonical full-world order, so a lookup's own order always agrees with
[`cell_to_index`](@ref).

The cheap structural invariants are enforced on *every* construction path
(this one, the raw positional form, `DD.rebuild` on a slice): sorted and unique
ids, `resolution` in range, and — since a Z7 id encodes its own resolution in
its digit slots — the first and last id checked against the declared one. That
endpoint check decodes through [`get_resolution`](@ref), so a *structurally*
invalid endpoint raises [`InvalidZ7Error`](@ref) rather than constructing.

`validate=true` additionally checks that every id is a structurally valid
cell *at* `resolution` — that is what rejects a mixed-resolution *interior*,
which the two endpoints cannot see. It is off by default because it is O(n)
decodes rather than two.

The second form is the globe-complete dimension:
`IGeo7Lookup(DGGSGlobeIds(IGEO7DGGS(), 8))` holds every res-8 cell in two words,
defaulting `resolution` to the ids' own level. There `validate` is a documented
no-op — a [`DGGSGlobeIds`](@ref) enumerates `ordinal_to_cell`, so its ids are
valid at `resolution` by construction, and running the per-id pass would be an
O(n) trap that establishes nothing.
"""
struct IGeo7Lookup{A<:AbstractVector{UInt64},M} <: DGG.AbstractDGGSLookup{UInt64}
    data::A
    resolution::Int
    metadata::M
    # The keyword form below used to be the only door these checks stood at,
    # with the default positional constructor open beside it — so
    # `IGeo7Lookup(res2_ids, 3, md)` built a lookup whose every `Contains`
    # selector answers "not stored" for cells that are stored. With the ids
    # sorted, checking the two endpoints means a whole vector at the wrong
    # resolution cannot hide behind them.
    function IGeo7Lookup(data::AbstractVector{UInt64}, resolution::Integer, metadata)
        _strictly_increasing(data) ||
            throw(ArgumentError("IGEO7 cell ids must be sorted ascending and unique"))
        res = Int(resolution)
        0 <= res <= IGeo7.MAX_RESOLUTION ||
            throw(ArgumentError("IGEO7 resolution must be in 0:$(IGeo7.MAX_RESOLUTION)"))
        isempty(data) || (IGeo7.get_resolution(first(data)) == res &&
                          IGeo7.get_resolution(last(data)) == res) ||
            throw(ArgumentError("IGEO7 cell ids do not all match resolution $res"))
        return new{typeof(data),typeof(metadata)}(data, res, metadata)
    end
end

function IGeo7Lookup(cell_ids; resolution::Integer, metadata=Dict{String,Any}(), validate::Bool=false)
    ids = [_to_id(id) for id in cell_ids]
    res = Int(resolution)
    user_metadata = metadata === Lookups.NoMetadata ? Dict{String,Any}() : metadata
    md = merge(
        Dict{String,Any}(
            "grid_name" => "igeo7",
            "resolution" => res,
            "indexing_scheme" => "z7-u64",
            "external_indexing_schemes" => ["z7-string"],
            "projection" => "snyder-isea (standard ISEA placement)",
        ),
        user_metadata,
    )
    # The cheap structural checks run first, in the inner constructor.
    lookup = IGeo7Lookup(ids, res, md)
    if validate
        all(id -> IGeo7.is_valid_cell(id) && IGeo7.get_resolution(id) == res, ids) ||
            throw(ArgumentError("IGEO7 cell ids do not all match resolution $res"))
    end
    return lookup
end

# The globe fast path. Both O(n) steps of the constructor above have to be
# skipped, not merely made faster: the comprehension would materialize every id
# and `validate`'s per-id pass would re-establish what the ordinal contract
# already guarantees. The inner constructor's own O(n) step is skipped by
# `Helpers.strictly_increasing(::DGGSGlobeIds)` in `core/globe_ids.jl`; the
# endpoint resolution check still runs, and still rejects a level that disagrees
# with `resolution`, at O(1).
#
# This is load-bearing rather than an optimization: `DD.rebuild` defaults
# `data=l.data` and routes through the keyword constructor, so without this
# method *any* rebuild — `format`, `set`, a metadata change — would silently
# materialize the whole globe.
function IGeo7Lookup(cell_ids::DGG.DGGSGlobeIds; resolution::Integer=cell_ids.level,
        metadata=Dict{String,Any}(), validate::Bool=false)
    # Rejected here rather than left to the inner constructor's endpoint
    # check, which reads a foreign id's bits and could conceivably agree by
    # accident — and because falling through to the generic method is the
    # materialization this whole path exists to prevent.
    cell_ids.system isa DGG.IGEO7DGGS || throw(ArgumentError(
        "an IGeo7Lookup cannot be built from a $(DGG.system_name(cell_ids.system)) globe"))
    res = Int(resolution)
    user_metadata = metadata === Lookups.NoMetadata ? Dict{String,Any}() : metadata
    md = merge(
        Dict{String,Any}(
            "grid_name" => "igeo7",
            "resolution" => res,
            "indexing_scheme" => "z7-u64",
            "external_indexing_schemes" => ["z7-string"],
            "projection" => "snyder-isea (standard ISEA placement)",
        ),
        user_metadata,
    )
    return IGeo7Lookup(cell_ids, res, md)
end

IGeo7Lookup(l::IGeo7Lookup; data=l.data, resolution=l.resolution, metadata=l.metadata) =
    IGeo7Lookup(data; resolution, metadata)

Lookups.metadata(l::IGeo7Lookup) = l.metadata
DD.parent(l::IGeo7Lookup) = l.data
DD.order(::IGeo7Lookup) = Lookups.ForwardOrdered()
Base.size(l::IGeo7Lookup) = size(l.data)
Base.length(l::IGeo7Lookup) = length(l.data)
Base.getindex(l::IGeo7Lookup, i::Int) = l.data[i]
Base.firstindex(l::IGeo7Lookup) = firstindex(l.data)
Base.lastindex(l::IGeo7Lookup) = lastindex(l.data)
Base.iterate(l::IGeo7Lookup, state...) = iterate(l.data, state...)

# The dimension tuple is part of an AbstractDimArray's type, so this dispatch
# cannot affect an ordinary vector or a raster with Cartesian dimensions.
function Base.eachindex(
        A::DD.AbstractDimArray{T,1,D}
    ) where {T,D<:Tuple{<:DD.Dimension{<:IGeo7Lookup}}}
    lookup = DD.val(only(DD.dims(A)))
    return IGeo7.IGEO7Indices(DD.parent(lookup), lookup.resolution)
end

function IGeo7.cell_to_position(
        A::DD.AbstractDimArray{T,1,D},
        index::IGeo7.IGEO7Index,
    ) where {T,D<:Tuple{<:DD.Dimension{<:IGeo7Lookup}}}
    lookup = DD.val(only(DD.dims(A)))
    position = DGG.cell_position(DD.parent(lookup), index.id)
    position === nothing && throw(BoundsError(A, index))
    return position
end

function IGeo7.position_to_cell(
        A::DD.AbstractDimArray{T,1,D},
        position::Integer,
    ) where {T,D<:Tuple{<:DD.Dimension{<:IGeo7Lookup}}}
    checkbounds(A, position)
    lookup = DD.val(only(DD.dims(A)))
    return IGeo7.IGEO7Index(@inbounds DD.parent(lookup)[position])
end

function Base.getindex(
        A::DD.AbstractDimArray{T,1,D},
        index::IGeo7.IGEO7Index,
    ) where {T,D<:Tuple{<:DD.Dimension{<:IGeo7Lookup}}}
    return A[IGeo7.cell_to_position(A, index)]
end

function Base.setindex!(
        A::DD.AbstractDimArray{T,1,D},
        value,
        index::IGeo7.IGEO7Index,
    ) where {T,D<:Tuple{<:DD.Dimension{<:IGeo7Lookup}}}
    return setindex!(A, value, IGeo7.cell_to_position(A, index))
end

function Base.checkbounds(
        ::Type{Bool},
        A::DD.AbstractDimArray{T,1,D},
        index::IGeo7.IGEO7Index,
    ) where {T,D<:Tuple{<:DD.Dimension{<:IGeo7Lookup}}}
    lookup = DD.val(only(DD.dims(A)))
    return DGG.cell_position(DD.parent(lookup), index.id) !== nothing
end

function IGeo7.neighbors(
        A::DD.AbstractDimArray{T,1,D},
        index::IGeo7.IGEO7Index,
    ) where {T,D<:Tuple{<:DD.Dimension{<:IGeo7Lookup}}}
    checkbounds(Bool, A, index) || throw(BoundsError(A, index))
    out = SmallVector{6,IGeo7.IGEO7Index}()
    for neighbor in IGeo7.neighbors(DGG.IGEO7DGGS(), index)
        checkbounds(Bool, A, neighbor) ||
            continue
        out = SmallCollections.push(out, neighbor)
    end
    return out
end

function IGeo7.cellbearing(
        A::DD.AbstractDimArray{T,1,D},
        from::IGeo7.IGEO7Index,
        to::IGeo7.IGEO7Index,
    ) where {T,D<:Tuple{<:DD.Dimension{<:IGeo7Lookup}}}
    checkbounds(Bool, A, from) || throw(BoundsError(A, from))
    checkbounds(Bool, A, to) || throw(BoundsError(A, to))
    from == to && return 0.0

    from_lon, from_lat = IGeo7.cell_center(from.id)
    to_lon, to_lat = IGeo7.cell_center(to.id)
    phi_from, phi_to = deg2rad(from_lat), deg2rad(to_lat)
    delta_lon = deg2rad(to_lon - from_lon)
    east = sin(delta_lon) * cos(phi_to)
    north = cos(phi_from) * sin(phi_to) -
            sin(phi_from) * cos(phi_to) * cos(delta_lon)
    return mod(rad2deg(atan(east, north)), 360.0)
end

function IGeo7.celldistance(
        A::DD.AbstractDimArray{T,1,D},
        from::IGeo7.IGEO7Index,
        to::IGeo7.IGEO7Index,
    ) where {T,D<:Tuple{<:DD.Dimension{<:IGeo7Lookup}}}
    checkbounds(Bool, A, from) || throw(BoundsError(A, from))
    checkbounds(Bool, A, to) || throw(BoundsError(A, to))
    from == to && return 0.0
    lookup = DD.val(only(DD.dims(A)))
    pfrom = DGG.cell_center(DGG.IGEO7DGGS(), lookup.resolution, from.id)
    pto = DGG.cell_center(DGG.IGEO7DGGS(), lookup.resolution, to.id)
    angle = GeometryOps.UnitSpherical.spherical_distance(pfrom, pto)
    return Float64(angle * IGeo7.R_AUTHALIC)
end

function IGeo7.cellarea(
        A::DD.AbstractDimArray{T,1,D},
        index::IGeo7.IGEO7Index,
    ) where {T,D<:Tuple{<:DD.Dimension{<:IGeo7Lookup}}}
    checkbounds(Bool, A, index) || throw(BoundsError(A, index))
    return IGeo7.cell_area(index.id)
end

function IGeo7.edges(
        A::DD.AbstractDimArray{T,1,D},
    ) where {T,D<:Tuple{<:DD.Dimension{<:IGeo7Lookup}}}
    return IGeo7.edges(eachindex(A))
end

function DD.rebuild(l::IGeo7Lookup; data=l.data, metadata=l.metadata, kw...)
    return IGeo7Lookup(data; resolution=l.resolution, metadata, kw...)
end

DD.format(dim::DD.Dimension{<:IGeo7Lookup}, ::AbstractRange) = dim
DD.Dimensions.format(l::IGeo7Lookup, ::Type, values, axis::AbstractRange) = l
Lookups.reducelookup(l::IGeo7Lookup) = Lookups.NoLookup(Base.OneTo(1))
Lookups.show_properties(io::IO, l::IGeo7Lookup) =
    print(io, " resolution: ", l.resolution, " ncells: ", length(l.data))
Lookups.show_properties(io::IO, mime, l::IGeo7Lookup) = Lookups.show_properties(io, l)

"""
    Touching(geometry)

Selector for every cell whose (spherical) polygon intersects `geometry`
(`Contains` selects by cell *center* instead). Answered by the package-level
spherical tree descent (`src/core/lookup_ops.jl`): `geometry` may cross the
antimeridian or enclose a pole, and its ring edges are great-circle arcs —
see [`DiscreteGlobalGrids.zonal`](@ref) for the full geometry conventions.
"""
struct Touching{G} <: Lookups.ArraySelector{G}
    val::G
end

Lookups.val(sel::Touching) = sel.val

# The lookup position holding a cell id, or `nothing`. One line here because
# the two branches — binary search over stored ids, `cell_to_ordinal` over a
# globe — are one generic method pair in core (`core/lookups.jl`), chosen by
# the type of the id vector and so invisible at every call site below.
_cell_position(l::IGeo7Lookup, id::UInt64) = DGG.cell_position(l.data, id)

Lookups.selectindices(l::IGeo7Lookup, sel::Lookups.StandardIndices) = sel

function Lookups.selectindices(l::IGeo7Lookup, sel::Lookups.At)
    id = _to_id(Lookups.val(sel))
    i = _cell_position(l, id)
    isnothing(i) && throw(ArgumentError("IGEO7 cell $(string(id; base=16)) is not present in the lookup"))
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

function Lookups.selectindices(l::IGeo7Lookup, sel::Lookups.Contains)
    value = Lookups.val(sel)
    point = _lonlat(value)
    if !isnothing(point)
        id = IGeo7.lonlat_to_cell(point[1], point[2], l.resolution)
        i = _cell_position(l, id)
        return isnothing(i) ? Int[] : [i]
    end

    return DGG._query_positions(l, value, :center)
end

Lookups.selectindices(l::IGeo7Lookup, sel::Touching) =
    DGG._query_positions(l, Lookups.val(sel), :touches)

# --- native API delegation -------------------------------------------------
# One method per function the sibling forwards to its native module. `kw...`
# carries the clean core's `orientation` keyword (design Section 5) through
# without changing any sibling-compatible call site.
#
# These are forwarders, not new API: semantics, argument meaning and error
# behavior are the core's, and the documentation for each is the docstring on
# the same-named `IGeo7` function (`src/IGeo7.jl`'s export block).
# They are deliberately left undocumented individually so this file stays a
# line-for-line port of the untainted sibling it came from.

lonlat_to_cell(lon::Real, lat::Real, resolution::Integer; kw...) =
    IGeo7.lonlat_to_cell(lon, lat, resolution; kw...)
lonlat_to_index(lon::Real, lat::Real, resolution::Integer; kw...) =
    IGeo7.lonlat_to_index(lon, lat, resolution; kw...)
lonlat_to_z7(lon::Real, lat::Real, resolution::Integer; kw...) =
    IGeo7.lonlat_to_z7(lon, lat, resolution; kw...)

num_cells(resolution::Integer) = IGeo7.num_cells(resolution)
res0_cells() = IGeo7.res0_cells()
cell_area(id) = IGeo7.cell_area(_to_id(id))
cell_boundary(id; closed_ring::Bool=true, kw...) =
    IGeo7.cell_boundary(_to_id(id); closed_ring, kw...)
cell_center(id; kw...) = IGeo7.cell_center(_to_id(id); kw...)
cell_to_index(id) = IGeo7.cell_to_index(_to_id(id))
cell_to_z7(id) = IGeo7.cell_to_z7(_to_id(id))
cell_to_children(id) = IGeo7.cell_to_children(_to_id(id))
cell_to_children(id, resolution::Union{Nothing,Integer}) =
    isnothing(resolution) ? IGeo7.cell_to_children(_to_id(id)) :
    IGeo7.cell_to_children(_to_id(id), resolution)
cell_to_parent(id) = IGeo7.cell_to_parent(_to_id(id))
cell_to_parent(id, resolution::Integer) =
    IGeo7.cell_to_parent(_to_id(id), resolution)
index_to_cell(index::Integer, resolution::Integer) =
    IGeo7.index_to_cell(index, resolution)

z7_base_cell(z7::Unsigned) = IGeo7.z7_base_cell(z7)
z7_child(z7::Unsigned, digit::Integer) = IGeo7.z7_child(z7, digit)
z7_children(z7::Unsigned) = IGeo7.z7_children(z7)
z7_digit(z7::Unsigned, resolution::Integer) = IGeo7.z7_digit(z7, resolution)
z7_from_hex(text::AbstractString) = IGeo7.z7_from_hex(text)
z7_from_string(text::AbstractString) = IGeo7.z7_from_string(text)
z7_is_descendant(z7::Unsigned, ancestor::Unsigned) =
    IGeo7.z7_is_descendant(z7, ancestor)
z7_is_pentagon(z7::Unsigned) = IGeo7.z7_is_pentagon(z7)
z7_parent(z7::Unsigned, resolution::Integer) =
    IGeo7.z7_parent(z7, resolution)
z7_parent(z7::Unsigned) = IGeo7.z7_parent(z7)
z7_resolution(z7::Unsigned) = IGeo7.z7_resolution(z7)
z7_to_cell(z7::Unsigned) = IGeo7.z7_to_cell(z7)
z7_to_hex(z7::Unsigned; kw...) = IGeo7.z7_to_hex(z7; kw...)
z7_to_string(z7::Unsigned) = IGeo7.z7_to_string(z7)

# --- geometry for the lookup as a whole ------------------------------------

"`(lon, lat)` center of every cell in the lookup, in lookup order."
cell_centers(l::IGeo7Lookup) = [cell_center(id) for id in l.data]

"""
    cell_polygon(id) -> GI.Polygon

The cell's closed boundary ring as a GeoInterface polygon in `(lon, lat)`
degrees. Rings are counter-clockwise seen from outside the sphere.
"""
function cell_polygon(id)
    ring = cell_boundary(id; closed_ring=true)
    isempty(ring) && throw(ArgumentError("IGEO7 cell boundary is empty"))
    return GI.Polygon([GI.LinearRing(ring)])
end

cell_polygons(l::IGeo7Lookup) = [cell_polygon(id) for id in l.data]

end
