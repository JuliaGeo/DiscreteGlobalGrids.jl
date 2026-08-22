# The Zarr half of the ancestor-subzone layout (see `src/io/subzones.jl`, which
# owns the arithmetic and the attribute vocabulary): Zarr calls, the chunk plan,
# and the DiskArrays wrapper that makes the two-dimensional store look like the
# one-dimensional cell axis it stands for.
#
#   create   subzonestore(dest, sys, level; ancestor_level, layers, overviews)
#   fill     dggwrite!(store, ancestor, values)   |   dggwrite!(store, cube)
#   one shot dggwrite(dest, cube; layout = :subzones, ancestor_level)
#   read     dggread(store)  ->  a cell-axis DimStack over `SubzoneCellArray`s
#
# A submodule, like the write half, so its helper names cannot collide with the
# extension's shared namespace.

module DGGSZarrSubzones

import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: AbstractCellIndex, AbstractHierarchicalGridSystem,
    CellLookup, Cells, ChunkedCellLookup, DGGSFormatError, SubzoneLayout,
    ANCESTOR_COORDINATE, ANCESTOR_DIMENSION, SUBZONE_DIMENSION,
    cellindex, columnindex, columnlength, columnpositions, level, ncells,
    positionindex, rawid, subzone_attrs, subzone_cellvector, subzone_columns,
    subzone_layout, subzone_runs, system, with_store_context
import ..DiscreteGlobalGridsZarrExt
using ..DiscreteGlobalGridsZarrExt: ARRAY_DIMENSIONS, selectvars, storeidentifier
import DimensionalData as DD
import DiskArrays
import Zarr

const OVERVIEW_LEVELS_ATTRIBUTE = "overview_levels"
const OVERVIEWS_ATTRIBUTE = "overviews"
const OVERVIEW_LEVEL_ATTRIBUTE = "dggs_overview_level"
const OVERVIEW_METHOD_ATTRIBUTE = "dggs_overview_method"

_overviewname(name::AbstractString, lev::Integer) = string(name, "_ovr", lev)

# ===========================================================================
# Creating a store
# ===========================================================================

"""
    SubzoneStore

An open ancestor-subzone store: the Zarr group, the [`SubzoneLayout`](@ref) its
attributes declare, and the data arrays by name.

Handed back by [`subzonestore`](@ref) and filled by [`dggwrite!`](@ref). It owns
no buffer and holds no lock: a column is a chunk is a file, so two tasks writing
two columns share nothing, and the handle may be used from all of them.
"""
struct SubzoneOverview{L<:SubzoneLayout}
    layout::L
    method::Symbol
    arrays::Dict{String,Zarr.ZArray}
end

struct SubzoneStore{L<:SubzoneLayout}
    group::Zarr.ZGroup
    layout::L
    arrays::Dict{String,Zarr.ZArray}
    overviews::Vector{SubzoneOverview}
    identifier::String
end

DGG.system(s::SubzoneStore) = system(s.layout)
DGG.level(s::SubzoneStore) = level(s.layout)
Base.keys(s::SubzoneStore) = sort!(collect(keys(s.arrays)))

Base.show(io::IO, s::SubzoneStore) = print(io, "SubzoneStore(", s.identifier, ", ",
    s.layout, ", ", join(keys(s), ", "), ")")

"""
    subzonestore(dest, system, level; ancestor_level, layers, overviews = Int[],
                 capacity = nothing, fill_value = NaN, ancestor_coordinate = true,
                 overview_method = :center,
                 attrs = Dict{String,Any}(), compressor = Zarr.BloscCompressor())
        -> SubzoneStore
    subzonestore(dest) -> SubzoneStore

Create an ancestor-subzone store, or reopen one for writing.

The store is a Zarr v2 group of `(ancestor, subzone)` arrays — `(capacity,
ncolumns)` in Julia's own order, chunked `(capacity, 1)`, so one chunk is one
level-`ancestor_level` subtree. Nothing is written into them here: the columns
are filled afterwards by [`dggwrite!`](@ref), in any order and from any number
of tasks, and a column nobody writes stays a chunk that was never stored and
reads back as `fill_value`.

  - `layers` names the data variables and their element types: `"elevation" =>
    Float32`, an iterable of such pairs, or a `NamedTuple`/`Dict` of them.
  - `capacity` is the row extent, which defaults to the measured longest subtree
    ([`subzone_capacity`](@ref)). Pass it where it is known — one pass over a
    level-6 ancestor grid is 1.2 million `descendant_range` calls.
  - `fill_value` is what an unwritten cell reads back as, `NaN` by default,
    which needs a floating-point layer; a layer of another element type has to
    be given one.
  - `ancestor_coordinate` writes the level-`ancestor_level` ids as a
    one-dimensional coordinate array. The column axis is IMPLICIT — position is
    the ancestor — and this reader never consults it; it is there so that a
    generic xdggs reader sees the ancestor axis for the cell axis it is.
  - `attrs` are the producer's own group attributes, which the layout's are
    stamped over.
  - `overviews` is a strictly increasing list of levels from `ancestor_level`
    through `level - 1`. Each gets a `<layer>_ovr<level>` array in the same
    ancestor-column layout. `overview_method = :center` samples the first base
    descendant of every overview cell; it is the only method implemented yet.

Reopening reads the layout back out of the attributes and checks it against the
arrays, so a mistyped path is refused rather than half-written into. `dest` is a
local directory path or an open writeable `Zarr.ZGroup`; a remote URL is refused.
"""
function DGG.subzonestore(dest::AbstractString, sys::AbstractHierarchicalGridSystem,
    lev::Integer; kw...)
    path = String(dest)
    _reject_remote(path)
    return _create(path, (attrs) -> Zarr.zgroup(path; attrs=attrs), sys, lev; kw...)
end

function DGG.subzonestore(dest::Zarr.ZGroup, sys::AbstractHierarchicalGridSystem,
    lev::Integer; kw...)
    _writeable(dest)
    return _create(storeidentifier(dest), (attrs) -> _stamp(dest, attrs), sys, lev; kw...)
end

function _create(identifier, opengroup, sys, lev; ancestor_level::Integer,
    layers, capacity::Union{Integer,Nothing}=nothing, fill_value=NaN,
    ancestor_coordinate::Bool=true, attrs=Dict{String,Any}(),
    overviews=Int[], overview_method::Symbol=:center,
    compressor=Zarr.BloscCompressor())

    specs = _layerspecs(layers)
    layout = SubzoneLayout(sys, lev, ancestor_level; capacity=capacity)
    overview_method === :center || throw(ArgumentError(
        "the primitive overview writer implements `overview_method = :center` only, " *
        "not $(repr(overview_method))."))
    overview_levels = _overviewlevels(overviews, layout)
    overview_layouts = SubzoneLayout[
        SubzoneLayout(sys, lo, ancestor_level; gridname=layout.gridname)
        for lo in overview_levels]
    _checkoverviewnames(specs, overview_levels)
    coordinate = ancestor_coordinate ? ANCESTOR_COORDINATE : nothing

    return with_store_context(identifier) do
        groupattrs = _attrs(attrs)
        fills = Dict{String,Any}(name => _fillvalue(fill_value, T, name)
                                 for (name, T) in specs)
        # One `fill_value` spelling for the whole store where the layers agree on
        # it, which they do whenever they share an element type; the attribute is
        # documentation, and the arrays each carry their own for real.
        merge!(groupattrs, subzone_attrs(layout;
            variables=String[name for (name, _) in specs], coordinate=coordinate,
            fill_value=_fillspelling(fills)))
        _overviewattrs!(groupattrs, overview_layouts, specs, overview_method)

        group = opengroup(groupattrs)
        arrays = Dict{String,Zarr.ZArray}()
        for (name, T) in specs
            arrays[name] = Zarr.zcreate(T, group, name, layout.capacity, layout.ncolumns;
                chunks=(layout.capacity, 1), fill_value=fills[name],
                fill_as_missing=false, compressor=compressor,
                attrs=Dict{String,Any}(ARRAY_DIMENSIONS =>
                    [ANCESTOR_DIMENSION, SUBZONE_DIMENSION]))
        end
        overview_handles = SubzoneOverview[]
        for overview_layout in overview_layouts
            overview_arrays = Dict{String,Zarr.ZArray}()
            for (name, T) in specs
                stored = _overviewname(name, overview_layout.level)
                overview_arrays[name] = Zarr.zcreate(T, group, stored,
                    overview_layout.capacity, overview_layout.ncolumns;
                    chunks=(overview_layout.capacity, 1), fill_value=fills[name],
                    fill_as_missing=false, compressor=compressor,
                    attrs=Dict{String,Any}(
                        ARRAY_DIMENSIONS => [ANCESTOR_DIMENSION, SUBZONE_DIMENSION],
                        OVERVIEW_LEVEL_ATTRIBUTE => overview_layout.level,
                        OVERVIEW_METHOD_ATTRIBUTE => String(overview_method)))
            end
            push!(overview_handles,
                SubzoneOverview(overview_layout, overview_method, overview_arrays))
        end
        coordinate === nothing || _writecoordinate(group, layout, coordinate)
        Zarr.consolidate_metadata(group)
        return SubzoneStore(group, layout, arrays, overview_handles, String(identifier))
    end
end

function _overviewlevels(overviews, layout::SubzoneLayout)
    raw = try
        collect(overviews)
    catch
        throw(ArgumentError("`overviews` is a list of integer levels."))
    end
    all(x -> x isa Integer && !(x isa Bool), raw) || throw(ArgumentError(
        "`overviews` is a list of integer levels, not $(repr(raw))."))
    levels = Int[x for x in raw]
    issorted(levels) && allunique(levels) || throw(ArgumentError(
        "`overviews` must be sorted and unique (strictly increasing), not " *
        "$(repr(levels))."))
    for lo in levels
        layout.ancestor_level <= lo < layout.level || throw(ArgumentError(
            "overview level $lo is outside $(layout.ancestor_level):$(layout.level - 1): " *
            "an overview reuses the level-$(layout.ancestor_level) ancestor columns " *
            "and must be coarser than the level-$(layout.level) base data."))
    end
    return levels
end

function _checkoverviewnames(specs, levels)
    base = Set(first.(specs))
    generated = String[_overviewname(name, lo) for lo in levels for (name, _) in specs]
    collisions = sort!(collect(intersect(base, Set(generated))))
    length(generated) == length(unique(generated)) && isempty(collisions) && return nothing
    throw(ArgumentError("the `<layer>_ovr<level>` overview names collide with layer " *
        "names: " * join(collisions, ", ") * ". Rename those layers."))
end

function _overviewattrs!(attrs, layouts, specs, method)
    block = attrs["dggs"][DGG.SUBZONE_BLOCK]
    block[OVERVIEW_LEVELS_ATTRIBUTE] = Int[l.level for l in layouts]
    block[OVERVIEWS_ATTRIBUTE] = Any[
        Dict{String,Any}(
            "level" => l.level,
            "method" => String(method),
            "subzone_count" => l.capacity,
            "chunk_shape" => Any[1, l.capacity],
            "variables" => Dict{String,Any}(
                name => _overviewname(name, l.level) for (name, _) in specs))
        for l in layouts]
    return attrs
end

# The column ids, as the one array in the store that is metadata rather than
# data: it is written once, whole, and never read back by this package.
# `grid_name`/`level` are the xdggs spelling of what it is — a level-`La` dense
# cell axis — which is exactly true of it and costs nothing to say.
function _writecoordinate(group, layout::SubzoneLayout, name::AbstractString)
    agrid = layout.ancestorgrid
    ids = [rawid(cellindex(agrid, i)) for i in 1:layout.ncolumns]
    z = Zarr.zcreate(eltype(ids), group, String(name), length(ids);
        chunks=(min(length(ids), 1 << 16),),
        attrs=Dict{String,Any}(ARRAY_DIMENSIONS => [ANCESTOR_DIMENSION],
            "grid_name" => layout.gridname, "level" => layout.ancestor_level))
    z[:] = ids
    return z
end

function DGG.subzonestore(dest::AbstractString)
    path = String(dest)
    _reject_remote(path)
    return DGG.subzonestore(Zarr.zopen(path, "w"))
end

function DGG.subzonestore(group::Zarr.ZGroup)
    _writeable(group)
    identifier = String(storeidentifier(group))
    return with_store_context(identifier) do
        attrs = Dict{String,Any}(group.attrs)
        layout = subzone_layout(attrs; store=identifier)
        declarations = _overviewdeclarations(attrs, layout; store=identifier)
        overview_names = Set(stored for d in declarations for stored in values(d.variables))
        arrays = Dict{String,Zarr.ZArray}()
        for (name, z) in group.arrays
            _issubzonearray(z) || continue
            name in overview_names && continue
            _checkshape(z, layout, name, identifier)
            arrays[name] = z
        end
        isempty(arrays) && throw(DGGSFormatError(check=:no_data_variables,
            store=identifier, observed=sort!(collect(keys(group.arrays))),
            detail="this ancestor-subzone store holds no array on the " *
                   "($ANCESTOR_DIMENSION, $SUBZONE_DIMENSION) dimensions."))
        overview_handles = SubzoneOverview[]
        for declaration in declarations
            Set(keys(declaration.variables)) == Set(keys(arrays)) || throw(DGGSFormatError(
                check=:overview_variable_mismatch, store=identifier,
                declared=sort!(collect(keys(declaration.variables))),
                observed=sort!(collect(keys(arrays))),
                detail="the level-$(declaration.layout.level) overview must name one " *
                       "array for every base variable."))
            overview_arrays = Dict{String,Zarr.ZArray}()
            for (name, stored) in declaration.variables
                haskey(group.arrays, stored) || throw(DGGSFormatError(
                    check=:missing_overview_variable, store=identifier,
                    declared=stored, observed=sort!(collect(keys(group.arrays))),
                    detail="the level-$(declaration.layout.level) overview of `$name` " *
                           "is declared as array `$stored`, which is absent."))
                z = group.arrays[stored]
                _issubzonearray(z) || throw(DGGSFormatError(
                    check=:invalid_overview_dimensions, store=identifier,
                    declared=[ANCESTOR_DIMENSION, SUBZONE_DIMENSION],
                    observed=get(z.attrs, ARRAY_DIMENSIONS, nothing),
                    detail="overview array `$stored` does not use the subzone dimensions."))
                _checkshape(z, declaration.layout, stored, identifier)
                overview_arrays[name] = z
            end
            push!(overview_handles, SubzoneOverview(
                declaration.layout, declaration.method, overview_arrays))
        end
        return SubzoneStore(group, layout, arrays, overview_handles, identifier)
    end
end

struct OverviewDeclaration{L<:SubzoneLayout}
    layout::L
    method::Symbol
    variables::Dict{String,String}
end

function _overviewdeclarations(attrs, base::SubzoneLayout; store::AbstractString="")
    block = attrs["dggs"][DGG.SUBZONE_BLOCK]
    levels_raw = get(block, OVERVIEW_LEVELS_ATTRIBUTE, Any[])
    entries = get(block, OVERVIEWS_ATTRIBUTE, Any[])
    levels_raw isa AbstractVector || _overviewformaterror(store,
        "`$OVERVIEW_LEVELS_ATTRIBUTE` is not a list.", levels_raw)
    entries isa AbstractVector || _overviewformaterror(store,
        "`$OVERVIEWS_ATTRIBUTE` is not a list.", entries)
    levels = Int[_overviewint(x, OVERVIEW_LEVELS_ATTRIBUTE, store) for x in levels_raw]
    issorted(levels) && allunique(levels) || _overviewformaterror(store,
        "overview levels must be strictly increasing.", levels)
    length(entries) == length(levels) || _overviewformaterror(store,
        "overview levels and declarations have different lengths.",
        (levels=length(levels), declarations=length(entries)))

    declarations = OverviewDeclaration[]
    for (k, entry) in pairs(entries)
        entry isa AbstractDict || _overviewformaterror(store,
            "overview declaration $k is not an object.", entry)
        lo = _overviewint(get(entry, "level", nothing), "overview level", store)
        lo == levels[k] || _overviewformaterror(store,
            "overview declaration $k names level $lo instead of $(levels[k]).", entry)
        base.ancestor_level <= lo < base.level || _overviewformaterror(store,
            "overview level $lo is not in $(base.ancestor_level):$(base.level - 1).", lo)
        method_raw = get(entry, "method", nothing)
        method_raw == "center" || throw(DGGSFormatError(
            check=:unsupported_overview_method, store=String(store),
            declared=method_raw, observed="center",
            detail="this primitive reader implements center-sampled overviews only."))
        capacity = _overviewint(get(entry, "subzone_count", nothing),
            "overview subzone_count", store)
        overview_layout = SubzoneLayout(base.system, lo, base.ancestor_level;
            gridname=base.gridname, capacity=capacity)
        get(entry, "chunk_shape", Any[1, capacity]) == [1, capacity] ||
            _overviewformaterror(store,
                "the level-$lo overview chunk shape is not one whole column.",
                get(entry, "chunk_shape", nothing))
        vars_raw = get(entry, "variables", nothing)
        vars_raw isa AbstractDict || _overviewformaterror(store,
            "the level-$lo overview variables are not an object.", vars_raw)
        variables = Dict{String,String}()
        for (name, stored) in vars_raw
            name isa AbstractString && stored isa AbstractString || _overviewformaterror(store,
                "the level-$lo overview variable map must contain string names.",
                name => stored)
            variables[String(name)] = String(stored)
        end
        push!(declarations, OverviewDeclaration(overview_layout, :center, variables))
    end
    return declarations
end

function _overviewint(x, field, store)
    x isa Integer && !(x isa Bool) && return Int(x)
    _overviewformaterror(store, "`$field` is not an integer.", x)
end

@noinline function _overviewformaterror(store, detail, observed)
    throw(DGGSFormatError(check=:invalid_overview_metadata, store=String(store),
        observed=observed, detail=detail))
end

# A store's data arrays are the ones on the layout's two dimensions, found by
# their dimension names rather than by the `variables` list in the attributes:
# the list is what the writer said, and the dimensions are what the array is.
_issubzonearray(z::Zarr.ZArray) =
    get(z.attrs, ARRAY_DIMENSIONS, nothing) ==
    [ANCESTOR_DIMENSION, SUBZONE_DIMENSION]

function _checkshape(z, layout::SubzoneLayout, name, identifier)
    size(z) == (layout.capacity, layout.ncolumns) && return nothing
    throw(DGGSFormatError(check=:subzone_shape_mismatch, store=String(identifier),
        declared=(layout.ncolumns, layout.capacity), observed=reverse(size(z)),
        detail="array `$name` is $(join(reverse(size(z)), "x")) where the layout " *
               "declares $(layout.ncolumns) columns of $(layout.capacity) subzones."))
end

_writeable(g::Zarr.ZGroup) = g.writeable || throw(ArgumentError(
    "a subzone store is written into a writeable group; this one was opened read-only."))

const REMOTE_SCHEMES = ("gs://", "s3://", "http://", "https://", "az://", "abfs://")

@noinline function _reject_remote(path)
    for scheme in REMOTE_SCHEMES
        startswith(path, scheme) || continue
        throw(ArgumentError("subzonestore writes local directory stores only; " *
                            "$(repr(path)) names a $(rstrip(scheme, ['/', ':'])) store. " *
                            "Write locally and upload, or open the remote group " *
                            "yourself and pass the ZGroup."))
    end
    return nothing
end

function _stamp(g::Zarr.ZGroup, attrs)
    isempty(g.arrays) || throw(DGGSFormatError(check=:destination_not_empty,
        observed=sort!(collect(keys(g.arrays))),
        detail="this group already holds arrays, and a subzone store is created " *
               "into an empty group. Reopen it with `subzonestore(group)` to add " *
               "columns to it, or write to a new group."))
    merge!(g.attrs, attrs)
    Zarr.writeattrs(g.zarr_format, g.storage, g.path, g.attrs)
    return g
end

_attrs(x) = Dict{String,Any}()
_attrs(md::AbstractDict) = Dict{String,Any}(String(k) => deepcopy(v) for (k, v) in md)
_attrs(md::NamedTuple) = Dict{String,Any}(String(k) => deepcopy(v) for (k, v) in pairs(md))

# --- the layer specification ------------------------------------------------

_layerspecs(p::Pair) = _layerspecs((p,))
_layerspecs(nt::NamedTuple) = _layerspecs(collect(pairs(nt)))

function _layerspecs(layers)
    out = Pair{String,DataType}[]
    for spec in layers
        spec isa Pair || throw(ArgumentError(
            "a layer is a `name => eltype` pair, and $(repr(spec)) is not one; " *
            "`layers = (\"elevation\" => Float32,)` names one layer."))
        name, T = spec
        T isa Type || throw(ArgumentError(
            "layer $(repr(String(name)))'s element type is $(repr(T)), which is not a type."))
        push!(out, String(name) => T)
    end
    isempty(out) && throw(ArgumentError("a subzone store has at least one layer."))
    sort!(out; by=first)
    allunique(first.(out)) || throw(ArgumentError(
        "two layers of this store would share a name: " *
        join(sort!(first.(out)), ", ") * "."))
    return out
end

# NaN is the default and needs a float to be one. Anything else the caller
# named is converted to the layer's own type, so an integer layer's fill is an
# integer and the store's declared fill is what the array really holds.
function _fillvalue(fill_value, ::Type{T}, name) where {T}
    fill_value === nothing && return nothing
    if fill_value isa AbstractFloat && isnan(fill_value)
        T <: AbstractFloat && return T(NaN)
        throw(ArgumentError(
            "layer $(repr(name)) has element type $T, which has no NaN; pass " *
            "`fill_value` for what an unwritten cell should read back as, or " *
            "`fill_value = nothing` to leave it at $T's zero."))
    end
    return convert(T, fill_value)
end

# What the attributes SAY the padding reads back as. A single spelling where
# the layers agree and a per-layer object where they do not — the arrays are
# the authority either way.
function _fillspelling(fills)
    values = unique(collect(Base.values(fills)))
    length(values) == 1 || return Dict{String,Any}(k => _spell(v) for (k, v) in fills)
    return _spell(only(values))
end

_spell(x) = x isa AbstractFloat && isnan(x) ? "NaN" : x
_spell(::Nothing) = nothing

# ===========================================================================
# Writing columns
# ===========================================================================

"""
    dggwrite!(store, ancestor, values; var = the store's only layer) -> store
    dggwrite!(store, cube) -> store

Fill columns of an open [`SubzoneStore`](@ref).

The first form writes ONE column: `ancestor` is a level-`ancestor_level` cell —
or its column index — and `values` is its subtree in ascending cell id, exactly
as long as that subtree really is. A `NamedTuple` or `Dict` of vectors writes
several layers of the one column. The rest of the column, where the ancestor is
one of the twelve pentagons, stays fill.

The second form writes a whole cube: every layer of a `DimStack`, or a
`DimArray`, over a cell axis at the store's own level. Its coverage has to be
whole columns — see [`subzone_runs`](@ref) — and it need not be all of them, or
contiguous.

Only the chunks of the columns named are touched, so two tasks writing disjoint
columns are safe on a directory store, and nothing shared is rewritten.
"""
function DGG.dggwrite!(store::SubzoneStore, ancestor, values::AbstractVector;
    var::Union{Symbol,AbstractString,Nothing}=nothing)
    name = _onlyvar(store, var)
    col = _column(store.layout, ancestor)
    h = columnlength(store.layout, col)
    length(values) == h || throw(ArgumentError(
        "column $col ($(DGG.columncell(store.layout, col))) holds $h cells and " *
        "$(length(values)) values were given. A pentagon's subtree is shorter " *
        "than a hexagon's; `columnlength(layout, column)` is its length."))
    store.arrays[name][1:h, col] = values
    for overview in store.overviews
        _writeoverview!(overview, store.layout, name, col, values)
    end
    return store
end

function _writeoverview!(overview::SubzoneOverview, base::SubzoneLayout,
    name::String, col::Int, values)
    sampled = _overviewvalues(overview.method, base, overview.layout, col, values)
    overview.arrays[name][1:length(sampled), col] = sampled
    return nothing
end

# Aggregation grows here: a new stored method gets its own implementation while
# the overview arrays keep exactly the same ancestor-subzone shape.
function _overviewvalues(method::Symbol, base::SubzoneLayout,
    overview::SubzoneLayout, col::Int, values)
    method === :center || throw(ArgumentError(
        "overview method $(repr(method)) is not implemented."))
    positions = columnpositions(overview, col)
    sampled = Vector{eltype(values)}(undef, length(positions))
    for (j, position) in enumerate(positions)
        overview_cell = cellindex(overview.grid, position)
        center_position = first(DGG.descendant_range(
            system(base), overview_cell, base.level))
        base_column, base_row = positionindex(base, center_position)
        base_column == col || error(
            "center descendant escaped ancestor column $col into $base_column")
        sampled[j] = values[base_row]
    end
    return sampled
end

function DGG.dggwrite!(store::SubzoneStore, ancestor, values::Union{NamedTuple,AbstractDict})
    for (name, v) in pairs(values)
        DGG.dggwrite!(store, ancestor, v; var=name)
    end
    return store
end

function DGG.dggwrite!(store::SubzoneStore, src::Union{DD.AbstractDimArray,DD.AbstractDimStack})
    return with_store_context(store.identifier) do
        celldim, lookup = _cellaxis(src)
        layers = _cubelayers(src, celldim)
        for (name, _) in layers
            haskey(store.arrays, name) || throw(ArgumentError(
                "this store has no layer `$name`; it holds " *
                join(keys(store), ", ") * "."))
        end
        runs = subzone_runs(store.layout, lookup)
        for (name, A) in layers
            data = DD.data(A)
            for run in runs
                # One column at a time, and by `getindex` rather than by a view:
                # the source may be lazy itself — a cube straight out of
                # `dggread` or a lazy regrid — and a view of one would be read
                # element by element on the way into the chunk. A run is one
                # chunk's worth, which is the block size a lazy source wants.
                DGG.dggwrite!(store, run.column, data[run.axis]; var=name)
            end
        end
        return store
    end
end

_column(layout, a::AbstractCellIndex) = columnindex(layout, a)
_column(layout, i::Integer) = (1 <= i <= layout.ncolumns ? Int(i) : throw(ArgumentError(
    "column $i is outside the store's 1:$(layout.ncolumns) columns.")))

function _onlyvar(store::SubzoneStore, var)
    var === nothing || return _knownvar(store, String(var))
    length(store.arrays) == 1 && return first(keys(store.arrays))
    throw(ArgumentError("this store holds " * join(keys(store), ", ") *
                        "; name the one to write with `var`."))
end

_knownvar(store::SubzoneStore, name) = haskey(store.arrays, name) ? name :
                                       throw(ArgumentError("this store has no layer `$name`; it holds " *
                                                           join(keys(store), ", ") * "."))

# The cube's cell dimension and its LOOKUP — not its ids. A `CellLookup` keeps
# position windows, and `subzone_runs` walks those, so a land-only cube of tens
# of millions of cells is planned without one id being materialized. This is
# deliberately not the one-dimensional writer's `_cellaxis`, which reads the
# whole axis out as raw ids on its way to a dense coordinate.
function _cellaxis(src)
    for d in DD.dims(src)
        lk = DD.val(d)
        lk isa Union{CellLookup,ChunkedCellLookup} || continue
        return d, lk
    end
    throw(ArgumentError("a subzone write needs a cube with a cell dimension: " *
                        "none of " *
                        join(map(d -> string(DD.name(d)), DD.dims(src)), ", ") *
                        " carries a CellLookup or a ChunkedCellLookup."))
end

_cubelayers(A::DD.AbstractDimArray, celldim) =
    _checklayers([(_layername(A), A)], celldim)
_cubelayers(s::DD.AbstractDimStack, celldim) =
    _checklayers([(String(k), s[k]) for k in keys(DD.layers(s))], celldim)

function _checklayers(layers, celldim)
    isempty(layers) && throw(ArgumentError("this cube has no layers to write."))
    for (name, A) in layers
        DD.hasdim(A, celldim) || throw(ArgumentError(
            "layer $(repr(name)) has no cell dimension; every layer of a store " *
            "shares one cell axis."))
        ndims(A) == 1 || throw(ArgumentError(
            "layer $(repr(name)) has $(ndims(A)) dimensions. The ancestor-subzone " *
            "layout spends both of its own on the cell axis, so a layer over it " *
            "is one-dimensional; write a cube with a time or band axis in the " *
            "one-dimensional cell layout, or one store per step."))
    end
    return sort!(layers; by=first)
end

function _layername(A)
    n = DD.name(A)
    n isa Symbol && n !== Symbol("") && return String(n)
    return "layer"
end

# ===========================================================================
# The one-shot write
# ===========================================================================

"""
    write_subzones(dest, opengroup, src; ancestor_level, kw...) -> ZGroup

`dggwrite(dest, cube; layout = :subzones, ancestor_level = k)`: create the store
the cube's own system and level imply, and fill every column it covers.

The layers and their element types come from the cube, so this is
[`subzonestore`](@ref) followed by [`dggwrite!`](@ref) and nothing else — the
incremental path is not a second implementation of the one-shot path, it IS the
one-shot path.
"""
function write_subzones(dest, src; ancestor_level::Union{Integer,Nothing}=nothing,
    capacity=nothing, fill_value=NaN, ancestor_coordinate::Bool=true,
    overviews=Int[], overview_method::Symbol=:center,
    compressor=Zarr.BloscCompressor(), kw...)
    isempty(kw) || throw(ArgumentError(
        "`layout = :subzones` takes `ancestor_level`, `capacity`, `fill_value`, " *
        "`ancestor_coordinate`, `overviews`, `overview_method` and `compressor`; " *
        join(sort!(String[string(k) for k in keys(kw)]), ", ") *
        " belongs to the one-dimensional cell layout, whose chunking and cell " *
        "coordinate this one has neither of."))
    ancestor_level === nothing && throw(ArgumentError(
        "`layout = :subzones` needs an `ancestor_level`: it is the level whose " *
        "subtrees become the store's columns, and there is no sensible default " *
        "for how coarse a chunk should be."))
    celldim, lookup = _cellaxis(src)
    layers = _cubelayers(src, celldim)
    sys = system(lookup)
    lev = level(lookup)
    store = DGG.subzonestore(dest, sys, lev; ancestor_level=ancestor_level,
        layers=[name => eltype(A) for (name, A) in layers], capacity=capacity,
        fill_value=fill_value, ancestor_coordinate=ancestor_coordinate,
        overviews=overviews, overview_method=overview_method,
        compressor=compressor, attrs=_groupattrs(src))
    DGG.dggwrite!(store, src)
    return store.group
end

_groupattrs(src) = _attrs(get(_attrs(DD.metadata(src)), "attrs", nothing))

_attrs(md::DD.Metadata) = _attrs(DD.val(md))

# ===========================================================================
# Reading: which columns the view spans
# ===========================================================================

# The cube position -> (column, row) map. `LevelColumns` is the whole store,
# needing no tables: a position is the level grid's own, so column and row are
# one `positionindex` call. `SelectedColumns` is a subset, whose positions are
# the selected subtrees concatenated, with offsets for the binary search.
struct LevelColumns end

struct SelectedColumns
    columns::Vector{Int}
    offsets::Vector{Int}      # length n+1, offsets[1] = 0
end

function SelectedColumns(layout::SubzoneLayout, columns::AbstractVector{<:Integer})
    cols = Int[Int(i) for i in columns]
    offsets = Vector{Int}(undef, length(cols) + 1)
    offsets[1] = 0
    for (k, i) in pairs(cols)
        offsets[k+1] = offsets[k] + columnlength(layout, i)
    end
    return SelectedColumns(cols, offsets)
end

axislength(::LevelColumns, layout) = Int(ncells(layout.grid))
axislength(s::SelectedColumns, layout) = s.offsets[end]

# The column one view position falls in, its row inside that column, and the
# last view position that lands in the same column — which is what turns a
# read of a range into one read per column touched.
@inline function locate(::LevelColumns, layout, p::Int)
    i, row = positionindex(layout, p)
    return i, row, p + (columnlength(layout, i) - row)
end

@inline function locate(s::SelectedColumns, layout, p::Int)
    k = searchsortedlast(s.offsets, p - 1)
    return s.columns[k], p - s.offsets[k], s.offsets[k+1]
end

chunksizes(::LevelColumns, layout) =
    Int[columnlength(layout, i) for i in 1:layout.ncolumns]

chunksizes(s::SelectedColumns, layout) = diff(s.offsets)

# ===========================================================================
# Reading: the lazy cell-axis vector
# ===========================================================================

"""
    SubzoneCellArray(z, layout, index)

One layer of an ancestor-subzone store as the ONE-dimensional cell-axis vector
it stands for: position `k` is the `k`th cell of the view's cell axis, and the
two-dimensional store behind it is not visible.

A `DiskArrays.AbstractDiskArray` whose chunks are the store's own columns,
published as `IrregularChunks`, so anything reading by chunk reads whole
subtrees with the pentagon padding dropped. A read inside one column is one
chunk read; a read spanning columns is one per column, in order. Nothing is
cached.

Read-only: a store is written through [`dggwrite!`](@ref), which writes whole
columns and can therefore keep the padding rule.
"""
struct SubzoneCellArray{T,L<:SubzoneLayout,I,Z} <: DiskArrays.AbstractDiskArray{T,1}
    z::Z
    layout::L
    index::I
    len::Int
    chunks::Base.RefValue{Any}
end

function SubzoneCellArray(z, layout::SubzoneLayout, index)
    n = axislength(index, layout)
    return SubzoneCellArray{eltype(z),typeof(layout),typeof(index),typeof(z)}(
        z, layout, index, n, Ref{Any}(nothing))
end

Base.size(A::SubzoneCellArray) = (A.len,)

DiskArrays.haschunks(::SubzoneCellArray) = DiskArrays.Chunked()

# Built on demand and memoized: the chunk sizes are one integer per column, and
# a whole-level view over a level-6 ancestor grid has 1 176 494 of them —
# nothing to build for a read that never asks how the store is chunked.
function DiskArrays.eachchunk(A::SubzoneCellArray)
    chunks = A.chunks[]
    chunks === nothing || return chunks
    sizes = chunksizes(A.index, A.layout)
    built = DiskArrays.GridChunks(DiskArrays.IrregularChunks(; chunksizes=sizes))
    A.chunks[] = built
    return built
end

function DiskArrays.readblock!(A::SubzoneCellArray, out, r::AbstractUnitRange)
    p = first(r)
    while p <= last(r)
        col, row, colstop = locate(A.index, A.layout, p)
        stop = min(last(r), colstop)
        n = stop - p + 1
        o = p - first(r) + 1
        out[o:(o+n-1)] = A.z[row:(row+n-1), col]
        p = stop + 1
    end
    return out
end

@noinline DiskArrays.writeblock!(A::SubzoneCellArray, _, ::AbstractUnitRange) =
    throw(ArgumentError(
        "an ancestor-subzone store is not written through its cell-axis view: a " *
        "column is one chunk and is written whole, padding rule included. Use " *
        "`dggwrite!(subzonestore(path), ancestor, values)`."))

# ===========================================================================
# Reading: the cube
# ===========================================================================

"""
    assemble(group, snapshot, identifier, vars, lazy, ancestors, requested_level) -> DimStack

`dggread` on an ancestor-subzone store: one `Cells` dimension carrying the
[`CellLookup`](@ref) the written columns spell, and one lazy
[`SubzoneCellArray`](@ref) per layer over it.

`ancestors` is `nothing` for the whole store — every column of the level, with
the ones nobody wrote reading back as fill — or the ancestor cells (or column
indices) to restrict the axis to.

`requested_level` is `nothing` (or the base level) for the base arrays, or one
of the configured overview levels. Overview layers keep their base variable
names on the returned stack even though their Zarr arrays are named
`<variable>_ovr<level>`.

The stack's metadata carries the source, the group attributes verbatim, and the
[`SubzoneLayout`](@ref) under `"layout"`, which is what a caller does column
arithmetic with. It carries no `StoreDescription`: that vocabulary describes a
one-dimensional axis and has no field this layout would fill.
"""
function assemble(group, snap, identifier, vars, lazy::Bool, ancestors, requested_level)
    base = subzone_layout(snap.attrs; store=identifier)
    declarations = _overviewdeclarations(snap.attrs, base; store=identifier)
    layout, method, stored_names = _readlevel(base, declarations, requested_level,
        snap.attrs, identifier)
    available = sort!(collect(keys(stored_names)))
    selected = selectvars(available, vars)
    columns = ancestors === nothing ? nothing : subzone_columns(layout, ancestors)
    index = columns === nothing ? LevelColumns() : SelectedColumns(layout, columns)
    lookup = CellLookup(subzone_cellvector(layout, columns))

    layers = map(selected) do name
        stored = stored_names[name]
        z = group[stored]
        _checkshape(z, layout, stored, identifier)
        A = SubzoneCellArray(z, layout, index)
        entry = DGG.getarray(snap, stored)
        DD.DimArray(lazy ? A : Array(A), (Cells(lookup),);
            name=Symbol(name), metadata=deepcopy(entry.attrs))
    end

    metadata = Dict{String,Any}(
        "source" => identifier,
        "layout" => layout,
        "attrs" => deepcopy(snap.attrs))
    if method !== nothing
        metadata["overview_level"] = layout.level
        metadata["overview_method"] = method
        metadata["base_level"] = base.level
    end
    return DD.DimStack(NamedTuple{Tuple(Symbol.(selected))}(Tuple(layers)); metadata)
end

function _readlevel(base, declarations, requested_level, attrs, identifier)
    block = attrs["dggs"][DGG.SUBZONE_BLOCK]
    base_variables = get(block, "variables", Any[])
    base_variables isa AbstractVector || _overviewformaterror(identifier,
        "base `variables` is not a list.", base_variables)
    all(x -> x isa AbstractString, base_variables) || _overviewformaterror(identifier,
        "base `variables` contains a non-string name.", base_variables)
    names = String[String(x) for x in base_variables]
    allunique(names) || _overviewformaterror(identifier,
        "base `variables` contains duplicate names.", names)

    if requested_level === nothing || requested_level == base.level
        return base, nothing, Dict{String,String}(name => name for name in names)
    end
    requested_level isa Integer && !(requested_level isa Bool) || throw(ArgumentError(
        "`level` is the base level $(base.level) or one of its integer overview " *
        "levels, not $(repr(requested_level))."))
    lo = Int(requested_level)
    for declaration in declarations
        declaration.layout.level == lo || continue
        return declaration.layout, declaration.method, declaration.variables
    end
    configured = Int[d.layout.level for d in declarations]
    throw(ArgumentError("level $lo is not stored: the base is level $(base.level) " *
        "and the configured overview levels are $(repr(configured))."))
end

end # module DGGSZarrSubzones
