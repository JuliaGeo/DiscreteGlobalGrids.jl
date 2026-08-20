# The Zarr half of the ancestor-subzone layout: the store, the incremental
# writer, and the lazy vector that makes a two-dimensional store look like the
# one-dimensional cell axis it stands for.
#
#   create   subzonestore(dest, sys, level; ancestor_level, layers)
#   fill     dggwrite!(store, ancestor, values)   |   dggwrite!(store, cube)
#   one shot dggwrite(dest, cube; layout = :subzones, ancestor_level)
#   read     dggread(store)  ->  a cell-axis DimStack over `SubzoneCellArray`s
#
# Everything format-semantic — the column arithmetic, the completeness rule, the
# attribute vocabulary — is `src/io/subzones.jl`'s. What is here is the Zarr
# calls, the chunk plan (which is not a plan: one column is one chunk), and the
# DiskArrays wrapper.
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
struct SubzoneStore{L<:SubzoneLayout}
    group::Zarr.ZGroup
    layout::L
    arrays::Dict{String,Zarr.ZArray}
    identifier::String
end

DGG.system(s::SubzoneStore) = system(s.layout)
DGG.level(s::SubzoneStore) = level(s.layout)
Base.keys(s::SubzoneStore) = sort!(collect(keys(s.arrays)))

Base.show(io::IO, s::SubzoneStore) = print(io, "SubzoneStore(", s.identifier, ", ",
    s.layout, ", ", join(keys(s), ", "), ")")

"""
    subzonestore(dest, system, level; ancestor_level, layers,
                 capacity = nothing, fill_value = NaN, ancestor_coordinate = true,
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
    compressor=Zarr.BloscCompressor())

    specs = _layerspecs(layers)
    layout = SubzoneLayout(sys, lev, ancestor_level; capacity=capacity)
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

        group = opengroup(groupattrs)
        arrays = Dict{String,Zarr.ZArray}()
        for (name, T) in specs
            arrays[name] = Zarr.zcreate(T, group, name, layout.capacity, layout.ncolumns;
                chunks=(layout.capacity, 1), fill_value=fills[name],
                fill_as_missing=false, compressor=compressor,
                attrs=Dict{String,Any}(ARRAY_DIMENSIONS =>
                    [ANCESTOR_DIMENSION, SUBZONE_DIMENSION]))
        end
        coordinate === nothing || _writecoordinate(group, layout, coordinate)
        Zarr.consolidate_metadata(group)
        return SubzoneStore(group, layout, arrays, String(identifier))
    end
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
        arrays = Dict{String,Zarr.ZArray}()
        for (name, z) in group.arrays
            _issubzonearray(z) || continue
            _checkshape(z, layout, name, identifier)
            arrays[name] = z
        end
        isempty(arrays) && throw(DGGSFormatError(check=:no_data_variables,
            store=identifier, observed=sort!(collect(keys(group.arrays))),
            detail="this ancestor-subzone store holds no array on the " *
                   "($ANCESTOR_DIMENSION, $SUBZONE_DIMENSION) dimensions."))
        return SubzoneStore(group, layout, arrays, identifier)
    end
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
    return store
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
            z = store.arrays[name]
            data = DD.data(A)
            for run in runs
                # One column at a time, and by `getindex` rather than by a view:
                # the source may be lazy itself — a cube straight out of
                # `dggread` or a lazy regrid — and a view of one would be read
                # element by element on the way into the chunk. A run is one
                # chunk's worth, which is the block size a lazy source wants.
                z[run.rows, run.column] = data[run.axis]
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
    compressor=Zarr.BloscCompressor(), kw...)
    isempty(kw) || throw(ArgumentError(
        "`layout = :subzones` takes `ancestor_level`, `capacity`, `fill_value`, " *
        "`ancestor_coordinate` and `compressor`; " *
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
        compressor=compressor, attrs=_groupattrs(src))
    DGG.dggwrite!(store, src)
    return store.group
end

_groupattrs(src) = _attrs(get(_attrs(DD.metadata(src)), "attrs", nothing))

_attrs(md::DD.Metadata) = _attrs(DD.val(md))

# ===========================================================================
# Reading: which columns the view spans
# ===========================================================================

# The cube position -> (column, row) map, in the two shapes it takes.
#
# `LevelColumns` is the whole store, and the reason it is its own type is that
# it needs NO tables: the view's positions are the level grid's own positions,
# so a column and a row are one `positionindex` call, and a store of 10^12 cells
# is opened in constant time and constant memory.
#
# `SelectedColumns` is a subset, whose positions are the selected subtrees
# concatenated; it carries the offsets that make that a binary search.
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

A `DiskArrays.AbstractDiskArray`, so it slices, iterates and materializes like
any lazy array — and its chunks are the store's own columns, published as
IRREGULAR chunks. That is the point of the layout: a hexagon subtree is `7^d`
cells and a pentagon's is `(5*7^d + 1)/6`, which Zarr's uniform chunk grid
cannot express and `DiskArrays.IrregularChunks` can. Anything that reads by
chunk — `eachchunk`, a lazy regrid, a copy into another array — therefore reads
whole subtrees, one chunk file per chunk, with the pentagon padding already
dropped.

A read inside one column is one chunk read; a read spanning columns is one per
column, in order. Nothing is cached: the chunk cache, where one is wanted, is
`DiskArrays.cache`'s business and not this wrapper's.

Writing is not supported here — a store is written through [`dggwrite!`](@ref),
which writes whole columns and can therefore keep the padding rule.
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
    assemble(group, snapshot, identifier, vars, lazy, ancestors) -> DimStack

`dggread` on an ancestor-subzone store: one `Cells` dimension carrying the
[`CellLookup`](@ref) the written columns spell, and one lazy
[`SubzoneCellArray`](@ref) per layer over it.

`ancestors` is `nothing` for the whole store — every column of the level, with
the ones nobody wrote reading back as fill — or the ancestor cells (or column
indices) to restrict the axis to.

The stack's metadata carries the source, the group attributes verbatim, and the
[`SubzoneLayout`](@ref) under `"layout"`, which is what a caller does column
arithmetic with. It carries no `StoreDescription`: that vocabulary describes a
one-dimensional axis and has no field this layout would fill.
"""
function assemble(group, snap, identifier, vars, lazy::Bool, ancestors)
    layout = subzone_layout(snap.attrs; store=identifier)
    available = String[a.name for a in snap.arrays
                       if a.dims == [ANCESTOR_DIMENSION, SUBZONE_DIMENSION]]
    selected = selectvars(available, vars)
    columns = ancestors === nothing ? nothing : subzone_columns(layout, ancestors)
    index = columns === nothing ? LevelColumns() : SelectedColumns(layout, columns)
    lookup = CellLookup(subzone_cellvector(layout, columns))

    layers = map(selected) do name
        z = group[name]
        _checkshape(z, layout, name, identifier)
        A = SubzoneCellArray(z, layout, index)
        entry = DGG.getarray(snap, name)
        DD.DimArray(lazy ? A : Array(A), (Cells(lookup),);
            name=Symbol(name), metadata=deepcopy(entry.attrs))
    end

    metadata = Dict{String,Any}(
        "source" => identifier,
        "layout" => layout,
        "attrs" => deepcopy(snap.attrs))
    return DD.DimStack(NamedTuple{Tuple(Symbol.(selected))}(Tuple(layers)); metadata)
end

end # module DGGSZarrSubzones
