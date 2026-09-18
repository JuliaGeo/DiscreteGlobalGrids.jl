# The Zarr half of the pyramid layout (see `src/io/pyramid.jl`, which owns the
# arithmetic and the attribute vocabulary): Zarr calls, the chunk-at-a-time
# write, and the DiskArrays wrapper that makes a slot-space array look like the
# cell axis it stands for.
#
#   write  dggwrite(dest, cube; layout = :pyramid)
#   read   dggread(store)  ->  StorePyramid
#
# The write is chunk-aligned throughout, which is what makes it sparse: a buffer
# is prefilled with the array's fill value, the cells that landed in that chunk
# are scattered into it, and it is written as one assignment covering exactly
# one chunk. A chunk no cell landed in is never touched and Zarr never stores
# it; a chunk whose every value came out fill is skipped for the same reason.
#
# A submodule, like the other two layouts, so its helper names cannot collide
# with the extension's shared namespace.

module DGGSZarrPyramid

import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: AbstractCellLookup, CellLookup, CellVector, Cells,
    DGGSFormatError, PyramidLayout, StorePyramid,
    DEFAULT_PYRAMID_CHUNK_EXPONENT, PYRAMID_LEVEL_PREFIX,
    ancestor, cellindex, chunkindices, chunklength, chunkrootlevel, chunkslots,
    descendant_range, globalindex, level, levelfromname, levelgrid, levelname,
    levels, levelslots, ncells, pyramid_attrs, pyramid_layout,
    pyramid_level_attrs, pyramid_dimension, pyramid_reduce, pyramid_runs,
    slotindex, system, with_store_context
import ..DiscreteGlobalGridsZarrExt
using ..DiscreteGlobalGridsZarrExt: ARRAY_DIMENSIONS, selectvars, storeidentifier
import DimensionalData as DD
import DiskArrays
import Zarr
using Statistics: mean

const REMOTE_SCHEMES = ("gs://", "s3://", "http://", "https://", "az://", "abfs://")

# ===========================================================================
# Writing
# ===========================================================================

"""
    write_pyramid(dest, src; aggregate, fill_value, chunks, overwrite,
                  compressor) -> ZGroup

`dggwrite(dest, cube; layout = :pyramid)`: write a cube's own level, then every
coarser level above it, into one group.

Each variable becomes a subgroup of `level0 … level{L}` arrays. The cube's level
is written first, chunk by chunk; each chunk is reduced to its parents as it is
written, so the level above is built from what was just stored rather than from
a second pass over the store, and no level of the pyramid is ever materialized
whole.

  - `aggregate` reduces the children that carry a value to their parent, and
    defaults to `mean` for a floating-point cube. A cube of another element type
    has to name one, because every level carries the cube's own element type and
    a mean does not land back in an integer.
  - `fill_value` is what an absent cell reads back as — `NaN` by default, which
    needs a floating-point layer. It is the layout's presence marker as well as
    its padding, so a value the aggregate can produce is refused.
  - `chunks` is `:auto` or a chunk length, which must be a power of the system's
    aperture.
  - `overwrite` clears what is in the way: the whole store where `dest` is a
    path, and the variables' own subgroups where it is an open group.
"""
function write_pyramid(dest, src; aggregate=nothing, fill_value=nothing,
    chunks=:auto, overwrite::Bool=false, compressor=Zarr.BloscCompressor(), kw...)

    isempty(kw) || throw(ArgumentError(
        "`layout = :pyramid` takes `aggregate`, `fill_value`, `chunks`, " *
        "`overwrite` and `compressor`; " *
        join(sort!(String[string(k) for k in keys(kw)]), ", ") *
        " belongs to the one-dimensional cell layout. This layout stores no cell " *
        "coordinate, so it has no encoding to choose and no runs to merge, and " *
        "its chunks are a power of the aperture rather than an element target."))

    celldim, lookup = _cellaxis(src)
    layers = _cubelayers(src, celldim)
    sys = system(lookup)
    L = level(lookup)
    layout = PyramidLayout(sys, L; chunkexponent=_chunkexponent(chunks, sys, L))

    # Both keywords are resolved BEFORE anything is created, so a cube whose
    # element type has no default fill leaves no half-written store behind.
    specs = [(name, A, _fillvalue(fill_value, eltype(A), name),
        _aggregate(aggregate, eltype(A), name)) for (name, A) in layers]

    group = _opengroup(dest, layers, overwrite)
    identifier = storeidentifier(group)
    return with_store_context(identifier) do
        attrs = _attrs(get(_attrs(DD.metadata(src)), "attrs", nothing))
        merge!(attrs, pyramid_attrs(layout;
            variables=String[name for (name, _, _, _) in specs],
            fill_value=_fillspelling(specs),
            aggregate=_aggregatename(aggregate)))
        _stamp(group, attrs)
        for (name, A, fillvalue, reduce) in specs
            _writevariable!(group, layout, name, A, lookup, fillvalue, reduce,
                compressor)
        end
        Zarr.consolidate_metadata(group)
        return group
    end
end

# --- the destination --------------------------------------------------------

# A path names a whole store, so `overwrite` clears the whole store; a variable
# added to a store that already exists goes through the `ZGroup` method, which
# clears only that variable. `zgroup` refuses a directory holding anything else.
function _opengroup(dest::AbstractString, layers, overwrite::Bool)
    path = String(dest)
    _reject_remote(path)
    overwrite && ispath(path) && rm(path; recursive=true)
    return Zarr.zgroup(path; attrs=Dict{String,Any}())
end

function _opengroup(dest::Zarr.ZGroup, layers, overwrite::Bool)
    dest.writeable || throw(ArgumentError(
        "dggwrite needs a writeable group; this one was opened read-only."))
    overwrite && _dropvariables(dest, layers)
    return _checkempty(dest, layers)
end

function _stamp(g::Zarr.ZGroup, attrs)
    merge!(g.attrs, attrs)
    Zarr.writeattrs(g.zarr_format, g.storage, g.path, g.attrs)
    return g
end

# What the attributes SAY an absent cell reads back as: one spelling where the
# variables agree on it, a per-variable object where they do not. The arrays
# each carry their own for real, and that is what the reader trusts.
function _fillspelling(specs)
    fills = Dict{String,Any}(name => _spell(f) for (name, _, f, _) in specs)
    values = unique(collect(Base.values(fills)))
    return length(values) == 1 ? only(values) : fills
end

_spell(x) = x isa AbstractFloat && isnan(x) ? "NaN" : x

# Provenance, and only where it is worth recording: an anonymous reduction has
# no name a later reader could act on, so the key is left out rather than filled
# with a gensym.
_aggregatename(::Nothing) = "mean"
function _aggregatename(f)
    s = string(f)
    return startswith(s, "#") ? nothing : s
end

function _checkempty(g::Zarr.ZGroup, layers)
    taken = sort!(String[name for (name, _) in layers if haskey(g.groups, name)])
    isempty(taken) || throw(DGGSFormatError(check=:destination_not_empty,
        declared=taken, observed=sort!(collect(keys(g.groups))),
        detail="this group already holds " * join(taken, ", ") *
               ", and dggwrite does not overwrite a variable. Pass " *
               "`overwrite = true`, write to a new group, or delete these first."))
    return g
end

# A variable is a subgroup, and a subgroup is a directory. Removing one is
# exact where the store is a directory and has no store-agnostic spelling
# anywhere else, so that is the only place a group-level `overwrite` works.
function _dropvariables(g::Zarr.ZGroup, layers)
    _dropvariables(g.storage, g.path, layers)
    for (name, _) in layers
        delete!(g.groups, String(name))
    end
    return nothing
end

function _dropvariables(store::Zarr.DirectoryStore, path::AbstractString, layers)
    for (name, _) in layers
        dir = joinpath(store.folder, String(path), String(name))
        isdir(dir) && rm(dir; recursive=true)
    end
    return nothing
end

@noinline _dropvariables(store::Zarr.AbstractStore, ::AbstractString, _) = throw(ArgumentError(
    "`overwrite = true` deletes the variable's subgroup, which this package " *
    "can do on a directory store and on no other: $(nameof(typeof(store))) has " *
    "no deletion this package may spell. Delete the variable yourself, or write " *
    "to a new group."))

@noinline function _reject_remote(path)
    for scheme in REMOTE_SCHEMES
        startswith(path, scheme) || continue
        throw(ArgumentError("dggwrite writes local directory stores only; " *
                            "$(repr(path)) names a $(rstrip(scheme, ['/', ':'])) store. " *
                            "Write locally and upload, or open the remote group yourself " *
                            "and pass the ZGroup."))
    end
    return nothing
end

# --- the chunk plan ---------------------------------------------------------

function _chunkexponent(chunks, sys, L::Int)
    chunks === :auto && return DEFAULT_PYRAMID_CHUNK_EXPONENT
    chunks isa Integer || throw(ArgumentError(
        "`chunks` is `:auto` or a chunk length, and $(repr(chunks)) is neither."))
    n = Int(chunks)
    a = DGG.slotcount(sys, 1) ÷ DGG.slotcount(sys, 0)
    k = 0
    p = 1
    while p < n
        p *= a
        k += 1
    end
    p == n || throw(ArgumentError(
        "the pyramid layout chunks at a power of the aperture, so that a chunk " *
        "is one subtree: $n is not a power of $a. The nearest are $(p ÷ a) and $p."))
    k <= L || throw(ArgumentError(
        "a chunk of $n slots is $a^$k, which is deeper than the cube's own " *
        "level $L; the largest chunk a level-$L pyramid can have is $(a^L)."))
    return k
end

# --- one variable -----------------------------------------------------------

function _writevariable!(group, layout::PyramidLayout, name, A, lookup, fillvalue,
    reduce, compressor)

    T = eltype(A)
    vargroup = Zarr.zgroup(group, String(name); attrs=_layerattrs(A))
    arrays = Dict{Int,Zarr.ZArray}()
    for J in levels(layout)
        arrays[J] = Zarr.zcreate(T, vargroup, levelname(layout, J),
            levelslots(layout, J); chunks=(chunklength(layout, J),),
            fill_value=fillvalue, fill_as_missing=false, compressor=compressor,
            attrs=pyramid_level_attrs(layout, J, pyramid_dimension(layout, J)))
    end

    L = layout.level
    indices, values = _writelevel!(arrays[L], layout, L,
        pyramid_runs(layout, L, lookup), parent(A), fillvalue, reduce, true)
    for J in (L-1):-1:first(levels(layout))
        runs = DGG.pyramid_index_runs(layout, J, indices)
        indices, values = _writelevel!(arrays[J], layout, J, runs, values, fillvalue,
            reduce, J > first(levels(layout)))
    end
    return nothing
end

# Write every chunk the runs touch, and — while each chunk's buffer is still in
# hand — reduce it to its parents. A parent's children never straddle two
# chunks, because a chunk is a subtree at least one level above, so this is the
# whole of the level above and needs no second pass.
function _writelevel!(z, layout::PyramidLayout, J::Int, runs, values, fillvalue,
    reduce, buildparents::Bool)

    T = eltype(z)
    w = chunklength(layout, J)
    buffer = Vector{T}(undef, w)
    parentindices = Int[]
    parentvalues = T[]

    i = 1
    n = length(runs)
    while i <= n
        chunk = runs[i].chunk
        j = i
        while j <= n && runs[j].chunk == chunk
            j += 1
        end
        slots = chunkslots(layout, J, chunk)
        fill!(buffer, fillvalue)
        touched = false
        for t in i:(j-1)
            run = runs[t]
            # By `getindex` rather than by a view: the source may be lazy, and a
            # view of a lazy array is read element by element on the way in. A
            # run is at most one chunk's worth, which is the block a lazy source
            # wants anyway.
            block = values[run.axis]
            for (u, row) in enumerate(run.rows)
                v = block[u]
                isequal(v, fillvalue) && continue
                buffer[slots[row]] = v
                touched = true
            end
        end
        if touched
            base = (chunk - 1) * w
            z[(base+1):(base+w)] = buffer
            if buildparents
                cells = chunkindices(layout, J, chunk)
                a, b = pyramid_reduce(layout, J, cells,
                    T[@inbounds buffer[s] for s in slots], reduce, fillvalue)
                append!(parentindices, a)
                append!(parentvalues, b)
            end
        end
        i = j
    end
    return parentindices, parentvalues
end

# --- keyword resolution -----------------------------------------------------

# NaN is the default and needs a float to be one. The fill value is this
# layout's presence marker, so unlike the subzone layout there is no "leave it
# at the element type's zero": a store whose absent cells read back as zero
# cannot tell a written zero from an unwritten cell.
function _fillvalue(fill_value, ::Type{T}, name) where {T}
    if fill_value === nothing
        T <: AbstractFloat && return T(NaN)
        throw(ArgumentError(
            "variable $(repr(String(name))) has element type $T, which has no " *
            "NaN to default to. A pyramid's fill value is what says `nothing " *
            "here`, and the coarse levels are the chunk index, so it cannot be " *
            "left to the element type's zero: pass `fill_value` for a $T that " *
            "the data never takes."))
    end
    return convert(T, fill_value)
end

function _aggregate(aggregate, ::Type{T}, name) where {T}
    aggregate === nothing || return aggregate
    T <: AbstractFloat && return mean
    throw(ArgumentError(
        "variable $(repr(String(name))) has element type $T, and `aggregate` " *
        "defaults to `mean` only for a floating-point one: every level of a " *
        "pyramid carries the cube's own element type, and a mean does not land " *
        "back in $T. Name the reduction the coarse levels should hold — `mode` " *
        "for a class, `vs -> round($T, sum(vs) / length(vs))` for a rounded " *
        "mean — or write the cube as a floating-point type."))
end

# --- the cube ---------------------------------------------------------------

# The cube's cell dimension and its LOOKUP, not its ids: a lookup keeps index
# windows, and `pyramid_runs` walks those, so a land-only cube of tens of
# millions of cells is planned without one id being materialized.
function _cellaxis(src)
    for d in DD.dims(src)
        lk = DD.val(d)
        lk isa AbstractCellLookup || continue
        return d, lk
    end
    throw(ArgumentError("a pyramid write needs a cube with a cell dimension: " *
                        "none of " *
                        join(map(d -> string(DD.name(d)), DD.dims(src)), ", ") *
                        " carries a cell lookup."))
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
            "layer $(repr(name)) has $(ndims(A)) dimensions. A pyramid store " *
            "spends its one dimension on the cell axis, so a layer over it is " *
            "one-dimensional; write a cube with a time or band axis in the " *
            "one-dimensional cell layout, or one store per step."))
        levelfromname(name) === nothing || throw(ArgumentError(
            "layer $(repr(name)) is named like one of the store's own level " *
            "arrays. A variable is a subgroup holding `$(PYRAMID_LEVEL_PREFIX)0`, " *
            "`$(PYRAMID_LEVEL_PREFIX)1`, … arrays, so it cannot be called that " *
            "itself; rename the layer."))
    end
    return sort!(layers; by=first)
end

function _layername(A)
    n = DD.name(A)
    n isa Symbol && n !== Symbol("") && return String(n)
    throw(ArgumentError(
        "a pyramid store keeps each variable in a subgroup of its own, and this " *
        "array has no name to call it. Pass a named `DimArray`, or a `DimStack`."))
end

_layerattrs(A) = _attrs(DD.metadata(A))
_attrs(x) = Dict{String,Any}()
_attrs(md::AbstractDict) = Dict{String,Any}(String(k) => deepcopy(v) for (k, v) in md)
_attrs(md::NamedTuple) = Dict{String,Any}(String(k) => deepcopy(v) for (k, v) in pairs(md))
_attrs(md::DD.Metadata) = _attrs(DD.val(md))

# ===========================================================================
# Reading: one level as the cell axis it stands for
# ===========================================================================

"""
    PyramidChunks(layout, J)

Level `J`'s chunk grid, as the DENSE cell-axis ranges it cuts the level into.

One entry per chunk that holds a cell, which is one per cell of the chunk root
level, so entry `k` is the subtree of the `k`th cell of that level and the
chunks under absent addresses are simply not there. Computed, never built: a
level-13 store has eight million of these and a reader that materialized them
would pay for the whole level to read one chunk of it.
"""
struct PyramidChunks{L<:PyramidLayout,G} <: DiskArrays.ChunkVector
    layout::L
    level::Int
    rootgrid::G
    n::Int
end

function PyramidChunks(layout::PyramidLayout, J::Integer)
    rl = chunkrootlevel(layout, J)
    rootgrid = levelgrid(system(layout), rl)
    return PyramidChunks(layout, Int(J), rootgrid, Int(ncells(system(layout), rl)))
end

Base.size(c::PyramidChunks) = (c.n,)

Base.@propagate_inbounds function Base.getindex(c::PyramidChunks, i::Int)
    @boundscheck checkbounds(c, i)
    return descendant_range(system(c.layout), cellindex(c.rootgrid, i), c.level)
end

function DiskArrays.findchunk(c::PyramidChunks, i::Int)
    1 <= i <= ncells(system(c.layout), c.level) || throw(BoundsError(c, i))
    sys = system(c.layout)
    grid = levelgrid(sys, c.level)
    a = ancestor(sys, cellindex(grid, i), level(c.rootgrid))
    return globalindex(c.rootgrid, a)::Int
end

# Bounded work: only the chunks the subset touches are looked at, where the
# generic fallback would size a vector by the whole chunk count.
function DiskArrays.subsetchunks(c::PyramidChunks, subsets::AbstractUnitRange)
    isempty(subsets) && return DiskArrays.RegularChunks(1, 0, 0)
    first(subsets) >= 1 || throw(BoundsError(c, first(subsets)))
    lo = DiskArrays.findchunk(c, first(subsets))
    hi = DiskArrays.findchunk(c, last(subsets))
    sizes = Int[length(intersect(c[k], subsets)) for k in lo:hi]
    return DiskArrays.IrregularChunks(; chunksizes=sizes)
end

"""
    PyramidLevelArray(z, layout, J)

One level of a pyramid store as the cell-axis vector it stands for: index `k` is
the `k`th cell of the complete level in canonical order, and the slot space the
array is really laid out on is not visible.

A `DiskArrays.AbstractDiskArray` whose chunks are the store's own — published
through [`PyramidChunks`](@ref) — so anything reading by chunk reads whole
subtrees with the absent addresses dropped. A read inside one chunk is one chunk
read; a read spanning chunks is one per chunk, in order. Nothing is cached; wrap
it in `DiskArrays.cache` for that.

Read-only: a pyramid is written through
[`dggwrite`](@ref DiscreteGlobalGrids.dggwrite), which writes whole chunks and
can therefore keep every level consistent with the one above it.
"""
struct PyramidLevelArray{T,L<:PyramidLayout,Z,G,C} <: DiskArrays.AbstractDiskArray{T,1}
    z::Z
    layout::L
    level::Int
    grid::G
    len::Int
    chunks::C
end

function PyramidLevelArray(z, layout::PyramidLayout, J::Integer)
    j = Int(J)
    sys = system(layout)
    grid = levelgrid(sys, j)
    chunks = DiskArrays.GridChunks(PyramidChunks(layout, j))
    return PyramidLevelArray{eltype(z),typeof(layout),typeof(z),typeof(grid),
        typeof(chunks)}(z, layout, j, grid, Int(ncells(sys, j)), chunks)
end

Base.size(A::PyramidLevelArray) = (A.len,)

DiskArrays.haschunks(::PyramidLevelArray) = DiskArrays.Chunked()
DiskArrays.eachchunk(A::PyramidLevelArray) = A.chunks

function DiskArrays.readblock!(A::PyramidLevelArray, out, r::AbstractUnitRange)
    isempty(r) && return out
    layout = A.layout
    sys = system(layout)
    rl = chunkrootlevel(layout, A.level)
    w = chunklength(layout, A.level)
    p = first(r)
    while p <= last(r)
        root = ancestor(sys, cellindex(A.grid, p), rl)
        chunk = slotindex(sys, root)
        cells = descendant_range(sys, root, A.level)
        stop = min(last(r), last(cells))
        row = p - first(cells) + 1
        n = stop - p + 1
        o = p - first(r) + 1
        base = (chunk - 1) * w
        slots = chunkslots(layout, A.level, chunk)
        if slots isa AbstractUnitRange
            # A complete subtree: cell position is slot offset, so the run is
            # one contiguous read of the store.
            out[o:(o+n-1)] = A.z[(base+row):(base+row+n-1)]
        else
            # A pentagon's subtree, where the absent addresses sit between the
            # cells. One read of the whole chunk, then a gather: the store would
            # read the whole chunk for any part of it anyway.
            whole = A.z[(base+1):(base+w)]
            for t in 1:n
                out[o+t-1] = whole[slots[row+t-1]]
            end
        end
        p = stop + 1
    end
    return out
end

@noinline DiskArrays.writeblock!(A::PyramidLevelArray, _, ::AbstractUnitRange) =
    throw(ArgumentError(
        "a pyramid store is not written through its cell-axis view: every level " *
        "has to stay consistent with the one above it, which a write to one " *
        "level cannot keep. Write the cube with " *
        "`dggwrite(dest, cube; layout = :pyramid)`."))

# ===========================================================================
# Reading: the pyramid
# ===========================================================================

"""
    assemble(group, snapshot, identifier, vars, lazy, cache) -> StorePyramid

`dggread` on a pyramid store: one lazy cube per level per variable, over that
level's complete cell axis, collected into a [`StorePyramid`](@ref DiscreteGlobalGrids.StorePyramid).

The whole pyramid is opened, not one level of it: the coarse levels ARE the
index to the fine ones, so a reader that wanted level 13 would have to read
level 7 to find out which of its chunks exist. They cost a sixth of the store
between them.

`cache` wraps each level in `DiskArrays.cache` — see
[`dggread`](@ref DiscreteGlobalGrids.dggread) for what that is worth and what
it costs.
"""
function assemble(group, snap, identifier, vars, lazy::Bool, cache)
    layout = pyramid_layout(snap.attrs; store=identifier)
    available = sort!(String[name for name in keys(group.groups)])
    isempty(available) && throw(DGGSFormatError(check=:no_data_variables,
        store=String(identifier), observed=sort!(collect(keys(group.arrays))),
        detail="this pyramid store holds no variable subgroup."))
    selected = selectvars(available, vars)

    stacks = Vector{Any}(undef, length(levels(layout)))
    fills = Dict{Symbol,Any}()
    arrays = Dict{String,Dict{Int,Zarr.ZArray}}()
    for name in selected
        arrays[name] = _variablelevels(group.groups[name], layout, name, identifier)
        fills[Symbol(name)] = _storedfill(arrays[name][layout.level], name, identifier)
    end

    for (k, J) in pairs(levels(layout))
        lookup = CellLookup(CellVector(levelgrid(system(layout), J)))
        layers = map(selected) do name
            z = arrays[name][J]
            A = PyramidLevelArray(z, layout, J)
            data = lazy ? _cached(A, cache) : Array(A)
            DD.DimArray(data, (Cells(lookup),);
                name=Symbol(name), metadata=Dict{String,Any}(z.attrs))
        end
        stacks[k] = DD.DimStack(
            NamedTuple{Tuple(Symbol.(selected))}(Tuple(layers));
            metadata=Dict{String,Any}("source" => identifier,
                "layout" => layout, "refinement_level" => J))
    end

    names = Tuple(Symbol.(selected))
    return StorePyramid(layout, [s for s in stacks],
        NamedTuple{names}(Tuple(fills[n] for n in names)), String(identifier))
end

# A cell axis read run by run visits one chunk many times over -- the tile in
# the hydrology tutorial is 9887 runs living in 76 chunks -- and without a cache
# each visit decompresses it again. The budget is in megabytes, and `false`
# leaves the array stateless, which is what a reader sharing it between tasks
# wants.
_cached(A, cache::Bool) = cache ? DiskArrays.cache(A) : A

function _cached(A, cache::Real)
    cache > 0 || throw(ArgumentError(
        "`cache` is a budget in megabytes and is positive, or `false` for none; " *
        "$cache is neither."))
    return DiskArrays.cache(A; maxsize=Int(cache))
end

@noinline _cached(A, cache) = throw(ArgumentError(
    "`cache` is `false`, `true`, or a budget in megabytes, not $(repr(cache))."))

function _variablelevels(vargroup, layout::PyramidLayout, name, identifier)
    out = Dict{Int,Zarr.ZArray}()
    for J in levels(layout)
        key = levelname(layout, J)
        haskey(vargroup.arrays, key) || throw(DGGSFormatError(
            check=:missing_pyramid_level, store=String(identifier), declared=key,
            observed=sort!(collect(keys(vargroup.arrays))),
            detail="variable `$name` is missing level $J. A pyramid store holds " *
                   "every level from $(first(levels(layout))) to " *
                   "$(layout.level): the coarse ones are how a reader finds out " *
                   "which chunks of the fine ones exist."))
        z = vargroup.arrays[key]
        _checkshape(z, layout, J, name, identifier)
        out[J] = z
    end
    return out
end

function _checkshape(z, layout::PyramidLayout, J::Integer, name, identifier)
    want = levelslots(layout, J)
    size(z) == (want,) || throw(DGGSFormatError(check=:pyramid_shape_mismatch,
        store=String(identifier), declared=want, observed=size(z),
        detail="`$name/$(levelname(layout, J))` is $(join(size(z), "x")) where " *
               "level $J's slot space is $want."))
    wantchunk = chunklength(layout, J)
    only(z.metadata.chunks) == wantchunk || throw(DGGSFormatError(
        check=:pyramid_chunk_mismatch, store=String(identifier),
        declared=wantchunk, observed=only(z.metadata.chunks),
        detail="`$name/$(levelname(layout, J))` is chunked at " *
               "$(only(z.metadata.chunks)) where the layout's exponent " *
               "$(layout.chunkexponent) makes level $J $wantchunk. A chunk has " *
               "to be one subtree, or the level above does not index this one."))
    return nothing
end

function _storedfill(z, name, identifier)
    f = z.metadata.fill_value
    f === nothing && throw(DGGSFormatError(check=:missing_fill_value,
        store=String(identifier), declared=name,
        detail="variable `$name` declares no `fill_value`. A pyramid's fill " *
               "value is what says `nothing here`, and its coarse levels index " *
               "its fine ones by it, so a store without one cannot be read."))
    return f
end

end # module DGGSZarrPyramid
