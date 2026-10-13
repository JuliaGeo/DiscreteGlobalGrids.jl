# The pyramid layout: every level of a hierarchy, each one array, each laid out
# on its own SLOT space rather than on its cells.
#
#     variable/level0      slotcount(sys, 0) values, chunked aᵐⁱⁿ⁽ᵏ'⁰⁾
#     variable/level1      slotcount(sys, 1) values, chunked aᵐⁱⁿ⁽ᵏ'¹⁾
#     ...
#     variable/level{L}    slotcount(sys, L) values, chunked aᵏ
#
# Two things follow from putting the data on slots and chunking at a power of
# the aperture, and the layout exists for both.
#
#   * A chunk IS a subtree. Slots `(c-1)aᵏ+1 : c·aᵏ` of level `J` are exactly the
#     descendants of slot `c` of level `J-k`, cells and absent addresses alike,
#     so a region written into the store touches only the chunks its subtrees
#     fall in and Zarr never stores the rest. Half an earth costs half an earth.
#   * The coarse level IS the chunk index. `level{J-k}` has one value per chunk
#     of `level{J}`, at the same position, and a parent is fill exactly when
#     every cell under it is — so a reader learns which chunks of a deep level
#     exist by reading a shallow one, with no LIST and no per-chunk probe. That
#     is why the overviews are not optional here.
#
# The price is the addresses that name no cell — a sixth of IGEO7's slots, the
# twelve pentagons' deleted branches. They fall in aligned blocks of `aᵏ⁻¹ … a⁰`
# per base, so all but the smallest are whole chunks that are never written, and
# what a full store really carries is `2·(aᵏ - 1)` fill values per level.
#
# Arithmetic, values and vocabulary only: no Zarr and no arrays. The store half
# is `ext/DiscreteGlobalGridsZarrExt/pyramid.jl`.
#
# Include order: after `conventions.jl`, whose `GRID_REFERENCE` and attribute
# helpers this reads, and before `api.jl`, whose stubs document the verbs the
# Zarr extension implements.

# ===========================================================================
# The on-disk vocabulary
# ===========================================================================

"""
    PYRAMID_LAYOUT

The `layout` string a pyramid store carries in its `dggs` attributes, and the
value [`ispyramidstore`](@ref) recognizes it by.
"""
const PYRAMID_LAYOUT = "implicit_pyramid"

"The layout revision this reader and writer understand."
const PYRAMID_LAYOUT_VERSION = 1

"The only `writer` whose pyramid stores this reader opens without a description."
const PYRAMID_WRITER = "DiscreteGlobalGrids.jl"

"The nested attribute object the layout is described in, inside `dggs`."
const PYRAMID_BLOCK = "pyramid_layout"

"""
    PYRAMID_LEVEL_PREFIX

What a level array is called inside a variable's group: this, followed by the
level number in decimal. `level0`, `level7`, `level13`.
"""
const PYRAMID_LEVEL_PREFIX = "level"

"""
    PYRAMID_ORDER

How the values of one level array are ordered: by ascending slot, which is the
level's canonical cell order with the absent addresses left in.
"""
const PYRAMID_ORDER = "ascending_slot"

"""
    PYRAMID_PADDING

What a slot naming no cell holds: the array's own `fill_value`, the same value
an unwritten chunk reads back as. The two are deliberately one thing — a reader
asking "is there data here" gets one answer, and the coarse level can be the
chunk index.
"""
const PYRAMID_PADDING = "absent_slot_fill"

"""
    DEFAULT_PYRAMID_CHUNK_EXPONENT

The `k` in `aᵏ` that `chunks = :auto` picks: 6, which on IGEO7 is 117 649 cells
per chunk — a megabyte at eight bytes a cell, and the same order as the
one-dimensional layout's element target.
"""
const DEFAULT_PYRAMID_CHUNK_EXPONENT = 6

# ===========================================================================
# The layout
# ===========================================================================

"""
    PyramidLayout(system, level; gridname, chunkexponent)

The shape of a pyramid store: every level from the system's root level to
`level`, each on its own slot space, each chunked at a power of the aperture.

  - `level` is the finest level stored, and the level the data was written at.
  - `chunkexponent` is the `k` of `aᵏ`; a level shallower than `k` is one chunk,
    so the exponent a level really uses is `min(k, J)` and its chunk count is
    `slotcount(system, J - min(k, J))`.

The mapping in both directions is arithmetic on one slot number:

```julia
slotchunk(layout, J, slot)      # (chunk, offset within it)
chunkslots(layout, J, chunk)    # where the chunk's real cells sit in it
chunkroot(layout, J, chunk)     # the coarse cell the chunk is the subtree of
```

Only systems with [`has_sorted_subtrees`](@ref) and the slot verbs
([`slotcount`](@ref DiscreteGlobalGrids.slotcount),
[`slotindex`](@ref DiscreteGlobalGrids.slotindex),
[`slotcell`](@ref DiscreteGlobalGrids.slotcell)) can be laid out this way.
"""
struct PyramidLayout{S}
    system::S
    gridname::String
    level::Int
    rootlevel::Int
    aperture::Int
    chunkexponent::Int
end

function PyramidLayout(sys::AbstractHierarchicalGridSystem, level::Integer;
    gridname::Union{AbstractString,Nothing}=nothing,
    chunkexponent::Integer=DEFAULT_PYRAMID_CHUNK_EXPONENT)

    has_sorted_subtrees(sys) || throw(ArgumentError(
        "$(nameof(typeof(sys))) has no descendant ranges, so a run of slots is " *
        "not a subtree and cannot be a chunk: the pyramid layout needs " *
        "`has_sorted_subtrees`."))
    root = first(levels(sys))
    L = Int(level)
    root <= L <= maxlevel(sys) || throw(ArgumentError(
        "a pyramid store holds levels $root through its own finest level, and " *
        "$L is outside $(levels(sys))."))
    k = Int(chunkexponent)
    k >= 0 || throw(ArgumentError(
        "the chunk exponent is the `k` of `aᵏ` and is at least zero, not $k."))

    a = _aperture(sys, root)
    _checkslotspace(sys, root, L, a)
    name = gridname === nothing ? gridnamefor(sys) : String(gridname)
    return PyramidLayout{typeof(sys)}(sys, name, L, root, a, k)
end

# The aperture, read off the slot space rather than asked for: a level has `a`
# times the slots of the one above it, which is what "aperture" means to this
# layout and the only thing it uses the number for.
function _aperture(sys, root::Int)
    n0 = slotcount(sys, root)
    n0 >= 1 || throw(ArgumentError(
        "$(nameof(typeof(sys)))'s root level has $n0 slots."))
    n1 = slotcount(sys, root + 1)
    a, r = divrem(n1, n0)
    (r == 0 && a >= 2) || throw(ArgumentError(
        "the pyramid layout chunks at a power of the aperture, and " *
        "$(nameof(typeof(sys))) has $n1 slots one level under $n0, which is " *
        "not a whole-number refinement."))
    return a
end

# A slot space that is not `slotcount(root)·a^d` has no power-of-`a` chunk grid
# that divides it, and every chunk calculation below would be off by a
# remainder. Checked once, over the levels the store will actually hold.
function _checkslotspace(sys, root::Int, L::Int, a::Int)
    n = slotcount(sys, root)
    for J in root:L
        want = n * a^(J - root)
        slotcount(sys, J) == want || throw(ArgumentError(
            "the pyramid layout needs each level's slot space to be the one " *
            "above it times the aperture: $(nameof(typeof(sys))) level $J has " *
            "$(slotcount(sys, J)) slots where $want was implied by the aperture " *
            "$a and level $root's $n."))
    end
    return nothing
end

function Base.:(==)(a::PyramidLayout, b::PyramidLayout)
    return a.system == b.system && a.gridname == b.gridname && a.level == b.level &&
           a.rootlevel == b.rootlevel && a.chunkexponent == b.chunkexponent
end

function Base.show(io::IO, l::PyramidLayout)
    print(io, "PyramidLayout(", l.gridname, ", levels ", l.rootlevel, ":", l.level,
        ", chunks ", l.aperture, "^", l.chunkexponent, ")")
end

system(l::PyramidLayout) = l.system
level(l::PyramidLayout) = l.level

"""
    levels(layout::PyramidLayout) -> UnitRange{Int}

Every level the store holds, coarsest first.
"""
levels(l::PyramidLayout) = l.rootlevel:l.level

"""
    levelname(layout, J) -> String

The array name level `J` is stored under inside a variable's group: `"level7"`.
"""
levelname(::PyramidLayout, J::Integer) = string(PYRAMID_LEVEL_PREFIX, Int(J))

"""
    levelfromname(name) -> Union{Int,Nothing}

The level `name` names, or `nothing` where it is not a level array's name at
all. `"level7"` is 7; `"level07"`, `"levels"` and `"cell_ids"` are `nothing`.
"""
function levelfromname(name::AbstractString)
    s = String(name)
    startswith(s, PYRAMID_LEVEL_PREFIX) || return nothing
    digits = SubString(s, ncodeunits(PYRAMID_LEVEL_PREFIX) + 1)
    isempty(digits) && return nothing
    all(isdigit, digits) || return nothing
    # A leading zero is a second spelling of a level that already has one, and
    # two spellings in one group is a store no reader can resolve.
    (length(digits) > 1 && first(digits) == '0') && return nothing
    return parse(Int, digits)
end

"""
    levelslots(layout, J) -> Int

The length of level `J`'s array: its slot space, cells and absent addresses
together.
"""
levelslots(l::PyramidLayout, J::Integer) = slotcount(l.system, _checklevel(l, J))

"""
    chunkexponent(layout, J) -> Int

The exponent level `J` is really chunked at, `min(k, J - rootlevel)`: a level
with fewer than `aᵏ` slots per root cell is one chunk per root cell rather than
a chunk that would not divide it.
"""
chunkexponent(l::PyramidLayout, J::Integer) =
    min(l.chunkexponent, _checklevel(l, J) - l.rootlevel)

"""
    chunklength(layout, J) -> Int

How many slots one chunk of level `J` holds: `aᵏ` at
[`chunkexponent`](@ref DiscreteGlobalGrids.chunkexponent)`(layout, J)`.
"""
chunklength(l::PyramidLayout, J::Integer) = l.aperture^chunkexponent(l, J)

"""
    chunkrootlevel(layout, J) -> Int

The level whose slots index level `J`'s chunks: `J - chunkexponent(layout, J)`.
Reading that level tells a reader which of level `J`'s chunks were written,
because a slot there is fill exactly when the chunk under it is.
"""
chunkrootlevel(l::PyramidLayout, J::Integer) = _checklevel(l, J) - chunkexponent(l, J)

"""
    chunkcount(layout, J) -> Int

How many chunks level `J` is cut into, which is
[`levelslots`](@ref DiscreteGlobalGrids.levelslots) of its
[`chunkrootlevel`](@ref DiscreteGlobalGrids.chunkrootlevel).
"""
chunkcount(l::PyramidLayout, J::Integer) = levelslots(l, chunkrootlevel(l, J))

function _checklevel(l::PyramidLayout, J::Integer)
    j = Int(J)
    l.rootlevel <= j <= l.level || throw(ArgumentError(
        "this store holds levels $(l.rootlevel):$(l.level), and $j is not one of them."))
    return j
end

# ===========================================================================
# Slots, chunks and subtrees
# ===========================================================================

"""
    slotchunk(layout, J, slot) -> (chunk, offset)

The chunk of level `J` that `slot` falls in, and its one-based offset inside it.
"""
function slotchunk(l::PyramidLayout, J::Integer, slot::Integer)
    n = levelslots(l, J)
    1 <= slot <= n || throw(BoundsError(l, Int(slot)))
    w = chunklength(l, J)
    q, r = divrem(Int(slot) - 1, w)
    return q + 1, r + 1
end

"""
    chunkroot(layout, J, chunk) -> Union{AbstractCellIndex,Nothing}

The cell whose subtree chunk `chunk` of level `J` is — at level
[`chunkrootlevel`](@ref DiscreteGlobalGrids.chunkrootlevel) and slot `chunk` —
or `nothing` where that slot names no cell, in which case the chunk holds no
cell either and is never written.
"""
function chunkroot(l::PyramidLayout, J::Integer, chunk::Integer)
    rl = chunkrootlevel(l, J)
    1 <= chunk <= chunkcount(l, J) || throw(BoundsError(l, Int(chunk)))
    return slotcell(l.system, rl, Int(chunk))
end

"""
    chunkindices(layout, J, chunk) -> UnitRange{Int}

The indices of the complete level-`J` grid whose cells live in `chunk` — the
chunk root's [`descendant_range`](@ref DiscreteGlobalGrids.descendant_range), and
empty where the chunk holds no cell.

This is what makes a chunk one contiguous piece of the cell axis, and the reason
a read of a range of cells is one read per chunk it spans.
"""
function chunkindices(l::PyramidLayout, J::Integer, chunk::Integer)
    a = chunkroot(l, J, chunk)
    a === nothing && return 1:0
    return descendant_range(l.system, a, Int(J))
end

"""
    chunkslots(layout, J, chunk) -> AbstractVector{Int}

Where chunk `chunk`'s real cells sit inside it: one-based offsets, ascending,
one per cell of [`chunkindices`](@ref DiscreteGlobalGrids.chunkindices) and in the
same order.

`1:chunklength(layout, J)` wherever the subtree is complete, which is every
chunk but the twelve rooted on a pentagon; those are walked once, per chunk, and
are the only place the layout pays for its absent addresses.
"""
function chunkslots(l::PyramidLayout, J::Integer, chunk::Integer)
    w = chunklength(l, J)
    r = chunkindices(l, J, chunk)
    isempty(r) && return 1:0
    length(r) == w && return 1:w
    grid = levelgrid(l.system, Int(J))
    base = (Int(chunk) - 1) * w
    out = Vector{Int}(undef, length(r))
    for (t, p) in enumerate(r)
        out[t] = slotindex(l.system, cellindex(grid, p)) - base
    end
    return out
end

"""
    cellslot(layout, cell) -> Int

The slot `cell` occupies in its own level's array — [`slotindex`](@ref DiscreteGlobalGrids.slotindex),
checked against the levels this store holds.
"""
function cellslot(l::PyramidLayout, c::AbstractCellIndex)
    _checklevel(l, level(c))
    return slotindex(l.system, c)
end

# ===========================================================================
# Cell axis -> chunks
# ===========================================================================

"""
    PyramidRun(chunk, rows, axis)

One contiguous piece of a cube's cell axis that lands inside one chunk: `axis`
indices of the cube map onto `rows` of chunk `chunk`, in order.

`rows` counts CELLS of the chunk, not slots — position `i` of
[`chunkslots`](@ref DiscreteGlobalGrids.chunkslots) — so a run says nothing about
the absent addresses between them and a writer looks the slots up once per chunk
rather than once per cell.
"""
struct PyramidRun
    chunk::Int
    rows::UnitRange{Int}
    axis::UnitRange{Int}
end

Base.:(==)(a::PyramidRun, b::PyramidRun) =
    a.chunk == b.chunk && a.rows == b.rows && a.axis == b.axis

Base.show(io::IO, r::PyramidRun) =
    print(io, "PyramidRun(", r.chunk, ", ", r.rows, ", ", r.axis, ")")

"""
    pyramid_runs(layout, J, cells) -> Vector{PyramidRun}

The chunks of level `J` a cube's cell axis covers, and which slice of the cube
goes into each.

`cells` is a [`CellLookup`](@ref DiscreteGlobalGrids.CellLookups.CellLookup), a
[`CellVector`](@ref), or any strictly ascending vector of level-`J` cells. The
first two are walked through their INDEX WINDOWS, so a cube of tens of millions
of cells is planned without one id being materialized.

Unlike the ancestor-subzone layout, a chunk need NOT be covered whole: a chunk
is written from a buffer prefilled with the array's fill value, so the cells the
axis does not name simply stay fill and read back as absent.
"""
pyramid_runs(l::PyramidLayout, J::Integer, lk::AbstractCellLookup) =
    pyramid_runs(l, J, parent(lk))

pyramid_runs(l::PyramidLayout, J::Integer, cv::AbstractCellVector) =
    pyramid_runs(l, J, region(cv))

function pyramid_runs(l::PyramidLayout, J::Integer, cv::CellVector)
    level(cv) == J || throw(ArgumentError(
        "the cube's cells are at level $(level(cv)) and this is the plan for level $J."))
    system(cv) == l.system || throw(ArgumentError(
        "the cube's cells are $(nameof(typeof(system(cv)))) cells and the store " *
        "holds $(nameof(typeof(l.system))) cells."))
    return _pyramid_runs(l, Int(J), Engine.intervals(cv.windows))
end

pyramid_runs(l::PyramidLayout, J::Integer, cells::AbstractVector) =
    _pyramid_runs(l, Int(J), _pyramid_intervals(l, Int(J), cells))

"""
    pyramid_index_runs(layout, J, indices) -> Vector{PyramidRun}

[`pyramid_runs`](@ref) from INDICES of the complete level-`J` grid rather than
from cells, strictly ascending.

This is what a level built by [`pyramid_reduce`](@ref) is written from: the
reduction names its parents by index, and turning those back into cells only to
turn them into indices again would be the same arithmetic twice.
"""
function pyramid_index_runs(l::PyramidLayout, J::Integer,
    indices::AbstractVector{<:Integer})
    _checklevel(l, J)
    ivs = Tuple{Int,Int}[]
    previous = 0
    for x in indices
        p = Int(x)
        p > previous || throw(ArgumentError(
            "the indices must be strictly ascending; $p follows $previous."))
        if !isempty(ivs) && p == previous + 1
            ivs[end] = (ivs[end][1], p)
        else
            push!(ivs, (p, p))
        end
        previous = p
    end
    return _pyramid_runs(l, Int(J), ivs)
end

# Ascending, disjoint index intervals from an explicit cell vector. The fallback
# path: a `CellVector` already keeps these and hands them over for free.
function _pyramid_intervals(l::PyramidLayout, J::Int, cells::AbstractVector)
    grid = levelgrid(l.system, J)
    ivs = Tuple{Int,Int}[]
    previous = 0
    for c in cells
        p = _gridindex(l, grid, J, c)
        p > previous || throw(ArgumentError(
            "the cell axis must be strictly ascending; index $p follows $previous."))
        if !isempty(ivs) && p == previous + 1
            ivs[end] = (ivs[end][1], p)
        else
            push!(ivs, (p, p))
        end
        previous = p
    end
    return ivs
end

function _gridindex(l::PyramidLayout, grid, J::Int, c::AbstractCellIndex)
    level(c) == J || throw(ArgumentError(
        "the cell axis holds a level-$(level(c)) cell and this is the plan for level $J."))
    p = globalindex(grid, c)
    p === nothing && throw(ArgumentError("$c names no cell of level $J."))
    return p
end

_gridindex(l::PyramidLayout, grid, J::Int, x::Integer) =
    _gridindex(l, grid, J, idcell(grid, x))

function _pyramid_runs(l::PyramidLayout, J::Int, intervals)
    runs = PyramidRun[]
    grid = levelgrid(l.system, J)
    rl = chunkrootlevel(l, J)
    offset = 0
    for (lo, hi) in intervals
        p = lo
        while p <= hi
            a = ancestor(l.system, cellindex(grid, p), rl)
            chunk = slotindex(l.system, a)
            r = descendant_range(l.system, a, J)
            stop = min(hi, last(r))
            n = stop - p + 1
            row = p - first(r) + 1
            push!(runs, PyramidRun(chunk, row:(row+n-1), (offset+1):(offset+n)))
            offset += n
            p = stop + 1
        end
    end
    return runs
end

# ===========================================================================
# The overviews
# ===========================================================================

"""
    pyramid_reduce(layout, J, indices, values, aggregate, fill) -> (indices, values)

Level `J-1` from level `J`: the parents of the cells `indices` names, each
reduced from the children that carry a value.

`indices` are indices of the complete level-`J` grid, strictly ascending, and
`values` are theirs. A child equal to `fill` is ABSENT — it is what an unwritten
cell reads back as, so a store cannot tell one from the other and neither does
this — and a parent with no present child is left out of the result entirely
rather than emitted as fill.

That is the invariant the layout's chunk index rests on: a stored value means
there is something under it, all the way up. `aggregate` returning the fill
value would break it, and raises.
"""
function pyramid_reduce(l::PyramidLayout, J::Integer, indices::AbstractVector{<:Integer},
    values::AbstractVector, aggregate, fillvalue)

    j = _checklevel(l, J)
    j > l.rootlevel || throw(ArgumentError(
        "level $j is this store's root level and has no parents to reduce to."))
    length(indices) == length(values) || throw(ArgumentError(
        "there are $(length(indices)) indices and $(length(values)) values."))

    sys = l.system
    grid = levelgrid(sys, j)
    parentgrid = levelgrid(sys, j - 1)
    T = eltype(values)

    outidx = Int[]
    outval = T[]
    buffer = T[]
    i = firstindex(indices)
    last_i = lastindex(indices)
    while i <= last_i
        a = ancestor(sys, cellindex(grid, Int(indices[i])), j - 1)
        r = descendant_range(sys, a, j)
        empty!(buffer)
        while i <= last_i && Int(indices[i]) <= last(r)
            v = values[i]
            isequal(v, fillvalue) || push!(buffer, v)
            i += 1
        end
        isempty(buffer) && continue
        push!(outidx, _parentindex(parentgrid, a))
        push!(outval, _reduced(T, aggregate, buffer, fillvalue, a))
    end
    return outidx, outval
end

function _parentindex(parentgrid, a)
    p = globalindex(parentgrid, a)
    p === nothing && throw(ArgumentError(
        "$a names no cell of level $(level(parentgrid)); the cell axis reached " *
        "outside its own hierarchy."))
    return p
end

function _reduced(::Type{T}, aggregate, buffer, fillvalue, cell) where {T}
    raw = aggregate(buffer)
    v = _aggregateconvert(T, raw, aggregate, cell)
    isequal(v, fillvalue) && throw(ArgumentError(
        "aggregating the $(length(buffer)) values under $cell gave $(repr(v)), " *
        "which is this store's fill value and therefore reads back as `no data " *
        "here`. A pyramid's coarse levels are its chunk index: a parent carries " *
        "a value exactly when something under it does. Choose a fill value the " *
        "aggregate cannot produce."))
    return v
end

function _aggregateconvert(::Type{T}, raw, aggregate, cell) where {T}
    raw isa T && return raw
    try
        return convert(T, raw)
    catch err
        err isa InexactError || err isa MethodError || rethrow()
        throw(ArgumentError(
            "aggregating under $cell gave $(repr(raw))::$(typeof(raw)), which is " *
            "not a $T. Every level of a pyramid store carries the element type " *
            "the data was written at, so `aggregate` has to land back in it — " *
            "`$(aggregate)` does not for $T. Pass one that does, such as " *
            "`vs -> round($T, sum(vs) / length(vs))`, or write the cube as $T's " *
            "floating-point counterpart."))
    end
end

# ===========================================================================
# The attributes
# ===========================================================================

"""
    pyramid_attrs(layout; variables = String[], fill_value = "NaN",
                  aggregate = nothing) -> Dict{String,Any}

The group attributes a pyramid store carries: a `dggs` object naming the grid
and the finest level, with everything the layout adds nested under
`pyramid_layout`.

No `zarr_conventions` declaration and no cell coordinate are written. The layout
carries no array of ids at all — that is the point of it — so there is nothing
for the one-dimensional convention to describe and nothing for a generic xdggs
reader to fingerprint. A store in this layout is read by a reader that knows the
block below, and by no other.
"""
function pyramid_attrs(l::PyramidLayout; variables=String[], fill_value="NaN",
    aggregate=nothing)
    block = Dict{String,Any}(
        "layout" => PYRAMID_LAYOUT,
        "version" => PYRAMID_LAYOUT_VERSION,
        "writer" => PYRAMID_WRITER,
        "levels" => Int[J for J in levels(l)],
        "level_prefix" => PYRAMID_LEVEL_PREFIX,
        "aperture" => l.aperture,
        "chunk_exponent" => l.chunkexponent,
        "slot_order" => PYRAMID_ORDER,
        "padding" => PYRAMID_PADDING,
        "padding_fill_value" => fill_value,
        "variables" => String[String(v) for v in variables])
    aggregate === nothing || (block["aggregate"] = String(aggregate))
    dggs = Dict{String,Any}(
        "name" => l.gridname,
        "refinement_level" => l.level,
        PYRAMID_BLOCK => block)
    ref = get(GRID_REFERENCE, l.gridname, nothing)
    ref === nothing || (dggs["indexing_scheme"] = String(ref.idscheme))
    return Dict{String,Any}("dggs" => dggs)
end

"""
    pyramid_level_attrs(layout, J, dimension) -> Dict{String,Any}

The attributes one level array carries: its dimension name, the level it holds,
and how it is cut up. Documentation for a reader looking at the store by hand —
the group's block is what this package reads.
"""
pyramid_level_attrs(l::PyramidLayout, J::Integer, dimension::AbstractString) =
    Dict{String,Any}(
        "_ARRAY_DIMENSIONS" => [String(dimension)],
        "refinement_level" => Int(J),
        "slot_count" => levelslots(l, J),
        "chunk_length" => chunklength(l, J))

"""
    pyramid_dimension(layout, J) -> String

The name of level `J`'s dimension. Each level has its own length, so each needs
its own dimension name for a reader that indexes dimensions by name.
"""
pyramid_dimension(l::PyramidLayout, J::Integer) =
    string("slots_", levelname(l, J))

"""
    ispyramidstore(attrs) -> Bool

Whether a group's attributes describe a pyramid store. Never throws: a store
this reader does not recognize is read the ordinary way, and a malformed one
fails in [`pyramid_layout`](@ref) with its own message.
"""
function ispyramidstore(attrs)
    attrs isa AbstractDict || return false
    dggs = get(attrs, "dggs", nothing)
    dggs isa AbstractDict || return false
    block = get(dggs, PYRAMID_BLOCK, nothing)
    block isa AbstractDict || return false
    return get(block, "layout", nothing) == PYRAMID_LAYOUT
end

"""
    pyramid_layout(attrs; store = "") -> PyramidLayout

The layout a store's group attributes describe, checked against what this
package can read: the layout name, a version it understands, a registered grid
name, and a level list and aperture that are the ones this grid really has.

The declared `aperture` and `levels` are CHECKED rather than believed. They are
the store's claim about the grid, and a store whose claim disagrees is one
written against another grid definition — which would read every chunk at the
wrong offset, silently.
"""
function pyramid_layout(attrs; store::AbstractString="")
    ispyramidstore(attrs) || throw(DGGSFormatError(check=:not_a_pyramid_store,
        store=String(store),
        detail="this group carries no `dggs.$PYRAMID_BLOCK` object naming the " *
               "`$PYRAMID_LAYOUT` layout."))
    dggs = attrs["dggs"]
    block = dggs[PYRAMID_BLOCK]
    convs = String[]

    version = get(block, "version", nothing)
    version == PYRAMID_LAYOUT_VERSION || throw(DGGSFormatError(
        check=:unsupported_layout_version, store=String(store),
        declared=version, observed=PYRAMID_LAYOUT_VERSION,
        detail="this reader implements version $PYRAMID_LAYOUT_VERSION of the " *
               "$PYRAMID_LAYOUT layout."))

    name = lowercase(strip(String(_require(dggs, "name", store, convs, "the dggs object"))))
    ref = gridreference(name; store=store)
    lev = _asint(_require(dggs, "refinement_level", store, convs, "the dggs object"),
        "refinement_level", store, convs)
    k = _asint(_require(block, "chunk_exponent", store, convs, "the pyramid layout"),
        "chunk_exponent", store, convs)

    order = get(block, "slot_order", PYRAMID_ORDER)
    order == PYRAMID_ORDER || throw(DGGSFormatError(check=:unsupported_slot_order,
        store=String(store), declared=order, observed=PYRAMID_ORDER,
        detail="this reader reads a level array in ascending slot order only."))
    pad = get(block, "padding", PYRAMID_PADDING)
    pad == PYRAMID_PADDING || throw(DGGSFormatError(check=:unsupported_pyramid_padding,
        store=String(store), declared=pad, observed=PYRAMID_PADDING,
        detail="this reader reads a slot naming no cell as the array's own fill " *
               "value, which is also what an unwritten chunk reads back as."))

    l = PyramidLayout(ref.system, lev; gridname=name, chunkexponent=k)

    declared = get(block, "aperture", nothing)
    declared === nothing || declared == l.aperture || throw(DGGSFormatError(
        check=:aperture_mismatch, store=String(store), declared=declared,
        observed=l.aperture,
        detail="the store declares aperture $declared and this grid refines by " *
               "$(l.aperture); the store was written against another grid."))
    want = Int[J for J in levels(l)]
    got = get(block, "levels", nothing)
    got === nothing || Int[_asint(x, "levels", store, convs) for x in got] == want ||
        throw(DGGSFormatError(check=:level_list_mismatch, store=String(store),
            declared=got, observed=want,
            detail="a pyramid store holds every level from the system's root to " *
                   "its own finest level; this one lists something else."))
    prefix = get(block, "level_prefix", PYRAMID_LEVEL_PREFIX)
    prefix == PYRAMID_LEVEL_PREFIX || throw(DGGSFormatError(
        check=:unsupported_level_prefix, store=String(store), declared=prefix,
        observed=PYRAMID_LEVEL_PREFIX,
        detail="this reader finds a level array by the name " *
               "`$PYRAMID_LEVEL_PREFIX` followed by the level number."))
    return l
end

# ===========================================================================
# The cube a pyramid store reads back as
# ===========================================================================

"""
    StorePyramid(layout, stacks, fills, identifier)

Every level of a pyramid store, lazily, as one [`AbstractPyramid`](@ref).

    pyr[3]            # level 3, as a `DimStack` over that level's `Cells` axis
    pyr[:elevation]   # one variable, at every level, as a `StorePyramid` again
    keys(pyr)         # the variables
    levels(pyr)       # the levels
    pyr.layout        # the `PyramidLayout` the store declares

A level is an ordinary cube over the COMPLETE level in canonical order: the slot
space the store is laid out on stops at the store, and a cell axis holds cells.
Values the store never had read back as `missing` through
[`cellvalues`](@ref DiscreteGlobalGrids.cellvalues), and as the array's own fill
value through the cube.

[`holdsdata`](@ref DiscreteGlobalGrids.holdsdata) and
[`cellvalues`](@ref DiscreteGlobalGrids.cellvalues) answer for ONE variable, so
they are asked of `pyr[:elevation]` rather than of a pyramid holding several.
That is the only thing the single-variable view is for.

Built by [`dggread`](@ref DiscreteGlobalGrids.dggread) on a store whose
attributes [`ispyramidstore`](@ref) recognizes.
"""
struct StorePyramid{L<:PyramidLayout,S,F} <: AbstractPyramid
    layout::L
    stacks::S
    fills::F
    identifier::String
end

system(p::StorePyramid) = system(p.layout)
levels(p::StorePyramid) = levels(p.layout)
Base.keys(p::StorePyramid) = keys(p.fills)

function Base.getindex(p::StorePyramid, J::Integer)
    j = Int(J)
    j in levels(p) || throw(ArgumentError(
        "this pyramid holds levels $(levels(p)), and $j is not one of them."))
    return p.stacks[j-first(levels(p))+1]
end

function Base.getindex(p::StorePyramid, var::Symbol)
    haskey(p.fills, var) || throw(ArgumentError(
        "this pyramid has no variable `$var`; it holds " *
        join(string.(keys(p.fills)), ", ") * "."))
    stacks = [s[(var,)] for s in p.stacks]
    return StorePyramid(p.layout, stacks, NamedTuple{(var,)}((p.fills[var],)),
        p.identifier)
end

Base.show(io::IO, p::StorePyramid) = print(io, "StorePyramid(", p.identifier, ", ",
    p.layout, ", ", join(string.(keys(p)), ", "), ")")

# The one variable the reading verbs answer for. A pyramid over several is a
# container; `holdsdata` and `cellvalues` are about values, and there is no
# honest single answer over more than one.
function _onlyvar(p::StorePyramid)
    ks = keys(p.fills)
    length(ks) == 1 && return only(ks)
    throw(ArgumentError(
        "this pyramid holds " * join(string.(ks), ", ") * ", and `holdsdata` " *
        "and `cellvalues` answer for one variable. Select one first: " *
        "`pyr[:$(first(ks))]`."))
end

function holdsdata(p::StorePyramid, c::AbstractCellIndex)
    var = _onlyvar(p)
    J = level(c)
    J in levels(p) || throw(ArgumentError(
        "$c is at level $J and this pyramid holds levels $(levels(p))."))
    i = globalindex(levelgrid(system(p), J), c)
    # A well-formed id naming no cell holds nothing, which is the honest answer
    # rather than an error: a descent meets these where a subtree was deleted.
    i === nothing && return false
    v = only(_leveldata(p, J, var)[i:i])
    return !isequal(v, p.fills[var])
end

function cellvalues(p::StorePyramid, J::Integer, cells)
    var = _onlyvar(p)
    j = Int(J)
    j in levels(p) || throw(ArgumentError(
        "this pyramid holds levels $(levels(p)), and $j is not one of them."))
    sys = system(p)
    grid = levelgrid(sys, j)
    fillvalue = p.fills[var]
    data = _leveldata(p, j, var)
    T = eltype(data)

    n = length(cells)
    out = Vector{Union{T,Missing}}(undef, n)
    fill!(out, missing)
    n == 0 && return out

    indices = Vector{Int}(undef, n)
    for (i, c) in enumerate(cells)
        indices[i] = _valueindex(grid, j, c)
    end

    # One read per CHUNK touched, which is one store read: the cells of a frame
    # are scattered over the level but clustered in the tree, so a per-cell read
    # would fetch the same chunk hundreds of times.
    order = sortperm(indices)
    rl = chunkrootlevel(p.layout, j)
    i = 1
    while i <= n
        start = indices[order[i]]
        if start == 0
            i += 1
            continue
        end
        r = descendant_range(sys, ancestor(sys, cellindex(grid, start), rl), j)
        k = i
        while k <= n && indices[order[k]] <= last(r)
            k += 1
        end
        stop = indices[order[k-1]]
        block = data[start:stop]
        for t in i:(k-1)
            slot = order[t]
            v = block[indices[slot]-start+1]
            isequal(v, fillvalue) || (out[slot] = v)
        end
        i = k
    end
    return out
end

# Zero for a cell this level does not have, which `cellvalues` reports as
# `missing` rather than refusing: a descent over a deleted branch asks.
function _valueindex(grid, j::Int, c::AbstractCellIndex)
    level(c) == j || throw(ArgumentError(
        "`cellvalues` was asked for level $j and given the level-$(level(c)) cell $c."))
    i = globalindex(grid, c)
    return i === nothing ? 0 : i
end

# The array behind one variable at one level: a lazy store-backed vector, or a
# plain one where the pyramid was materialized. `parent` of a `DimArray` is
# whichever it is, and the reading verbs need nothing else from the cube.
_leveldata(p::StorePyramid, J::Integer, var::Symbol) = parent(p[J][var])
