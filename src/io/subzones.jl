# The ancestor-subzone layout: one column per coarse cell, one row per subzone.
#
# INTERIM. Zarr v2/v3 chunks are uniform, so a chunk grid that follows the tree
# exactly is not expressible on a one-dimensional cell axis: a pentagon subtree
# holds `p(d) = (5*7^d + 1)/6` cells where a hexagon's holds `7^d`. This layout
# buys tree-aligned chunking by spending a dimension on it, and goes away once
# Zarr supports variable chunk sizes.
#
#     dim 1 (fastest)   subzone position 1:capacity within one ancestor's subtree
#     dim 2             the ancestor, at position `i` of the complete level-La grid
#     chunks            (capacity, 1) — one chunk per ancestor column
#
# So a chunk is a subtree, an unwritten ancestor is a chunk that was never
# stored (Zarr reads it back as fill), and the twelve pentagon columns carry
# `p(d)` real values followed by fill. Position within a column is OGC API-DGGS
# SUB-ZONE ORDER, here ascending cell id.
#
# Arithmetic and vocabulary only: no Zarr, no arrays, and nothing that
# materializes a level-L id vector.
#
# Include order: after `conventions.jl`, whose `GRID_REFERENCE` this reads, and
# before `api.jl`, whose stubs document the verbs the Zarr extension implements.

# ===========================================================================
# The on-disk vocabulary
# ===========================================================================

"""
    SUBZONE_LAYOUT

The `layout` string an ancestor-subzone store carries in its `dggs` attributes,
and the value [`issubzonestore`](@ref) recognizes it by.
"""
const SUBZONE_LAYOUT = "ancestor_subzone"

"The layout revision this reader and writer understand."
const SUBZONE_LAYOUT_VERSION = 1

"The only `writer` whose subzone stores this reader opens without a description."
const SUBZONE_WRITER = "DiscreteGlobalGrids.jl"

"The nested attribute object the layout is described in, inside `dggs`."
const SUBZONE_BLOCK = "subzone_layout"

"The name of the fastest-varying (Julia dimension 1) axis: position in a subtree."
const SUBZONE_DIMENSION = "subzone"

"The name of the slowest-varying axis: the level-`ancestor_level` cell."
const ANCESTOR_DIMENSION = "ancestor"

"The optional coordinate array of level-`ancestor_level` ids."
const ANCESTOR_COORDINATE = "ancestor_cell_ids"

"""
    SUBZONE_ORDER

How the values inside one column are ordered: ascending cell id, which on a
system with sorted subtrees is ascending position within the ancestor's
[`descendant_range`](@ref).
"""
const SUBZONE_ORDER = "ascending_id"

"""
    SUBZONE_PADDING

What a column shorter than the array's row extent holds after its cells: the
array's own `fill_value`, at the END of the column. A pentagon's subtree is
`p(d) = (5*7^d + 1)/6` cells against a hexagon's `7^d`, so twelve columns of an
IGEO7 store are short; which twelve is derivable from the grid and is not
recorded.
"""
const SUBZONE_PADDING = "trailing_fill"

# ===========================================================================
# The layout
# ===========================================================================

"""
    SubzoneLayout(system, level, ancestor_level; gridname, capacity)

The shape of an ancestor-subzone store: a level-`level` cell axis cut into the
subtrees of the complete level-`ancestor_level` grid, one subtree per column.

  - `ncolumns` is `ncells(levelgrid(system, ancestor_level))`, and column `i` is
    the `i`th level-`ancestor_level` cell in canonical order. The column axis is
    IMPLICIT: position is the ancestor, and a store may carry the ids as a
    coordinate array for interop but is not read through it.
  - `capacity` is the number of ROWS the store needs, which is the longest
    column: `7^d` for IGEO7, `4^d` for a quad-face system. Shorter columns —
    the twelve pentagon-rooted ones — hold their cells first and fill after
    ([`SUBZONE_PADDING`](@ref)).

Both directions of the mapping are O(level) arithmetic on one cell id:

```julia
subzoneindex(layout, cell)      # (column, row)
columnpositions(layout, i)      # the level grid positions column i holds
columncell(layout, i)           # the ancestor cell itself
```

Only systems with [`has_sorted_subtrees`](@ref) can be laid out this way: a
column is a contiguous run of the level grid, and a system whose descendants are
scattered has no such run.

The default `capacity` is measured by one pass over the ancestor grid; pass it
where it is already known (`subzone_count` in a store's attributes).
"""
struct SubzoneLayout{S,G,A}
    system::S
    gridname::String
    level::Int
    ancestor_level::Int
    grid::G
    ancestorgrid::A
    capacity::Int
    ncolumns::Int
end

function SubzoneLayout(sys::AbstractHierarchicalGridSystem, level::Integer,
    ancestor_level::Integer; gridname::Union{AbstractString,Nothing}=nothing,
    capacity::Union{Integer,Nothing}=nothing)
    # The trait first: a system whose subtrees are scattered has no columns
    # whatever it is called, and hearing about its store spelling instead would
    # send the caller to register a name that would not help.
    has_sorted_subtrees(sys) || throw(ArgumentError(
        "$(nameof(typeof(sys))) has no descendant ranges, so a cell's subtree is " *
        "not one run of its level and cannot be a column: the ancestor-subzone " *
        "layout needs `has_sorted_subtrees`."))
    L, La = Int(level), Int(ancestor_level)
    0 <= La <= L || throw(ArgumentError(
        "an ancestor level is between 0 and the store's own level, and $La is " *
        "not in 0:$L."))
    grid = levelgrid(sys, L)
    agrid = levelgrid(sys, La)
    nc = Int(ncells(agrid))
    cap = capacity === nothing ? subzone_capacity(sys, La, L) : _checkedcapacity(capacity)
    name = gridname === nothing ? gridnamefor(sys) : String(gridname)
    return SubzoneLayout{typeof(sys),typeof(grid),typeof(agrid)}(
        sys, name, L, La, grid, agrid, cap, nc)
end

_checkedcapacity(capacity::Integer) = capacity >= 1 ? Int(capacity) : throw(ArgumentError(
    "a column holds at least one cell, so the row extent is at least one, not $capacity"))

function Base.:(==)(a::SubzoneLayout, b::SubzoneLayout)
    return a.system == b.system && a.gridname == b.gridname && a.level == b.level &&
           a.ancestor_level == b.ancestor_level && a.capacity == b.capacity &&
           a.ncolumns == b.ncolumns
end

function Base.show(io::IO, l::SubzoneLayout)
    print(io, "SubzoneLayout(", l.gridname, ", level ", l.level, " under level ",
        l.ancestor_level, ", ", l.ncolumns, " columns of ", l.capacity, ")")
end

system(l::SubzoneLayout) = l.system
level(l::SubzoneLayout) = l.level

"""
    subzone_depth(layout) -> Int

`d = level - ancestor_level`, the depth of one column's subtree.
"""
subzone_depth(l::SubzoneLayout) = l.level - l.ancestor_level

"""
    subzone_capacity(system, ancestor_level, level) -> Int

The row extent an ancestor-subzone store needs: the size of the LONGEST
level-`ancestor_level` subtree at `level`.

Measured rather than assumed, because the closed form is the grid's business and
not this layer's: aperture 7 gives `7^d` for a hexagon and `(5*7^d + 1)/6` for a
pentagon, a quad-face system gives `4^d` throughout, and a system yet to be
written gives whatever [`descendant_range`](@ref) says. The pass is one O(level)
range per level-`ancestor_level` cell — a fifth of a second for the 1 176 494
columns of an IGEO7 level-6 ancestor grid, and it happens once, when a store is
created. Every read takes the number out of the store's attributes instead.
"""
function subzone_capacity(sys::AbstractHierarchicalGridSystem, ancestor_level::Integer,
    level::Integer)
    agrid = levelgrid(sys, ancestor_level)
    cap = 0
    for i in 1:Int(ncells(agrid))
        cap = max(cap, length(descendant_range(sys, cellindex(agrid, i), level)))
    end
    return cap
end

"""
    gridnamefor(system) -> String

The canonical store spelling of `system`, read backwards out of
[`GRID_REFERENCE`](@ref).

A grid name pins the id packing, so a system with no registered name has no
store spelling either and this raises rather than inventing one.
"""
function gridnamefor(sys)
    for (name, ref) in GRID_REFERENCE
        ref.system == sys && return name
    end
    throw(DGGSFormatError(check=:unknown_grid_name, observed=nameof(typeof(sys)),
        detail="$(nameof(typeof(sys))) has no canonical store name; registered " *
               "names are " * join(sort!(collect(keys(GRID_REFERENCE))), ", ") *
               ". Add one with `register_grid!(name, GridReference(...))` before " *
               "writing."))
end

# ===========================================================================
# The mapping, both ways
# ===========================================================================

"""
    columncell(layout, i) -> AbstractCellIndex

The level-`ancestor_level` cell column `i` holds the subtree of.
"""
function columncell(l::SubzoneLayout, i::Integer)
    1 <= i <= l.ncolumns || throw(BoundsError(l, Int(i)))
    return cellindex(l.ancestorgrid, Int(i))
end

"""
    columnindex(layout, ancestor) -> Int

The column an ancestor cell occupies: its position in the complete
level-`ancestor_level` grid. Throws for a cell of another level.
"""
function columnindex(l::SubzoneLayout, a::AbstractCellIndex)
    level(a) == l.ancestor_level || throw(ArgumentError(
        "$a is a level-$(level(a)) cell; the columns of this store are the " *
        "level-$(l.ancestor_level) cells."))
    p = cellposition(l.ancestorgrid, a)
    p === nothing && throw(ArgumentError(
        "$a names no cell of levelgrid($(nameof(typeof(l.system))), $(l.ancestor_level))."))
    return p
end

"""
    columnpositions(layout, i) -> UnitRange{Int}

The positions of the complete level grid that column `i` holds — the ancestor's
[`descendant_range`](@ref), which is what makes a column one contiguous piece of
the cell axis.
"""
columnpositions(l::SubzoneLayout, i::Integer) =
    descendant_range(l.system, columncell(l, i), l.level)

"""
    columnlength(layout, i) -> Int

How many cells column `i` really holds: `capacity` for a full subtree and less
for a pentagon-rooted one, whose remaining rows are fill.
"""
columnlength(l::SubzoneLayout, i::Integer) = length(columnpositions(l, i))

"""
    subzoneindex(layout, cell) -> (column, row)

Where one level-`level` cell sits in the store: its ancestor's column, and its
one-based position in the ancestor's subtree in ascending id order.

O(level) digit arithmetic on the id, so a whole store is addressed without any
level-`level` id vector ever existing.
"""
function subzoneindex(l::SubzoneLayout, c::AbstractCellIndex)
    level(c) == l.level || throw(ArgumentError(
        "$c is a level-$(level(c)) cell; this store holds level $(l.level)."))
    p = cellposition(l.grid, c)
    p === nothing && throw(ArgumentError(
        "$c names no cell of levelgrid($(nameof(typeof(l.system))), $(l.level))."))
    return positionindex(l, p)
end

"""
    positionindex(layout, p) -> (column, row)

[`subzoneindex`](@ref) from a POSITION of the complete level grid rather than
from a cell id — what a run of the cell axis is walked with.
"""
function positionindex(l::SubzoneLayout, p::Integer)
    c = cellindex(l.grid, Int(p))
    a = ancestor(l.system, c, l.ancestor_level)
    i = cellposition(l.ancestorgrid, a)
    i === nothing && throw(ArgumentError(
        "the level-$(l.ancestor_level) ancestor $a of position $p names no cell " *
        "of its own level grid."))
    r = descendant_range(l.system, a, l.level)
    return i, Int(p) - first(r) + 1
end

# ===========================================================================
# Cell axis -> columns
# ===========================================================================

"""
    SubzoneRun(column, rows, axis)

One contiguous piece of a cube's cell axis that lands inside one column:
`axis` positions of the cube map onto `rows` of column `column`, in order.

`length(rows) == length(axis)` always; a run whose `rows` is not the whole
column is a partially covered subtree, which [`subzone_runs`](@ref) refuses.
"""
struct SubzoneRun
    column::Int
    rows::UnitRange{Int}
    axis::UnitRange{Int}
end

Base.:(==)(a::SubzoneRun, b::SubzoneRun) =
    a.column == b.column && a.rows == b.rows && a.axis == b.axis

"""
    subzone_runs(layout, cells; complete = true) -> Vector{SubzoneRun}

The columns a cube's cell axis covers, and which slice of the cube goes into
each.

`cells` is a [`CellLookup`](@ref), a [`CellVector`](@ref), or any strictly
ascending vector of level-`level` cell ids. The first two are walked through
their POSITION WINDOWS — one step per column touched, and no id is ever
materialized, which is what lets a land-only cube of tens of millions of cells
be planned in microseconds. A plain vector costs one `cellposition` per cell.

`complete = true` refuses a column the axis covers only part of. That is the
layout's central restriction and not an implementation limit: a column is one
chunk, a chunk is written whole, and a store whose columns were written from
partial coverage would read back with data cells indistinguishable from fill.
Ancestor-snapped coverage — the way [`covering`](@ref) and a multi-order query
name a region — satisfies it by construction.
"""
function subzone_runs(l::SubzoneLayout, lk::CellLookup; complete::Bool=true)
    return subzone_runs(l, parent(lk); complete=complete)
end

function subzone_runs(l::SubzoneLayout, cv::CellVector; complete::Bool=true)
    level(cv) == l.level || throw(ArgumentError(
        "the cube's cells are at level $(level(cv)) and the store holds level $(l.level)."))
    system(cv) == l.system || throw(ArgumentError(
        "the cube's cells are $(nameof(typeof(system(cv)))) cells and the store " *
        "holds $(nameof(typeof(l.system))) cells."))
    return _subzone_runs(l, Engine.intervals(cv.windows), complete)
end

function subzone_runs(l::SubzoneLayout, cells::AbstractVector; complete::Bool=true)
    return _subzone_runs(l, _position_intervals(l, cells), complete)
end

# Ascending, disjoint position intervals from an explicit id vector. The
# fallback path: `CellVector` already keeps these and hands them over for free.
function _position_intervals(l::SubzoneLayout, cells::AbstractVector)
    ivs = Tuple{Int,Int}[]
    previous = 0
    for c in cells
        p = _cellposition_checked(l, c)
        p > previous || throw(ArgumentError(
            "the cell axis must be strictly ascending; position $p follows $previous."))
        if !isempty(ivs) && p == previous + 1
            ivs[end] = (ivs[end][1], p)
        else
            push!(ivs, (p, p))
        end
        previous = p
    end
    return ivs
end

function _cellposition_checked(l::SubzoneLayout, c::AbstractCellIndex)
    level(c) == l.level || throw(ArgumentError(
        "the cell axis holds a level-$(level(c)) cell and the store holds level $(l.level)."))
    p = cellposition(l.grid, c)
    p === nothing && throw(ArgumentError("$c names no cell of level $(l.level)."))
    return p
end

_cellposition_checked(l::SubzoneLayout, x::Integer) =
    _cellposition_checked(l, idcell(l.grid, x))

function _subzone_runs(l::SubzoneLayout, intervals, complete::Bool)
    runs = SubzoneRun[]
    offset = 0
    for (lo, hi) in intervals
        p = lo
        while p <= hi
            i, row = positionindex(l, p)
            r = descendant_range(l.system, columncell(l, i), l.level)
            stop = min(hi, last(r))
            n = stop - p + 1
            push!(runs, SubzoneRun(i, row:(row+n-1), (offset+1):(offset+n)))
            offset += n
            p = stop + 1
        end
    end
    complete && _checkcomplete(l, runs)
    return runs
end

# A column is written whole or not at all. The check is per RUN — a column
# covered in two disconnected pieces produces two runs, neither of which is the
# whole column, and is caught by the same comparison.
@noinline function _checkcomplete(l::SubzoneLayout, runs)
    for run in runs
        h = columnlength(l, run.column)
        (first(run.rows) == 1 && length(run.rows) == h) || throw(DGGSFormatError(
            check=:incomplete_subtree, declared=h, observed=length(run.rows),
            detail="the ancestor-subzone layout writes a column whole: the cell " *
                   "axis covers rows $(run.rows) of column $(run.column) " *
                   "($(columncell(l, run.column))), which holds $h cells. Snap the " *
                   "cube's coverage to level-$(l.ancestor_level) subtrees — " *
                   "`covering` on a level-$(l.ancestor_level) cell vector expanded " *
                   "to level $(l.level) is one way — or write it in the " *
                   "one-dimensional cell layout instead."))
    end
    return nothing
end

# ===========================================================================
# Columns -> cell axis
# ===========================================================================

"""
    subzone_cellvector(layout, columns) -> CellVector

The cell axis a set of columns spells: their subtrees concatenated, in ascending
column order. `columns` is `nothing` for the whole store — every column of the
level, which is the complete level grid and one window.

This is the read side of [`subzone_runs`](@ref) and its exact inverse: writing a
cube's columns and reading them back gives the same axis.
"""
subzone_cellvector(l::SubzoneLayout, ::Nothing) = CellVector(l.grid)

function subzone_cellvector(l::SubzoneLayout, columns::AbstractVector{<:Integer})
    ranges = UnitRange{Int}[]
    previous = 0
    for i in columns
        i > previous || throw(ArgumentError(
            "columns are read in ascending order and each at most once; column " *
            "$i follows $previous."))
        previous = Int(i)
        r = columnpositions(l, i)
        # Windows are kept MAXIMAL — two `CellVector`s over the same cells must
        # compare equal, and `RangeWindows` compares boundaries — so adjacent
        # columns, which the whole-level case is nothing but, merge here.
        if !isempty(ranges) && first(r) == last(ranges[end]) + 1
            ranges[end] = first(ranges[end]):last(r)
        else
            push!(ranges, r)
        end
    end
    return Engine.CellVector(Engine._range_windows(ranges), l.grid, nothing, l.level)
end

subzone_cellvector(l::SubzoneLayout, columns) =
    subzone_cellvector(l, Int[Int(i) for i in columns])

"""
    subzone_columns(layout, ancestors) -> Vector{Int}

The column indices a set of ancestor cells names, ascending and unique.
`ancestors` may hold cells or column indices; `nothing` means every column.
"""
subzone_columns(l::SubzoneLayout, ::Nothing) = collect(1:l.ncolumns)

function subzone_columns(l::SubzoneLayout, ancestors)
    out = Int[_column(l, a) for a in ancestors]
    sort!(out)
    unique!(out)
    return out
end

_column(l::SubzoneLayout, a::AbstractCellIndex) = columnindex(l, a)

function _column(l::SubzoneLayout, i::Integer)
    1 <= i <= l.ncolumns || throw(ArgumentError(
        "column $i is outside the store's 1:$(l.ncolumns) columns."))
    return Int(i)
end

# ===========================================================================
# The attributes
# ===========================================================================

"""
    subzone_attrs(layout; variables = String[], coordinate = nothing,
                  fill_value = "NaN") -> Dict{String,Any}

The group attributes an ancestor-subzone store carries: a `dggs` object naming
the grid and the level, with everything the layout adds nested under
`subzone_layout`. No `zarr_conventions` declaration is written — this is not the
convention's one-dimensional layout.

`coordinate` names the array of level-`ancestor_level` ids where the store
carries one. The column axis is implicit either way and is never read through
it, so this is interop and provenance rather than structure.
"""
function subzone_attrs(l::SubzoneLayout; variables=String[],
    coordinate::Union{AbstractString,Nothing}=nothing, fill_value="NaN")
    block = Dict{String,Any}(
        "layout" => SUBZONE_LAYOUT,
        "version" => SUBZONE_LAYOUT_VERSION,
        "writer" => SUBZONE_WRITER,
        "ancestor_level" => l.ancestor_level,
        "ancestor_dimension" => ANCESTOR_DIMENSION,
        "ancestor_count" => l.ncolumns,
        "ancestor_order" => SUBZONE_ORDER,
        "subzone_dimension" => SUBZONE_DIMENSION,
        "subzone_count" => l.capacity,
        "subzone_order" => SUBZONE_ORDER,
        "padding" => SUBZONE_PADDING,
        "padding_fill_value" => fill_value,
        # Zarr's own order, outermost dimension first, which is the order the
        # `.zarray` files hold: one chunk is one whole column.
        "chunk_shape" => Any[1, l.capacity],
        "variables" => String[String(v) for v in variables])
    coordinate === nothing || (block["ancestor_coordinate"] = String(coordinate))
    dggs = Dict{String,Any}(
        "name" => l.gridname,
        "refinement_level" => l.level,
        SUBZONE_BLOCK => block)
    ref = get(GRID_REFERENCE, l.gridname, nothing)
    ref === nothing || (dggs["indexing_scheme"] = String(ref.idscheme))
    return Dict{String,Any}("dggs" => dggs)
end

"""
    issubzonestore(attrs) -> Bool

Whether a group's attributes describe an ancestor-subzone store. Never throws:
a store this reader does not recognize is read the ordinary way, and a malformed
one fails in [`subzone_layout`](@ref) with its own message.
"""
function issubzonestore(attrs)
    attrs isa AbstractDict || return false
    dggs = get(attrs, "dggs", nothing)
    dggs isa AbstractDict || return false
    block = get(dggs, SUBZONE_BLOCK, nothing)
    block isa AbstractDict || return false
    return get(block, "layout", nothing) == SUBZONE_LAYOUT
end

"""
    subzone_layout(attrs; store = "") -> SubzoneLayout

The layout a store's group attributes describe, checked against what this
package can read: the layout name, a version it understands, a registered grid
name, and a level pair whose column count and row extent are the ones the store
declares.

The declared `ancestor_count` and `subzone_count` are CHECKED against the
arithmetic rather than believed. They are the store's claim about the grid, and
a store whose claim disagrees is one written against another grid definition —
which would read every column at the wrong offset, silently.
"""
function subzone_layout(attrs; store::AbstractString="")
    issubzonestore(attrs) || throw(DGGSFormatError(check=:not_a_subzone_store,
        store=String(store),
        detail="this group carries no `dggs.$SUBZONE_BLOCK` object naming the " *
               "`$SUBZONE_LAYOUT` layout."))
    dggs = attrs["dggs"]
    block = dggs[SUBZONE_BLOCK]
    convs = String[]
    version = get(block, "version", nothing)
    version == SUBZONE_LAYOUT_VERSION || throw(DGGSFormatError(
        check=:unsupported_layout_version, store=String(store),
        declared=version, observed=SUBZONE_LAYOUT_VERSION,
        detail="this reader implements version $SUBZONE_LAYOUT_VERSION of the " *
               "$SUBZONE_LAYOUT layout."))
    name = lowercase(strip(String(_require(dggs, "name", store, convs, "the dggs object"))))
    ref = gridreference(name; store=store)
    lev = _asint(_require(dggs, "refinement_level", store, convs, "the dggs object"),
        "refinement_level", store, convs)
    anc = _asint(_require(block, "ancestor_level", store, convs, "the subzone layout"),
        "ancestor_level", store, convs)
    cap = _asint(_require(block, "subzone_count", store, convs, "the subzone layout"),
        "subzone_count", store, convs)
    l = SubzoneLayout(ref.system, lev, anc; gridname=name, capacity=cap)
    ncol = _asint(_require(block, "ancestor_count", store, convs, "the subzone layout"),
        "ancestor_count", store, convs)
    ncol == l.ncolumns || throw(DGGSFormatError(check=:ancestor_count_mismatch,
        store=String(store), declared=ncol, observed=l.ncolumns,
        detail="the store declares $ncol level-$anc columns and this grid has " *
               "$(l.ncolumns); the store was written against another grid."))
    order = get(block, "subzone_order", SUBZONE_ORDER)
    order == SUBZONE_ORDER || throw(DGGSFormatError(check=:unsupported_subzone_order,
        store=String(store), declared=order, observed=SUBZONE_ORDER,
        detail="this reader reads columns in ascending cell id only."))
    pad = get(block, "padding", SUBZONE_PADDING)
    pad == SUBZONE_PADDING || throw(DGGSFormatError(check=:unsupported_subzone_padding,
        store=String(store), declared=pad, observed=SUBZONE_PADDING,
        detail="this reader drops padding from the END of a short column only."))
    return l
end

"""
    subzone_coordinate(attrs) -> Union{String,Nothing}

The name of the ancestor-id coordinate array a subzone store carries, where it
carries one.
"""
function subzone_coordinate(attrs)
    issubzonestore(attrs) || return nothing
    name = get(attrs["dggs"][SUBZONE_BLOCK], "ancestor_coordinate", nothing)
    return name isa AbstractString ? String(name) : nothing
end
