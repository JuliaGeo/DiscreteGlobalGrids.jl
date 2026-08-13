# ---------------------------------------------------------------------------
# `H3System` and its level grids
#
# The canonical dense order, which everything else here is a consequence of:
#
#     cells are numbered base cell by base cell (0:121, the order
#     `getRes0Cells` already returns), and within a base cell by H3's own child
#     position — the digit path read as a number with pentagon gaps skipped.
#
# That is exactly raw-`UInt64` order restricted to one level (see `cell.jl`), so
# `cellindex` of a level grid comes out sorted, and a subtree occupies a
# contiguous block of positions. One cumulative table per resolution turns a
# position into a base cell with a single binary search, and libh3's
# `cellToChildPos` / `childPosToCell` do the O(1) within-base-cell half — pentagon
# holes included, which is the whole reason the ordinals are hole-free.
# ---------------------------------------------------------------------------

"""
    H3System() <: AbstractHierarchicalGridSystem

Uber's [H3](https://h3geo.org) grid: an aperture-7 hexagonal hierarchy on an
icosahedron, with twelve pentagons, over resolutions `0:15`.

Canonical id: [`H3Cell`](@ref). Geometry, location and adjacency all come from
libh3 itself through `H3_jll`, so this system agrees with every other H3
implementation cell for cell rather than approximating one.

# Conventions this system documents

  - **Boundaries** are libh3's `cellToBoundary` ring, unchanged: counter-clockwise
    seen from outside, implicitly closed, and 5 to 10 vertices — a cell that
    crosses an icosahedron face edge carries extra *distortion* vertices there,
    and they are part of the exact boundary rather than noise to be cleaned off.
  - **`cellat` ties** are libh3's own: the point is handed to `latLngToCell`,
    whose answer on a shared edge is deterministic and is the same answer every
    other H3 binding gives.
  - **Neighbour order** is rotational: [`neighbors`](@ref) is rings `1..k`
    concatenated outward, each ring counter-clockwise seen from outside, so
    [`ring`](@ref) is the tail block of [`neighbors`](@ref). The walk's starting
    direction is libh3's own and is deterministic per cell rather than a uniform
    compass bearing.
  - **Centroids** are libh3's `cellToLatLng`, the centre the hierarchy is
    actually built around.

`has_sorted_subtrees` is `true`: see the module docstring on the canonical order.
"""
struct H3System <: AbstractHierarchicalGridSystem end

"""
    H3Grid(level) <: AbstractGrid

The complete H3 grid at one resolution — build one with
[`levelgrid(H3System(), l)`](@ref levelgrid) rather than by calling this.

A lightweight descriptor: it stores the resolution and nothing else, so
constructing the res-15 grid (569,707,381,193,162 cells) is free.
"""
struct H3Grid <: AbstractGrid
    level::Int
end

system(::H3Grid) = H3System()
level(g::H3Grid) = g.level

Base.show(io::IO, g::H3Grid) = print(io, "H3Grid(res ", g.level, ")")
Base.show(io::IO, ::H3System) = print(io, "H3System()")

# ===========================================================================
# The base tessellation, and the per-resolution prefix sums
#
# Both tables are pure Julia, deliberately: they are computed at precompile
# time, and reaching into a JLL there is the kind of thing that works until it
# does not. Every entry is checked against libh3 in `test/systems/H3/`, which is
# where the oracle belongs.
# ===========================================================================

# The twelve pentagon base cells. Fixed by H3's icosahedron orientation; the
# test suite asserts this tuple against `getPentagons(0)`.
const PENTAGON_BASE_CELLS = (4, 14, 24, 38, 49, 58, 63, 72, 83, 97, 107, 117)

# A res-0 index: mode 1, resolution 0, base cell `b`, all fifteen digits 7.
_res0_index(b::Integer) = (UInt64(1) << 59) | (UInt64(b) << 45) | UInt64(0x1FFFFFFFFFFF)

const _H3_ROOT_IDS = ntuple(i -> _res0_index(i - 1), 122)

# Descendants of one base cell at resolution `r`. A hexagon's subtree is a
# clean `7^r`; a pentagon deletes one child at every level, which telescopes to
# `1 + 5(7^r - 1)/6`.
_base_cell_descendants(b::Int, r::Int) =
    b in PENTAGON_BASE_CELLS ? 1 + 5 * (7^r - 1) ÷ 6 : 7^r

# `_H3_ROOT_ENDS[r + 1][j]` is the number of res-`r` cells in base cells
# `0:(j-1)` — the cumulative table a position is binary-searched against.
const _H3_ROOT_ENDS = ntuple(MAX_RESOLUTION + 1) do i
    r = i - 1
    return cumsum([_base_cell_descendants(b, r) for b in 0:121])
end

# ===========================================================================
# Required system interface
# ===========================================================================

cellindextype(::H3System) = H3Cell

levels(::H3System) = 0:MAX_RESOLUTION

function levelgrid(::H3System, l::Integer)
    lvl = Int(l)
    0 <= lvl <= MAX_RESOLUTION || throw(ArgumentError(
        "H3 resolution $lvl is outside levels(H3System()) = 0:$MAX_RESOLUTION"))
    return H3Grid(lvl)
end

"""
    rootcells(::H3System)

The 122 base cells, ascending. These are H3's res-0 tessellation: 110 hexagons
and 12 pentagons.
"""
rootcells(::H3System) = [H3Cell(id) for id in _H3_ROOT_IDS]

"""
    parent(::H3System, c::H3Cell) -> H3Cell

The cell one resolution coarser, by `cellToParent` — digit arithmetic, no
lookup. Throws an `ArgumentError` on a res-0 cell, which has no parent.
"""
function Base.parent(::H3System, c::H3Cell)
    l = level(c)
    l > 0 || throw(ArgumentError(
        "H3 cell $c is a root cell (resolution 0) and has no parent"))
    return H3Cell(H3Native.cell_to_parent(c.id, l - 1))
end

"""
    children(::H3System, c::H3Cell) -> SmallVector{7,H3Cell}

The seven children of a hexagon, or the six of a pentagon, ascending.

The container is fixed-capacity and the call **does not allocate**: libh3
writes the children into a stack buffer
([`H3Native.cell_to_children_7`](@ref)), which matters because tree descent
calls this once per node. A pentagon simply returns six of the seven slots,
which is the whole pentagon special case — generic code that reads `length`
rather than assuming the aperture needs nothing else.
"""
function children(::H3System, c::H3Cell)
    l = level(c)
    l < MAX_RESOLUTION || throw(ArgumentError(
        "H3 cell $c is at max_level $MAX_RESOLUTION and has no children"))
    out = SmallVector{7,H3Cell}()
    # `cellToChildren` emits ascending indices (digit 0 first, deleted digits
    # skipped) into the leading slots, so this preserves the canonical order
    # rather than imposing one. A pentagon leaves the last slot `0`, which is
    # never a valid H3 index.
    for id in H3Native.cell_to_children_7(c.id, l + 1)
        id == 0 && continue
        out = SmallCollections.push(out, H3Cell(id))
    end
    return out
end

# ===========================================================================
# Traits
# ===========================================================================

# The canonical order is base-cell-major then child-position, and a subtree is
# a contiguous run of child positions -- see the header of this file.
has_sorted_subtrees(::H3System) = true

"""
    cap_inflation(::H3System) -> Float64

`1.2`, the generic default, kept because it is measured to be sound here.

Under aperture 7 the children of a hexagon overhang their parent, so a cell's
own cap is not a legal [`node_extent`](@ref). The overhang converges: sampling
every base cell and descending nine levels along the outermost branch puts the
worst ratio of descendant-vertex distance to the cell's own cap radius at
**1.0522**, with per-level increments already down to ~1e-5 and shrinking
geometrically. `1.2` clears that by 14%.

`test/systems/H3/` re-measures the overhang rather than trusting this comment,
but deliberately re-measures it *weakly* — a beam of 60 over 6 levels, asserted
under 1.10 — so the suite stays fast. That 1.10 is a test threshold, not the
converged figure; the offline 1.0522 above is. The committed check would catch
a factor that had become unsound by a wide margin, and the separate covering-law
testset walks a chain to max_level and checks containment directly, which is
the property that actually matters.
"""
cap_inflation(::H3System) = 1.2

"""
    max_neighbors(::H3System, connectivity) -> Int

`6`, for either connectivity.

H3's cells are hexagons and pentagons, where sharing a vertex and sharing an
edge are the same relation — three cells meet at every vertex and any two of
them already share an edge — so [`Vertex()`](@ref Vertex) and [`Edge()`](@ref Edge)
coincide, and the bound is the hexagon's six. The twelve pentagons have five.
"""
max_neighbors(::H3System, ::Connectivity=Vertex()) = 6

# ===========================================================================
# The dense order: positions <-> ids
# ===========================================================================

function ncells(grid::H3Grid)
    return @inbounds _H3_ROOT_ENDS[grid.level+1][122]
end

"""
    cellindex(grid::H3Grid, i::Int) -> H3Cell

The id at position `i`: one binary search of the base-cell prefix sums, then
libh3's `childPosToCell` for the position within that base cell.
"""
function cellindex(grid::H3Grid, i::Int)
    n = ncells(grid)
    1 <= i <= n || throw(BoundsError(grid, i))
    r = grid.level
    ends = @inbounds _H3_ROOT_ENDS[r+1]
    # First base cell whose cumulative count reaches `i`; counts are strictly
    # increasing (every base cell has at least one descendant), so there are no
    # ties to break.
    b = searchsortedfirst(ends, i)
    previous = b == 1 ? 0 : @inbounds ends[b-1]
    root = @inbounds _H3_ROOT_IDS[b]
    return H3Cell(H3Native.child_pos_to_cell(i - previous - 1, root, r))
end

"""
    cellposition(grid::H3Grid, c::H3Cell) -> Union{Int,Nothing}

The position of `c`, or `nothing` when `c` is not a cell of this grid — a
different resolution, or not a valid index at all.

`nothing` rather than an error is the contract, and it is what makes asking
"is this cell here?" the normal way to intersect an id set with a grid. The
validity check is not paranoia: libh3 will happily compute a child position for
a malformed index, and returning a confident wrong position for one is worse
than returning nothing.
"""
function cellposition(grid::H3Grid, c::H3Cell)
    level(c) == grid.level || return nothing
    H3Native.is_valid_cell(c.id) || return nothing
    b = H3Native.get_base_cell(c.id)
    ends = @inbounds _H3_ROOT_ENDS[grid.level+1]
    previous = b == 0 ? 0 : @inbounds ends[b]
    return Int(previous + H3Native.cell_to_child_pos(c.id, 0) + 1)
end

# ===========================================================================
# Derived hierarchy: closed-form overrides of the generic walks
# ===========================================================================

"""
    ancestor(::H3System, c::H3Cell, l::Integer) -> H3Cell

The ancestor at resolution `l`, in one `cellToParent` call rather than
`level(c) - l` of them.
"""
function ancestor(::H3System, c::H3Cell, l::Integer)
    target = Int(l)
    lc = level(c)
    target <= lc || throw(ArgumentError(
        "ancestor level $target is deeper than the cell's own level $lc"))
    target >= 0 || throw(ArgumentError(
        "ancestor level $target is above the root level 0"))
    target == lc && return c
    return H3Cell(H3Native.cell_to_parent(c.id, target))
end

"""
    descendants(::H3System, c::H3Cell, l::Integer) -> Vector{H3Cell}

Every descendant at resolution `l`, ascending. `cellToChildren` spans any
number of levels in one call and already emits ascending indices, so this never
expands level by level.

O(subtree) and materialising, as the contract says — reach for
[`descendant_range`](@ref) instead wherever positions will do.
"""
function descendants(::H3System, c::H3Cell, l::Integer)
    target = Int(l)
    lc = level(c)
    target >= lc || throw(ArgumentError(
        "descendant level $target is above the cell's own level $lc"))
    target <= MAX_RESOLUTION || throw(ArgumentError(
        "descendant level $target is past max_level $MAX_RESOLUTION"))
    target == lc && return H3Cell[c]
    return [H3Cell(id) for id in H3Native.cell_to_children(c.id, target)]
end

"""
    descendant_range(::H3System, c::H3Cell, l::Integer) -> UnitRange{Int}

The contiguous interval of **positions** in `levelgrid(H3System(), l)` that the
descendants of `c` at resolution `l` occupy.

Two calls into libh3 and no enumeration: the subtree's first descendant is
child position 0 under `c`, and `cellToChildrenSize` counts the rest in closed
form with the pentagon gaps already taken out. Because the canonical order is
base-cell-major and child positions are contiguous within a base cell, the
descendants of `c` are exactly the positions in `[first, first + count)` — no
holes, including across pentagons, where the missing digit paths are absent
from the numbering rather than skipped over inside it.
"""
function descendant_range(sys::H3System, c::H3Cell, l::Integer)
    target = Int(l)
    lc = level(c)
    target >= lc || throw(ArgumentError(
        "descendant level $target is above the cell's own level $lc"))
    target <= MAX_RESOLUTION || throw(ArgumentError(
        "descendant level $target is past max_level $MAX_RESOLUTION"))
    grid = levelgrid(sys, target)
    # The `l == level(c)` case is the cell's own one-element range: window
    # descent in `HierarchicalGridCursor` asks for it at the level above the
    # leaves and must not get an exception.
    if target == lc
        p = cellposition(grid, c)
        p === nothing && throw(ArgumentError("$c is not a valid H3 cell"))
        return p:p
    end
    first_child = H3Cell(H3Native.child_pos_to_cell(0, c.id, target))
    p = cellposition(grid, first_child)
    p === nothing && throw(ArgumentError("$c is not a valid H3 cell"))
    count = Int(H3Native.cell_to_children_size(c.id, target))
    return p:(p+count-1)
end
