# Canonical order is base-cell-major, then H3 child position with pentagon gaps
# omitted. Per-resolution prefix tables locate the base cell.

"""
    H3System() <: AbstractHierarchicalGridSystem

H3's aperture-7 hexagonal hierarchy with twelve pentagons at resolutions
`0:15`. Libh3 supplies geometry, location, and adjacency. Boundaries are
implicitly closed counter-clockwise rings and include distortion vertices.
Canonical ordering makes `has_sorted_subtrees` true.
"""
struct H3System <: AbstractHierarchicalGridSystem end

# Grid descriptor for all cells at one H3 resolution.
const LevelGrid = HierarchicalLevelGrid{H3System}

Base.show(io::IO, ::H3System) = print(io, "H3System()")

# ===========================================================================
# The base tessellation, and the per-resolution prefix sums
#
# These tables permit precompilation without loading libh3.
# ===========================================================================

# The twelve pentagon base cells in H3's icosahedron orientation.
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
# `0:(j-1)` — the cumulative table an index is binary-searched against.
const _H3_ROOT_ENDS = ntuple(MAX_RESOLUTION + 1) do i
    r = i - 1
    return cumsum([_base_cell_descendants(b, r) for b in 0:121])
end

# ===========================================================================
# Required system interface
# ===========================================================================

cellindextype(::H3System) = H3Cell

levels(::H3System) = 0:MAX_RESOLUTION

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

The seven children of a hexagon, or the six of a pentagon, in ascending id
order. The result uses a fixed-capacity vector.
"""
function children(::H3System, c::H3Cell)
    l = level(c)
    l < MAX_RESOLUTION || throw(ArgumentError(
        "H3 cell $c is at maxlevel $MAX_RESOLUTION and has no children"))
    out = SmallVector{7,H3Cell}()
    # `cellToChildren` emits ascending ids with deleted digits skipped. A
    # pentagon leaves the final slot as the invalid sentinel `0`.
    for id in H3Native.cell_to_children_7(c.id, l + 1)
        id == 0 && continue
        out = SmallCollections.push(out, H3Cell(id))
    end
    return out
end

# ===========================================================================
# Traits
# ===========================================================================

# Base-cell-major child-position order makes each subtree contiguous.
has_sorted_subtrees(::H3System) = true

# `latLngToCell` is libh3's own inverse, so location reads the point.
has_direct_location(::H3System) = true

"""
    cap_inflation(::H3System) -> Float64

`1.2`, the generic default.

Children overhang their parents, so [`node_extent`](@ref) must be inflated. The
measured maximum ratio of a descendant *boundary vertex*'s distance from an
ancestor's cell-cap centre to that cap's radius is `1.0522`; `1.2` preserves the
covering invariant. Descendant caps are not the quantity bounded and may exceed
it.
"""
cap_inflation(::H3System) = 1.2

"""
    maxneighbors(::H3System, connectivity) -> Int

`6`, for either connectivity.

H3's cells are hexagons and pentagons, where sharing a vertex and sharing an
edge are the same relation — three cells meet at every vertex and any two of
them already share an edge — so [`Vertex()`](@ref Vertex) and [`Edge()`](@ref Edge)
coincide, and the bound is the hexagon's six. The twelve pentagons have five.
"""
maxneighbors(::H3System, ::Connectivity=Vertex()) = 6

DGG.winding(::H3System, ::Connectivity = Vertex()) = DGG.CounterClockwise()

"""
    maxring(::H3System, k, connectivity) -> Int

`6k`, as on any hexagonal system: tight at a hexagon, an over-bound at the
twelve pentagons, whose rings hold `5k`.
"""
DGG.maxring(::H3System, k::Integer, ::Connectivity = Vertex()) =
    (steps = Int(k); steps == 0 ? 1 : 6 * steps)

# ===========================================================================
# The dense order: indices <-> ids
# ===========================================================================

function ncells(::H3System, l::Integer)
    return @inbounds _H3_ROOT_ENDS[Int(l)+1][122]
end

"""
    cellindex(::H3System, l::Integer, i::Int) -> H3Cell

The id at index `i` of resolution `l`, computed from the base-cell prefix
sums and libh3's `childPosToCell`. The grid must bounds-check `i` first.
"""
function cellindex(::H3System, l::Integer, i::Int)
    r = Int(l)
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
    globalindex(::H3System, c::H3Cell) -> Union{Int,Nothing}

The index of `c` in its resolution's dense order, or `nothing` when `c` is
not a valid index. The grid must reject cells from another resolution first.

Malformed ids return `nothing`; libh3 child-position arithmetic is not itself a
validity check.
"""
function globalindex(::H3System, c::H3Cell)
    H3Native.is_valid_cell(c.id) || return nothing
    b = H3Native.get_base_cell(c.id)
    ends = @inbounds _H3_ROOT_ENDS[level(c)+1]
    previous = b == 0 ? 0 : @inbounds ends[b]
    return Int(previous + H3Native.cell_to_child_pos(c.id, 0) + 1)
end

# ===========================================================================
# Derived hierarchy: closed-form overrides of the generic walks
# ===========================================================================

"""
    ancestor(::H3System, c::H3Cell, l::Integer) -> H3Cell

The ancestor at resolution `l`, computed by one `cellToParent` call.
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

Every descendant at resolution `l`, in ascending id order. `cellToChildren`
handles any depth in one call.

O(subtree) and materialising, as the contract says — reach for
[`descendant_range`](@ref) instead wherever indices will do.
"""
function descendants(::H3System, c::H3Cell, l::Integer)
    target = Int(l)
    lc = level(c)
    target >= lc || throw(ArgumentError(
        "descendant level $target is above the cell's own level $lc"))
    target <= MAX_RESOLUTION || throw(ArgumentError(
        "descendant level $target is past maxlevel $MAX_RESOLUTION"))
    target == lc && return H3Cell[c]
    return [H3Cell(id) for id in H3Native.cell_to_children(c.id, target)]
end

"""
    descendant_range(::H3System, c::H3Cell, l::Integer) -> UnitRange{Int}

The contiguous interval of **indices** in `levelgrid(H3System(), l)` that the
descendants of `c` at resolution `l` occupy.

Computed without enumeration from child position zero and
`cellToChildrenSize`. Pentagon gaps are absent from H3 child positions, so the
range is contiguous and hole-free.
"""
function descendant_range(sys::H3System, c::H3Cell, l::Integer)
    target = Int(l)
    lc = level(c)
    target >= lc || throw(ArgumentError(
        "descendant level $target is above the cell's own level $lc"))
    target <= MAX_RESOLUTION || throw(ArgumentError(
        "descendant level $target is past maxlevel $MAX_RESOLUTION"))
    grid = levelgrid(sys, target)
    # The `l == level(c)` case is the cell's own one-element range: window
    # descent in `HierarchicalGridCursor` asks for it at the level above the
    # leaves and must not get an exception.
    if target == lc
        p = globalindex(grid, c)
        p === nothing && throw(ArgumentError("$c is not a valid H3 cell"))
        return p:p
    end
    first_child = H3Cell(H3Native.child_pos_to_cell(0, c.id, target))
    p = globalindex(grid, first_child)
    p === nothing && throw(ArgumentError("$c is not a valid H3 cell"))
    count = Int(H3Native.cell_to_children_size(c.id, target))
    return p:(p+count-1)
end
