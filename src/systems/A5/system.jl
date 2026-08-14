# Canonical positions are quintant-major, then Hilbert-state-major:
# `quintant * 4^(level-1) + S + 1`. Level 0 uses origin order.

"""
    A5System() <: AbstractHierarchicalGridSystem

The [A5](https://a5geo.org) dodecahedral pentagonal grid at resolutions `0:29`.
The canonical id is [`A5Cell`](@ref); projection, hierarchy, geometry, and
adjacency follow upstream a5. A5 is equal-area on the *ellipsoid* and its
coordinates are geodetic, so unit-sphere [`cell_area`](@ref) carries the
authalic conversion and varies about 1% peak to peak within a level, still
summing to 4π.

The encoding also represents a level `-1` world cell and 42 of 60 level-30
quintants; neither belongs to a complete system level.
"""
struct A5System <: AbstractHierarchicalGridSystem end

# `levelgrid(A5System(), l)` returns the package's `HierarchicalLevelGrid`, a
# lightweight resolution descriptor. A5's fast paths dispatch on this alias;
# the five primitives it forwards to are the `(sys, ...)` methods below.
#
# Because A5 lacks sorted subtrees, treeifying a complete grid materializes all
# root positions. Use a `PartialGrid` for deep queries.
const LevelGrid = HierarchicalLevelGrid{A5System}

Base.show(io::IO, ::A5System) = print(io, "A5System()")

# `MAX_LEVEL` (29) is defined beside the encoding it is a fact about, in
# `cell.jl`.

# ===========================================================================
# Required system interface
# ===========================================================================

cellindextype(::A5System) = A5Cell

levels(::A5System) = 0:MAX_LEVEL

"""
    rootcells(::A5System)

The twelve dodecahedron faces, ascending. These are A5's res-0 tessellation, and
the only level whose cells are pentagons of the *base solid* rather than of the
lattice.
"""
rootcells(::A5System) = [A5Cell(id) for id in A5Native.res0_cells()]

"""
    parent(::A5System, c::A5Cell) -> A5Cell

The valid cell one resolution coarser, computed by id arithmetic. Throws
`ArgumentError` for a level-0 cell or invalid id.

The three regimes meet here: dropping to level 1 keeps the quintant and clears
the Hilbert state, dropping to level 0 divides the quintant by five, and every
step below level 1 drops two Hilbert bits.

Validation occurs before truncation so malformed padding cannot become a valid
parent id.
"""
function Base.parent(::A5System, c::A5Cell)
    l = level(c)
    l == 0 && throw(ArgumentError(
        "A5 cell $c is a root cell (resolution 0) and has no parent"))
    1 <= l <= MAX_LEVEL || throw(ArgumentError(
        "A5 cell $c is at resolution $l, outside levels(A5System()) = 0:$MAX_LEVEL"))
    isvalid(c) || throw(ArgumentError("A5 cell $c is not a valid cell"))
    return A5Cell(A5Native.cell_to_parent(c.id, l - 1))
end

"""
    children(::A5System, c::A5Cell) -> SmallVector{5,A5Cell}

The five quintants of a res-0 face, or the four Hilbert children of anything
deeper, ascending.

Returns a non-allocating fixed-capacity vector in ascending id order. Level-0
cells have five children; deeper cells have four.
"""
function children(::A5System, c::A5Cell)
    l = level(c)
    l < MAX_LEVEL || throw(ArgumentError(
        "A5 cell $c is at max_level $MAX_LEVEL and has no children"))
    cell = _decode(c.id)
    cell === nothing && throw(ArgumentError("A5 cell $c is not a valid cell"))
    out = SmallVector{5,A5Cell}()
    if l == 0
        # `serialize` keys on `segment_n = mod(segment - first_quintant, 5)`, so
        # walking `n` and mapping back to a segment is walking the ids in order.
        for n in 0:4
            segment = mod(n + cell.origin.first_quintant, 5)
            out = SmallCollections.push(out,
                A5Cell(A5Native.serialize(NativeCell(cell.origin, segment, UInt64(0), 1))))
        end
    else
        shifted = cell.S << 2
        for i in UInt64(0):UInt64(3)
            out = SmallCollections.push(out,
                A5Cell(A5Native.serialize(
                    NativeCell(cell.origin, cell.segment, shifted + i, l + 1))))
        end
    end
    return out
end

# ===========================================================================
# Traits
# ===========================================================================

"""
    has_sorted_subtrees(::A5System) -> Bool

`false`: no [`descendant_range`](@ref) contract is asserted across the level-0
to level-1 quintant fan-out.

  - [`treeify`](@ref) uses selection mode and materializes root positions;
    prefer a [`PartialGrid`](@ref) for deep grids — a complete one is O(cells)
    in memory and not viable past about level 12.
  - `MultiOrderCellSet` orders by `(level, position)` rather than by curve
    interval, and `level_ranges` on one raises an `ArgumentError`.
  - [`descendants`](@ref) is overridden to avoid level-by-level expansion.
"""
has_sorted_subtrees(::A5System) = false

"""
    cap_inflation(::A5System) -> Float64

`1.75`.

A5's four Hilbert children cover their parent's area but can extend beyond its
footprint. The measured descendant-to-cell-cap ratio reaches `1.45363`, with an
extrapolated bound of `1.47078`; `1.75` preserves the [`node_extent`](@ref)
covering invariant.
"""
cap_inflation(::A5System) = 1.75

"""
    max_neighbors(::A5System, connectivity) -> Int

`11` under [`Vertex()`](@ref Vertex) and `5` under [`Edge()`](@ref Edge).

A5 has corner-only neighbours, so `Vertex()` and `Edge()` differ.

| resolution | `Edge()` | `Vertex()` |
|---|---|---|
| 0 | 5 | 5 |
| 1 | 3 | 11 |
| ≥ 2 | 5 | 6, 7 or 8 |

The global bounds are therefore `5` and `11`; the latter occurs at level 1.
"""
max_neighbors(::A5System, ::Vertex) = 11
max_neighbors(::A5System, ::Edge) = 5

# ===========================================================================
# The dense order: positions <-> ids
# ===========================================================================

# `4^(level - 1)`, the number of cells per quintant, as an `Int64`. At level 29
# this is 7.2e16 and the whole level is 4.3e18 — inside `Int64` with a factor of
# two to spare, which is exactly why `levels` stops where it does.
_quintant_span(l::Int) = Int64(4)^(l - 1)

function ncells(::A5System, l::Integer)
    Int(l) == 0 && return 12
    return Int(60 * _quintant_span(Int(l)))
end

"""
    cellindex(::A5System, l::Integer, i::Int) -> A5Cell

The id at position `i` of resolution `l`: one `divrem` into
`(quintant, Hilbert state)` and one `serialize`. No table, no search, and O(1)
at every level. The grid has already bounds-checked `i`.
"""
function cellindex(::A5System, lvl::Integer, i::Int)
    l = Int(lvl)
    l == 0 && return A5Cell(@inbounds A5Native.res0_cells()[i])
    quintant, S = divrem(Int64(i) - 1, _quintant_span(l))
    origin = @inbounds A5Native.ORIGINS[Int(quintant ÷ 5)+1]
    segment = mod(Int(quintant % 5) + origin.first_quintant, 5)
    return A5Cell(A5Native.serialize(NativeCell(origin, segment, UInt64(S), l)))
end

"""
    cellposition(::A5System, c::A5Cell) -> Union{Int,Nothing}

The position of `c` in its own resolution's dense order, or `nothing` when `c`
is not a cell at all — the world cell, a res-30 id, or an id that is malformed
in any of the ways [`isvalid`](@ref) rejects. The grid has already rejected a
cell from another resolution.

Malformed ids, including ids with nonzero padding, return `nothing`.
"""
function cellposition(::A5System, c::A5Cell)
    l = level(c)
    cell = _decode(c.id)
    cell === nothing && return nothing
    quintant = Int64(c.id >> A5Native.HILBERT_START_BIT)
    l == 0 && return Int(quintant) + 1
    return Int(quintant * _quintant_span(l) + Int64(cell.S) + 1)
end

# ===========================================================================
# Derived hierarchy: closed-form overrides of the generic walks
# ===========================================================================

"""
    ancestor(::A5System, c::A5Cell, l::Integer) -> A5Cell

The ancestor at resolution `l`, in one `cell_to_parent` call rather than
`level(c) - l` of them: the bit arithmetic truncates to any coarser level
directly.

Validation precedes the identity case: invalid cells have no ancestors,
including themselves.
"""
function ancestor(::A5System, c::A5Cell, l::Integer)
    target = Int(l)
    lc = level(c)
    target <= lc || throw(ArgumentError(
        "ancestor level $target is deeper than the cell's own level $lc"))
    target >= 0 || throw(ArgumentError(
        "ancestor level $target is above the root level 0"))
    isvalid(c) || throw(ArgumentError("A5 cell $c is not a valid cell"))
    target == lc && return c
    return A5Cell(A5Native.cell_to_parent(c.id, target))
end

"""
    descendants(::A5System, c::A5Cell, l::Integer) -> Vector{A5Cell}

Every descendant at resolution `l`, ascending. `cell_to_children` spans any
number of levels in one call, so this never expands level by level.

The result is checked and sorted only if the level-0 segment rotation requires
it.

O(subtree) and materialising, as the contract says: level-`l` descendants of a
res-0 cell number `5·4^(l-1)`.
"""
function descendants(::A5System, c::A5Cell, l::Integer)
    target = Int(l)
    lc = level(c)
    target >= lc || throw(ArgumentError(
        "descendant level $target is above the cell's own level $lc"))
    target <= MAX_LEVEL || throw(ArgumentError(
        "descendant level $target is past max_level $MAX_LEVEL"))
    # Ahead of the identity case, so that the rule is the same one `ancestor`
    # states: an invalid cell has no relatives, itself included.
    isvalid(c) || throw(ArgumentError("A5 cell $c is not a valid cell"))
    target == lc && return A5Cell[c]
    ids = collect(UInt64, A5Native.cell_to_children(c.id, target))
    issorted(ids) || sort!(ids)
    return [A5Cell(id) for id in ids]
end
