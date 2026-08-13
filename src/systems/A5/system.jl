# ---------------------------------------------------------------------------
# `A5System` and its level grids
#
# The canonical dense order, which everything else here is a consequence of:
#
#     cells are numbered quintant by quintant (`5·origin + segment`, 0:59), and
#     within a quintant by the Hilbert state `S` (0:4^(level-1)-1), so
#
#         position = quintant · 4^(level-1) + S + 1
#
#     Level 0 is the twelve dodecahedron faces in origin order.
#
# That is exactly raw-`UInt64` order restricted to one level (see `cell.jl`), so
# `cellindex` of a level grid comes out sorted, and both directions of the
# position <-> id map are closed-form arithmetic with no table and no search.
#
# What it is NOT is a curve order across the hierarchy. The res-0 -> res-1
# fan-out cuts a face into five quintants that are five *separate* top-level
# keys, so a res-0 cell's descendants are five disjoint runs rather than one —
# hence `has_sorted_subtrees == false` and no `descendant_range`.
# ---------------------------------------------------------------------------

"""
    A5System() <: AbstractHierarchicalGridSystem

The [A5](https://a5geo.org) grid: an equal-area pentagonal tiling built on a
dodecahedron, over resolutions `0:29`.

Canonical id: [`A5Cell`](@ref). The projection, the codec, the hierarchy and the
adjacency are upstream a5's own (ported in [`A5Native`](@ref)), so this system
agrees with other a5 implementations cell for cell rather than approximating
one.

# Conventions this system documents

  - **Boundaries** are a5's own ring, **densified**: a cell's edges are straight
    in the dodecahedral face plane, not on the sphere, so each is split into
    `2^(6-level)` great-circle segments down to level 6 (one above it) before
    being inverse-projected. Counter-clockwise seen from outside, implicitly
    closed. A level-0 ring has 5 corners, a level-1 ring 3 (a quintant is a
    triangular slice of a face) and every deeper ring 5 — times that
    subdivision.
  - **Centroids** are a5's `cell_to_lonlat`: the face-plane centroid of the
    cell's polygon pulled back through the equal-area projection, which is the
    centre the lattice is built around rather than the mean of the ring.
  - **`cellat` ties** are a5's own `lonlat_to_cell`: project to the nearest
    face, estimate the lattice cell, and accept it when the face-plane polygon
    contains the point; otherwise widen over a fixed spiral of 24 samples and
    then over the estimates' neighbours, and fall back to the nearest miss.
    Deterministic per platform, and the same cell any other a5 implementation
    names for that point.
  - **Neighbour order** is rotational: [`neighbors`](@ref) is rings `1..k`
    concatenated outward, each ring counter-clockwise seen from outside, so
    [`ring`](@ref) is the tail block of [`neighbors`](@ref). Ring 1 **starts at
    the neighbour with the smallest [`A5Cell`](@ref) id**; see `neighbors.jl`.
  - **`Vertex()` and `Edge()` genuinely differ here** — A5's pentagons are not
    edge-to-edge — and both are implemented. See [`max_neighbors`](@ref).
  - **Areas** are the published ring's spherical area, and the ring *is* the
    cell. A5 is equal-area on the **ellipsoid**, and its coordinates are
    geodetic, so on the unit sphere a level's areas carry the authalic → geodetic
    latitude conversion and vary about 1% peak to peak; they still sum to 4π.

# The two levels that are not levels

`levels(A5System())` is `0:29`, not `0:30`. The encoding reaches resolution 30,
but only 42 of the 60 quintants still fit 64 bits there, so there is no complete
res-30 grid — and adjacency at res 30 is worse than incomplete, since a
neighbour in an unsupported quintant comes back as a *res-29* id. Below the
range, `A5Cell(0)` is upstream's whole-sphere "resolution -1" world cell.
Neither is a cell of any grid; see [`A5Cell`](@ref) and `isvalid`.
"""
struct A5System <: AbstractHierarchicalGridSystem end

"""
    A5Grid(level) <: AbstractGrid

The complete A5 grid at one resolution — build one with
[`levelgrid(A5System(), l)`](@ref levelgrid) rather than by calling this.

A lightweight descriptor: it stores the resolution and nothing else, so
constructing the res-29 grid (4,323,455,642,275,676,160 cells) is free.

!!! warning "Do not `treeify` a deep complete grid"
    A5 has no [`has_sorted_subtrees`](@ref), so [`treeify`](@ref) builds a
    *selection-mode* cursor whose root materialises `1:ncells(grid)`. That is
    fine at the shallow levels and impossible past about level 12. Everything
    that treeifies — the generic [`cellat`](@ref), [`query`](@ref) on a grid —
    inherits the limit; `cellat` is overridden here precisely so it does not,
    and a deep query should be run as `query(sys, pred; level)` over a
    [`PartialGrid`](@ref) of the region instead.
"""
struct A5Grid <: AbstractGrid
    level::Int
end

system(::A5Grid) = A5System()
level(g::A5Grid) = g.level

Base.show(io::IO, g::A5Grid) = print(io, "A5Grid(res ", g.level, ")")
Base.show(io::IO, ::A5System) = print(io, "A5System()")

# `MAX_LEVEL` (29) is defined beside the encoding it is a fact about, in
# `cell.jl`.

# ===========================================================================
# Required system interface
# ===========================================================================

cellindextype(::A5System) = A5Cell

levels(::A5System) = 0:MAX_LEVEL

function levelgrid(::A5System, l::Integer)
    lvl = Int(l)
    0 <= lvl <= MAX_LEVEL || throw(ArgumentError(
        "A5 resolution $lvl is outside levels(A5System()) = 0:$MAX_LEVEL"))
    return A5Grid(lvl)
end

"""
    rootcells(::A5System)

The twelve dodecahedron faces, ascending. These are A5's res-0 tessellation, and
the only level whose cells are pentagons of the *base solid* rather than of the
lattice.
"""
rootcells(::A5System) = [A5Cell(id) for id in A5Native.res0_cells()]

"""
    parent(::A5System, c::A5Cell) -> A5Cell

The cell one resolution coarser, by `cell_to_parent` — bit arithmetic on the id,
no lookup. Throws an `ArgumentError` on a res-0 cell, which has no parent.

The three regimes meet here: dropping to level 1 keeps the quintant and clears
the Hilbert state, dropping to level 0 divides the quintant by five, and every
step below level 1 drops two Hilbert bits.
"""
function Base.parent(::A5System, c::A5Cell)
    l = level(c)
    l == 0 && throw(ArgumentError(
        "A5 cell $c is a root cell (resolution 0) and has no parent"))
    1 <= l <= MAX_LEVEL || throw(ArgumentError(
        "A5 cell $c is at resolution $l, outside levels(A5System()) = 0:$MAX_LEVEL"))
    return A5Cell(A5Native.cell_to_parent(c.id, l - 1))
end

"""
    children(::A5System, c::A5Cell) -> SmallVector{5,A5Cell}

The five quintants of a res-0 face, or the four Hilbert children of anything
deeper, ascending.

The container is fixed-capacity and the call **does not allocate**, which
matters because tree descent calls this once per node. Both branches emit
ascending ids *by construction* rather than by sorting: at level 0 the loop runs
over the quintant number `5·origin + n` and converts to a5's segment numbering
(a rotation of it) inside the loop, and below level 0 the four children are
`S<<2 .+ (0:3)` in one quintant, which is already the order their raw ids take.

Generic code must never assume a fixed child count — five here, four there — and
reading `length` is all it takes not to.
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

`false` — the conservative default, kept **deliberately**, and the reason A5 is
the one system in this package whose grids take the generic slow paths.

The arithmetic suggests the trait could hold. A cell at level `m ≥ 1` in
quintant `q` with Hilbert state `S` owns level-`l` states
`S·4^(l-m) .. (S+1)·4^(l-m) - 1`, which by the position formula at the top of
this file is one contiguous run; and a res-0 face owns quintants `5·origin`
through `5·origin + 4`, which is another. Nothing here has *verified* the
two-sided contract of [`descendant_range`](@ref) — every descendant in the
range **and** every position in the range a descendant — at the res-0 → res-1
regime change, where `serialize` keys the quintant on
`mod(segment - first_quintant, 5)`, a rotation of a5's own segment walk. A trait
declared falsely produces silently wrong subtree answers with nothing
downstream able to detect them, so it stays undeclared until that work is done,
and no [`descendant_range`](@ref) method exists.

What that costs, concretely, and what to do about it:

  - [`treeify`](@ref) returns a **selection-mode** `HierarchicalGridCursor`.
    Every node materialises the positions it owns, and the root materialises
    `1:ncells(grid)` — so treeifying a *complete* level grid is O(cells) in
    memory and is not viable past about level 12. Treeify a
    [`PartialGrid`](@ref) of the region instead; that is what the mode is for.
  - `MultiOrderCellSet` orders by `(level, position)` rather than by curve
    interval, and `level_ranges` on one raises an `ArgumentError`.
  - The generic [`descendants`](@ref) would expand level by level; this system
    overrides it, so that cost is not paid.

A5 is the first real system on both of those paths, which is why
`test/systems/A5/` drives treeify, query and multi-order coverage through them
rather than leaving the substrate's mocks as their only coverage.
"""
has_sorted_subtrees(::A5System) = false

"""
    cap_inflation(::A5System) -> Float64

`1.75`, well above the package default of `1.2`, and the number is measured.

The cause is structural. An A5 cell is a fixed pentagon placed at a lattice
point and scaled by `2^-resolution`, and that pentagon is **not a rep-4 tile**:
the four cells a Hilbert digit names as children cover their parent's *area* but
not its *footprint*, so a descendant's vertices reach well outside the parent's
own bounding cap. Sampling random points, about 37% of them land in a res-`r`
cell whose parent is not the res-`r-1` cell containing them (0% for 0 → 1, where
the quintant cut is exact).

The CAP-VALIDATION sweep in `test/systems/A5/` measures the union ratio — the
farthest descendant vertex from a cell's cap centre, over that cell's own
radius before inflation — exhaustively at res 0-2 and on ordinal samples at res
3, 5 and 8, out to depth 8, and extrapolates the geometric tail. Worst measured
**1.45363** (res 8, depth 6); the increments halve cleanly from depth 4 on, so
the tail puts the supremum at **1.47078**. H3 measures 1.052 and IGeo7 1.048 by
comparison. `1.75` clears the extrapolated figure by 19%, at the price of caps
about 2.4× the default's area, and no descendant in the sweep reaches beyond
83% of the inflated radius.

Setting it lower is a correctness bug — see the covering law in
[`node_extent`](@ref). The suite also walks a chain of children from every
res-0 cell all the way to `max_level` and checks containment directly, which is
the property that actually matters; the worst point on those twelve chains sits
0.46 rad *inside* the root's extent.
"""
cap_inflation(::A5System) = 1.75

"""
    max_neighbors(::A5System, connectivity) -> Int

`11` under [`Vertex()`](@ref Vertex) and `5` under [`Edge()`](@ref Edge).

**The two connectivities genuinely differ on A5**, which makes it the first
non-quadrilateral system in this package where they do. The interface's
`Vertex()` docstring says hexagonal and pentagonal grids coincide; that is true
of the *icosahedral* hex-plus-12-pentagon systems, where exactly three cells
meet at every vertex, and it is false here. A5's pentagons tile in the manner of
a Cairo tiling: four cells meet at some corners, so a cell has neighbours it
shares one corner with and no edge. Measured over complete levels 1-3 against
the ring geometry itself, every cell in a5's `edge_only = false` set shares
**exactly one** corner with the subject where an `edge_only = true` neighbour
shares **two**, with no cell in one set and not the other — so upstream's two
modes are exactly Moore and von Neumann, and both are wired.

The degrees are constant per regime, not merely bounded, which is what makes
them safe to extrapolate past the sampled levels:

| resolution | `Edge()` | `Vertex()` |
|---|---|---|
| 0 | 5 | 5 |
| 1 | 3 | 11 |
| ≥ 2 | 5 | 6, 7 or 8 |

Res 0 is the one level where the two coincide: a dodecahedron face shares a
vertex with no face it does not also share an edge with. Res 1 is the outlier
in the other direction — a quintant is a *triangle* with three edges, while the
ten cells it touches at a corner include the four other quintants of its own
face (they all meet at the face centre) and six across two dodecahedron
vertices. So `11` is the bound, attained only at level 1, and `8` is the most
any cell of a level below it reports.

Sweeping complete levels 0-4 (5,112 cells) and 200-cell ordinal samples at res
5, 9, 15, 22 and 29 found no cell outside those figures, and the neighbour
relation symmetric in both directions and closed within its level throughout.
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

function ncells(grid::A5Grid)
    grid.level == 0 && return 12
    return Int(60 * _quintant_span(grid.level))
end

"""
    cellindex(grid::A5Grid, i::Int) -> A5Cell

The id at position `i`: one `divrem` into `(quintant, Hilbert state)` and one
`serialize`. No table, no search, and O(1) at every level.
"""
function cellindex(grid::A5Grid, i::Int)
    n = ncells(grid)
    1 <= i <= n || throw(BoundsError(grid, i))
    l = grid.level
    l == 0 && return A5Cell(@inbounds A5Native.res0_cells()[i])
    quintant, S = divrem(Int64(i) - 1, _quintant_span(l))
    origin = @inbounds A5Native.ORIGINS[Int(quintant ÷ 5)+1]
    segment = mod(Int(quintant % 5) + origin.first_quintant, 5)
    return A5Cell(A5Native.serialize(NativeCell(origin, segment, UInt64(S), l)))
end

"""
    cellposition(grid::A5Grid, c::A5Cell) -> Union{Int,Nothing}

The position of `c`, or `nothing` when `c` is not a cell of this grid — a
different resolution, the world cell, a res-30 id, or an id that is malformed in
any of the ways [`isvalid`](@ref) rejects.

`nothing` rather than an error is the contract, and it is what makes asking "is
this cell here?" the normal way to intersect an id set with a grid. The validity
check is not paranoia: the a5 arithmetic will happily decode an index with junk
in its padding bits into a neighbouring cell's `(quintant, S)`, and returning
that cell's position confidently is worse than returning nothing.
"""
function cellposition(grid::A5Grid, c::A5Cell)
    l = grid.level
    level(c) == l || return nothing
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
"""
function ancestor(::A5System, c::A5Cell, l::Integer)
    target = Int(l)
    lc = level(c)
    target <= lc || throw(ArgumentError(
        "ancestor level $target is deeper than the cell's own level $lc"))
    target >= 0 || throw(ArgumentError(
        "ancestor level $target is above the root level 0"))
    target == lc && return c
    return A5Cell(A5Native.cell_to_parent(c.id, target))
end

"""
    descendants(::A5System, c::A5Cell, l::Integer) -> Vector{A5Cell}

Every descendant at resolution `l`, ascending. `cell_to_children` spans any
number of levels in one call, so this never expands level by level.

The result is sorted rather than trusted: below level 0 a5 walks the Hilbert
state, which *is* id order, but the res-0 fan-out walks a5's segments while the
ids key on `mod(segment - first_quintant, 5)` — a rotation of that walk. One
`issorted` pass covers that case for O(n) rather than the O(n log n) of sorting
unconditionally.

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
    target == lc && return A5Cell[c]
    isvalid(c) || throw(ArgumentError("A5 cell $c is not a valid cell"))
    ids = collect(UInt64, A5Native.cell_to_children(c.id, target))
    issorted(ids) || sort!(ids)
    return [A5Cell(id) for id in ids]
end
