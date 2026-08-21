# ---------------------------------------------------------------------------
# Stencils on subsets.
#
# `neighbors` / `ring` on a SUBSET of a complete level are the complete level's
# answer clipped to membership — `ring(sub, c, k) == filter(in(sub), ring(complete,
# c, k))`, stated in `src/interface/grid.jl`. Two consequences shape this file:
#
#   * the per-system automata are reached through the complete level grid, so a
#     subset gets the fast path instead of the geometric tree walk; and
#   * distance is the SYSTEM's, so nothing here does a breadth-first walk of its
#     own and a hole in the subset cannot lengthen a path around itself.
#
# The position forms are the same answers read as indices into a data vector.
# `adjacency` is those for a whole region at once.
# ---------------------------------------------------------------------------

# ===========================================================================
# Membership, and the error for a cell that has none
# ===========================================================================

@noinline function _not_a_member(sub, c)
    throw(ArgumentError(
        "$c is not a cell of $sub. On a subset `neighbors` and `ring` are the " *
        "complete level's answer clipped to membership, so a cell the subset " *
        "does not hold has no neighbourhood here — asking the complete level " *
        "would answer about cells this collection does not contain"))
end

@inline function _member_or_throw(sub, c::AbstractCellIndex)
    p = cellposition(sub, c)
    p === nothing && _not_a_member(sub, c)
    return p
end

# ===========================================================================
# The clipped ids
# ===========================================================================

# One shape for both verbs and both containers: run the complete level's answer
# and keep what the subset holds, in the order it came in. That order IS the
# rotational contract, so clipping preserves it for free.
#
# The membership check is the CALLER's, not this function's: an argument is
# evaluated before the call, so checking here would compute a whole `k == 3`
# disc for a cell the subset does not hold and throw afterwards.
# The clip keeps the container it was handed. A clipped ring is a subset of the
# unclipped one, so `k == 1` stays in the `SmallVector` the system answered in
# rather than moving to the heap for the membership test; `k >= 2` has no static
# bound, arrives in a `Vector`, and leaves in one.
function _clip(sub, cells::SmallCollections.SmallVector{N,T}) where {N,T}
    out = SmallCollections.SmallVector{N,T}()
    for nb in cells
        cellposition(sub, nb) === nothing ||
            (out = SmallCollections.push(out, nb))
    end
    return out
end

function _clip(sub, cells)
    out = Vector{eltype(cells)}()
    for nb in cells
        cellposition(sub, nb) === nothing || push!(out, nb)
    end
    return out
end

"""
    neighbors(pg::PartialGrid, c, k = 1; connectivity = Vertex())
    ring(pg::PartialGrid, c, k; connectivity = Vertex())

The complete level's neighbourhood of `c`, clipped to the subset — which is
what [`neighbors`](@ref) means on a subset. The complete grid is
`levelgrid(system(pg), level(pg))`, so this reaches the system's own automaton
and never the geometric tree walk; the clip is one `cellposition` per candidate.

`c` outside `pg` throws an `ArgumentError` rather than answering about cells
`pg` does not hold.
"""
function neighbors(pg::PartialGrid, c::AbstractCellIndex, k::Integer = 1;
        connectivity::Connectivity = Vertex())
    _member_or_throw(pg, c)
    return _clip(pg, neighbors(pg.complete, c, Int(k); connectivity))
end

function ring(pg::PartialGrid, c::AbstractCellIndex, k::Integer;
        connectivity::Connectivity = Vertex())
    _member_or_throw(pg, c)
    return _clip(pg, ring(pg.complete, c, Int(k); connectivity))
end

"""
    neighbors(cv::CellVector, c, k = 1; connectivity = Vertex())
    ring(cv::CellVector, c, k; connectivity = Vertex())

[`neighbors`](@ref) and [`ring`](@ref) on the compressed collection, with the
same clipped-to-membership meaning a [`PartialGrid`](@ref) has. Membership is
the window search — `O(log #windows)` — so the clip costs nothing the vector
was not already able to answer.
"""
function neighbors(cv::CellVector, c::AbstractCellIndex, k::Integer = 1;
        connectivity::Connectivity = Vertex())
    _member_or_throw(cv, c)
    return _clip(cv, neighbors(cv.grid, c, Int(k); connectivity))
end

function ring(cv::CellVector, c::AbstractCellIndex, k::Integer;
        connectivity::Connectivity = Vertex())
    _member_or_throw(cv, c)
    return _clip(cv, ring(cv.grid, c, Int(k); connectivity))
end

# ===========================================================================
# Neighbour counts default to the length of the clipped one-ring. Systems with
# structural degree information may specialize this method.
# ===========================================================================

neighborcount(grid::AbstractGrid, c::AbstractCellIndex;
        connectivity::Connectivity = Vertex()) =
    length(neighbors(grid, c, 1; connectivity))

neighborcount(cv::CellVector, c::AbstractCellIndex;
        connectivity::Connectivity = Vertex()) =
    length(neighbors(cv, c, 1; connectivity))

# ===========================================================================
# The position forms
#
# A bare `Int` is a position, so these are the same two verbs read as indices
# into a data vector laid out against the collection — the id answer mapped
# through `cellposition`, element for element. The ROTATIONAL order is carried
# across: slot `i` of a row names the same direction as slot `i` of the ring,
# which is the whole content of the order contract and the thing an oriented
# stencil reads.
# ===========================================================================

# The generic route, for a grid with no faster membership test than its own
# `cellposition`. A complete level grid lands here too, where the clip is a
# no-op.
function _positions(sub, cells)
    out = Int[]
    for nb in cells
        p = cellposition(sub, nb)
        p === nothing || push!(out, p)
    end
    return out
end

"""
    neighbors(grid::AbstractGrid, p::Int, k = 1; connectivity = Vertex()) -> Vector{Int}
    ring(grid::AbstractGrid, p::Int, k; connectivity = Vertex()) -> Vector{Int}

The neighbourhood of the cell at **position** `p`, as in-set positions in the
same counter-clockwise order the id form answers in. `p` outside
`1:ncells(grid)` is a `BoundsError`, from [`cellindex`](@ref).

See [`neighbors`](@ref) for the contract these keep, and [`adjacency`](@ref)
for a whole region's worth at once.
"""
neighbors(grid::AbstractGrid, p::Int, k::Integer = 1;
        connectivity::Connectivity = Vertex()) =
    _positions(grid, neighbors(grid, cellindex(grid, p), Int(k); connectivity))

ring(grid::AbstractGrid, p::Int, k::Integer;
        connectivity::Connectivity = Vertex()) =
    _positions(grid, ring(grid, cellindex(grid, p), Int(k); connectivity))

"""
    neighbors(sub::PartialGrid, p::Int, k = 1; connectivity = Vertex()) -> Vector{Int}
    ring(sub::PartialGrid, p::Int, k; connectivity = Vertex()) -> Vector{Int}
    neighbors(cv::CellVector, p::Int, k = 1; connectivity = Vertex()) -> Vector{Int}
    ring(cv::CellVector, p::Int, k; connectivity = Vertex()) -> Vector{Int}

The position forms on the two subset collections: positions in the subset, in
the ring order [`neighbors`](@ref) states.

The complete level's candidates are clipped ONCE here rather than by going
through the id form and testing membership again on the way back — a position
already names a cell of the subset, so there is nothing for the id form's
out-of-set check to decide.
"""
neighbors(pg::PartialGrid, p::Int, k::Integer = 1;
        connectivity::Connectivity = Vertex()) =
    _positions(pg, neighbors(pg.complete, cellindex(pg, p), Int(k); connectivity))

ring(pg::PartialGrid, p::Int, k::Integer;
        connectivity::Connectivity = Vertex()) =
    _positions(pg, ring(pg.complete, cellindex(pg, p), Int(k); connectivity))

neighbors(cv::CellVector, p::Int, k::Integer = 1;
        connectivity::Connectivity = Vertex()) =
    _positions(cv, neighbors(cv.grid, cv[p], Int(k); connectivity))

ring(cv::CellVector, p::Int, k::Integer;
        connectivity::Connectivity = Vertex()) =
    _positions(cv, ring(cv.grid, cv[p], Int(k); connectivity))

# The subtree's position block, or `nothing` when this grid is not one. The
# count decides it: the ids are validated ascending and inside the root's range,
# so holding as many of them as the range is long means holding all of it.
#
# This is the one predicate `halo`, `border`, `interior` and `adjacency` all
# branch on, read from one function so that "is a subtree" cannot come to mean
# two things.
function _whole_subtree_range(pg::PartialGrid)
    _is_rooted(pg) || return nothing
    has_sorted_subtrees(pg.system) || return nothing
    r = descendant_range(pg.system, pg.root_id, pg.level)
    return length(r) == ncells(pg) ? r : nothing
end

# ===========================================================================
# Cross-level adjacency on a multi-order set
# ===========================================================================

"""
    member_neighbors(set::MultiOrderCellSet, c; connectivity = Vertex()) -> Vector

The members of `set` that share a boundary with the member `c` — its neighbours
in a set whose cells sit at different levels, so a neighbour may be COARSER than
`c` (an ancestor of one of its system-neighbours) or FINER (a cell deep inside
one), and this is the one call that finds both.

`c` must be a member; anything else is an `ArgumentError`. `c` itself is never
in the result. The order is the package's one order — counter-clockwise about
`cell_centroid(system(set), c)` seen from outside, as [`neighbors`](@ref)
states — measured on each neighbour's own centroid, whatever level it sits at.
The start is the smallest-`(level, position)` neighbour, and exact azimuth ties
break the same way, so the answer is deterministic and independent of how the
walk found them.

# The algorithm

A member's subtree is a region, and two regions of the hierarchy touch exactly
when two of their cells touch at a common depth. The set already names one:
its REFERENCE LEVEL `L`, the depth its covering guarantee is stated at and no
shallower than any member. So the question becomes a question about level `L`,
and it is answered without expanding anything:

 1. walk `c`'s subtree BORDER at `L` — [`EdgeCellIterator`](@ref), `O(border)` time
    and `O(depth)` memory, and the border is where every contact with the outside
    lives by its own definition;
 2. take each border cell's one-ring at `L`, under the requested connectivity;
 3. map each of those cells to the member that contains it, if any. Members are
    disjoint subtrees, so there is at most one, and on a sorted-subtree system
    it is found by binary search over the set's own curve keys: the keys are the
    members' reference-level range starts, ascending and disjoint, so
    `searchsortedlast` plus one range test decides it. A neighbour inside `c`'s
    own subtree maps back to `c` and is dropped.

`O(|border(c, L)| · degree · log|set|)` time and `O(|answer|)` memory beyond the
border walk's own `O(depth)`. The border is the square root of the subtree, so a
member many levels above `L` costs the perimeter of its block and never its
area.

# Exactness, and where it is the hierarchy's answer rather than geometry's

Two members are neighbours here when their level-`L` cells are, which is
*geometric* boundary sharing exactly where the system's refinement is congruent
— HEALPix, S2, ISEA4R, ISEA4T, rHEALPix/AusPIX and the IVEA/RTEA rhombi, whose
children tile their parent, so a member's footprint is the union of its level-`L`
descendants'. Under `Edge()` the same equivalence holds for shared edges, so
vertex-only contact is excluded rather than approximated away.

Where children do not tile their parent — IGEO7 and H3 (aperture 7), A5, and the
central-place ISEA3H and ISEA4H — a member's footprint is NOT its descendants'
union, and no level-`L` statement can be a statement about the drawn polygons;
see [`MultiOrderCoverage`](@ref) for the size of that gap. The answer there is
the hierarchy's, which is the same relation [`border`](@ref) is defined by, and
it is consistent with every other subtree verb in this package. That is a
carve-out about the SYSTEMS, not about this walk.

!!! note "A5 and the IVEA/RTEA rhombi pay for their missing primitives here too"
    Without [`has_sorted_subtrees`](@ref) there are no curve keys to binary
    search, so the member lookup is a set built per call, `O(|set|)`; and the
    border iterator materialises the subtree rather than walking it. The answer
    is unchanged.
"""
function member_neighbors(set::MultiOrderCellSet, c::AbstractCellIndex;
        connectivity::Connectivity = Vertex())
    sys = set.system
    L = set.reference_level
    grid = levelgrid(sys, L)
    lookup = _member_lookup(set)
    home = _home_index(lookup, set, L, c)
    home === nothing && throw(ArgumentError(
        "$c is not a member of this MultiOrderCellSet, so it has no neighbours " *
        "in it; the verb is about members, not about cells of a level"))
    hits = Int[]
    seen = Set{Int}((home,))
    for d in EdgeCellIterator(sys, c, L; connectivity)
        for nb in neighbors(grid, d, 1; connectivity)
            i = _member_of(lookup, set, grid, nb)
            (i === nothing || i in seen) && continue
            push!(seen, i)
            push!(hits, i)
        end
    end
    length(hits) <= 1 && return [set.cells[i] for i in hits]
    # `(level, key)` IS ascending `(level, position)`: within one level the keys
    # are that level's own order, whichever branch built them. It is the tie
    # break and the start anchor, not the order — the order is the ring's.
    rank(i) = (level(set.cells[i]), set.keys[i])
    sort!(hits; by = rank)
    # Each member is read on its own level's grid, because a coarse neighbour
    # has no centroid at `L`.
    centroid(x) = cell_centroid(levelgrid(sys, level(x)), x)
    centre = centroid(c)
    anchor = centroid(set.cells[first(hits)])
    e1, e2 = _tangent_basis(centre, anchor)
    spoke = _azimuth(centre, e1, e2, anchor)
    sort!(hits; by = i -> (_phase(_azimuth(centre, e1, e2,
            centroid(set.cells[i])) - spoke), rank(i)))
    return [set.cells[i] for i in hits]
end

# The curve keys are the members' reference-level range starts — ascending and
# disjoint — so no index has to be built at all.
struct CurveKeyLookup end

# A5 has no ranges to search, so the members go into a dictionary and the walk
# up from a level-`L` cell asks it once per generation.
struct AncestorLookup{ID}
    index::Dict{ID,Int}
end

_member_lookup(set::MultiOrderCellSet) =
    has_sorted_subtrees(set.system) ? CurveKeyLookup() :
    AncestorLookup(Dict(c => i for (i, c) in enumerate(set.cells)))

# Is `c` a member, and which one — the same two structures answering about a
# cell at its OWN level rather than at the reference level. Binary search again
# rather than a scan, so the membership check is not the expensive step.
function _home_index(::CurveKeyLookup, set::MultiOrderCellSet, L::Int,
        c::AbstractCellIndex)
    level(c) <= L || return nothing
    i = searchsortedfirst(set.keys, first(descendant_range(set.system, c, L)))
    (i <= length(set.cells) && set.cells[i] == c) || return nothing
    return i
end

_home_index(lk::AncestorLookup, ::MultiOrderCellSet, ::Int, c::AbstractCellIndex) =
    get(lk.index, c, nothing)

function _member_of(::CurveKeyLookup, set::MultiOrderCellSet, grid::AbstractGrid,
        nb::AbstractCellIndex)
    p = cellposition(grid, nb)
    p === nothing && return nothing
    i = searchsortedlast(set.keys, p)
    i >= 1 || return nothing
    return p <= last(descendant_range(set.system, set.cells[i], level(grid))) ?
           i : nothing
end

function _member_of(lk::AncestorLookup, set::MultiOrderCellSet, grid::AbstractGrid,
        nb::AbstractCellIndex)
    top = first(levels(set.system))
    a = nb
    while true
        i = get(lk.index, a, nothing)
        i === nothing || return i
        level(a) <= top && return nothing
        a = Base.parent(set.system, a)
    end
end
