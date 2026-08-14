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
# The position forms are the same answers read as indices into a data vector,
# and `halo_table` is those for a whole subset at once. On a rooted subtree it
# takes T20's split: an interior cell has no neighbour outside the subtree by
# definition, so its row needs no membership test at all.
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
# The position forms
#
# A bare `Int` is a position, so these are the same two verbs read as indices
# into a data vector laid out against the collection. ASCENDING rather than
# rotational: an index list is consumed by membership, and sorting it here is
# what lets a stencil table be built without a round trip through ids.
# ===========================================================================

# The generic route, for a grid with no faster membership test than its own
# `cellposition`. A complete level grid lands here too, where the clip is a
# no-op and the sort is the only work.
function _positions(sub, cells)
    out = Int[]
    for nb in cells
        p = cellposition(sub, nb)
        p === nothing || push!(out, p)
    end
    return sort!(out)
end

"""
    neighbors(grid::AbstractGrid, p::Int, k = 1; connectivity = Vertex()) -> Vector{Int}
    ring(grid::AbstractGrid, p::Int, k; connectivity = Vertex()) -> Vector{Int}

The neighbourhood of the cell at **position** `p`, as in-set positions in
ascending order. `p` outside `1:ncells(grid)` is a `BoundsError`, from
[`cellindex`](@ref).

See [`neighbors`](@ref) for the contract these keep, and [`halo_table`](@ref)
for the whole grid's worth at once.
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

The position forms on the two subset collections: positions in the subset,
ascending.

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

# ===========================================================================
# The halo table
# ===========================================================================

halo_table(grid::AbstractGrid, k::Integer = 1;
        connectivity::Connectivity = Vertex()) =
    [neighbors(grid, p, Int(k); connectivity) for p in 1:ncells(grid)]

halo_table(cv::CellVector, k::Integer = 1;
        connectivity::Connectivity = Vertex()) =
    [neighbors(cv, p, Int(k); connectivity) for p in 1:length(cv)]

# The rooted-subtree fast path. A subtree's INTERIOR is exactly the descendants
# with no neighbour outside it, so an interior cell's whole row is in-set by
# definition and its positions are the complete grid's minus the block offset —
# no membership test, no search. Only the rim, which is O(sqrt(subtree)) of the
# cells, has anything to decide.
#
# Three conditions, all load-bearing:
#
#   * the grid must hold the WHOLE subtree, not a subset of one, or "interior"
#     stops meaning "every neighbour is present";
#   * the iterators must be built with the SAME connectivity the halo is asked
#     for, since `Edge()`'s interior is strictly larger than `Vertex()`'s and
#     using the wrong one would admit a neighbour that is not there; and
#   * `k` must be 1. The T20 split proves a statement about the ONE-ring only:
#     a cell two steps inside the rim can still reach outside at `k == 2`, and
#     there is no "k-interior" iterator to ask.
function halo_table(pg::PartialGrid, k::Integer = 1;
        connectivity::Connectivity = Vertex())
    steps = Int(k)
    r = steps == 1 ? _whole_subtree_range(pg) : nothing
    r === nothing && return [neighbors(pg, p, steps; connectivity) for p in 1:ncells(pg)]
    return _rooted_halo(pg, r, connectivity)
end

# The subtree's position block, or `nothing` when this grid is not one. The
# count decides it: the ids are validated ascending and inside the root's range,
# so holding as many of them as the range is long means holding all of it.
function _whole_subtree_range(pg::PartialGrid)
    _is_rooted(pg) || return nothing
    has_sorted_subtrees(pg.system) || return nothing
    r = descendant_range(pg.system, pg.root_id, pg.level)
    return length(r) == ncells(pg) ? r : nothing
end

# Two passes rather than one merged walk. The two iterators are each ascending
# and together the descendants, so interleaving them by head comparison would
# make the row index a counter and save one `cellposition` per cell — but the
# two engines have DIFFERENT state types, so the interleaved loop is a
# union-typed `iterate` per step and measured 12x slower than the two walks
# separately. The lookup is the cheaper of the two costs.
#
# One loop per engine keeps each monomorphic, which is the whole reason they are
# written out twice instead of being parameterised over a membership predicate.
function _rooted_halo(pg::PartialGrid, r::UnitRange{Int}, connectivity::Connectivity)
    sys, root, complete = pg.system, pg.root_id, pg.complete
    l = pg.level
    out = Vector{Vector{Int}}(undef, ncells(pg))
    _inner_rows!(out, complete, InnerCellIterator(sys, root, l; connectivity),
        first(r), last(r), connectivity)
    _rim_rows!(out, complete, EdgeCellIterator(sys, root, l; connectivity),
        first(r), last(r), connectivity)
    return out
end

# `cellposition(complete, ·)` cannot be `nothing` for either the walked cell or
# its neighbours: both are cells OF the complete level, which is the grid being
# asked. The `::Int` says so where a `nothing` would otherwise be silent.
#
# The interior pass takes "every neighbour is inside the block" on the
# iterator's word, and that word is the ONLY thing keeping it from writing a
# wrong row: a position outside the block still lands somewhere in `out` after
# the offset. So it is checked rather than assumed — an interior cell with an
# outside neighbour means the split is broken, and this says so at the cell.
# Both passes check the ROW index the same way, for the same reason: what they
# are told to write is only as trustworthy as the walk that named it.
@noinline _escaped(d, nb) = error(
    "interior cell $d has neighbour $nb outside its own subtree: the " *
    "interior/rim split is wrong, and the halo table would be silently bad")

@noinline _outside(d) = error(
    "walked cell $d is outside the subtree it was walked from: the " *
    "interior/rim split is wrong, and the halo table would be silently bad")

function _inner_rows!(out, complete, inner, lo::Int, hi::Int,
        connectivity::Connectivity)
    for d in inner
        row = Int[]
        for nb in neighbors(complete, d, 1; connectivity)
            q = cellposition(complete, nb)::Int
            @boundscheck (lo <= q <= hi) || _escaped(d, nb)
            push!(row, q - lo + 1)
        end
        p = cellposition(complete, d)::Int
        @boundscheck (lo <= p <= hi) || _outside(d)
        out[p-lo+1] = sort!(row)
    end
    return out
end

function _rim_rows!(out, complete, rim, lo::Int, hi::Int, connectivity::Connectivity)
    for d in rim
        row = Int[]
        for nb in neighbors(complete, d, 1; connectivity)
            q = cellposition(complete, nb)::Int
            (q < lo || q > hi) && continue
            push!(row, q - lo + 1)
        end
        p = cellposition(complete, d)::Int
        @boundscheck (lo <= p <= hi) || _outside(d)
        out[p-lo+1] = sort!(row)
    end
    return out
end

# ===========================================================================
# The subset halo: the same boundary, seen from outside
# ===========================================================================

"""
    halo(sub; connectivity = Vertex())

The **out-of-set** cells adjacent to the subset: every cell of the complete
level that `sub` does not hold but that has a `connectivity`-neighbor it does,
each exactly once, in the complete level's ascending order.

On a rooted [`PartialGrid`](@ref) holding a WHOLE subtree the answer is
[`NeighborCellIterator`](@ref)`(system, root, level; connectivity)` itself —
lazy, and no membership machinery at all, because a subtree's outside ring is a
fact of the hierarchy. Everything else — a non-rooted `PartialGrid`, a
[`CellVector`](@ref), a rooted grid missing cells, any grid of A5, whose order
establishes no [`descendant_range`](@ref) to test wholeness against — takes
the eager fallback:
scan the members' one-rings, keep what [`cellposition`](@ref) cannot find,
deduplicate, sort. Built two ways over the same cells, the two paths agree
element for element.

They do return different types — an iterator and a `Vector` — so `collect` the
result where a vector is needed. See [`subtree_halo`](@ref) for the pure-cell
form and [`halo_table`](@ref) for the in-set stencil this is the complement of.
"""
function halo(pg::PartialGrid; connectivity::Connectivity = Vertex())
    r = _whole_subtree_range(pg)
    r === nothing || return NeighborCellIterator(pg.system, pg.root_id, pg.level;
        connectivity)
    return _clipped_out_ring(pg, pg.complete,
        (cellindex(pg, p) for p in 1:ncells(pg)), connectivity)
end

halo(cv::CellVector; connectivity::Connectivity = Vertex()) =
    _clipped_out_ring(cv, cv.grid, cv, connectivity)

# The fallback: `O(members · degree)` candidates with one membership test each —
# the same `cellposition === nothing` clip the tables above run, with the
# polarity reversed. Members arrive as an argument because the two containers
# spell iteration differently; the complete grid is where the one-rings and the
# final order both come from.
function _clipped_out_ring(sub, complete, members, connectivity::Connectivity)
    out = cellindextype(system(sub))[]
    for c in members
        for nb in neighbors(complete, c, 1; connectivity)
            cellposition(sub, nb) === nothing && push!(out, nb)
        end
    end
    return unique!(sort!(out))
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
in the result. The order is ascending `(level, position)`, deterministic and
independent of how the walk found them.

# The algorithm

A member's subtree is a region, and two regions of the hierarchy touch exactly
when two of their cells touch at a common depth. The set already names one:
its REFERENCE LEVEL `L`, the depth its covering guarantee is stated at and no
shallower than any member. So the question becomes a question about level `L`,
and it is answered without expanding anything:

 1. walk `c`'s subtree RIM at `L` — [`EdgeCellIterator`](@ref), `O(rim)` time
    and `O(depth)` memory, and the rim is where every contact with the outside
    lives by its own definition;
 2. take each rim cell's one-ring at `L`, under the requested connectivity;
 3. map each of those cells to the member that contains it, if any. Members are
    disjoint subtrees, so there is at most one, and on a sorted-subtree system
    it is found by binary search over the set's own curve keys: the keys are the
    members' reference-level range starts, ascending and disjoint, so
    `searchsortedlast` plus one range test decides it. A neighbour inside `c`'s
    own subtree maps back to `c` and is dropped.

`O(|rim(c, L)| · degree · log|set|)` time and `O(|answer|)` memory beyond the
rim walk's own `O(depth)`. The rim is the square root of the subtree, so a
member many levels above `L` costs the perimeter of its block and never its
area.

# Exactness, and where it is the hierarchy's answer rather than geometry's

Two members are neighbours here when their level-`L` cells are, which is
*geometric* boundary sharing exactly where the system's refinement is congruent
— HEALPix, S2 and ISEA4R, whose four children tile their parent, so a member's
footprint is the union of its level-`L` descendants'. Under `Edge()` the same
equivalence holds for shared edges, so vertex-only contact is excluded rather
than approximated away.

Where children do not tile their parent — IGEO7 and H3 (aperture 7) and A5 —
a member's footprint is NOT its descendants' union, and no level-`L` statement
can be a statement about the drawn polygons; see [`MultiOrderCoverage`](@ref)
for the size of that gap. The answer there is the hierarchy's, which is the
same relation [`subtree_border`](@ref) is defined by, and it is consistent with
every other subtree verb in this package. That is a carve-out about the
SYSTEMS, not about this walk.

!!! note "A5 pays for its missing primitives here too"
    Without [`has_sorted_subtrees`](@ref) there are no curve keys to binary
    search, so the member lookup is a set built per call, `O(|set|)`; and the
    rim iterator materialises the subtree rather than walking it, as T20
    documents. The answer is the same one.
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
    # `(level, key)` IS ascending `(level, position)`: within one level the keys
    # are that level's own order, whichever branch built them.
    sort!(hits; by = i -> (level(set.cells[i]), set.keys[i]))
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
