# The engines behind `halo`: lazy walks over the cells outside a subtree or a
# subset that touch it. Systems specialize `halo_engine`; conservative candidate
# bands must be filtered by the requested adjacency before yielding.

"""
    SubtreeHaloIterator(sys, c, l; connectivity = Vertex())

The halo of `c`'s subtree at level `l`, lazily: every level-`l` cell that is
**not** a descendant of `c` but has a neighbour that is, in ascending canonical
order, each cell exactly once.

What [`halo`](@ref) returns for a region that is a whole rooted subtree, with
`cells = true`. `l == level(c)` is `c`'s own one-ring, sorted. `l < level(c)`
and `l > maxlevel(sys)` throw an `ArgumentError`.

`cellposition(levelgrid(sys, l), x)` is strictly increasing over the walk, which
is what `halo`'s position form reads. This differs from the rotational ordering
of [`neighbors`](@ref).

Construction does not materialize the halo. The iterator holds `O(depth)` walk
state and bounded neighbour containers.

[`Base.IteratorSize`](@ref) is `HasLength()` only when an engine derives an
exact count; otherwise it is `SizeUnknown()`, and
[`sizehint`](@ref DiscreteGlobalGrids.sizehint) is the inexact estimate.
"""
struct SubtreeHaloIterator{S<:AbstractHierarchicalGridSystem,C<:AbstractCellIndex,
        K<:Connectivity,E}
    system::S
    cell::C
    level::Int
    connectivity::K
    engine::E
end

function SubtreeHaloIterator(sys::AbstractHierarchicalGridSystem,
        c::AbstractCellIndex, l::Integer; connectivity::Connectivity = Vertex())
    target = Int(l)
    return SubtreeHaloIterator(sys, c, target, connectivity,
        halo_engine(sys, c, target, connectivity))
end

Base.iterate(it::SubtreeHaloIterator) = iterate(it.engine)
Base.iterate(it::SubtreeHaloIterator, state) = iterate(it.engine, state)

# Derive `eltype` from the input cell type so engine dispatch unions do not
# widen collected results to `Vector{Any}`.
Base.eltype(::Type{<:SubtreeHaloIterator{S,C}}) where {S,C} = C
Base.IteratorSize(::Type{<:SubtreeHaloIterator{S,C,K,E}}) where {S,C,K,E} =
    Base.IteratorSize(E)

# Engines without a constant-time count intentionally provide no `length`.
Base.length(it::SubtreeHaloIterator) = length(it.engine)

Base.show(io::IO, it::SubtreeHaloIterator) = print(io, "SubtreeHaloIterator(",
    it.system, ", ", it.cell, ", ", it.level, "; connectivity = ",
    it.connectivity, ")")

# ===========================================================================
# The one-ring engine: `l == level(c)`
# ===========================================================================

# Emit the rotationally ordered native ring in ascending order without copying
# its fixed-capacity container.
"""
    RingHaloEngine(ring)

`c`'s own one-ring, ascending, by selection emit. `O(degree^2)` time with
`degree <= maxneighbors(sys, connectivity)`, no allocation, isbits state.

`length` equals `length(ring)`, requiring native one-rings to contain no
duplicates. `collect_subtree` reports a count mismatch if this invariant fails.
"""
struct RingHaloEngine{V,C}
    ring::V
end

RingHaloEngine(ring) = RingHaloEngine{typeof(ring),eltype(ring)}(ring)

Base.eltype(::Type{<:RingHaloEngine{V,C}}) where {V,C} = C
Base.IteratorSize(::Type{<:RingHaloEngine}) = Base.HasLength()
Base.length(e::RingHaloEngine) = length(e.ring)

struct RingHaloState{C}
    emitted::Int
    last::C
end

function Base.iterate(e::RingHaloEngine{V,C}) where {V,C}
    isempty(e.ring) && return nothing
    best = first(e.ring)
    for x in e.ring
        x < best && (best = x)
    end
    return (best, RingHaloState(1, best))
end

function Base.iterate(e::RingHaloEngine{V,C}, s::RingHaloState{C}) where {V,C}
    s.emitted >= length(e.ring) && return nothing
    found = false
    best = s.last
    for x in e.ring
        x <= s.last && continue
        if !found || x < best
            best = x
            found = true
        end
    end
    found || return nothing
    return (best, RingHaloState(s.emitted + 1, best))
end

# ===========================================================================
# Generic engine construction
# ===========================================================================

# Shared level validation keeps generic and specialized engines consistent.
function check_halo_level(sys::AbstractHierarchicalGridSystem,
        c::AbstractCellIndex, target::Int)
    lc = level(c)
    target >= lc || throw(ArgumentError(
        "halo: level $target is above the cell's own level $lc"))
    target <= maxlevel(sys) || throw(ArgumentError(
        "halo: level $target is past maxlevel $(maxlevel(sys))"))
    return nothing
end

"""
    generic_halo_engine(sys, c, target, connectivity)

Return the generic halo engine after level validation: [`RingHaloEngine`](@ref)
at depth zero, [`OutsideWalkEngine`](@ref) for sorted subtrees, or
[`ScanHaloEngine`](@ref) when no descendant range is available. Specialized
engines call this method when their preconditions fail.
"""
function generic_halo_engine(sys::AbstractHierarchicalGridSystem,
        c::AbstractCellIndex, target::Int, connectivity::Connectivity)
    check_halo_level(sys, c, target)
    lc = level(c)
    grid = levelgrid(sys, target)
    target == lc && return RingHaloEngine(neighbors(grid, c, 1; connectivity))
    return outside_walk_engine(sys, c, target, connectivity, IndexedNeighbors())
end

halo_engine(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex,
    target::Int, connectivity::Connectivity) =
    generic_halo_engine(sys, c, target, connectivity)

# Authalic wrapping does not change the hierarchy or adjacency.
halo_engine(sys::AuthalicSystem, c::AbstractCellIndex, target::Int,
    connectivity::Connectivity) = halo_engine(sys.system, c, target, connectivity)


"""
    SubsetHaloIterator(subset, connectivity, engine)

What [`halo`](@ref) returns, with `cells = true`, for a region that is not a
rooted complete subtree. Construction is O(1); iteration uses an O(depth)
frame stack and prunes with [`subset_span`](@ref). Its size is unknown.
"""
struct SubsetHaloIterator{S,K<:Connectivity,E}
    subset::S
    connectivity::K
    engine::E
end

Base.iterate(it::SubsetHaloIterator) = iterate(it.engine)
Base.iterate(it::SubsetHaloIterator, state) = iterate(it.engine, state)

# The subset engine union is narrow enough for inference to resolve `eltype(E)`.
Base.eltype(::Type{<:SubsetHaloIterator{S,K,E}}) where {S,K,E} = eltype(E)
Base.IteratorSize(::Type{<:SubsetHaloIterator{S,K,E}}) where {S,K,E} =
    Base.IteratorSize(E)

# Subset engines do not provide a constant-time `length`.
Base.length(it::SubsetHaloIterator) = length(it.engine)

Base.show(io::IO, it::SubsetHaloIterator) = print(io, "SubsetHaloIterator(",
    it.subset, "; connectivity = ", it.connectivity, ")")

# ===========================================================================
# The same walk in complete-grid position space, without first materializing ids.
# ===========================================================================

"""
    HaloPositionIterator(halo, grid)

A halo walk read as `cellposition`s on `grid`, lazily — what [`halo`](@ref)
returns by default, and what [`halo_positions`](@ref) wraps an id walk in.

Yields `Int`, strictly increasing, one per cell of the underlying walk and in
the same order. Everything else is the wrapped iterator's:
[`Base.IteratorSize`](@ref), the `length` that exists on exactly two engines and
on no others, resumability, and the `O(depth)` state.
"""
struct HaloPositionIterator{I,G}
    halo::I
    grid::G
end

# Subtree positions refer to the target level; subset positions refer to the
# complete grid from which the subset was drawn.
_halo_grid(it::SubtreeHaloIterator) = levelgrid(it.system, it.level)
_halo_grid(it::SubsetHaloIterator) = it.engine.grid

"""
    halo_positions(it) -> HaloPositionIterator

An id halo walk read as POSITIONS on the grid it was cut from: strictly
increasing, lazily, with the walk's own `O(depth)` state. `halo(region)` already
answers in positions; this is the wrapper it uses, for a walk obtained with
`cells = true`.
"""
halo_positions(it::Union{SubtreeHaloIterator,SubsetHaloIterator}) =
    HaloPositionIterator(it, _halo_grid(it))

# Engines may override the default `cellposition` conversion when state already
# contains the position.
@inline _halo_position(::Any, grid, x, state) = cellposition(grid, x)::Int

Base.eltype(::Type{<:HaloPositionIterator}) = Int
Base.IteratorSize(::Type{<:HaloPositionIterator{I,G}}) where {I,G} =
    Base.IteratorSize(I)

# Position conversion does not change the wrapped iterator's countability.
Base.length(it::HaloPositionIterator) = length(it.halo)

Base.show(io::IO, it::HaloPositionIterator) =
    print(io, "halo_positions(", it.halo, ")")

function Base.iterate(it::HaloPositionIterator)
    r = iterate(it.halo)
    r === nothing && return nothing
    x, s = r
    return (_halo_position(it.halo.engine, it.grid, x, s), s)
end

function Base.iterate(it::HaloPositionIterator, state)
    r = iterate(it.halo, state)
    r === nothing && return nothing
    x, s = r
    return (_halo_position(it.halo.engine, it.grid, x, s), s)
end

# ===========================================================================
# An approximate size, which is deliberately not a `length`
# ===========================================================================

# The per-engine half of `sizehint`, which is where the estimates live. A new
# engine has none unless it defines one explicitly.
_halo_sizehint(::Any) = nothing

# A one-ring has an exact declared length.
_halo_sizehint(e::RingHaloEngine) = length(e)

# ===========================================================================
# The adjacency providers
#
# Candidate enumeration is separate from adjacency testing. The indexed
# provider uses native one-rings; the geometry provider compares boundaries.
# ===========================================================================

"""
    IndexedNeighbors()

Adjacency by the system's native one-ring: `x` touches the subtree when one of
its neighbours has `root` as its ancestor. `O(degree · depth)` per candidate and
allocation-free wherever `neighbors` is.
"""
struct IndexedNeighbors end

"""
    ForcedGeometry()

Adjacency by unit-sphere boundary comparison: `x` touches the subtree when its
boundary shares a vertex ([`Vertex`](@ref)) or two ([`Edge`](@ref)) with the
boundary of some target-level descendant of `root`.

Descendants are visited through a second pruned hierarchy walk using `O(depth)`
memory. This provider is independent of [`IndexedNeighbors`](@ref)'s index
arithmetic and can validate specialized engines.

Comparison is against **target-level descendant** boundaries, not `root`'s own
polygon: H3, IGeo7 and A5 descendants can overhang their parent.
"""
struct ForcedGeometry end

"""
    SubsetMembership(subset, complete)

Adjacency to a SUBSET rather than to a subtree: `x` is a halo cell when the
subset does not hold it and does hold one of its neighbours,

    cellposition(subset, x) === nothing &&
        any(nb -> cellposition(subset, nb) !== nothing,
            neighbors(complete, x, 1; connectivity))

This predicate tests both non-membership and contact because an arbitrary subset
has no root ancestor or descendant range. Consequently, an absent interior cell
is emitted when it touches a member.

`complete` is `levelgrid(system, level)`, because subset adjacency is clipped to
membership. Membership costs O(log(number of windows)) for [`CellVector`](@ref)
and O(log(number of cells)) for [`PartialGrid`](@ref).
"""
struct SubsetMembership{S,G}
    subset::S
    complete::G
end

# ---------------------------------------------------------------------------
# The subset's shape, read as position spans
# ---------------------------------------------------------------------------

const _SPAN_NONE = 0
const _SPAN_SOME = 1
const _SPAN_ALL = 2

"""
    subset_span(subset, lo::Int, hi::Int) -> Int

How much of the position block `lo:hi` of the subset's own complete level the
subset holds: `_SPAN_NONE`, `_SPAN_SOME` or `_SPAN_ALL`.

The block is a node's [`descendant_range`](@ref). Classification costs
O(log(number of windows)) for [`CellVector`](@ref) and O(log(number of cells))
for [`PartialGrid`](@ref), without scanning the block.

Both containers store cells in ascending position order, so full containment is
decided by the endpoints plus a count: the
stored entries between them number `hi - lo + 1` only if none is missing.
"""
function subset_span end

# --- the indexed test -------------------------------------------------------

# A candidate touches the subtree when its native one-ring contains a descendant
# of the root.
@inline function _touches_root(sys, grid, root, rootlevel::Int, x,
        connectivity::Connectivity)
    for nb in neighbors(grid, x, 1; connectivity)
        ancestor(sys, nb, rootlevel) == root && return true
    end
    return false
end

@inline _touches_subtree(::IndexedNeighbors, e, x) =
    _touches_root(e.system, e.grid, e.root, e.rootlevel, x, e.connectivity)

# --- the subset test --------------------------------------------------------

# For subsets, test non-membership and adjacency to a member.
@inline function _touches_subtree(p::SubsetMembership, e, x)
    cellposition(p.subset, x) === nothing || return false
    for nb in neighbors(p.complete, x, 1; connectivity = e.connectivity)
        cellposition(p.subset, nb) === nothing || return true
    end
    return false
end

# --- the geometry test ------------------------------------------------------

# The same rule the generic geometric `neighbors` fallback uses
# (`src/fallbacks/locate.jl`): shared corners within a tolerance scaled from the
# candidate's own shortest edge. One shared vertex is contact, two is an edge.
@inline _needed_contacts(::Vertex) = 1
@inline _needed_contacts(::Edge) = 2

function _touches_subtree(::ForcedGeometry, e, x)
    xb = cell_boundary(e.grid, x)
    tol = _match_tolerance(xb)
    needed = _needed_contacts(e.connectivity)
    xcap = cell_cap(e.grid, x)
    # At depth zero, compare directly with the root; a descendant cursor has no
    # target-level child to expand.
    if e.rootlevel == e.target
        intersects_cap(cell_cap(e.grid, e.root), xcap) || return false
        return _shared_vertices(xb, cell_boundary(e.grid, e.root), tol) >= needed
    end
    return _descendant_touches(e, xb, tol, needed, xcap)
end

# Walk target-level descendants lazily and prune nodes whose extents miss the
# candidate cap. Uses O(depth) memory.
function _descendant_touches(e, xb, tol::Float64, needed::Int, xcap)
    sys = e.system
    st = _empty_walk_stack(e.root)
    st = Helpers.small_push(st, HaloFrame(e.root, 0x1))
    while !isempty(st)
        k = length(st)
        f = @inbounds st[k]
        # Treat target-level frames as leaves; `children` may reject max-level
        # cells.
        if level(f.cell) == e.target
            st = Helpers.small_pop(st)
            continue
        end
        kids = children(sys, f.cell)
        if f.next > length(kids)
            st = Helpers.small_pop(st)
            continue
        end
        d = @inbounds kids[f.next]
        st = Helpers.small_setlast(st, HaloFrame(f.cell, f.next + 0x1))
        if level(d) == e.target
            intersects_cap(cell_cap(e.grid, d), xcap) || continue
            _shared_vertices(xb, cell_boundary(e.grid, d), tol) >= needed &&
                return true
            continue
        end
        intersects_cap(node_extent(sys, d), xcap) || continue
        st = Helpers.small_push(st, HaloFrame(d, 0x1))
    end
    return false
end

# ===========================================================================
# The outside-first walk
# ===========================================================================

# An outside-first walk considers each candidate once in canonical order, so it
# needs neither a seen set nor a final sort.
#
# Cap pruning is sound when a halo cell `x` shares a point with some
# descendant `d` of `root`. `node_extent(sys, root)` provably contains every
# boundary point of every descendant of `root`, so that shared point lies in it;
# `node_extent(sys, n)` likewise contains `x`'s boundary for any ancestor `n` of
# `x`. So the two caps intersect, and a node whose cap misses the root cap can
# contain no halo cell.
#
# This requires native neighbours to share a boundary point. A system with
# topological adjacency between geometrically disjoint cells must provide its
# own `halo_engine` or widen the pruning cap.
#
# Unit-sphere caps avoid longitude/latitude seam and pole degeneracies.
#
# ---------------------------------------------------------------------------
# Subset pruning uses position spans rather than geometry.
#
# `SubsetMembership` has no root, so there is no root cap to compare against and
# no covering law to lean on. What it has instead is the subset's own position
# spans, and the prune those support is the COARSE-CONTAINMENT LAW:
#
#     for every pair of cells `x`, `y` that the system calls VERTEX-adjacent at
#     level `l`, `parent(y)` is `parent(x)` or a vertex-neighbour of it.
#
# This implies that a
# neighbour of a level-`t` cell has its level-`lc` ancestor inside the CLOSED
# one-ring of that cell's own level-`lc` ancestor. If neither a node nor any
# level-`lc` neighbour of it holds a member, no target-level descendant of it can
# have a member neighbour, and none of them is a halo cell. That is what
# `_near_subset` tests, so traversal follows the subset boundary.
#
# A system whose refinement violates this coarse-containment law must provide a
# custom subset engine.
#
# Probe with `Vertex()` because an edge halo is a subset of the vertex halo and
# the coarse-containment law is required only for vertex adjacency.

# One frame per level strictly above the target, so a full-depth walk from the
# root generation pushes at most `maxlevel` of them — 30 on S2, the deepest
# registered system. 34 is that plus four spare.
const _HALO_STACK_CAP = 34

struct HaloFrame{C}
    cell::C
    next::UInt8
end

const HaloStack{C} = Helpers.SmallList{_HALO_STACK_CAP,HaloFrame{C}}

@inline _empty_walk_stack(c::C) where {C} =
    Helpers.empty_small_list(Val(_HALO_STACK_CAP), HaloFrame(c, 0x0))

"""
    OutsideWalkEngine(system, grid, root, rootlevel, target, lo, hi, rootcap,
                      roots, provider, connectivity)

The correctness fallback: outside cells in canonical hierarchy order, the
subject subtree skipped whole, nodes pruned by cap, survivors tested by
`provider`.

Each outside cell is considered at most once, so there is no seen-set and no
final sort; the walk state is one `O(depth)` frame stack of `(cell, next child)`
pairs, and the geometry provider's descendant cursor is a second one. This path
may do more work than an indexed specialization — it is the correctness
fallback, and with [`ForcedGeometry`](@ref) it is the independent oracle those
specializations are checked against.

Requires more than [`has_sorted_subtrees`](@ref), which promises only that a
subtree's target-level descendants are CONTIGUOUS in position. This walk emits
in the order it meets cells, with no sort to repair it, so it additionally
requires that `children(sys, c)` and `rootcells(sys)` are each ordered by their
elements' TARGET-LEVEL descendant ranges — sibling `i` before sibling `j`
exactly when `first(descendant_range(sys, kids[i], target)) <
first(descendant_range(sys, kids[j], target))`, at every level and for every
target. Contiguity without that ordering still produces contiguous blocks, just
visited out of order, and the walk would emit a mis-sorted halo with no error
raised anywhere. Every bundled system satisfies it, and
`test/systems/crosssystem/subtree_halos.jl` is what says so: its law compares
this walk element for element against an ascending-POSITION scan of the target
level. A system that does not satisfy it must supply its own `halo_engine`
rather than inherit this one.
"""
struct OutsideWalkEngine{S,G,C,R,P,K}
    system::S
    grid::G
    root::C
    rootlevel::Int
    target::Int
    lo::Int
    hi::Int
    rootcap::Cap
    roots::R
    provider::P
    connectivity::K
end

Base.eltype(::Type{<:OutsideWalkEngine{S,G,C,R,P,K}}) where {S,G,C,R,P,K} = C
Base.IteratorSize(::Type{<:OutsideWalkEngine}) = Base.SizeUnknown()

struct OutsideWalk{C}
    top::Int
    stack::HaloStack{C}
end

const _HALO_SKIP = 0
const _HALO_EMIT = 1
const _HALO_DESCEND = 2

# Only the forced-geometry provider applies a target-level cap prune. Its
# adjacency check walks the root subtree, so the cap comparison can avoid much
# more work. Native indexed adjacency and subset membership are cheaper than
# constructing the cap and proceed directly to their checks.
#
# Sound for the same reason the internal-node prune is: a halo cell shares a
# boundary point with a descendant of `root`, that point is inside `rootcap` by
# the covering law, and it is one of the candidate's own corners, hence inside
# the candidate's cap.
@inline _target_prune(::IndexedNeighbors, e, c) = true
@inline _target_prune(::ForcedGeometry, e, c) =
    intersects_cap(cell_cap(e.grid, c), e.rootcap)

# Subtree providers prune by ancestry and geometry; `SubsetMembership` prunes
# by position spans. Provider dispatch keeps each admission path monomorphic.
@inline _admit(e::OutsideWalkEngine, c) = _admit(e.provider, e, c)

# --- the subtree arm --------------------------------------------------------

# Three questions, cheapest first. A node whose whole descendant range sits
# inside the subject's is the subject subtree itself (ranges nest or are
# disjoint), and integer comparison retires it without touching geometry.
#
# Containment is checked only at the root level. Root descendants are
# unreachable because the root is skipped whole; every other deeper subtree
# has a disjoint descendant range.
@inline function _admit(p::Union{IndexedNeighbors,ForcedGeometry},
        e::OutsideWalkEngine, c)
    lc = level(c)
    if lc == e.rootlevel
        r = descendant_range(e.system, c, e.target)
        (first(r) >= e.lo && last(r) <= e.hi) && return _HALO_SKIP
    end
    if lc == e.target
        _target_prune(p, e, c) || return _HALO_SKIP
        return _touches_subtree(p, e, c) ? _HALO_EMIT : _HALO_SKIP
    end
    intersects_cap(node_extent(e.system, c), e.rootcap) || return _HALO_SKIP
    return _HALO_DESCEND
end

# --- the subset arm ---------------------------------------------------------

# Two integer questions and no geometry at all — see `subset_halo_engine` for
# why the cap is gone and what replaces it.
#
#   * a node the subset holds ENTIRELY holds no halo cell, because a halo cell
#     is by definition one the subset does not hold. One `subset_span` retires
#     the whole block, which is what makes a chunk's own interior free;
#   * a node the subset touches at all must be descended, because the boundary
#     between held and unheld cells is inside it;
#   * a node the subset does not touch is descended only if a level-`lc`
#     NEIGHBOUR of it does — the coarse-containment law below.
@inline function _admit(p::SubsetMembership, e::OutsideWalkEngine, c)
    lc = level(c)
    if lc == e.target
        return _touches_subtree(p, e, c) ? _HALO_EMIT : _HALO_SKIP
    end
    r = descendant_range(e.system, c, e.target)
    s = subset_span(p.subset, first(r), last(r))
    s == _SPAN_ALL && return _HALO_SKIP
    s == _SPAN_SOME && return _HALO_DESCEND
    return _near_subset(p, e, c, lc) ? _HALO_DESCEND : _HALO_SKIP
end

# The one-ring of `c` at `c`'s OWN level, asked of the subset's spans. The probe
# is `Vertex()` whatever connectivity was requested, for the seam band's reason:
# the `Edge()` halo is a subset of the `Vertex()` one, so one conservative
# superset serves both, and the requested connectivity is still what
# `_touches_subtree` filters by at the target level.
@inline function _near_subset(p::SubsetMembership, e::OutsideWalkEngine, c, lc::Int)
    coarse = levelgrid(e.system, lc)
    for nb in neighbors(coarse, c, 1; connectivity = Vertex())
        r = descendant_range(e.system, nb, e.target)
        subset_span(p.subset, first(r), last(r)) == _SPAN_NONE || return true
    end
    return false
end

Base.iterate(e::OutsideWalkEngine{S,G,C}) where {S,G,C} =
    iterate(e, OutsideWalk(0, _empty_walk_stack(e.root)))

function Base.iterate(e::OutsideWalkEngine{S,G,C}, w::OutsideWalk{C}) where {S,G,C}
    top = w.top
    st = w.stack
    while true
        if isempty(st)
            top += 1
            top > length(e.roots) && return nothing
            c = @inbounds e.roots[top]
            v = _admit(e, c)
            v == _HALO_SKIP && continue
            v == _HALO_EMIT && return (c, OutsideWalk(top, st))
            st = Helpers.small_push(st, HaloFrame(c, 0x1))
            continue
        end
        k = length(st)
        f = @inbounds st[k]
        kids = children(e.system, f.cell)
        if f.next > length(kids)
            st = Helpers.small_pop(st)
            continue
        end
        c = @inbounds kids[f.next]
        st = Helpers.small_setlast(st, HaloFrame(f.cell, f.next + 0x1))
        v = _admit(e, c)
        v == _HALO_SKIP && continue
        v == _HALO_EMIT && return (c, OutsideWalk(top, st))
        st = Helpers.small_push(st, HaloFrame(c, 0x1))
    end
end

# ===========================================================================
# The linear scan: canonical order without a descendant range
#
# A5 uses the canonical-order scan because it does not provide
# `descendant_range`. The scan costs O(ncells) time and O(1) memory. A faster
# engine would require a validated boundary traversal and canonical ordering of
# descendant streams.
# ===========================================================================

"""
    ScanHaloEngine(system, grid, root, rootlevel, target, provider, connectivity)

Every cell of the target level in position order, the descendants skipped and
the rest tested. `O(1)` memory and canonical by construction, but `O(ncells)`
time — the price of a system with no [`descendant_range`](@ref) to prune by, and
A5 is the only one. See the comment above this type for what a dedicated A5
engine would have to prove first.
"""
struct ScanHaloEngine{S,G,C,P,K}
    system::S
    grid::G
    root::C
    rootlevel::Int
    target::Int
    provider::P
    connectivity::K
end

Base.eltype(::Type{<:ScanHaloEngine{S,G,C,P,K}}) where {S,G,C,P,K} = C
Base.IteratorSize(::Type{<:ScanHaloEngine}) = Base.SizeUnknown()

# Outside-ness, stated because the scan has nothing structural to lean on: it
# meets every cell of the level, the subject's own descendants included, where
# `_admit` would already have retired the subject subtree by integer range
# containment. `SubsetMembership` is the one provider whose `_touches_subtree`
# decides outside-ness itself, so it answers `true` here rather than paying for
# the membership search twice — and it must, since a subset has no ancestor to
# compare and a `rootlevel` below the shallowest level would throw.
@inline _scan_outside(::Union{IndexedNeighbors,ForcedGeometry}, e, x) =
    ancestor(e.system, x, e.rootlevel) != e.root
@inline _scan_outside(::SubsetMembership, e, x) = true

# The scan state already contains the position following the emitted cell.
@inline _halo_position(::ScanHaloEngine, ::Any, ::Any, state::Int) = state - 1

Base.iterate(e::ScanHaloEngine) = iterate(e, 1)
function Base.iterate(e::ScanHaloEngine, p::Int)
    n = ncells(e.grid)
    while p <= n
        x = cellindex(e.grid, p)
        if _scan_outside(e.provider, e, x) && _touches_subtree(e.provider, e, x)
            return (x, p + 1)
        end
        p += 1
    end
    return nothing
end

# ===========================================================================
# Building either one
# ===========================================================================

function outside_walk_engine(sys::AbstractHierarchicalGridSystem,
        c::AbstractCellIndex, target::Int, connectivity::Connectivity, provider)
    grid = levelgrid(sys, target)
    lc = level(c)
    has_sorted_subtrees(sys) ||
        return ScanHaloEngine(sys, grid, c, lc, target, provider, connectivity)
    r = descendant_range(sys, c, target)
    return OutsideWalkEngine(sys, grid, c, lc, target, first(r), last(r),
        node_extent(sys, c), rootcells(sys), provider, connectivity)
end

"""
    geometry_halo_engine(sys, c, target, connectivity)

The generic walk forced onto [`ForcedGeometry`](@ref), whatever fast path the
system would otherwise take. Not a public verb and not reachable from
[`SubtreeHaloIterator`](@ref)'s keyword constructor: it exists so a test can
build the oracle explicitly and hand it to the positional constructor.
"""
geometry_halo_engine(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex,
        target::Int, connectivity::Connectivity) =
    outside_walk_engine(sys, c, target, connectivity, ForcedGeometry())

"""
    subset_halo_engine(sys, subset, complete, target, connectivity)

Return the outside-first engine for an arbitrary subset. It uses
[`SubsetMembership`](@ref) and subset position spans instead of a subtree root
range or cap. The subtree-skip fields are inert so an absent cell inside the
subset's span can still be emitted as a halo cell. Nodes are descended only when
the subset touches them or one of their same-level neighbours.

`O(depth)` memory beyond the subset's own storage. A system with no
[`descendant_range`](@ref) gets [`ScanHaloEngine`](@ref), for the same reason it
does at a subtree.

This is not a system extension point because its subject is an arbitrary
membership predicate rather than a root cell.
"""
function subset_halo_engine(sys::AbstractHierarchicalGridSystem, subset,
        complete::AbstractGrid, target::Int, connectivity::Connectivity)
    provider = SubsetMembership(subset, complete)
    seed = first(rootcells(sys))
    nolevel = first(levels(sys)) - 1
    has_sorted_subtrees(sys) ||
        return ScanHaloEngine(sys, complete, seed, nolevel, target, provider,
            connectivity)
    return OutsideWalkEngine(sys, complete, seed, nolevel, target, 1, 0,
        full_sphere_cap(), rootcells(sys), provider, connectivity)
end

# ===========================================================================
# The exterior-perimeter walk, shared by the quad-face family
# ===========================================================================

# A subtree of HEALPix, S2 or ISEA4R is an aligned `side x side` block in one
# face's (diamond's) square lattice, and its halo is the width-1 band around
# it. Two cases, one walk.
#
# Away from a face edge the band uses the in-face 3×3 lattice and equals the halo,
# with the four diagonal-contact corners dropped under `Edge()`. Nothing needs
# checking and nothing needs a seam table.
#
# Across a seam, build one conservative candidate rectangle per touched face;
# filter every candidate through the native one-ring before yielding.
#
# Each face occupies the contiguous id range
# `[face·n², (face+1)·n²)` and faces are numbered ascending, so global canonical
# order is (face ascending, in-face curve code ascending). Walking the rectangle
# list in face order is therefore already a canonical merge of the face-local
# streams: no heap, no seen-set, no sort. Merging same-face rectangles into
# their bounding rectangle is what makes it one stream per face, which is what
# makes it impossible to emit a cell twice.
#
# The band is not one curve interval. The engine
# interval of any one subtree, so there is no `base + offset` to run. The engine
# descends each rectangle's whole FACE in curve order and prunes every quadrant
# that misses the rectangle: curve order at each level is ascending id by
# construction, so survivors arrive ascending. A quadrant of side `h` survives
# only if it meets the rectangle, and at most `O(perimeter/h + 1)` of each size
# do, so the walk visits `O(halo + depth)` nodes in `O(depth)` memory — the
# descent tracks the running curve code and lattice origin and restores both on
# pop by re-reading the parent's last step.
#
# At depth zero, `side == 1` routes to `RingHaloEngine`: a
# one-cell band is the plain eight-cell lattice neighbourhood, which is wrong at
# every seam (a cube-corner cell has seven neighbours, a HEALPix degree-3 vertex
# seven, an ISEA4R cell at icosahedral vertex 0 or 11 nine), and the native
# one-ring answers it exactly in one call.

# `mask` would be dead weight here — pruning is by lattice overlap, not by
# flush sides — so the frame is the two bytes the descent actually needs.
# `HaloFrame` above is the outside walk's frame and carries a cell; this one
# carries a curve state, so the two cannot share a name.
struct SquareBandFrame
    orientation::UInt8
    next::UInt8
end

const SquareBandStack = Helpers.SmallList{_SQUARE_CAP,SquareBandFrame}

@inline _empty_band_stack() =
    Helpers.empty_small_list(Val(_SQUARE_CAP), SquareBandFrame(0x0, 0x0))

"""
    FaceRect(face, orientation, x0, y0, x1, y1)

One face's candidate rectangle: the inclusive lattice box `[x0, x1] x [y0, y1]`
on 0-based `face`, to be descended under curve state `orientation` (the state
that face's ROOT is read under, from [`face_orientation`](@ref)).

The rectangles of a [`SquareBandEngine`](@ref) are one per face and sorted by
`face`, which is what makes walking them a canonical merge.

`Int32` BOUNDS BIND AT LEVEL 32, NOT AT `maxlevel`. A level-`l` lattice
coordinate runs to `2^l - 1`, so `Int32` holds one through level 31
(`2^31 - 1 == typemax(Int32)`) and overflows at level 32. S2's `maxlevel` of 30
is the deepest registered system, so there is exactly ONE level of headroom, and
the quantity to compare a future `maxlevel` bump against is 31 — not 30, and
not `_SQUARE_CAP`. Past it the failure is an `InexactError` raised by this
constructor from inside `square_halo_engine`, i.e. from iterator construction,
which is loud but says nothing about the cause; widen these six fields to
`Int64` (they are `Int32` only to keep `_BAND_RECT_CAP` rectangles inline and
cheap to copy) rather than clamping. `test/systems/crosssystem/subtree_halos.jl`
walks a `maxlevel` block on all three systems, so the level-31 boundary is
approached from one level below on every run.
"""
struct FaceRect
    face::Int32
    orientation::UInt8
    x0::Int32
    y0::Int32
    x1::Int32
    y1::Int32
end

# Twelve is HEALPix's face count and the largest of the three (ISEA4R has ten
# diamonds, S2 six faces); one rectangle per face is the hard ceiling, because
# same-face rectangles are merged. Measured worst case is seven.
const _BAND_RECT_CAP = 12
const BandRects = Helpers.SmallList{_BAND_RECT_CAP,FaceRect}

@inline _empty_band_rects() = Helpers.empty_small_list(Val(_BAND_RECT_CAP),
    FaceRect(0, 0x0, 0, 0, 0, 0))

# Merge by face, so each face appears once and no cell can be reached twice.
#
# `r.orientation` is DISCARDED and the incumbent's kept, which is correct because
# the branch is only taken when `q.face == r.face` and both orientations came
# from `face_orientation(sys, face)` — a pure function of the face alone, with no
# dependence on the cell, the level or the rectangle. So the two values are
# necessarily equal and the merge has no choice to make. That is an invariant of
# the `face_orientation` contract (`src/interface/system.jl`): a system whose
# orientation varied within a face would break the whole descent, not only this
# line, because `SquareBandEngine` seeds each rectangle's descent at that face's
# ROOT and reads no per-cell state at all.
function _merge_rect(rects::BandRects, r::FaceRect)
    for i in 1:length(rects)
        q = @inbounds rects[i]
        q.face == r.face || continue
        return Helpers.small_setindex(rects, FaceRect(q.face, q.orientation,
            min(q.x0, r.x0), min(q.y0, r.y0),
            max(q.x1, r.x1), max(q.y1, r.y1)), i)
    end
    return Helpers.small_push(rects, r)
end

# ---------------------------------------------------------------------------
# Checked and unchecked, as a type rather than a flag
#
# The in-face band is EXACT — proved by the interval guard in
# `square_halo_engine`, and the reason no `neighbors` call appears on that path.
# The seam band is a conservative superset and every candidate owes the native
# test. Carrying the difference as the engine's second type parameter keeps both
# paths monomorphic: the exact path's emit rule inlines to the corner test and
# the checked path's to the one-ring, with no branch on a field in either.
# ---------------------------------------------------------------------------

"""
    NoCheck()

The emit rule of an exact band: a surviving leaf IS a halo cell, subject only to
`Edge()` dropping the four diagonal corners. Zero-size, so an engine carrying it
costs nothing for the distinction.
"""
struct NoCheck end

"""
    NativeCheck(system, grid, root, rootlevel, connectivity)

The emit rule of a conservative band: a candidate is a halo cell only if the
system's own one-ring puts a descendant of `root` next to it,

    any(nb -> ancestor(sys, nb, rootlevel) == root,
        neighbors(grid, x, 1; connectivity))

which is the exact definition, applied to every candidate before it is yielded.
`O(degree · depth)` per candidate.
"""
struct NativeCheck{S,G,C,K}
    system::S
    grid::G
    root::C
    rootlevel::Int
    connectivity::K
end

"""
    SquareBandEngine(curve, check, level, faceside, homeface, x0, y0, side, corners, rects)

The halo of the `side x side` block at lattice origin `(x0, y0)` of `homeface`,
in ascending id, by one pruned quadtree descent per rectangle in `rects` — which
are one per face and ascending by face, so the concatenation is already the
canonical merge.

`faceside` is a face's full lattice side at `level`. Yields [`LevelIndex`](@ref)
on [`SquareBorderEngine`](@ref)'s reasoning and takes the same
[`quadrant_step`](@ref) curves. `O(candidates + depth)` time, `O(depth)` memory.

`check` decides both the emit rule and the count contract:

  - [`NoCheck`](@ref) — the block is nowhere flush with its face edge, the one
    rectangle is the width-1 band, and band and halo are the same set. The count
    is closed form, `4·side + 4` or `4·side` under `Edge()`, so
    [`Base.IteratorSize`](@ref) is `HasLength()`.
  - [`NativeCheck`](@ref) — the block is flush somewhere, the rectangles are a
    conservative superset, and every candidate is tested before it is yielded.
    No perimeter formula survives a seam (a cube corner is three cells, not
    four; an ISEA4R icosahedral vertex is five), so `IteratorSize` is
    `SizeUnknown()` and there is **no `length` method at all**.

SIZE, WHICH IS NOT FREE EVEN WHERE TIME IS. `rects` is a fixed
`_BAND_RECT_CAP`-slot inline list, so an engine is 352 bytes on the in-face path
(`SquareBandEngine{MortonCurve,NoCheck}`), 296 of them the rectangle list with
exactly one slot used. Nothing is heap-allocated and no measured time is spent on
the unused slots — the descent reads `length(rects)`, never the capacity — but
the in-face path is "free" in TIME only, and an engine returned by value costs
that copy. Shrinking it would mean a second engine type for the one-rectangle
case, which is not worth two monomorphic walks.
"""
struct SquareBandEngine{V,K}
    curve::V
    check::K
    level::Int
    faceside::Int64
    homeface::Int32
    x0::Int64
    y0::Int64
    side::Int64
    # `Vertex()`, i.e. keep the four diagonal-contact corners. Read only by
    # `_band_emit(::NoCheck, ...)`: under `NativeCheck` the one-ring is the
    # filter and this field is inert, since the connectivity the check applies
    # is the one stored on the `NativeCheck` itself.
    corners::Bool
    rects::BandRects
end

"""
    SquareBandWalk(rect, stack, code, x, y)

The walk state: which rectangle is being descended, the frame stack, and the
within-face curve code and lattice origin of the sub-square on top of it. Frame
`k`'s node has side `faceside >> (k - 1)`, so none of the last three needs
storing per frame — all are restored on pop by replaying the parent's last
[`quadrant_step`](@ref).
"""
struct SquareBandWalk
    rect::Int
    stack::SquareBandStack
    code::Int64
    x::Int64
    y::Int64
end

@inline _rect_overlaps(r::FaceRect, x::Int64, y::Int64, h::Int64) =
    x <= r.x1 && x + h - 1 >= r.x0 && y <= r.y1 && y + h - 1 >= r.y0

# The block is not its own halo. On the home face this prunes it whole — one
# integer test retires a quadrant of any size — and it is needed on BOTH paths:
# a block cell's neighbours include other block cells, so the native check would
# happily admit one.
@inline _inside_block(e::SquareBandEngine, x::Int64, y::Int64, h::Int64) =
    x >= e.x0 && x + h - 1 <= e.x0 + e.side - 1 &&
    y >= e.y0 && y + h - 1 <= e.y0 + e.side - 1

# The four cells that touch the block at a vertex only, halo under `Vertex()`
# and not under `Edge()`.
@inline _band_corner(e::SquareBandEngine, x::Int64, y::Int64) =
    (x == e.x0 - 1 || x == e.x0 + e.side) &&
    (y == e.y0 - 1 || y == e.y0 + e.side)

@inline _band_emit(::NoCheck, e::SquareBandEngine, cell, cx::Int64, cy::Int64) =
    e.corners || !_band_corner(e, cx, cy)

# The lattice position is the in-face band's business, not the one-ring's: the
# check is the halo's own definition and needs only the cell. Hence the three
# unnamed argument types — the signature exists to match `_band_emit`'s shape,
# and naming arguments it ignores would suggest it consults them.
@inline _band_emit(chk::NativeCheck, ::SquareBandEngine, cell, ::Int64, ::Int64) =
    _touches_root(chk.system, chk.grid, chk.root, chk.rootlevel, cell,
        chk.connectivity)

Base.eltype(::Type{<:SquareBandEngine}) = LevelIndex

Base.IteratorSize(::Type{<:SquareBandEngine{V,NoCheck}}) where {V} =
    Base.HasLength()
Base.IteratorSize(::Type{<:SquareBandEngine{V,K}}) where {V,K<:NativeCheck} =
    Base.SizeUnknown()

# Declared only for the exact band. A seam engine has no closed-form count, so
# the `MethodError` from here is the contract being kept — see the type's
# docstring.
Base.length(e::SquareBandEngine{V,NoCheck}) where {V} =
    e.corners ? Int(4 * e.side + 4) : Int(4 * e.side)

# The exact band's hint is its own count. The seam band's is the band plus two
# cells per corner: a seam corner can contribute a second cell where more than
# three faces meet, which is ISEA4R at icosahedral vertex 0 or 11 and nowhere
# else measured. See `sizehint` for the sweep those two sentences come
# from, and note that a hint three cells generous of the worst case measured is
# a `sizehint!` and not a `length` — the count contract above is untouched.
_halo_sizehint(e::SquareBandEngine{V,NoCheck}) where {V} = length(e)
_halo_sizehint(e::SquareBandEngine) = Int(4 * e.side + 8)

Base.iterate(e::SquareBandEngine) =
    iterate(e, SquareBandWalk(0, _empty_band_stack(),
        Int64(0), Int64(0), Int64(0)))

function Base.iterate(e::SquareBandEngine, w::SquareBandWalk)
    r = w.rect
    st = w.stack
    code = w.code
    x = w.x
    y = w.y
    while true
        if isempty(st)
            # One rectangle finished; start the next face's descent at its root.
            r += 1
            r > length(e.rects) && return nothing
            st = Helpers.small_push(_empty_band_stack(),
                SquareBandFrame((@inbounds e.rects[r]).orientation, 0x0))
            code = Int64(0)
            x = Int64(0)
            y = Int64(0)
            continue
        end
        rect = @inbounds e.rects[r]
        k = length(st)
        f = @inbounds st[k]
        if f.next > 0x3
            st = Helpers.small_pop(st)
            if !isempty(st)
                # Undo the push that entered the frame just dropped: its parent
                # is back on top and has advanced past the position it used.
                pk = length(st)
                pf = @inbounds st[pk]
                p = Int(pf.next) - 1
                i, j, _ = quadrant_step(e.curve, pf.orientation, p)
                half = e.faceside >> pk
                code -= p * half * half
                x -= i * half
                y -= j * half
            end
            continue
        end
        p = Int(f.next)
        st = Helpers.small_setlast(st, SquareBandFrame(f.orientation, f.next + 0x1))
        i, j, o = quadrant_step(e.curve, f.orientation, p)
        half = e.faceside >> k
        cx = x + i * half
        cy = y + j * half
        _rect_overlaps(rect, cx, cy, half) || continue
        rect.face == e.homeface && _inside_block(e, cx, cy, half) && continue
        child = code + p * half * half
        if half == 1
            cell = LevelIndex(e.level,
                Int64(rect.face) * e.faceside * e.faceside + child)
            _band_emit(e.check, e, cell, cx, cy) || continue
            # Emitted without descending, so `code`/`x`/`y` still describe the
            # frame on top and need no restoring on the way back in.
            return (cell, SquareBandWalk(r, st, code, x, y))
        end
        st = Helpers.small_push(st, SquareBandFrame(o, 0x0))
        code = child
        x = cx
        y = cy
    end
end

# ===========================================================================
# Building one: the interval guard, and the seam derivation
# ===========================================================================

"""
    square_halo_engine(sys, curve, c, target, connectivity, x0, y0, side, face, n)

The halo engine for the `side x side` block at lattice origin `(x0, y0)` of
0-based `face`, on a face of side `n` at level `target`. The quad-face family's
[`halo_engine`](@ref border_engine) is this call plus [`lattice_decode`](@ref).

Away from the face edge it is the exact width-1 band, unchecked and counted.
Flush with it, `_seam_band_engine` takes over. `side == 1` never reaches here.
"""
function square_halo_engine(sys::AbstractHierarchicalGridSystem, curve,
        c::AbstractCellIndex, target::Int, connectivity::Connectivity,
        x0::Int64, y0::Int64, side::Int64, face::Int64, n::Int64)
    home = FaceRect(face, face_orientation(sys, face),
        max(Int64(0), x0 - 1), max(Int64(0), y0 - 1),
        min(n - 1, x0 + side), min(n - 1, y0 + side))
    # The width-1 band fits inside the face: every halo cell is in-face, the
    # band IS the halo, and nothing below needs a neighbour query.
    if 1 <= x0 && x0 + side <= n - 1 && 1 <= y0 && y0 + side <= n - 1
        return SquareBandEngine(curve, NoCheck(), target, n, Int32(face),
            x0, y0, side, connectivity isa Vertex,
            Helpers.small_push(_empty_band_rects(), home))
    end
    return _seam_band_engine(sys, curve, c, target, connectivity,
        x0, y0, side, face, n, home)
end

# The quad-face family's wiring. The halo — the outside face of the subtree
# boundary — is the width-1 band around the block, walked lazily by the
# face-quadtree descent. Away from the face edge that band is entirely in-face,
# where adjacency is the plain 3x3 lattice, so the band IS the halo (minus its
# four corners under `Edge()`). Flush with the edge it crosses the seam onto
# other faces, and `square_halo_engine` derives those candidates by asking
# `neighbors` about a few border cells, then filters every one of them with the
# native one-ring. No seam table is read here.
#
# The block's origin comes from the PARENT's `(ix, iy)` shifted left by `d`, not
# from decoding the block's first id: min-Morton is the min corner, but a Hilbert
# block's first position is whichever corner the curve enters by, so decoding it
# would name a different corner per orientation.
#
# `d == 0` is depth zero, which the generic engine answers with the cell's own
# one-ring — exact at the irregular vertices, where a band of one is not.
function halo_engine(sys::AbstractQuadFaceGridSystem, c::LevelIndex, target::Int,
        connectivity::Connectivity)
    check_halo_level(sys, c, target)
    checked_id(sys, c)
    d = target - level(c)
    d == 0 && return generic_halo_engine(sys, c, target, connectivity)
    ix, iy, face = lattice_decode(sys, c)
    return square_halo_engine(sys, subtree_curve(sys), c, target, connectivity,
        Int64(ix) << d, Int64(iy) << d, Int64(1) << d, Int64(face), nside(target))
end

# ---------------------------------------------------------------------------
# Deriving candidate rectangles without a local seam table
#
# Only block sides flush with a face edge can have off-face halo cells. For each
# flush side, the two endpoint cells provide the extreme foreign neighbours.
# The S2, ISEA4R, and HEALPix seam maps are monotone along an edge, so all
# interior images form a contiguous run bounded by those endpoint images.
#
# A block corner that meets a face corner is already an endpoint probe for both
# incident sides. Its vertex-neighbour query also includes any cells reached on
# a third face, so no separate corner rule is required. Probes always use
# `Vertex()` to cover the `Edge()` subset; `NativeCheck` applies the requested
# connectivity when candidates are emitted.
#
# The derived rectangles are exact bounding boxes of probe images. If their
# count exceeds `_BAND_RECT_CAP`, the implementation uses the generic halo
# engine rather than overflowing the fixed-capacity list.
# ---------------------------------------------------------------------------

# The distinct probe positions of a block, at most four — a corner of a block
# flush on two sides is an endpoint of both. Eight is the number of `(side,
# endpoint)` pairs, so the list never has to reject a push.
const _PROBE_CAP = 8
const ProbeList = Helpers.SmallList{_PROBE_CAP,NTuple{2,Int64}}

@inline _empty_probe_list() =
    Helpers.empty_small_list(Val(_PROBE_CAP), (Int64(0), Int64(0)))

@inline function _add_probe(probes::ProbeList, sx::Int64, sy::Int64)
    for i in 1:length(probes)
        p = @inbounds probes[i]
        p[1] == sx && p[2] == sy && return probes
    end
    return Helpers.small_push(probes, (sx, sy))
end

# One probe: everything the native one-ring of the border cell at `(sx, sy)` can
# see, bucketed by face. In-face neighbours already inside the home band box are
# dropped rather than merged, so the home rectangle stays the tight band;
# anything else on the home face — a face adjacent to ITSELF across a seam,
# which none of the three has — would widen it, and correctly.
function _seam_probe(sys, grid, rects::BandRects, target::Int, face::Int64,
        sx::Int64, sy::Int64, home::FaceRect)
    for nb in neighbors(grid, lattice_cell(sys, target, sx, sy, face), 1;
            connectivity = Vertex())
        ix, iy, g = lattice_decode(sys, nb)
        gx = Int64(ix)
        gy = Int64(iy)
        g == face && home.x0 <= gx <= home.x1 && home.y0 <= gy <= home.y1 &&
            continue
        rects = _merge_rect(rects, FaceRect(g, face_orientation(sys, g),
            gx, gy, gx, gy))
    end
    return rects
end

function _seam_band_engine(sys::AbstractHierarchicalGridSystem, curve,
        c::AbstractCellIndex, target::Int, connectivity::Connectivity,
        x0::Int64, y0::Int64, side::Int64, face::Int64, n::Int64,
        home::FaceRect)
    # The faces ARE the root cells, and the root level is `first(levels(sys))` —
    # not 0, which is only the systems bundled here.
    nfaces = ncells(levelgrid(sys, first(levels(sys))))
    nfaces <= _BAND_RECT_CAP ||
        return generic_halo_engine(sys, c, target, connectivity)
    grid = levelgrid(sys, target)
    rects = Helpers.small_push(_empty_band_rects(), home)
    # The two extreme border cells of every flush side — eight positions naming at
    # most four distinct cells, because a flush CORNER is an endpoint of both of
    # its sides and a whole-face block names each of its four corners twice.
    # `_merge_rect` already makes a repeat idempotent, so the deduplication is
    # for cost, not correctness: `_seam_probe` is a `neighbors` call, which
    # allocates a `Vector` on S2, and a whole-face block would otherwise pay for
    # eight of them to learn what four say. Still O(1) in the halo's size either
    # way — construction is a constant-time act however deep the target.
    probes = _empty_probe_list()
    if x0 == 0
        probes = _add_probe(probes, x0, y0)
        probes = _add_probe(probes, x0, y0 + side - 1)
    end
    if x0 + side == n
        probes = _add_probe(probes, x0 + side - 1, y0)
        probes = _add_probe(probes, x0 + side - 1, y0 + side - 1)
    end
    if y0 == 0
        probes = _add_probe(probes, x0, y0)
        probes = _add_probe(probes, x0 + side - 1, y0)
    end
    if y0 + side == n
        probes = _add_probe(probes, x0, y0 + side - 1)
        probes = _add_probe(probes, x0 + side - 1, y0 + side - 1)
    end
    for i in 1:length(probes)
        p = @inbounds probes[i]
        rects = _seam_probe(sys, grid, rects, target, face, p[1], p[2], home)
    end
    # Face order is canonical order, so the list is sorted once, here, by
    # walking the faces rather than the rectangles. Faces are `0:nfaces-1` and
    # there are at most twelve of them, so this is a fixed 144-comparison pass
    # and not a sort in any sense that could grow.
    ordered = _empty_band_rects()
    for g in 0:(nfaces - 1)
        for i in 1:length(rects)
            q = @inbounds rects[i]
            if q.face == g
                ordered = Helpers.small_push(ordered, q)
                break
            end
        end
    end
    # A rectangle whose face is not in `0:nfaces-1` never gets picked up, and a
    # short `ordered` is a SHORT HALO — the one way this walk could answer wrong
    # rather than fall back. `_merge_rect` keeps at most one rectangle per face,
    # so the lengths agree exactly when every face was seen. If they do not, the
    # system's face numbering is not the one the id law assumes, and this file's
    # rule is to fall back, never to approximate — the same answer the
    # `_BAND_RECT_CAP` guard above gives to the other assumption it cannot check.
    length(ordered) == length(rects) ||
        return generic_halo_engine(sys, c, target, connectivity)
    return SquareBandEngine(curve,
        NativeCheck(sys, grid, c, level(c), connectivity),
        target, n, Int32(face), x0, y0, side, connectivity isa Vertex, ordered)
end

# ===========================================================================
# The calibrated directed walk, shared by the two aperture-7 systems
# ===========================================================================

# H3 and IGeo7 subtrees have hexagonal spiral borders rather than rectangular
# lattice bounds. Their border automata expose an arc `(L, s)` of lattice
# directions, and the halo is reached by walking the arcs of neighbouring
# subtrees that face `root`.
#
# Nested adjacency guarantees that every target-level halo cell descends from
# a same-level neighbour of `root`: if two cells are adjacent, their parents
# are equal or adjacent. `_hex_calibrate` determines the exposed arc from the
# neighbour's children that touch `root`, avoiding parity-specific seed tables
# for the two systems.
#
# Same-level neighbour subtrees are disjoint. Sorting neighbours by the first
# position of their target-level descendant range therefore produces ascending,
# non-overlapping candidate blocks, while each border automaton emits its own block
# in ascending id order.
#
# Candidates cannot belong to `root` because a one-ring never contains its
# subject and every candidate descends from a one-ring neighbour. The native
# adjacency check still filters every candidate before emission. At depth one,
# `HexChildHaloEngine` emits the touching children directly without starting a
# border automaton.

# Six is a hexagon's neighbour count and five a pentagon's, so the list is never
# more than six long on either system. Eight is that plus slack, so a system with
# a wider ring reaches the guard in `hex_halo_engine` rather than a `BoundsError`
# from inside `small_push`.
const _HEX_RING_CAP = 8

"""
    HexNeighbour(cell, lo, arclen, start)

One of `root`'s same-level neighbours, with the arc `(arclen, start)` calibrated
to face `root` and `lo = first(descendant_range(sys, cell, target))`, the key the
ring is kept sorted by. `arclen == 0` marks an uncalibrated entry, which only the
depth-one engine holds.
"""
struct HexNeighbour{C}
    cell::C
    lo::Int
    arclen::Int8
    start::Int8
end

const HexRing{C} = Helpers.SmallList{_HEX_RING_CAP,HexNeighbour{C}}

@inline _empty_hex_ring(c::C) where {C} =
    Helpers.empty_small_list(Val(_HEX_RING_CAP), HexNeighbour(c, 0, Int8(0), Int8(0)))

# Insertion sort over at most six entries: the whole list fits in registers and a
# heap would be a heap-allocation. Shifts down from the end, so the list is
# sorted by `lo` at every point and the walk can read it straight through.
function _hex_insert(ring::HexRing{C}, e::HexNeighbour{C}) where {C}
    ring = Helpers.small_push(ring, e)
    i = length(ring)
    while i > 1
        prev = @inbounds ring[i - 1]
        prev.lo <= e.lo && break
        ring = Helpers.small_setindex(ring, prev, i)
        i -= 1
    end
    return Helpers.small_setindex(ring, e, i)
end

# ---------------------------------------------------------------------------
# Calibration
# ---------------------------------------------------------------------------

"""
    _minimal_arc(p, q) -> (arclen, start)

The shortest arc of the six-direction ring containing both `p` and `q`, or
`(0, 0)` when there is not exactly one such arc or its length is not 2 or 3.

The scan is the definition rather than the closed form (`arclen` is the ring
distance plus one, taken the short way round) so that the UNIQUENESS question is
answered by counting, not argued: two directions exactly opposite each other
admit two four-arcs and no shorter one, and this returns `(0, 0)` there instead
of picking one. Thirty-six iterations of integer work, once per neighbour.
"""
function _minimal_arc(p::Int, q::Int)
    bestlen = 7
    beststart = 0
    ties = 0
    for L in 1:6, s in 0:5
        (mod(p - s, 6) < L && mod(q - s, 6) < L) || continue
        if L < bestlen
            bestlen = L
            beststart = s
            ties = 1
        elseif L == bestlen
            ties += 1
        end
    end
    ties == 1 || return (0, 0)
    (bestlen == 2 || bestlen == 3) || return (0, 0)
    return (bestlen, beststart)
end

"""
    _hex_calibrate(sys, grid1, root, rootlevel, nb, connectivity) -> (arclen, start)

The arc of `nb` that faces `root`, derived by asking which of `nb`'s children are
halo cells of `root`'s subtree one level down, or `(0, 0)` if any of the guards
fails.

Measured over 52,182 `(root, neighbour)` pairs on both systems, root levels 0-11:
the number of touching children is ALWAYS exactly two, the minimal covering arc
is always unique, and its length is 2 in 98.6% of cases and 3 in the remaining
1.4% — never 1, never 4 or more. The arc-3 case is fully characterised: `nb` is a
pentagon and its deleted direction lies strictly between the two touching
children. Nothing else predicts it, which is the other reason this is measured
per call rather than looked up.

So none of the returns below was observed to fire. That is what a guard on a
structural claim should look like — it is here because the claim is evidence and
not yet a theorem, and a system whose one-ring changed would meet it rather than
meeting a wrong answer.
"""
function _hex_calibrate(sys, grid1, root, rootlevel::Int, nb,
        connectivity::Connectivity)
    p = -1
    q = -1
    n = 0
    for k in children(sys, nb)
        _touches_root(sys, grid1, root, rootlevel, k, connectivity) || continue
        n += 1
        n > 2 && return (0, 0)
        d = hex_child_direction(sys, k)
        # The centre child is enclosed by its six siblings and cannot reach out
        # of its parent at all, so it has no direction to put on the ring. If one
        # ever tested as touching, the ring model is wrong for this system and
        # the arc must not be guessed.
        d < 0 && return (0, 0)
        n == 1 ? (p = d) : (q = d)
    end
    n == 2 || return (0, 0)
    return _minimal_arc(p, q)
end

"""
    _hex_validate(sys, root, rootlevel, ring, connectivity) -> Bool

Return whether every depth-two halo cell under each neighbour occurs in that
neighbour's calibrated walk. `hex_halo_engine` runs this bounded validation for
depths of at least three. Both sequences are ascending, so containment uses a
two-pointer merge without allocation.
"""
function _hex_validate(sys, root, rootlevel::Int, ring::HexRing{C},
        connectivity::Connectivity) where {C}
    target = rootlevel + 2
    grid = levelgrid(sys, target)
    for i in 1:length(ring)
        e = @inbounds ring[i]
        arc = seeded_border_engine(sys, e.cell, target, Int(e.arclen), Int(e.start))
        r = iterate(arc)
        for k in children(sys, e.cell), x in children(sys, k)
            _touches_root(sys, grid, root, rootlevel, x, connectivity) || continue
            while r !== nothing && r[1] < x
                r = iterate(arc, r[2])
            end
            (r !== nothing && r[1] == x) || return false
            r = iterate(arc, r[2])
        end
    end
    return true
end

# ---------------------------------------------------------------------------
# Depth one: the calibration is the answer
# ---------------------------------------------------------------------------

"""
    HexChildHaloEngine(system, grid, root, rootlevel, target, connectivity, ring)

`target == rootlevel + 1`: each neighbour's children in turn, native-checked, the
neighbours in descendant-range order. No automaton and no calibration — at this
depth the children ARE the candidates.

`children` is re-read per step rather than stored, which keeps the walk state two
integers and costs a bounded-container rebuild of at most seven ids.
"""
struct HexChildHaloEngine{S,G,C,K}
    system::S
    grid::G
    root::C
    rootlevel::Int
    target::Int
    connectivity::K
    ring::HexRing{C}
end

Base.eltype(::Type{<:HexChildHaloEngine{S,G,C,K}}) where {S,G,C,K} = C
Base.IteratorSize(::Type{<:HexChildHaloEngine}) = Base.SizeUnknown()

"""
    HexChildWalk(slot, next)

Which neighbour of the ring is being read, and which of its children comes next.
"""
struct HexChildWalk
    slot::Int
    next::Int
end

Base.iterate(e::HexChildHaloEngine) = iterate(e, HexChildWalk(1, 1))

function Base.iterate(e::HexChildHaloEngine, w::HexChildWalk)
    slot = w.slot
    next = w.next
    while slot <= length(e.ring)
        kids = children(e.system, (@inbounds e.ring[slot]).cell)
        if next > length(kids)
            slot += 1
            next = 1
            continue
        end
        x = @inbounds kids[next]
        next += 1
        _touches_subtree(IndexedNeighbors(), e, x) || continue
        return (x, HexChildWalk(slot, next))
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Deeper: one seeded automaton per neighbour
# ---------------------------------------------------------------------------

"""
    HexArcHaloEngine(system, grid, root, rootlevel, target, connectivity, ring)

`target > rootlevel + 1`: the system's own border automaton seeded with each
neighbour's calibrated arc, walked to `target`, every leaf native-checked before
it is yielded. The ring is in descendant-range order and the blocks are disjoint,
so concatenating the neighbours' streams is already the canonical merge.

Memory is `O(depth)`: one seeded engine and frame stack plus the fixed ring.
[`Base.IteratorSize`](@ref) is `SizeUnknown()` and `length` is not defined. The
formula used by [`sizehint`](@ref DiscreteGlobalGrids.sizehint) has not been
derived for every seeded
transition and therefore is not an exact-length contract.
"""
struct HexArcHaloEngine{S,G,C,K}
    system::S
    grid::G
    root::C
    rootlevel::Int
    target::Int
    connectivity::K
    ring::HexRing{C}
end

Base.eltype(::Type{<:HexArcHaloEngine{S,G,C,K}}) where {S,G,C,K} = C
Base.IteratorSize(::Type{<:HexArcHaloEngine}) = Base.SizeUnknown()

# `3^(d+1) + 3` is exact around a hexagon and conservatively bounds the smaller
# pentagon census. Use it only as an allocation hint.
_halo_sizehint(e::HexChildHaloEngine) = _hex_sizehint(e.target - e.rootlevel)
_halo_sizehint(e::HexArcHaloEngine) = _hex_sizehint(e.target - e.rootlevel)

@inline _hex_sizehint(d::Int) = 3^(d + 1) + 3

"""
    HexArcWalk(slot, arc, state)

The walk state: which neighbour of the ring is being descended, that neighbour's
seeded automaton, and the automaton's own frame stack. All three are isbits, and
the engine is rebuilt rather than stored per slot, so resuming needs nothing that
was not returned.
"""
struct HexArcWalk{A,W}
    slot::Int
    arc::A
    state::W
end

@inline function _hex_arc_engine(e::HexArcHaloEngine, slot::Int)
    nb = @inbounds e.ring[slot]
    return seeded_border_engine(e.system, nb.cell, e.target, Int(nb.arclen),
        Int(nb.start))
end

# Open neighbour `slot` and every later one until a candidate passes the check.
# Written as a fresh-start scan rather than folded into the resume path so the
# walk state is only ever built from a narrowed, concrete iterate result — a
# `Union{Nothing,W}` state field would infect the engine's whole return type.
function _hex_arc_open(e::HexArcHaloEngine, slot::Int)
    while slot <= length(e.ring)
        arc = _hex_arc_engine(e, slot)
        r = iterate(arc)
        while r !== nothing
            x, w = r
            _touches_subtree(IndexedNeighbors(), e, x) &&
                return (x, HexArcWalk(slot, arc, w))
            r = iterate(arc, w)
        end
        slot += 1
    end
    return nothing
end

Base.iterate(e::HexArcHaloEngine) = _hex_arc_open(e, 1)

function Base.iterate(e::HexArcHaloEngine, s::HexArcWalk)
    r = iterate(s.arc, s.state)
    while r !== nothing
        x, w = r
        _touches_subtree(IndexedNeighbors(), e, x) &&
            return (x, HexArcWalk(s.slot, s.arc, w))
        r = iterate(s.arc, w)
    end
    return _hex_arc_open(e, s.slot + 1)
end

# ---------------------------------------------------------------------------
# Building one
# ---------------------------------------------------------------------------

"""
    hex_halo_engine(sys, c, target, connectivity)

Return the calibrated halo engine used by H3 and IGeo7. Fall back to
[`generic_halo_engine`](@ref) whenever a precondition fails.

The system must have sorted subtrees, the target must be deeper than `c`, the
neighbour ring must fit its fixed capacity, every neighbour must calibrate, and
depths of at least three must pass [`_hex_validate`](@ref).

The initial ring probe uses `Vertex()` because an edge halo is a subset of the
vertex halo. Calibration and validation use the requested connectivity. If that
connectivity produces an unsupported arc, the method returns the generic engine.
"""
function hex_halo_engine(sys::AbstractHierarchicalGridSystem,
        c::AbstractCellIndex, target::Int, connectivity::Connectivity)
    lc = level(c)
    (has_sorted_subtrees(sys) && lc < target <= maxlevel(sys)) ||
        return generic_halo_engine(sys, c, target, connectivity)
    nbs = neighbors(levelgrid(sys, lc), c, 1; connectivity = Vertex())
    length(nbs) <= _HEX_RING_CAP ||
        return generic_halo_engine(sys, c, target, connectivity)
    grid = levelgrid(sys, target)
    ring = _empty_hex_ring(c)
    if target == lc + 1
        for nb in nbs
            ring = _hex_insert(ring, HexNeighbour(nb,
                first(descendant_range(sys, nb, target)), Int8(0), Int8(0)))
        end
        return HexChildHaloEngine(sys, grid, c, lc, target, connectivity, ring)
    end
    grid1 = levelgrid(sys, lc + 1)
    for nb in nbs
        arclen, start = _hex_calibrate(sys, grid1, c, lc, nb, connectivity)
        arclen == 0 && return generic_halo_engine(sys, c, target, connectivity)
        ring = _hex_insert(ring, HexNeighbour(nb,
            first(descendant_range(sys, nb, target)), Int8(arclen), Int8(start)))
    end
    if target - lc >= 3 && !_hex_validate(sys, c, lc, ring, connectivity)
        return generic_halo_engine(sys, c, target, connectivity)
    end
    return HexArcHaloEngine(sys, grid, c, lc, target, connectivity, ring)
end
