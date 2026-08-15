# The outside face of a subtree boundary: the level-`l` cells that are NOT
# descendants of `c` but have a neighbour that is.
#
# One public type, one private ENGINE per system, exactly as the rim and
# interior walks are built (`subtree_iterators.jl`). The engine is chosen once,
# at construction, by `halo_engine` — private multiple dispatch, so an engine's
# whole protocol is monomorphic and a system adds a fast path by adding one
# method, not by setting a flag anyone can read.
#
# Every engine here is EXACT. A specialization may enumerate a conservative
# candidate band, but it must run the adjacency test on every candidate before
# yielding it: approximation never reaches this file's public surface.

"""
    SubtreeHaloIterator(sys, c, l; connectivity = Vertex())

The halo of `c`'s subtree at level `l`, lazily: every level-`l` cell that is
**not** a descendant of `c` but has a neighbour that is, in ascending canonical
order, each cell exactly once.

`collect` of this is [`subtree_halo`](@ref), element for element. `l == level(c)`
is `c`'s own one-ring, sorted. `l < level(c)` and `l > max_level(sys)` throw an
`ArgumentError`.

Construction never materialises the halo: the iterator holds `O(depth)` walk
state and bounded native neighbour containers, and nothing sized by the answer.
Taking the first few cells of a deep halo therefore costs what those cells cost,
not what the whole ring would.

[`Base.IteratorSize`](@ref) is `SizeUnknown()` wherever no count is proved —
face seams, poles and pentagons break perimeter formulas — and `HasLength()`
only where an engine derives the count in closed form. There is no `length` that
would silently cost a traversal.

See also [`halo`](@ref) for the same question about a subset, and
[`EdgeCellIterator`](@ref) for the inside face of the same boundary.
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
Base.eltype(::Type{<:SubtreeHaloIterator{S,C,K,E}}) where {S,C,K,E} = eltype(E)
Base.IteratorSize(::Type{<:SubtreeHaloIterator{S,C,K,E}}) where {S,C,K,E} =
    Base.IteratorSize(E)

# Deliberately delegated rather than defined: an engine that cannot count
# without walking defines no `length`, and the `MethodError` from here is the
# contract being kept.
Base.length(it::SubtreeHaloIterator) = length(it.engine)

# `collect` must BE the guarded path, not merely parallel to it: `collect`'s own
# `HasLength` route would size its vector from a miscounting engine and hand
# back an `undef` tail as cell ids. See `collect_subtree`.
Base.collect(it::SubtreeHaloIterator) = collect_subtree(it)

Base.show(io::IO, it::SubtreeHaloIterator) = print(io, "SubtreeHaloIterator(",
    it.system, ", ", it.cell, ", ", it.level, "; connectivity = ",
    it.connectivity, ")")

# ===========================================================================
# The one-ring engine: `l == level(c)`
# ===========================================================================

# The ring comes back in ROTATIONAL order and the halo owes ascending order, so
# the sort is a selection emit over a container of at most `max_neighbors`
# elements: the state is the count emitted plus the last value, both isbits, and
# nothing is copied or allocated. A `sort` here would have to leave the fixed
# capacity container and heap-allocate (see `systems/H3/neighbors.jl`).
"""
    RingHaloEngine(ring)

`c`'s own one-ring, ascending, by selection emit. `O(degree^2)` time with
`degree <= max_neighbors(sys, connectivity)`, no allocation, isbits state.

`length` is `length(ring)`, which is honest only because a native one-ring
never lists the same cell twice — the selection emit yields each DISTINCT value
once, so a repeated neighbour would make the engine yield fewer cells than it
claims. No bundled system repeats; if one ever did, `collect_subtree` is the
backstop that turns the miscount into an `error` rather than an `undef` tail.
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

# Named rather than folded into the method below so a system override has a
# fallback to CALL, not merely to shadow.
function generic_halo_engine(sys::AbstractHierarchicalGridSystem,
        c::AbstractCellIndex, target::Int, connectivity::Connectivity)
    lc = level(c)
    target >= lc || throw(ArgumentError(
        "subtree_halo: level $target is above the cell's own level $lc"))
    target <= max_level(sys) || throw(ArgumentError(
        "subtree_halo: level $target is below the system's deepest level " *
        "$(max_level(sys))"))
    grid = levelgrid(sys, target)
    target == lc && return RingHaloEngine(neighbors(grid, c, 1; connectivity))
    return outside_walk_engine(sys, c, target, connectivity, IndexedNeighbors())
end

halo_engine(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex,
    target::Int, connectivity::Connectivity) =
    generic_halo_engine(sys, c, target, connectivity)

# The ellipsoid wrapper is transparent: the halo of a subtree is a question about
# the discrete hierarchy, which the authalic transform does not touch.
halo_engine(sys::AuthalicSystem, c::AbstractCellIndex, target::Int,
    connectivity::Connectivity) = halo_engine(sys.system, c, target, connectivity)

"""
    subtree_halo(sys, c, l; connectivity = Vertex()) -> Vector

The level-`l` cells outside `c`'s subtree that touch it, ascending. The
explicitly materialising form of [`SubtreeHaloIterator`](@ref) — the halo can be
far larger than the subtree's rim, so collecting it is the caller's decision,
never the constructor's.

The motivating read is a chunk plus its stencil margin:
`descendant_range(sys, c, l)` is the chunk's contiguous position block, and
`cellposition.(Ref(levelgrid(sys, l)), subtree_halo(sys, c, l))` is the extra
fetch list a one-ring stencil needs — with no halo table built at all.
"""
subtree_halo(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex,
        l::Integer; connectivity::Connectivity = Vertex()) =
    collect_subtree(SubtreeHaloIterator(sys, c, l; connectivity))

"""
    halo(subset; connectivity = Vertex())

The cells immediately outside a same-level subset — the subset face of
[`SubtreeHaloIterator`](@ref).

**This function has no methods yet.** The name is declared and exported here so
the two halo verbs live in one file and the docs can cross-reference them; the
methods arrive with the subset containers (`PartialGrid`, `CellVector`,
`CellLookup`). Until then every call is a `MethodError` — deliberately, because
a stub returning `nothing` or an empty iterator would answer a question it
cannot yet answer.
"""
function halo end

# ===========================================================================
# The adjacency providers
#
# Candidate enumeration and adjacency testing are kept apart on purpose: the
# walk below decides WHICH cells to consider, a provider decides whether one of
# them touches the subtree. The indexed provider is the DGGS's own one-ring; the
# geometry provider compares unit-sphere boundaries and shares no topology with
# it, which is what makes it an oracle rather than a second opinion.
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

Descendants are visited through a second pruned hierarchy walk, so this costs
`O(depth)` memory rather than the subtree. It shares no index arithmetic with
[`IndexedNeighbors`](@ref) and is the oracle every specialized engine is tested
against — never `neighbors`, and never [`subtree_border`](@ref), which do share
it.

Comparison is against **target-level descendant** boundaries, not `root`'s own
polygon: H3, IGeo7 and A5 descendants can overhang their parent.
"""
struct ForcedGeometry end

# --- the indexed test -------------------------------------------------------

@inline function _touches_subtree(::IndexedNeighbors, e, x)
    for nb in neighbors(e.grid, x, 1; connectivity = e.connectivity)
        ancestor(e.system, nb, e.rootlevel) == e.root && return true
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
    # Depth zero has no descent to make: `root` IS its own only target-level
    # descendant, so the comparison is against the root's own boundary. Handled
    # here rather than left to the cursor because a cursor seeded at the target
    # level has no children to expand — it would descend past the target to
    # `max_level` and throw.
    if e.rootlevel == e.target
        intersects_cap(cell_cap(e.grid, e.root), xcap) || return false
        return _shared_vertices(xb, cell_boundary(e.grid, e.root), tol) >= needed
    end
    return _descendant_touches(e, xb, tol, needed, xcap)
end

# A lazy cursor over `root`'s target-level descendants, pruned by the candidate's
# own cap: a node whose extent misses `xcap` holds no descendant that can touch
# `x`. `O(depth)` memory, and the descendants are never materialised.
function _descendant_touches(e, xb, tol::Float64, needed::Int, xcap)
    sys = e.system
    st = _empty_walk_stack(e.root)
    st = Helpers.small_push(st, HaloFrame(e.root, 0x1))
    while !isempty(st)
        k = length(st)
        f = @inbounds st[k]
        # Defence in depth behind the depth-zero early return above: a frame at
        # the target level is a leaf of THIS walk, whatever put it there.
        # `children` of a `max_level` cell throws, so never ask.
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

# Why outside-first at all: a descendant-first algorithm (walk the rim, take
# each rim cell's neighbours, drop the ones still inside) produces duplicates
# out of order, and therefore needs a seen-set and a final sort — both sized by
# the halo. Reversing the search removes both. Each outside cell is considered
# at most once, so there is nothing to deduplicate and nothing to sort.
#
# Why the cap prune is sound: a halo cell `x` shares at least a point with some
# descendant `d` of `root`. `node_extent(sys, root)` provably contains every
# boundary point of every descendant of `root`, so that shared point lies in it;
# `node_extent(sys, n)` likewise contains `x`'s boundary for any ancestor `n` of
# `x`. So the two caps intersect, and a node whose cap misses the root cap can
# contain no halo cell.
#
# THE ASSUMPTION THAT COUPLES THE TWO HALVES OF THIS FILE. That argument is
# GEOMETRIC — it is the covering law over boundary POINTS — but the default
# provider's adjacency is TOPOLOGICAL, the system's native one-ring. The step
# "a halo cell shares at least a point with some descendant" is therefore an
# assumption about `neighbors`, not a theorem: it holds only while every pair of
# cells the system calls adjacent has boundaries that share a point. No bundled
# system violates it (checked by the forced-geometry agreement testset). A
# system that did — an adjacency defined by index arithmetic across a seam with
# no shared drawn corner — would keep the prune's arithmetic intact and its
# soundness not, and the walk would silently drop that halo cell. Such a system
# needs its own `halo_engine`, or a rootcap widened to cover the discrepancy.
#
# Why not lon/lat: longitude/latitude boxes are unusable at seams and poles, so
# everything here stays in unit-sphere XYZ.

# One frame per level strictly above the target, so a full-depth walk from the
# root generation pushes at most `max_level` of them — 30 on S2, the deepest
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

# How expensive a cap prune is BEFORE the adjacency test, per provider. The
# indexed provider is one native `neighbors` call, cheaper than building a cap
# from a boundary, so it declines: at the target level it goes straight to the
# test. The forced-geometry provider walks the ROOT's subtree once per
# candidate, so a single cap comparison that retires the candidate outright is
# worth many times its cost — and this is the oracle path, which the
# differential tests hammer. Dispatch on the provider so the indexed path pays
# nothing for the distinction, not a branch on a field.
#
# Sound for the same reason the internal-node prune is: a halo cell shares a
# boundary point with a descendant of `root`, that point is inside `rootcap` by
# the covering law, and it is one of the candidate's own corners, hence inside
# the candidate's cap.
@inline _target_prune(::IndexedNeighbors, e, c) = true
@inline _target_prune(::ForcedGeometry, e, c) =
    intersects_cap(cell_cap(e.grid, c), e.rootcap)

# Three questions, cheapest first. A node whose whole descendant range sits
# inside the subject's is the subject subtree itself (ranges nest or are
# disjoint), and integer comparison retires it without touching geometry.
#
# That containment can only hold AT the root's own level, which is why the guard
# is `==` and not `>=`. Deeper than the root there are two cases and neither can
# be contained: a node that is one of the root's own descendants is unreachable,
# because the root was already skipped whole and the walk never descended into
# it; and a node that is not has a descendant range disjoint from the root's, so
# the containment test is a `descendant_range` call that cannot succeed. On
# IGeo7 with a level-2 root at `l = 5` that is 630 of the 689 nodes the walk
# visits.
@inline function _admit(e::OutsideWalkEngine, c)
    lc = level(c)
    if lc == e.rootlevel
        r = descendant_range(e.system, c, e.target)
        (first(r) >= e.lo && last(r) <= e.hi) && return _HALO_SKIP
    end
    if lc == e.target
        _target_prune(e.provider, e, c) || return _HALO_SKIP
        return _touches_subtree(e.provider, e, c) ? _HALO_EMIT : _HALO_SKIP
    end
    intersects_cap(node_extent(e.system, c), e.rootcap) || return _HALO_SKIP
    return _HALO_DESCEND
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
# ===========================================================================

"""
    ScanHaloEngine(system, grid, root, rootlevel, target, provider, connectivity)

Every cell of the target level in position order, the descendants skipped and
the rest tested. `O(1)` memory and canonical by construction, but `O(ncells)`
time — the price of a system with no [`descendant_range`](@ref) to prune by, and
A5 is the only one.
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

Base.iterate(e::ScanHaloEngine) = iterate(e, 1)
function Base.iterate(e::ScanHaloEngine, p::Int)
    n = ncells(e.grid)
    while p <= n
        x = cellindex(e.grid, p)
        if ancestor(e.system, x, e.rootlevel) != e.root &&
           _touches_subtree(e.provider, e, x)
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
