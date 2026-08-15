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

# ===========================================================================
# The exterior-perimeter walk, shared by the three aperture-4 systems
# ===========================================================================

# A subtree of HEALPix, S2 or ISEA4R is an aligned `side x side` block in one
# face's (diamond's) square lattice, and its halo is the width-1 band around
# it. Two cases, one walk.
#
# AWAY FROM THE FACE EDGE the band is entirely in-face, where adjacency is the
# plain 3×3 lattice on all three systems, so band and halo are the SAME SET —
# with the four diagonal-contact corners dropped under `Edge()`. Nothing needs
# checking and nothing needs a seam table.
#
# ACROSS A SEAM the halo leaves the face: cells on up to six other faces, whose
# ids no in-face arithmetic can name. They are reached by generalising the walk
# from one band to a short, FACE-ORDERED list of candidate rectangles, one per
# face the halo can touch — see `_seam_band_engine` for how those rectangles are
# derived and why they cover. That list is a conservative SUPERSET, so every
# candidate is put through the native one-ring before it is yielded.
#
# WHY CONCATENATION IS A MERGE. Each face occupies the contiguous id range
# `[face·n², (face+1)·n²)` and faces are numbered ascending, so global canonical
# order is (face ascending, in-face curve code ascending). Walking the rectangle
# list in face order is therefore already a canonical merge of the face-local
# streams: no heap, no seen-set, no sort. Merging same-face rectangles into
# their bounding rectangle is what makes it one stream per face, which is what
# makes it impossible to emit a cell twice.
#
# WHY A QUADTREE DESCENT AND NOT AN OFFSET WALK. The band is not a curve
# interval of any one subtree, so there is no `base + offset` to run. The engine
# descends each rectangle's whole FACE in curve order and prunes every quadrant
# that misses the rectangle: curve order at each level is ascending id by
# construction, so survivors arrive ascending. A quadrant of side `h` survives
# only if it meets the rectangle, and at most `O(perimeter/h + 1)` of each size
# do, so the walk visits `O(halo + depth)` nodes in `O(depth)` memory — the
# descent tracks the running curve code and lattice origin and restores both on
# pop by re-reading the parent's last step.
#
# DEPTH ZERO IS NOT HERE. `side == 1` routes to `RingHaloEngine` instead: a
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
on [`SquareRimEngine`](@ref)'s reasoning and takes the same
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

@inline function _band_emit(chk::NativeCheck, e::SquareBandEngine, cell,
        cx::Int64, cy::Int64)
    for nb in neighbors(chk.grid, cell, 1; connectivity = chk.connectivity)
        ancestor(chk.system, nb, chk.rootlevel) == chk.root && return true
    end
    return false
end

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
0-based `face`, on a face of side `n` at level `target`. The three aperture-4
systems' `halo_engine` methods are this call plus their own lattice decode.

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

# ---------------------------------------------------------------------------
# Deriving the candidate rectangles, with no seam table of this file's own
#
# WHAT HAS TO BE COVERED. A halo cell is a neighbour of a block cell. Block
# cells with a neighbour off the face are exactly the rim cells of the sides
# that are FLUSH with the face edge, and their off-face neighbours are the
# images of the extended-lattice positions one step outside that edge. So for a
# block flush on, say, `x = 0`, the foreign candidates are the images of
# `(-1, y')` for `y'` running over `[y0 - 1, y0 + side]` clipped to the face,
# plus — where a corner of the block is a corner of the face — the cells that a
# DIAGONAL step off two edges at once reaches.
#
# WHY TWO PROBES PER SIDE SUFFICE. All three systems' seam maps are edge-to-edge
# affine with a sign: S2's `wrap_xyf` computes `k = (σ·b + n - 1) >> 1` from the
# centred along-edge coordinate, which is `y` or `n - 1 - y`; ISEA4R's
# `lattice_neighbors` reads the paired rim slot at `n - 1 - j`; HEALPix's
# `nested_neighbors` applies the `NB_SWAPARRAY` mirrors and transpose. Each is
# monotone along the edge AT EVERY `n`, so the images of an interval of `y'` are
# a contiguous run lying between the images of its ENDPOINTS.
#
# The two extreme rim cells of a flush side see both endpoints directly: the
# extreme cell at `(0, y0)` has `(-1, y0 - 1)` among its eight neighbours, and
# the one at `(0, y0 + side - 1)` has `(-1, y0 + side)`. So the bounding box of
# what the two probes see already contains the image of every interior rim cell
# of that side, and nothing needs widening. That is not an argument this file
# takes on trust — see the verification note below, and note that the tests
# would fail if a bound were moved one cell inward.
#
# WHERE A THIRD FACE APPEARS. At a flush CORNER the block's corner cell is a
# probe of both flush sides, and its own neighbour list already contains
# whatever the diagonal step reaches: nothing on S2 (`wrap_xyf` returns
# `nothing` at a cube corner) or at HEALPix's `-1` entries in `NB_FACEARRAY` —
# which occur only for the double-out `nbnum`, i.e. only at a face corner, never
# at a run endpoint — and the interior of `CORNER_FANS` on ISEA4R, which is the
# two extra diamonds meeting at icosahedral vertex 0 or 11. So corners need no
# separate rule; they need only that the corner cell is probed, which flushness
# already guarantees.
#
# THE PROBE ASKS FOR `Vertex()` WHATEVER WAS REQUESTED, so one coverage argument
# serves both connectivities: the `Edge()` halo is a subset of the `Vertex()`
# one, and a superset of the larger is a superset of the smaller. The requested
# connectivity is still what `NativeCheck` filters by.
#
# COVERAGE IS EXHAUSTIVELY VERIFIED, not merely argued: every flush block of
# every size at every origin on every face, levels 1 through 6, both
# connectivities, on all three systems — zero halo cells outside the derived
# rectangles, worst case seven rectangles. Because the seam maps are affine at
# every `n`, levels past 6 add no new structure, only longer runs. The
# differential tests in `test/systems/crosssystem/subtree_halos.jl` re-run the
# same claim through the forced-geometry oracle, which is the only oracle in
# that file that can see a candidate this derivation never proposed.
#
# WHAT FALLS BACK. One configuration only: a system with more faces than
# `_BAND_RECT_CAP`, which none of the three is. Everything else — every flush
# side, every face corner, every whole-face block, both connectivities — is
# walked. The guard is kept because it is the one assumption the derivation
# cannot check itself, and a `small_push` past capacity would be a `BoundsError`
# from inside an iterator rather than an honest fallback.
# ---------------------------------------------------------------------------

# One probe: everything the native one-ring of the rim cell at `(sx, sy)` can
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
    nfaces = ncells(levelgrid(sys, 0))
    nfaces <= _BAND_RECT_CAP ||
        return generic_halo_engine(sys, c, target, connectivity)
    grid = levelgrid(sys, target)
    rects = Helpers.small_push(_empty_band_rects(), home)
    # The two extreme rim cells of every flush side. At most eight probes, at
    # most four distinct cells, and O(1) in the halo's size — construction stays
    # a constant-time act however deep the target.
    if x0 == 0
        rects = _seam_probe(sys, grid, rects, target, face, x0, y0, home)
        rects = _seam_probe(sys, grid, rects, target, face, x0, y0 + side - 1, home)
    end
    if x0 + side == n
        rects = _seam_probe(sys, grid, rects, target, face, x0 + side - 1, y0, home)
        rects = _seam_probe(sys, grid, rects, target, face, x0 + side - 1,
            y0 + side - 1, home)
    end
    if y0 == 0
        rects = _seam_probe(sys, grid, rects, target, face, x0, y0, home)
        rects = _seam_probe(sys, grid, rects, target, face, x0 + side - 1, y0, home)
    end
    if y0 + side == n
        rects = _seam_probe(sys, grid, rects, target, face, x0, y0 + side - 1, home)
        rects = _seam_probe(sys, grid, rects, target, face, x0 + side - 1,
            y0 + side - 1, home)
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
    return SquareBandEngine(curve,
        NativeCheck(sys, grid, c, level(c), connectivity),
        target, n, Int32(face), x0, y0, side, connectivity isa Vertex, ordered)
end
