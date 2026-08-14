# Lazy rim and interior of a subtree; `subtree_border` / `subtree_interior` are
# `collect` of these.
#
# One public type per verb. The per-system algorithm lives in an ENGINE the type
# forwards the whole iteration protocol to, so `rim_engine` / `interior_engine`
# are the single extension point a system overrides — one method each, beside
# its `subtree_border`.
#
# Walk state travels in the value `iterate` threads, never in the engine: an
# iterator is restartable, and every state but the materialising fallback's is
# `isbits`, so a partial walk allocates nothing.

# ===========================================================================
# The two public types
# ===========================================================================

"""
    EdgeCellIterator(sys, c, l; connectivity = Vertex())

The rim of `c`'s subtree at level `l`, lazily: every level-`l` descendant of `c`
with a neighbour that is not one, in ascending canonical order.

`collect` of this is [`subtree_border`](@ref), element for element.
`l == level(c)` yields exactly `c` — a cell is its own rim. `l < level(c)` and
`l > max_level(sys)` throw an `ArgumentError`, as the eager verb does.

Memory is `O(depth)` and independent of the rim's size: IGeo7, H3, HEALPix,
ISEA4R and S2 walk a pruned subtree from an `isbits` stack and allocate nothing
per element. `eltype` is `cellindextype(sys)`.

[`Base.IteratorSize`](@ref) is `HasLength()` wherever the count is closed-form —
`3^(d+1)-3` / `5(3^d-1)/2` on the two hexagonal systems, `4·2^d-4` on the three
square ones, `d = l - level(c)` — and `SizeUnknown()` for the generic scan,
which has no count short of running it. There is no `length` that would.

See also [`InnerCellIterator`](@ref), the complement, and
[`descendant_range`](@ref).
"""
struct EdgeCellIterator{S<:AbstractHierarchicalGridSystem,C<:AbstractCellIndex,
        K<:Connectivity,E}
    system::S
    cell::C
    level::Int
    connectivity::K
    engine::E
end

function EdgeCellIterator(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex,
        l::Integer; connectivity::Connectivity = Vertex())
    target = Int(l)
    return EdgeCellIterator(sys, c, target, connectivity,
        rim_engine(sys, c, target, connectivity))
end

"""
    InnerCellIterator(sys, c, l; connectivity = Vertex())

The interior of `c`'s subtree at level `l`, lazily: the level-`l` descendants
that are **not** on the rim, in ascending canonical order.

`collect` of this is [`subtree_interior`](@ref), element for element, and
together with [`EdgeCellIterator`](@ref) it partitions
[`descendants`](@ref)`(sys, c, l)`. `l == level(c)` is empty.

The interior is generated from the rim walk's PRUNED branches — the automaton
prunes exactly where a branch goes wholly interior, so the interior is a
disjoint union of complete sub-subtrees, emitted in place. No membership set,
and the rim is never materialised to be subtracted.

Memory and `eltype` are [`EdgeCellIterator`](@ref)'s. The count is closed-form
on the five systems with an automaton (`(2^d-2)^2` on the square ones), so
`IteratorSize` is `HasLength()` there and `SizeUnknown()` for the generic scan.
"""
struct InnerCellIterator{S<:AbstractHierarchicalGridSystem,C<:AbstractCellIndex,
        K<:Connectivity,E}
    system::S
    cell::C
    level::Int
    connectivity::K
    engine::E
end

function InnerCellIterator(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex,
        l::Integer; connectivity::Connectivity = Vertex())
    target = Int(l)
    return InnerCellIterator(sys, c, target, connectivity,
        interior_engine(sys, c, target, connectivity))
end

"""
    NeighborCellIterator(sys, c, l; connectivity = Vertex())

The **halo** of `c`'s subtree at level `l`, lazily: the level-`l` cells that
are *not* descendants of `c` but have a neighbour that is, in ascending
canonical order, each exactly once. The outside face of the boundary whose
inside face is [`EdgeCellIterator`](@ref).

`collect` of this is [`subtree_halo`](@ref), element for element.
`l == level(c)` yields exactly the cell's own one-ring, sorted. `l < level(c)`
and `l > max_level(sys)` throw an `ArgumentError`, as the eager verb does.

The generic engine **materializes**: halo cells arrive rim-cell by rim-cell
with duplicates and out of order, so the global dedup-and-sort wants the whole
set in hand — `O(rim · degree)` time and memory, the same standing the T20
walkers give A5. HEALPix, S2 and ISEA4R walk the exterior perimeter of their
aligned block directly, in `O(depth)` memory, wherever the block is not flush
with a face edge; see [`SquareHaloEngine`](@ref).

[`Base.IteratorSize`](@ref) is `SizeUnknown` as a contract — face seams, poles
and pentagons break perimeter formulas, so no closed-form count is promised.
An engine that can prove one may declare `HasLength` (the square walk and the
materialised fallbacks do), and `collect` guards it against under-delivery.
"""
struct NeighborCellIterator{S<:AbstractHierarchicalGridSystem,C<:AbstractCellIndex,
        K<:Connectivity,E}
    system::S
    cell::C
    level::Int
    connectivity::K
    engine::E
end

function NeighborCellIterator(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex,
        l::Integer; connectivity::Connectivity = Vertex())
    target = Int(l)
    return NeighborCellIterator(sys, c, target, connectivity,
        neighbor_engine(sys, c, target, connectivity))
end

const SubtreeIterator{S,C,K,E} = Union{EdgeCellIterator{S,C,K,E},
                                       InnerCellIterator{S,C,K,E},
                                       NeighborCellIterator{S,C,K,E}}

Base.iterate(it::SubtreeIterator) = iterate(it.engine)
Base.iterate(it::SubtreeIterator, state) = iterate(it.engine, state)

Base.eltype(::Type{<:EdgeCellIterator{S,C,K,E}}) where {S,C,K,E} = eltype(E)
Base.eltype(::Type{<:InnerCellIterator{S,C,K,E}}) where {S,C,K,E} = eltype(E)
Base.eltype(::Type{<:NeighborCellIterator{S,C,K,E}}) where {S,C,K,E} = eltype(E)
Base.IteratorSize(::Type{<:EdgeCellIterator{S,C,K,E}}) where {S,C,K,E} =
    Base.IteratorSize(E)
Base.IteratorSize(::Type{<:InnerCellIterator{S,C,K,E}}) where {S,C,K,E} =
    Base.IteratorSize(E)
Base.IteratorSize(::Type{<:NeighborCellIterator{S,C,K,E}}) where {S,C,K,E} =
    Base.IteratorSize(E)

# Deliberately delegated rather than defined: an engine that cannot count
# without walking defines no `length`, and the `MethodError` from here is the
# contract ("no `length` that silently costs a full traversal") being kept.
Base.length(it::SubtreeIterator) = length(it.engine)

# The docstrings advertise `collect` as the eager verbs, so `collect` must BE
# the guarded path, not merely parallel to it: `collect`'s own `HasLength` route
# would size its vector from a miscounting automaton and hand back an `undef`
# tail as cell ids. See `collect_subtree`.
Base.collect(it::SubtreeIterator) = collect_subtree(it)

# Connectivity is shown even though it changes nothing on five of the six
# systems: on A5 it changes the answer, and that is exactly when someone is
# reading this.
Base.show(io::IO, it::EdgeCellIterator) = print(io, "EdgeCellIterator(",
    it.system, ", ", it.cell, ", ", it.level, "; connectivity = ",
    it.connectivity, ")")
Base.show(io::IO, it::InnerCellIterator) = print(io, "InnerCellIterator(",
    it.system, ", ", it.cell, ", ", it.level, "; connectivity = ",
    it.connectivity, ")")
Base.show(io::IO, it::NeighborCellIterator) = print(io, "NeighborCellIterator(",
    it.system, ", ", it.cell, ", ", it.level, "; connectivity = ",
    it.connectivity, ")")

# ===========================================================================
# The materialising engine
# ===========================================================================

"""
    EagerEngine(cells)

An engine over an already-built vector: the depth-0 answers, and the fallback
for a system that cannot walk lazily.
"""
struct EagerEngine{C}
    cells::Vector{C}
end

Base.iterate(e::EagerEngine) = iterate(e.cells)
Base.iterate(e::EagerEngine, state) = iterate(e.cells, state)
Base.eltype(::Type{EagerEngine{C}}) where {C} = C
Base.IteratorSize(::Type{<:EagerEngine}) = Base.HasLength()
Base.length(e::EagerEngine) = length(e.cells)

# ===========================================================================
# The generic scan: lazy over positions, one `ancestor` test per cell
# ===========================================================================

# A neighbour is outside the subtree exactly when its ancestor at the root's
# level is not the root — the same test the eager fallback used, and still
# cheaper than a membership set over the subtree.
@inline function _has_outside_neighbor(sys, grid, root, rootlevel, d, connectivity)
    for nb in neighbors(grid, d, 1; connectivity)
        ancestor(sys, nb, rootlevel) != root && return true
    end
    return false
end

"""
    ScanRimEngine(sys, grid, root, rootlevel, first, last, connectivity)

Walk the subtree's position range and keep the cells with a neighbour outside
it. `O(subtree · degree)` time, `O(1)` memory, and lazy — the rim is never
collected. Requires [`has_sorted_subtrees`](@ref).
"""
struct ScanRimEngine{S,G,C,K}
    system::S
    grid::G
    root::C
    rootlevel::Int
    first::Int
    last::Int
    connectivity::K
end

Base.eltype(::Type{<:ScanRimEngine{S,G,C,K}}) where {S,G,C,K} = C
Base.IteratorSize(::Type{<:ScanRimEngine}) = Base.SizeUnknown()

Base.iterate(e::ScanRimEngine) = iterate(e, e.first)
function Base.iterate(e::ScanRimEngine, i::Int)
    while i <= e.last
        d = cellindex(e.grid, i)
        _has_outside_neighbor(e.system, e.grid, e.root, e.rootlevel, d,
            e.connectivity) && return (d, i + 1)
        i += 1
    end
    return nothing
end

"""
    ScanInteriorEngine(sys, grid, root, rootlevel, first, last, connectivity)

[`ScanRimEngine`](@ref)'s complement, on the same terms.
"""
struct ScanInteriorEngine{S,G,C,K}
    system::S
    grid::G
    root::C
    rootlevel::Int
    first::Int
    last::Int
    connectivity::K
end

Base.eltype(::Type{<:ScanInteriorEngine{S,G,C,K}}) where {S,G,C,K} = C
Base.IteratorSize(::Type{<:ScanInteriorEngine}) = Base.SizeUnknown()

Base.iterate(e::ScanInteriorEngine) = iterate(e, e.first)
function Base.iterate(e::ScanInteriorEngine, i::Int)
    while i <= e.last
        d = cellindex(e.grid, i)
        _has_outside_neighbor(e.system, e.grid, e.root, e.rootlevel, d,
            e.connectivity) || return (d, i + 1)
        i += 1
    end
    return nothing
end

# ===========================================================================
# Generic engine construction
# ===========================================================================

function rim_engine(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex,
        target::Int, connectivity::Connectivity)
    lc = level(c)
    target >= lc || throw(ArgumentError(
        "subtree_border: level $target is above the cell's own level $lc"))
    C = cellindextype(sys)
    target == lc && return EagerEngine(C[c])
    has_sorted_subtrees(sys) || return EagerEngine(
        _eager_border(sys, c, target, connectivity))
    r = descendant_range(sys, c, target)
    return ScanRimEngine(sys, levelgrid(sys, target), c, lc,
        first(r), last(r), connectivity)
end

function interior_engine(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex,
        target::Int, connectivity::Connectivity)
    lc = level(c)
    target >= lc || throw(ArgumentError(
        "subtree_interior: level $target is above the cell's own level $lc"))
    C = cellindextype(sys)
    target == lc && return EagerEngine(C[])
    has_sorted_subtrees(sys) || return EagerEngine(
        _eager_interior(sys, c, target, connectivity))
    r = descendant_range(sys, c, target)
    return ScanInteriorEngine(sys, levelgrid(sys, target), c, lc,
        first(r), last(r), connectivity)
end

# The halo's generic engine materializes on EVERY system, not just A5: the
# candidates arrive rim-cell by rim-cell, duplicated wherever two rim cells
# share an outside neighbour and ordered by which rim cell found them, so the
# global dedup-and-sort needs the whole set in hand. `O(rim · degree)` time and
# memory — the documented fallback, not a hidden one. The aligned-block systems
# override with `SquareHaloEngine` below for non-flush blocks; IGeo7 and H3
# stay here deliberately, since an aperture-7 subtree boundary is a fractal
# with no perimeter arithmetic to walk, and A5 stays for T20's reason.
neighbor_engine(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex,
    target::Int, connectivity::Connectivity) =
    generic_neighbor_engine(sys, c, target, connectivity)

# Named rather than folded into the method above so a system override has a
# fallback to CALL, not just to shadow: the square systems reach it for the
# blocks their automaton must decline (flush with a face edge, or depth 0).
function generic_neighbor_engine(sys::AbstractHierarchicalGridSystem,
        c::AbstractCellIndex, target::Int, connectivity::Connectivity)
    lc = level(c)
    target >= lc || throw(ArgumentError(
        "subtree_halo: level $target is above the cell's own level $lc"))
    target == lc && return EagerEngine(
        sort!(cellindextype(sys)[nb for nb in
            neighbors(levelgrid(sys, target), c, 1; connectivity)]))
    return EagerEngine(_eager_halo(sys, c, target, connectivity))
end

# A5 alone lands here: its canonical order establishes no `descendant_range`, so
# there is no position interval to walk and no closed-form child adjacency to
# build an automaton from. The subtree is materialised internally — the iterator
# is still an iterator, but its memory is `O(subtree)`, not `O(depth)`.
function _eager_border(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex,
        target::Int, connectivity::Connectivity)
    lc = level(c)
    grid = levelgrid(sys, target)
    out = cellindextype(sys)[]
    for d in descendants(sys, c, target)
        _has_outside_neighbor(sys, grid, c, lc, d, connectivity) && push!(out, d)
    end
    return out
end

function _eager_interior(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex,
        target::Int, connectivity::Connectivity)
    lc = level(c)
    grid = levelgrid(sys, target)
    out = cellindextype(sys)[]
    for d in descendants(sys, c, target)
        _has_outside_neighbor(sys, grid, c, lc, d, connectivity) || push!(out, d)
    end
    return out
end

# The defining law run forwards: every outside neighbour of a rim cell is halo,
# and every halo cell is an outside neighbour of some rim cell (its in-subtree
# neighbour is rim by definition). The rim arrives through the system's own
# `EdgeCellIterator`, so a system with an `O(rim)` rim walk pays `O(rim)` here
# too — only the dedup-and-sort is generic. Level validation is the rim
# engine's, which is why `target > max_level` need not be re-checked above.
function _eager_halo(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex,
        target::Int, connectivity::Connectivity)
    lc = level(c)
    grid = levelgrid(sys, target)
    out = cellindextype(sys)[]
    for d in EdgeCellIterator(sys, c, target; connectivity)
        for nb in neighbors(grid, d, 1; connectivity)
            ancestor(sys, nb, lc) == c || push!(out, nb)
        end
    end
    return unique!(sort!(out))
end

# The wrapper is transparent to all three walks: the rim and halo of a subtree
# are questions about the discrete hierarchy, which the ellipsoid does not touch.
rim_engine(sys::AuthalicSystem, c::AbstractCellIndex, target::Int,
    connectivity::Connectivity) = rim_engine(sys.system, c, target, connectivity)
interior_engine(sys::AuthalicSystem, c::AbstractCellIndex, target::Int,
    connectivity::Connectivity) = interior_engine(sys.system, c, target, connectivity)
neighbor_engine(sys::AuthalicSystem, c::AbstractCellIndex, target::Int,
    connectivity::Connectivity) = neighbor_engine(sys.system, c, target, connectivity)

# ===========================================================================
# The square-block walk, shared by the three aperture-4 systems
# ===========================================================================

# A subtree of an aperture-4 system is an aligned `2^d x 2^d` lattice block whose
# descendants occupy the contiguous id block `[id * 4^d, (id+1) * 4^d)`; a
# descendant's offset in it is its position along the system's space-filling
# curve. The rim is the block's perimeter, so the walk is: descend the quadtree
# in curve order, carrying which of the four outer sides the sub-square is still
# flush with, and prune a quadrant that is flush with none. Visiting the four
# curve positions in order emits in ascending id by construction.
#
# The curve enters only through `quadrant_step`, which is what lets HEALPix and
# ISEA4R (Morton) and S2 (Hilbert) share one walker.

const _RIM_XMIN = 0x1
const _RIM_XMAX = 0x2
const _RIM_YMIN = 0x4
const _RIM_YMAX = 0x8
const _RIM_ALL = _RIM_XMIN | _RIM_XMAX | _RIM_YMIN | _RIM_YMAX

"""
    MortonCurve()

Z-order: curve position `p`'s low bit is the `x` half and its high bit the `y`
half, at every level and in every orientation.
"""
struct MortonCurve end

"""
    quadrant_step(curve, orientation, p) -> (i, j, child_orientation)

Which half of its parent square in `x` (`i`) and `y` (`j`) the sub-square at
curve position `p` is, and the orientation state its own children are read
under. Orientation is inert for [`MortonCurve`](@ref); S2's Hilbert curve
advances it.

Dispatches on the CURVE, not on the system — it is the one parameter of
[`SquareRimEngine`](@ref), not a system extension point, which is why it lives
here rather than beside `rim_engine` in `src/interface/`. A system supplies a
curve type and one method of this; three systems currently supply two curves.
"""
function quadrant_step end

@inline quadrant_step(::MortonCurve, orientation::UInt8, p::Int) =
    (p & 1, (p >> 1) & 1, orientation)

# A frame carries only what its parent cannot recompute: the sides it is flush
# with, its curve orientation, and the child position to try next (retired at
# 4). Its size is `side >> (depth - 1)` and its code is threaded alongside the
# stack, both restored on pop by one shift and one subtraction — which is what
# keeps the frame three bytes, and the whole stack a register-sized inline
# value. 32 frames is past every registered system's depth (S2's 30 is deepest).
struct SquareFrame
    mask::UInt8
    orientation::UInt8
    next::UInt8
end

const _SQUARE_CAP = 32
const SquareStack = Helpers.SmallList{_SQUARE_CAP,SquareFrame}

@inline _empty_square_stack() =
    Helpers.empty_small_list(Val(_SQUARE_CAP), SquareFrame(0x0, 0x0, 0x0))

@inline _bump(f::SquareFrame) =
    SquareFrame(f.mask, f.orientation, f.next + 0x1)

# Which outer sides a quadrant inherits: the low-`x` half can be exposed at the
# block's `x` minimum but never at its maximum, and so on. Zero means every side
# of it is interior, which is the pruning rule.
@inline _child_mask(mask::UInt8, i::Int, j::Int) =
    (i == 0 ? mask & _RIM_XMIN : mask & _RIM_XMAX) |
    (j == 0 ? mask & _RIM_YMIN : mask & _RIM_YMAX)

"""
    SquareRimEngine(curve, base, level, side, orientation)

The perimeter of the `side x side` block whose first cell has raw index `base`,
in ascending id. `4·side - 4` cells for `side > 1`, `O(depth)` memory, `O(rim)`
time.

Yields [`LevelIndex`](@ref) rather than a parameterized cell type, and not by
oversight: the walk's premise is that a descendant's id is `base` plus its
offset along the curve, which is a statement about a dense integer id, and
`LevelIndex` is what wraps one. All three aperture-4 systems use it. A fourth
that did not would need a cell-construction parameter here, not just an
`eltype`.
"""
struct SquareRimEngine{V}
    curve::V
    base::Int64
    level::Int
    side::Int64
    orientation::UInt8
end

"""
    SquareWalk(stack, code)

The walk state: the frame stack, and the within-subtree code of the sub-square
on top of it. Frame `k` has side `side >> (k - 1)`, so neither size nor code
needs storing per frame.
"""
struct SquareWalk
    stack::SquareStack
    code::Int64
end

@inline _root_walk(e) = SquareWalk(
    Helpers.small_push(_empty_square_stack(),
        SquareFrame(_RIM_ALL, e.orientation, 0x0)), Int64(0))

# Undo the push that entered the frame just dropped: its parent is now on top
# and has already advanced past the position it used, so that position is
# `next - 1` and the sub-square it entered had side `side >> k`.
@inline function _restore_code(e, st::SquareStack, code::Int64)
    k = length(st)
    k == 0 && return code
    p = Int((@inbounds st[k]).next) - 1
    half = e.side >> k
    return code - p * half * half
end

Base.eltype(::Type{<:SquareRimEngine}) = LevelIndex
Base.IteratorSize(::Type{<:SquareRimEngine}) = Base.HasLength()
Base.length(e::SquareRimEngine) = e.side == 1 ? 1 : Int(4 * e.side - 4)

function Base.iterate(e::SquareRimEngine)
    e.side == 1 &&
        return (LevelIndex(e.level, e.base), SquareWalk(_empty_square_stack(), Int64(0)))
    return iterate(e, _root_walk(e))
end

function Base.iterate(e::SquareRimEngine, w::SquareWalk)
    st = w.stack
    code = w.code
    while !isempty(st)
        k = length(st)
        f = @inbounds st[k]
        if f.next > 0x3
            st = Helpers.small_pop(st)
            code = _restore_code(e, st, code)
            continue
        end
        p = Int(f.next)
        st = Helpers.small_setlast(st, _bump(f))
        i, j, o = quadrant_step(e.curve, f.orientation, p)
        m = _child_mask(f.mask, i, j)
        m == 0x0 && continue
        half = e.side >> k
        child = code + p * half * half
        # A rim leaf is emitted without descending, so `code` still describes the
        # frame on top and needs no restoring on the way back in.
        half == 1 && return (LevelIndex(e.level, e.base + child), SquareWalk(st, code))
        st = Helpers.small_push(st, SquareFrame(m, o, 0x0))
        code = child
    end
    return nothing
end

"""
    SquareInteriorEngine(curve, base, level, side, orientation)

The block's interior, as the disjoint union of the sub-squares the rim walk
prunes: each is emitted as its contiguous run of ids, so no rim membership is
ever tested or stored. `(side - 2)^2` cells, `O(depth)` memory.
"""
struct SquareInteriorEngine{V}
    curve::V
    base::Int64
    level::Int
    side::Int64
    orientation::UInt8
end

# A pruned sub-square contributes the id run `[next, stop)`; the walk is where
# it resumes once that run is drained.
struct SquareInteriorState
    walk::SquareWalk
    next::Int64
    stop::Int64
end

Base.eltype(::Type{<:SquareInteriorEngine}) = LevelIndex
Base.IteratorSize(::Type{<:SquareInteriorEngine}) = Base.HasLength()
Base.length(e::SquareInteriorEngine) =
    e.side <= 2 ? 0 : Int((e.side - 2) * (e.side - 2))

function Base.iterate(e::SquareInteriorEngine)
    e.side == 1 && return nothing
    return iterate(e, SquareInteriorState(_root_walk(e), Int64(0), Int64(0)))
end

function Base.iterate(e::SquareInteriorEngine, s::SquareInteriorState)
    st = s.walk.stack
    code = s.walk.code
    nxt = s.next
    stop = s.stop
    while true
        if nxt < stop
            return (LevelIndex(e.level, e.base + nxt),
                SquareInteriorState(SquareWalk(st, code), nxt + Int64(1), stop))
        end
        isempty(st) && return nothing
        k = length(st)
        f = @inbounds st[k]
        if f.next > 0x3
            st = Helpers.small_pop(st)
            code = _restore_code(e, st, code)
            continue
        end
        p = Int(f.next)
        st = Helpers.small_setlast(st, _bump(f))
        i, j, o = quadrant_step(e.curve, f.orientation, p)
        half = e.side >> k
        child = code + p * half * half
        m = _child_mask(f.mask, i, j)
        if m == 0x0
            nxt = child                      # a whole interior sub-square
            stop = child + half * half
            continue
        end
        half == 1 && continue                # a rim cell
        st = Helpers.small_push(st, SquareFrame(m, o, 0x0))
        code = child
    end
end

# ===========================================================================
# The exterior-perimeter walk, shared by the same three systems
# ===========================================================================

# The halo of an aligned block that is nowhere flush with its face's edge is
# the width-1 band around it — every band cell is on the block's own face, and
# in-face adjacency is the plain 3×3 lattice on all three systems, so band and
# halo coincide, with the four diagonal-contact corners dropped under `Edge()`.
#
# Unlike the rim, the band is not a curve interval of any one subtree, so there
# is no `base + offset` walk to run. Instead the engine descends the whole
# FACE's quadtree in curve order and prunes every quadrant that misses the
# band: curve order at each level is ascending id by construction, so whatever
# survives arrives ascending with no sort and no dedup. A quadrant of side `h`
# survives only if it intersects the band, and at most `O(perimeter / h + 1)`
# of each size do, so the walk visits `O(halo + depth)` nodes in `O(depth)`
# memory — the descent tracks the running curve code and lattice origin, and
# restores both on pop by re-reading the parent's last step.
#
# A block flush with a face edge has halo cells across the seam, on other
# faces, whose ids no in-face arithmetic can name — those blocks take the
# materialising generic engine instead, which is why this type needs no seam
# machinery and no system tables, only the curve.

# `mask` would be dead weight here — pruning is by lattice overlap, not by
# flush sides — so the frame is the two bytes the descent actually needs.
struct HaloFrame
    orientation::UInt8
    next::UInt8
end

const HaloStack = Helpers.SmallList{_SQUARE_CAP,HaloFrame}

@inline _empty_halo_stack() =
    Helpers.empty_small_list(Val(_SQUARE_CAP), HaloFrame(0x0, 0x0))

"""
    SquareHaloEngine(curve, facebase, level, faceside, orientation, x0, y0, side, corners)

The halo of the `side x side` block at lattice origin `(x0, y0)` on the face
whose first cell has raw index `facebase`, in ascending id: `4·side + 4` band
cells, or `4·side` with `corners = false` (`Edge()` drops the four
diagonal-contact corners). Valid only for a block not flush with the face —
`1 <= x0` and `x0 + side <= faceside - 1`, likewise `y0` — which the
constructing system checks; `O(halo + depth)` time, `O(depth)` memory.

`faceside` is the face's full lattice side at `level` and `orientation` the
curve state its root is read under: the walk descends the face, not the block,
because the band cells are outside the block's own curve interval. Yields
[`LevelIndex`](@ref) on [`SquareRimEngine`](@ref)'s reasoning, and takes the
same [`quadrant_step`](@ref) curves.
"""
struct SquareHaloEngine{V}
    curve::V
    facebase::Int64
    level::Int
    faceside::Int64
    orientation::UInt8
    x0::Int64
    y0::Int64
    side::Int64
    corners::Bool
end

"""
    HaloWalk(stack, code, x, y)

The walk state: the frame stack, and the within-face curve code and lattice
origin of the sub-square on top of it. Frame `k`'s node has side
`faceside >> (k - 1)`, so none of the three needs storing per frame — all are
restored on pop by replaying the parent's last `quadrant_step`.
"""
struct HaloWalk
    stack::HaloStack
    code::Int64
    x::Int64
    y::Int64
end

# Does a node overlap the band's bounding square at all, and is it swallowed by
# the block? A surviving node overlaps the outer square without sitting inside
# the block, which forces it to overlap the band itself — the outer square
# minus the block IS the band.
@inline _halo_overlaps(e::SquareHaloEngine, x::Int64, y::Int64, h::Int64) =
    x <= e.x0 + e.side && x + h - 1 >= e.x0 - 1 &&
    y <= e.y0 + e.side && y + h - 1 >= e.y0 - 1

@inline _inside_block(e::SquareHaloEngine, x::Int64, y::Int64, h::Int64) =
    x >= e.x0 && x + h - 1 <= e.x0 + e.side - 1 &&
    y >= e.y0 && y + h - 1 <= e.y0 + e.side - 1

# The four cells that touch the block at a vertex only, halo under `Vertex()`
# and not under `Edge()`.
@inline _band_corner(e::SquareHaloEngine, x::Int64, y::Int64) =
    (x == e.x0 - 1 || x == e.x0 + e.side) &&
    (y == e.y0 - 1 || y == e.y0 + e.side)

Base.eltype(::Type{<:SquareHaloEngine}) = LevelIndex
Base.IteratorSize(::Type{<:SquareHaloEngine}) = Base.HasLength()
Base.length(e::SquareHaloEngine) =
    e.corners ? Int(4 * e.side + 4) : Int(4 * e.side)

Base.iterate(e::SquareHaloEngine) = iterate(e, HaloWalk(
    Helpers.small_push(_empty_halo_stack(), HaloFrame(e.orientation, 0x0)),
    Int64(0), Int64(0), Int64(0)))

function Base.iterate(e::SquareHaloEngine, w::HaloWalk)
    st = w.stack
    code = w.code
    x = w.x
    y = w.y
    while !isempty(st)
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
        st = Helpers.small_setlast(st, HaloFrame(f.orientation, f.next + 0x1))
        i, j, o = quadrant_step(e.curve, f.orientation, p)
        half = e.faceside >> k
        cx = x + i * half
        cy = y + j * half
        _halo_overlaps(e, cx, cy, half) || continue
        _inside_block(e, cx, cy, half) && continue
        child = code + p * half * half
        if half == 1
            # A surviving leaf is a band cell; only the corner rule is left.
            !e.corners && _band_corner(e, cx, cy) && continue
            # Emitted without descending, so `code`/`x`/`y` still describe the
            # frame on top and need no restoring on the way back in.
            return (LevelIndex(e.level, e.facebase + child), HaloWalk(st, code, x, y))
        end
        st = Helpers.small_push(st, HaloFrame(o, 0x0))
        code = child
        x = cx
        y = cy
    end
    return nothing
end
