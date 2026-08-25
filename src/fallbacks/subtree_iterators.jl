# Lazy subtree border and interior iterators, which `border` and `interior` read
# for a region that is a whole rooted subtree. Systems specialize `border_engine` and
# `interior_engine`; iteration state is returned by `iterate`, so engines remain
# restartable.

# ===========================================================================
# The two public types
# ===========================================================================

"""
    EdgeCellIterator(sys, c, l; connectivity = Vertex())

The border of `c`'s subtree at level `l`, lazily: every level-`l` descendant of `c`
with a neighbour that is not one, in ascending canonical order.

The walk `border(subtree(sys, c, l))` reads.
`l == level(c)` yields exactly `c` — a cell is its own border. `l < level(c)` and
`l > maxlevel(sys)` throw an `ArgumentError`, as the eager verb does.

Memory is `O(depth)` and independent of the border's size: IGeo7, H3, HEALPix,
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
        border_engine(sys, c, target, connectivity))
end

"""
    InnerCellIterator(sys, c, l; connectivity = Vertex())

The interior of `c`'s subtree at level `l`, lazily: the level-`l` descendants
that are **not** on the border, in ascending canonical order.

The walk `interior(subtree(sys, c, l))` reads, and
together with [`EdgeCellIterator`](@ref) it partitions
[`descendants`](@ref)`(sys, c, l)`. `l == level(c)` is empty.

The interior is generated from branches pruned by the border walk. The automaton
prunes exactly where a branch goes wholly interior, so the interior is a
disjoint union of complete sub-subtrees, emitted in place. No membership set,
and the border is never materialised to be subtracted.

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

const SubtreeIterator{S,C,K,E} = Union{EdgeCellIterator{S,C,K,E},
                                       InnerCellIterator{S,C,K,E}}

Base.iterate(it::SubtreeIterator) = iterate(it.engine)
Base.iterate(it::SubtreeIterator, state) = iterate(it.engine, state)

# Derive `eltype` from the input cell type so engine dispatch unions do not widen
# collected results to `Vector{Any}`.
Base.eltype(::Type{<:EdgeCellIterator{S,C}}) where {S,C} = C
Base.eltype(::Type{<:InnerCellIterator{S,C}}) where {S,C} = C
Base.IteratorSize(::Type{<:EdgeCellIterator{S,C,K,E}}) where {S,C,K,E} =
    Base.IteratorSize(E)
Base.IteratorSize(::Type{<:InnerCellIterator{S,C,K,E}}) where {S,C,K,E} =
    Base.IteratorSize(E)

# Engines without a constant-time count intentionally provide no `length`.
Base.length(it::SubtreeIterator) = length(it.engine)

# Guard declared lengths during collection; a wrong count must not expose an
# uninitialized tail.
Base.collect(it::SubtreeIterator) = collect_subtree(it)

# Include connectivity because it affects A5 results.
Base.show(io::IO, it::EdgeCellIterator) = print(io, "EdgeCellIterator(",
    it.system, ", ", it.cell, ", ", it.level, "; connectivity = ",
    it.connectivity, ")")
Base.show(io::IO, it::InnerCellIterator) = print(io, "InnerCellIterator(",
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
# The generic scan: lazy over indices, one `ancestor` test per cell
# ===========================================================================

# A neighbour is outside the subtree when its root-level ancestor differs from
# the subtree root.
@inline function _has_outside_neighbor(sys, grid, root, rootlevel, d, connectivity)
    for nb in neighbors(grid, d, 1; connectivity)
        ancestor(sys, nb, rootlevel) != root && return true
    end
    return false
end

"""
    ScanBorderEngine(sys, grid, root, rootlevel, first, last, connectivity)

Walk the subtree's index range and keep the cells with a neighbour outside
it. `O(subtree · degree)` time, `O(1)` memory, and lazy — the border is never
collected. Requires [`has_sorted_subtrees`](@ref).
"""
struct ScanBorderEngine{S,G,C,K}
    system::S
    grid::G
    root::C
    rootlevel::Int
    first::Int
    last::Int
    connectivity::K
end

Base.eltype(::Type{<:ScanBorderEngine{S,G,C,K}}) where {S,G,C,K} = C
Base.IteratorSize(::Type{<:ScanBorderEngine}) = Base.SizeUnknown()

Base.iterate(e::ScanBorderEngine) = iterate(e, e.first)
function Base.iterate(e::ScanBorderEngine, i::Int)
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

[`ScanBorderEngine`](@ref)'s complement, on the same terms.
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

function border_engine(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex,
        target::Int, connectivity::Connectivity)
    lc = level(c)
    target >= lc || throw(ArgumentError(
        "border: level $target is above the cell's own level $lc"))
    C = cellindextype(sys)
    target == lc && return EagerEngine(C[c])
    has_sorted_subtrees(sys) || return EagerEngine(
        _eager_border(sys, c, target, connectivity))
    r = descendant_range(sys, c, target)
    return ScanBorderEngine(sys, levelgrid(sys, target), c, lc,
        first(r), last(r), connectivity)
end

function interior_engine(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex,
        target::Int, connectivity::Connectivity)
    lc = level(c)
    target >= lc || throw(ArgumentError(
        "interior: level $target is above the cell's own level $lc"))
    C = cellindextype(sys)
    target == lc && return EagerEngine(C[])
    has_sorted_subtrees(sys) || return EagerEngine(
        _eager_interior(sys, c, target, connectivity))
    r = descendant_range(sys, c, target)
    return ScanInteriorEngine(sys, levelgrid(sys, target), c, lc,
        first(r), last(r), connectivity)
end

# Systems without `descendant_range` materialize descendants before iterating,
# using O(subtree) rather than O(depth) memory. A5 currently uses this path.
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

# Authalic wrapping does not change the hierarchy or subtree membership.
border_engine(sys::AuthalicSystem, c::AbstractCellIndex, target::Int,
    connectivity::Connectivity) = border_engine(sys.system, c, target, connectivity)
interior_engine(sys::AuthalicSystem, c::AbstractCellIndex, target::Int,
    connectivity::Connectivity) = interior_engine(sys.system, c, target, connectivity)

# ===========================================================================
# The square-block walk, shared by the quad-face family
# ===========================================================================

# A subtree of an aperture-4 system is an aligned `2^d x 2^d` lattice block whose
# descendants occupy the contiguous id block `[id * 4^d, (id+1) * 4^d)`; a
# descendant's offset in it is its position along the system's space-filling
# curve. The border is the block's perimeter, so the walk is: descend the quadtree
# in curve order, carrying which of the four outer sides the sub-square is still
# flush with, and prune a quadrant that is flush with none. Visiting the four
# curve positions in order emits in ascending id by construction.
#
# `quadrant_step` abstracts the Morton and Hilbert curve differences.

const _BORDER_XMIN = 0x1
const _BORDER_XMAX = 0x2
const _BORDER_YMIN = 0x4
const _BORDER_YMAX = 0x8
const _BORDER_ALL = _BORDER_XMIN | _BORDER_XMAX | _BORDER_YMIN | _BORDER_YMAX

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

Dispatches on the curve used by [`SquareBorderEngine`](@ref). A system supplies a
curve type and a method for its quadrant transition.
"""
function quadrant_step end

@inline quadrant_step(::MortonCurve, orientation::UInt8, p::Int) =
    (p & 1, (p >> 1) & 1, orientation)

# A frame stores exposed sides, curve orientation, and the next child. The
# square size and curve code are derived from stack depth. The 32-frame capacity
# exceeds the deepest registered system level (S2 level 30).
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
    (i == 0 ? mask & _BORDER_XMIN : mask & _BORDER_XMAX) |
    (j == 0 ? mask & _BORDER_YMIN : mask & _BORDER_YMAX)

"""
    SquareBorderEngine(curve, base, level, side, orientation)

The perimeter of the `side x side` block whose first cell has raw index `base`,
in ascending id. `4·side - 4` cells for `side > 1`, `O(depth)` memory, `O(border)`
time.

Yields [`LevelIndex`](@ref) because the traversal constructs each dense id as
`base` plus its curve offset. Systems using another id type require a different
engine or a cell-construction parameter.
"""
struct SquareBorderEngine{V}
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
        SquareFrame(_BORDER_ALL, e.orientation, 0x0)), Int64(0))

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

Base.eltype(::Type{<:SquareBorderEngine}) = LevelIndex
Base.IteratorSize(::Type{<:SquareBorderEngine}) = Base.HasLength()
Base.length(e::SquareBorderEngine) = e.side == 1 ? 1 : Int(4 * e.side - 4)

function Base.iterate(e::SquareBorderEngine)
    e.side == 1 &&
        return (LevelIndex(e.level, e.base), SquareWalk(_empty_square_stack(), Int64(0)))
    return iterate(e, _root_walk(e))
end

function Base.iterate(e::SquareBorderEngine, w::SquareWalk)
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
        # A border leaf is emitted without descending, so `code` still describes the
        # frame on top and needs no restoring on the way back in.
        half == 1 && return (LevelIndex(e.level, e.base + child), SquareWalk(st, code))
        st = Helpers.small_push(st, SquareFrame(m, o, 0x0))
        code = child
    end
    return nothing
end

"""
    SquareInteriorEngine(curve, base, level, side, orientation)

The block's interior, as the disjoint union of the sub-squares the border walk
prunes: each is emitted as its contiguous run of ids, so no border membership is
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
        half == 1 && continue                # a border cell
        st = Helpers.small_push(st, SquareFrame(m, o, 0x0))
        code = child
    end
end
