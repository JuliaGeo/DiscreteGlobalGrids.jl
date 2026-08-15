# The subtree rim is generated in O(3^d) without enumerating its O(7^d)
# interior. A rim cell is exposed along a contiguous arc of the six
# lattice directions, and the arc a child inherits is a function of the parent's
# arc and which digit the child is. The state is `(L, s)`: the arc of exposed
# directions `s, s+1, ..., s+L-1` (mod 6), in `_H3_DIGIT_DIR` positions.
# `L == 6` is the subtree root, fully exposed and the one state with no arc
# ends; `L == 0` means interior, and the walk prunes there. Digit 0 (the centre
# child) is never on the rim.
#
# All pentagon base cells delete digit 1 (the K axis). The pentagon flag drops
# after the root because a rim suffix contains no zero digit.
# The transition table is fitted and parity-dependent. Its parity roles are
# swapped relative to IGeo7: the branch H3 takes
# at an even child level is the one IGeo7 takes at an odd one.
#
# The walk is a resumable iterator: an explicit frame stack, one inline value, so
# a partial rim costs nothing beyond `O(depth)`. The interior comes off the same
# walk — where the automaton PRUNES, the whole branch below is interior, so
# `H3InteriorEngine` descends it in full instead of dropping it, and never needs
# a rim membership set.

# Extended, not shadowed: this file defines H3's method on the package generic,
# so `H3.subtree_border` and `DiscreteGlobalGrids.subtree_border` are one
# function and generic code reaches the automaton without knowing it exists.
import ..DiscreteGlobalGrids: subtree_border

"""
    _H3_DIGIT_DIR[digit] -> Int

Position `0:5` on the direction ring of the IJK digit `digit in 1:6`, i.e. the
inverse of the geometric cyclic order `K(1), JK(3), J(2), IJ(6), I(4), IK(5)`.
"""
const _H3_DIGIT_DIR = (0, 2, 1, 4, 5, 3)

"""
    _h3_border_step(state, digit, level) -> NTuple{2,Int}

The state a `digit` child at absolute resolution `level` inherits from a cell in
`state`, or `(0, 0)` when that child is interior to the subtree.
"""
@inline function _h3_border_step(state::NTuple{2,Int}, digit::Int, level::Int)
    L, s = state
    (L == 0 || digit == 0) && return (0, 0)
    t = @inbounds _H3_DIGIT_DIR[digit]
    o = mod(t - s, 6)
    o < L || return (0, 0)
    if iseven(level)                         # child level is Class II
        L < 6 && o == 0 && return (2, mod(t + 1, 6))
        L < 6 && o == L - 1 && return (3, mod(t - 1, 6))
        return (4, mod(t - 1, 6))
    else                                     # child level is Class III (mirror)
        L < 6 && o == 0 && return (3, mod(t - 1, 6))
        L < 6 && o == L - 1 && return (2, mod(t - 2, 6))
        return (4, mod(t - 2, 6))
    end
end

# `next` is the digit to try when this frame is next on top; digits run `1:6` on
# the rim walk and `0:6` below a pruned branch, so a frame retires at `next > 6`.
# `full` marks a branch the automaton pruned: everything under it is interior.
# Resolution is not stored — frame `k` sits at `root + k - 1`.
struct H3Frame
    L::Int8
    s::Int8
    next::Int8
    pentagon::Bool
    full::Bool
end

const _H3_STACK_CAP = MAX_RESOLUTION + 1
const H3Stack = DGG.Helpers.SmallList{_H3_STACK_CAP,H3Frame}

@inline _h3_empty_stack() = DGG.Helpers.empty_small_list(Val(_H3_STACK_CAP),
    H3Frame(Int8(0), Int8(0), Int8(0), false, false))

"""
    H3Walk(z, stack)

The walk state. `z` carries the digits of the current path — each descent
overwrites one digit slot, so the parent's id needs no restoring — and already
holds `target` in its resolution field, which is what makes every leaf come out
fully formed and in ascending order by construction.
"""
struct H3Walk
    z::UInt64
    stack::H3Stack
end

# The rim walk starts at digit 1: digit 0 is the centre child, which the
# automaton would prune anyway. The interior walk must start at 0, because that
# centre child is exactly the branch it is looking for.
@inline function _h3_root_walk(e, start::Int8)
    f = H3Frame(Int8(6), Int8(0), start, e.pentagon, false)
    return H3Walk(e.z, DGG.Helpers.small_push(_h3_empty_stack(), f))
end

@inline _h3_bump(f::H3Frame) =
    H3Frame(f.L, f.s, f.next + Int8(1), f.pentagon, f.full)

# Descendants of a cell `d` levels down: the all-zero chain plus five branches
# under a pentagon, a full `7^d` otherwise.
@inline _h3_subtree_count(pentagon::Bool, d::Int) =
    pentagon ? 1 + 5 * (7^d - 1) ÷ 6 : 7^d

@inline _h3_rim_count(pentagon::Bool, d::Int) =
    d == 0 ? 1 : (pentagon ? (5 * (3^d - 1)) ÷ 2 : 3 * 3^d - 3)

"""
    H3RimEngine(z, res, target, pentagon)

The rim automaton as a resumable iterator over `H3Cell`s. `z` is the subtree
root with its resolution field already moved to `target`.
"""
struct H3RimEngine
    z::UInt64
    res::Int
    target::Int
    pentagon::Bool
end

Base.eltype(::Type{H3RimEngine}) = H3Cell
Base.IteratorSize(::Type{H3RimEngine}) = Base.HasLength()
Base.length(e::H3RimEngine) = _h3_rim_count(e.pentagon, e.target - e.res)

function Base.iterate(e::H3RimEngine)
    e.target == e.res && return (H3Cell(e.z), H3Walk(e.z, _h3_empty_stack()))
    return iterate(e, _h3_root_walk(e, Int8(1)))
end

Base.iterate(e::H3RimEngine, w::H3Walk) = _h3_rim_advance(e.res, e.target, w)

# The walk itself, taking the two numbers it actually reads rather than an
# engine: the seeded arc engine below runs the same automaton from a different
# root frame, and this is the whole of what they share. Every frame past the root
# is built here, so a seeded entry differs from a rim entry in exactly one
# `H3Frame` and in nothing else.
function _h3_rim_advance(res0::Int, target::Int, w::H3Walk)
    z = w.z
    st = w.stack
    while !isempty(st)
        k = length(st)
        f = @inbounds st[k]
        if f.next > Int8(6)
            st = DGG.Helpers.small_pop(st)
            continue
        end
        digit = Int(f.next)
        st = DGG.Helpers.small_setlast(st, _h3_bump(f))
        f.pentagon && digit == 1 && continue   # the K axis, deleted under a pentagon
        res = res0 + k - 1
        child = _h3_border_step((Int(f.L), Int(f.s)), digit, res + 1)
        child[1] == 0 && continue
        shift = _h3_digit_shift(res + 1)
        z = (z & ~(UInt64(7) << shift)) | (UInt64(digit) << shift)
        res + 1 == target && return (H3Cell(z), H3Walk(z, st))
        st = DGG.Helpers.small_push(st,
            H3Frame(Int8(child[1]), Int8(child[2]), Int8(1), false, false))
    end
    return nothing
end

"""
    H3ArcEngine(z, res, target, pentagon, L, s)

The same automaton entered at an arbitrary arc: the level-`target` descendants of
the cell `z` reachable along the exposed directions `s, s+1, …, s+L-1 (mod 6)`,
ascending, in `O(depth)` memory. `H3RimEngine` is this with `(L, s) == (6, 0)`,
which is the one state with no arc ends and therefore the one the census
formulas describe — so this engine declares `SizeUnknown()` and no `length`.

`pentagon` is the SEED CELL's own flag, not a subtree root's. A calibrated arc is
seeded at a neighbour, which may itself be a pentagon, and dropping the flag
would walk the deleted K axis and yield ids for cells that do not exist. It drops
below the root frame for the usual reason: a rim suffix contains no zero digit,
so no descendant reached from here is a pentagon.
"""
struct H3ArcEngine
    z::UInt64
    res::Int
    target::Int
    pentagon::Bool
    L::Int8
    s::Int8
end

Base.eltype(::Type{H3ArcEngine}) = H3Cell
Base.IteratorSize(::Type{H3ArcEngine}) = Base.SizeUnknown()

Base.iterate(e::H3ArcEngine) = iterate(e, H3Walk(e.z,
    DGG.Helpers.small_push(_h3_empty_stack(),
        H3Frame(e.L, e.s, Int8(1), e.pentagon, false))))

Base.iterate(e::H3ArcEngine, w::H3Walk) = _h3_rim_advance(e.res, e.target, w)

"""
    H3InteriorEngine(z, res, target, pentagon)

The complement, off the same automaton: a branch the rim walk prunes is wholly
interior, so this one descends it in full (digits `0:6`) instead of dropping it.
No rim membership is ever tested or stored.
"""
struct H3InteriorEngine
    z::UInt64
    res::Int
    target::Int
    pentagon::Bool
end

Base.eltype(::Type{H3InteriorEngine}) = H3Cell
Base.IteratorSize(::Type{H3InteriorEngine}) = Base.HasLength()
function Base.length(e::H3InteriorEngine)
    d = e.target - e.res
    d == 0 && return 0
    return _h3_subtree_count(e.pentagon, d) - _h3_rim_count(e.pentagon, d)
end

function Base.iterate(e::H3InteriorEngine)
    e.target == e.res && return nothing
    return iterate(e, _h3_root_walk(e, Int8(0)))
end

function Base.iterate(e::H3InteriorEngine, w::H3Walk)
    z = w.z
    st = w.stack
    while !isempty(st)
        k = length(st)
        f = @inbounds st[k]
        if f.next > Int8(6)
            st = DGG.Helpers.small_pop(st)
            continue
        end
        digit = Int(f.next)
        st = DGG.Helpers.small_setlast(st, _h3_bump(f))
        f.pentagon && digit == 1 && continue
        res = e.res + k - 1
        # Below a branch the automaton pruned, every cell is interior; on the rim
        # path, the child it prunes is where the interior starts.
        child = f.full ? (0, 0) : _h3_border_step((Int(f.L), Int(f.s)), digit, res + 1)
        shift = _h3_digit_shift(res + 1)
        z = (z & ~(UInt64(7) << shift)) | (UInt64(digit) << shift)
        if child[1] != 0
            res + 1 == e.target && continue      # a rim cell: not ours
            st = DGG.Helpers.small_push(st,
                H3Frame(Int8(child[1]), Int8(child[2]), Int8(0), false, false))
            continue
        end
        res + 1 == e.target && return (H3Cell(z), H3Walk(z, st))
        # A cell stays a pentagon exactly while its digits are all zero, which is
        # reachable here — unlike on the rim walk — because digit 0 is taken.
        st = DGG.Helpers.small_push(st,
            H3Frame(Int8(0), Int8(0), Int8(0), f.pentagon && digit == 0, true))
    end
    return nothing
end

function DGG.rim_engine(sys::H3System, c::H3Cell, target::Int,
        connectivity::Connectivity)
    lvl, pentagon = _h3_border_checked(c, target)
    return H3RimEngine(_h3_with_resolution(c.id, target), lvl, target, pentagon)
end

function DGG.interior_engine(sys::H3System, c::H3Cell, target::Int,
        connectivity::Connectivity)
    lvl, pentagon = _h3_border_checked(c, target)
    return H3InteriorEngine(_h3_with_resolution(c.id, target), lvl, target, pentagon)
end

# ---------------------------------------------------------------------------
# The two lines the shared calibrated halo walk needs (`hex_halo_engine`)
# ---------------------------------------------------------------------------

# The ring position of the step from a cell's parent to the cell, read off the
# cell's own last digit through the same table `_h3_border_step` uses. Digit 0 is
# the centre child, which has no direction, and a base cell has no parent.
function DGG.hex_child_direction(::H3System, c::H3Cell)
    res = level(c)
    res == 0 && return -1
    digit = Int((c.id >> _h3_digit_shift(res)) & UInt64(7))
    digit == 0 && return -1
    return @inbounds _H3_DIGIT_DIR[digit]
end

# Unvalidated on purpose: `hex_halo_engine` owns the level guard and only ever
# passes cells that came out of `neighbors`. `_h3_border_checked` is the entry
# point for the public verbs, which do not know that.
DGG.seeded_rim_engine(::H3System, c::H3Cell, target::Int, arclen::Int,
        start::Int) =
    H3ArcEngine(_h3_with_resolution(c.id, target), level(c), target,
        H3Native.is_pentagon(c.id), Int8(arclen), Int8(start))

# The halo is approached from the neighbouring subtrees rather than from the
# root's own, because a subtree's halo is not an interval of anything H3 can
# name. See `hex_halo_engine` for the calibration, the containment argument, and
# the guards that send a case back to the generic walk.
DGG.halo_engine(sys::H3System, c::H3Cell, target::Int,
    connectivity::Connectivity) =
    DGG.hex_halo_engine(sys, c, target, connectivity)

# libh3 validates neither `cellToChildren` nor `cellToChildrenSize`, so a
# malformed index would otherwise come back as a confidently enumerated rim
# of cells that do not exist.
function _h3_border_checked(c::H3Cell, target::Int)
    lvl = level(c)
    target >= lvl || throw(ArgumentError(
        "border level $target is above the cell's own level $lvl"))
    target <= MAX_RESOLUTION || throw(ArgumentError(
        "border level $target is past max_level $MAX_RESOLUTION"))
    H3Native.is_valid_cell(c.id) || throw(ArgumentError(
        "H3 cell $c is not a valid cell"))
    return lvl, H3Native.is_pentagon(c.id)
end

"""
    subtree_border(sys::H3System, c::H3Cell, l::Integer; connectivity = Vertex()) -> Vector{H3Cell}

Return, in ascending order, level-`l` descendants with a neighbour outside
`c`'s subtree. Complexity is `O(3^depth)` rather than full-subtree
`O(7^depth)`. Counts are `3^(depth+1) - 3` for hexagons and
`5(3^depth - 1)/2` for pentagons; depth zero returns `[c]`. H3 vertex and edge
adjacency coincide, so `connectivity` does not affect the result.

`collect` of [`EdgeCellIterator`](@ref), which is the same automaton lazily.
"""
subtree_border(sys::H3System, c::H3Cell, l::Integer;
    connectivity::Connectivity=Vertex()) =
    DGG.collect_subtree(DGG.EdgeCellIterator(sys, c, l; connectivity))
