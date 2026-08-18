# The subtree border is generated in O(3^d) without enumerating its O(7^d)
# interior. A border cell is exposed along a contiguous arc of the six
# lattice directions, and the arc a child inherits is a function of the parent's
# arc and which digit the child is. The state is `(L, s)`: the arc of exposed
# directions `s, s+1, ..., s+L-1` (mod 6), in `_H3_DIGIT_DIR` positions.
# `L == 6` is the subtree root, fully exposed and the one state with no arc
# ends; `L == 0` means interior, and the walk prunes there. Digit 0 (the centre
# child) is never on the border.
#
# All pentagon base cells delete digit 1 (the K axis). The pentagon flag drops
# after the root because a border suffix contains no zero digit.
# The transition table depends on child-level parity. H3's even-level rule
# corresponds to IGeo7's odd-level rule.
#
# The resumable iterator stores an inline frame stack in `O(depth)` memory.
# `H3InteriorEngine` fully descends branches that the border iterator prunes.

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
# the border walk and `0:6` below a pruned branch, so a frame retires at `next > 6`.
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

Iterator state containing the current path id and inline frame stack. `z`
already contains the target resolution; descending overwrites one digit slot.
"""
struct H3Walk
    z::UInt64
    stack::H3Stack
end

# The border walk starts at digit 1: digit 0 is the centre child, which the
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

@inline _h3_border_count(pentagon::Bool, d::Int) =
    d == 0 ? 1 : (pentagon ? (5 * (3^d - 1)) ÷ 2 : 3 * 3^d - 3)

"""
    H3BorderEngine(z, res, target, pentagon)

The border automaton as a resumable iterator over `H3Cell`s. `z` is the subtree
root with its resolution field already moved to `target`.
"""
struct H3BorderEngine
    z::UInt64
    res::Int
    target::Int
    pentagon::Bool
end

Base.eltype(::Type{H3BorderEngine}) = H3Cell
Base.IteratorSize(::Type{H3BorderEngine}) = Base.HasLength()
Base.length(e::H3BorderEngine) = _h3_border_count(e.pentagon, e.target - e.res)

function Base.iterate(e::H3BorderEngine)
    e.target == e.res && return (H3Cell(e.z), H3Walk(e.z, _h3_empty_stack()))
    return iterate(e, _h3_root_walk(e, Int8(1)))
end

Base.iterate(e::H3BorderEngine, w::H3Walk) = _h3_border_advance(e.res, e.target, w)

# Advance either a full-border or seeded-arc walk from its current frame stack.
function _h3_border_advance(res0::Int, target::Int, w::H3Walk)
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

Iterate level-`target` descendants of `z` along exposed directions
`s:s+L-1 (mod 6)`, in ascending id order and `O(depth)` memory. Unlike
`H3BorderEngine`, arbitrary arcs have no closed-form length and report
`SizeUnknown()`.

`pentagon` describes the seed cell. It suppresses the deleted K-axis digit at
the root; descendants reached along a border suffix are not pentagons.
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

Base.iterate(e::H3ArcEngine, w::H3Walk) = _h3_border_advance(e.res, e.target, w)

"""
    H3InteriorEngine(z, res, target, pentagon)

Iterate interior descendants. Branches pruned by the border automaton are wholly
interior and are descended over digits `0:6` without storing a border set.
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
    return _h3_subtree_count(e.pentagon, d) - _h3_border_count(e.pentagon, d)
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
        # Below a branch the automaton pruned, every cell is interior; on the border
        # path, the child it prunes is where the interior starts.
        child = f.full ? (0, 0) : _h3_border_step((Int(f.L), Int(f.s)), digit, res + 1)
        shift = _h3_digit_shift(res + 1)
        z = (z & ~(UInt64(7) << shift)) | (UInt64(digit) << shift)
        if child[1] != 0
            res + 1 == e.target && continue      # a border cell: not ours
            st = DGG.Helpers.small_push(st,
                H3Frame(Int8(child[1]), Int8(child[2]), Int8(0), false, false))
            continue
        end
        res + 1 == e.target && return (H3Cell(z), H3Walk(z, st))
        # A cell stays a pentagon exactly while its digits are all zero, which is
        # reachable here — unlike on the border walk — because digit 0 is taken.
        st = DGG.Helpers.small_push(st,
            H3Frame(Int8(0), Int8(0), Int8(0), f.pentagon && digit == 0, true))
    end
    return nothing
end

function DGG.border_engine(sys::H3System, c::H3Cell, target::Int,
        connectivity::Connectivity)
    lvl, pentagon = _h3_border_checked(c, target)
    return H3BorderEngine(_h3_with_resolution(c.id, target), lvl, target, pentagon)
end

function DGG.interior_engine(sys::H3System, c::H3Cell, target::Int,
        connectivity::Connectivity)
    lvl, pentagon = _h3_border_checked(c, target)
    return H3InteriorEngine(_h3_with_resolution(c.id, target), lvl, target, pentagon)
end

# ---------------------------------------------------------------------------
# Hexagonal halo support
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

# `hex_halo_engine` validates the target level and supplies neighbouring cells.
DGG.seeded_border_engine(::H3System, c::H3Cell, target::Int, arclen::Int,
        start::Int) =
    H3ArcEngine(_h3_with_resolution(c.id, target), level(c), target,
        H3Native.is_pentagon(c.id), Int8(arclen), Int8(start))

# Generate a subtree halo from calibrated arcs of neighbouring subtrees.
DGG.halo_engine(sys::H3System, c::H3Cell, target::Int,
    connectivity::Connectivity) =
    DGG.hex_halo_engine(sys, c, target, connectivity)

# libh3 validates neither `cellToChildren` nor `cellToChildrenSize`, so a
# malformed index would otherwise come back as a confidently enumerated border
# of cells that do not exist.
function _h3_border_checked(c::H3Cell, target::Int)
    lvl = level(c)
    target >= lvl || throw(ArgumentError(
        "border level $target is above the cell's own level $lvl"))
    target <= MAX_RESOLUTION || throw(ArgumentError(
        "border level $target is past maxlevel $MAX_RESOLUTION"))
    H3Native.is_valid_cell(c.id) || throw(ArgumentError(
        "H3 cell $c is not a valid cell"))
    return lvl, H3Native.is_pentagon(c.id)
end
