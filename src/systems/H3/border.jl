# ---------------------------------------------------------------------------
# The subtree border automaton
#
# The rim of a subtree — the descendants at some deeper resolution that have a
# neighbour outside the subtree — without enumerating the subtree. A subtree of
# depth `d` has `7^d` cells but only `O(3^d)` of them on the rim, so walking the
# rim directly is the difference between a feasible query and an infeasible one.
#
# HOW IT WORKS. A cell on the rim is exposed along a contiguous *arc* of the six
# lattice directions, and the arc a child inherits is a function of the parent's
# arc and which digit the child is. The state is `(L, s)`: the arc of exposed
# directions `s, s+1, ..., s+L-1` (mod 6), in `_H3_DIGIT_DIR` positions.
# `L == 6` is the subtree root, fully exposed and the one state with no arc
# ends; `L == 0` means interior, and the walk prunes there. Digit 0 (the centre
# child) is never on the rim.
#
# PENTAGONS. All twelve pentagon base cells delete the same child, digit 1, the
# K axis — H3 has no per-cell deleted digit, so the arc table needs no pentagon
# variant, only the skip. The flag drops to `false` one level down for good:
# H3's pentagon descendants are exactly the all-digit-0 chain, and a rim suffix
# carries no zero digit, so only the subtree root can ever be a pentagon on a
# rim path.
#
# PROVENANCE — this table was *fitted from observation, not derived*, and is
# carried over from the pre-redesign kernel with its evidence intact:
# exhaustive agreement with the definition (enumerate the subtree, keep the
# cells with a `gridDisk` neighbour outside it) at depths 0-4 over 63 roots
# spanning both level parities, all 12 pentagon base cells, higher-resolution
# pentagons, base-cell-crossing roots and random roots at res 1-8; depths 5-6
# over a five-root subset; sampled soundness at depth 10; and an exhaustive
# search over all 720 digit-to-direction assignments that leaves exactly the six
# rotations of `_H3_DIGIT_DIR`, which the root state `(6, 0)` makes immaterial.
# `test/systems/H3/` re-runs the brute-force agreement rather than trusting this.
#
# Note that the parity roles are SWAPPED relative to IGeo7: the branch H3 takes
# at an even child level is the one IGeo7 takes at an odd one.
#
# STATUS. There is no generic subtree border/interior hook in `src/interface/`
# or `src/fallbacks/` as of T5, so this is exposed as an H3-module function
# rather than as a method on a package generic. Wiring it (and IGeo7's and
# HEALPix's) to a shared `subtree_border` generic is a T7 job.
# ---------------------------------------------------------------------------

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

"""
    _h3_fill_border!(out, z, res, target, state, pentagon) -> out

Digit-lexicographic DFS over the border automaton, appending every res-`target`
rim descendant of `z` to `out`. `z` already carries `target` in its resolution
field, so each step overwrites one digit slot and the leaves come out fully
formed — and in ascending order, by construction rather than by sorting.
"""
function _h3_fill_border!(out::Vector{UInt64}, z::UInt64, res::Int, target::Int,
        state::NTuple{2,Int}, pentagon::Bool)
    shift = _h3_digit_shift(res + 1)
    cleared = z & ~(UInt64(7) << shift)
    for digit in 1:6
        pentagon && digit == 1 && continue    # the K axis, deleted under a pentagon
        child_state = _h3_border_step(state, digit, res + 1)
        child_state[1] == 0 && continue
        child = cleared | (UInt64(digit) << shift)
        if res + 1 == target
            push!(out, child)
        else
            _h3_fill_border!(out, child, res + 1, target, child_state, false)
        end
    end
    return out
end

"""
    subtree_border(sys::H3System, c::H3Cell, l::Integer) -> Vector{H3Cell}

The **rim** of `c`'s subtree at resolution `l`: every descendant of `c` at `l`
that has a neighbour outside the subtree, ascending.

`O(3^depth)` rather than the `O(7^depth)` of enumerating the subtree and
testing each cell — the rim of a depth-12 subtree is 1,594,320 cells out of
13,841,287,201.

`subtree_border(sys, c, level(c))` is `[c]`: a depth-0 subtree is the cell
itself, and its whole neighbourhood is outside it.

The count is exactly `3^(depth+1) - 3` for a hexagon and `5(3^depth - 1)/2` for
a pentagon.

!!! note "Not yet a package generic"
    The interface has no subtree border/interior generic as of T5, so this is
    reached as `DiscreteGlobalGrids.H3.subtree_border`. T7 is where the three
    systems' automata get a shared name.
"""
function subtree_border(::H3System, c::H3Cell, l::Integer)
    target = Int(l)
    lvl = level(c)
    target >= lvl || throw(ArgumentError(
        "border level $target is above the cell's own level $lvl"))
    target <= MAX_RESOLUTION || throw(ArgumentError(
        "border level $target is past max_level $MAX_RESOLUTION"))
    # libh3 validates neither `cellToChildren` nor `cellToChildrenSize`, so a
    # malformed index would otherwise come back as a confidently enumerated rim
    # of cells that do not exist.
    H3Native.is_valid_cell(c.id) || throw(ArgumentError(
        "H3 cell $c is not a valid cell"))
    target == lvl && return H3Cell[c]
    # The resolution field moves to `target` once; the digit slots between
    # `lvl` and `target` are filled in on the way down, and the ones below
    # `target` keep the 7 padding they already carry.
    z = _h3_with_resolution(c.id, target)
    pentagon = H3Native.is_pentagon(c.id)
    depth = target - lvl
    p3 = 3^depth
    out = Vector{UInt64}()
    sizehint!(out, pentagon ? (5 * (p3 - 1)) ÷ 2 : 3 * p3 - 3)
    _h3_fill_border!(out, z, lvl, target, (6, 0), pentagon)
    return [H3Cell(id) for id in out]
end
