# ---------------------------------------------------------------------------
# Edge neighbors by GBT digit arithmetic — the global, geometry-free kernel
#
# PROVENANCE — NOT clean-room.
#
# Everything else in this module is clean-room (paper, first principles, or
# fitted to black-box oracle output; see `IGeo7`'s docstring). This file is
# not: the eight tables below and the carry/rotation procedure are ported from
# IGEO7.jl (https://github.com/allixender/IGEO7.jl, `src/IGEO7.jl`), whose own
# comments trace them to `cpp_source/library.h`. The port is used **with the
# author's permission**; IGEO7.jl is AGPL-3.0, so the license this file may be
# distributed under is a decision recorded outside the source — see the note
# in `README.md`. What changed in the port is shape, not arithmetic: static
# tuples for the tables, `SmallList` for the result, this module's
# `InvalidZ7Error` validation, and a single fused pass over the six directions
# instead of three array passes.
#
# Why it exists. `_cell_neighbors_geometric` (grid.jl) answers the same
# question by leaving the index entirely: decode to the Eisenstein lattice,
# step, project with `dev_to_xyz`, and re-decode through `_xyz_to_z7` — a
# floating-point round trip with a search in it, per neighbor. It is correct
# (it is the oracle-validated decoder) and it is ~30x too slow to sweep a
# continent with. Adjacency, though, is not a geometric fact about IGEO7; it
# is an arithmetic one, and it stays inside the digit string.
#
# The arithmetic. Z7's digits are Generalized Balanced Ternary: the seven
# digits at a level are the center plus the six unit directions of a hexagonal
# lattice, and *stepping one cell over* is GBT addition of a direction digit —
# a 7x7 table, because the hexagon's seven-element aggregate is closed under
# it. The addition is not carry-free: a step out of the parent's aggregate
# lands in a sibling, so the table returns a pair, `(carry, digit)`, and the
# carry is added into the digit one level coarser. Exactly a positional
# numeral's ripple, and it terminates the same way — most steps stop at the
# first level.
#
# Two complications, both at the icosahedron:
#
#   * Chirality alternates with level parity (the same `isodd(k)` alternation
#     `_encode_lattice` applies to the Eisenstein basis), so odd levels add
#     clockwise and even levels counter-clockwise — two table pairs, not one.
#
#   * A carry off level 1 leaves the base cell. The direction it left by picks
#     the neighboring base from `BASE_CELL_NEIGHBOURS`, and because adjacent
#     faces of an icosahedron do not share a digit frame, the whole digit
#     string is then rotated into the new face's frame — a multiply by a power
#     of 5 mod 7 (5 generates the multiplicative group mod 7, i.e. it rotates
#     the six directions by one step). `ROTATIONS` and `POLE_0_ROTATIONS`
#     carry how many such steps each crossing needs, and the pentagon's
#     missing digit adds the last correction, `EXCLUDE_NEIGHBOURS`.
#
# Pentagons need no dedup here, unlike the geometric path: the deleted digit
# is known per base, so the missing direction is skipped rather than computed
# and discovered to be a duplicate. Verified over 216,092 cells (every cell of
# res 0-5, plus random samples at res 7/10/14/19): the six raw results are
# always distinct and never a padding sentinel.
#
# Agreement with `_cell_neighbors_geometric` is exact — same ids, same order —
# on every cell of resolutions 0-5 (196,092 cells), on 3,000 random cells at
# each of resolutions 6, 8, 10, 12, 15 and 19, and on all twelve pentagon
# chains through resolution 8. That sweep is `test/IGeo7/test_gbt_neighbors.jl`.
# ---------------------------------------------------------------------------

# --- ported tables ---------------------------------------------------------

"""
    BASE_CELL_NEIGHBOURS

`BASE_CELL_NEIGHBOURS[base+1][carry]` — the base cell reached by carrying off
level 1 of `base` in direction `carry`, in space-filling-curve order from the
`k = 1` poking direction. Every base is a pentagon (the twelve icosahedron
vertices), so one of the six slots is not a real direction; which one is
[`EXCLUDE_NEIGHBOURS`](@ref), and the slot is filled with a repeat of its
neighbor rather than a sentinel.
"""
const BASE_CELL_NEIGHBOURS = (
    (5, 4, 4, 2, 1, 3),      # base 0  (north pole)
    (5, 0, 0, 6, 10, 2),     # base 1
    (1, 0, 0, 7, 6, 3),      # base 2
    (2, 0, 0, 8, 7, 4),      # base 3
    (3, 0, 0, 9, 8, 5),      # base 4
    (4, 0, 0, 10, 9, 1),     # base 5
    (10, 2, 1, 11, 11, 7),   # base 6
    (6, 3, 2, 11, 11, 8),    # base 7
    (7, 4, 3, 11, 11, 9),    # base 8
    (8, 5, 4, 11, 11, 10),   # base 9
    (9, 1, 5, 11, 11, 6),    # base 10
    (9, 6, 10, 8, 8, 7),     # base 11 (south pole)
)

"""
    EXCLUDE_NEIGHBOURS

`EXCLUDE_NEIGHBOURS[base+1]` — the one direction slot of
[`BASE_CELL_NEIGHBOURS`](@ref) that is not a real edge of the pentagon at
`base`, and equally the digit its pentagon chain deletes: `2` for the northern
bases `0:5`, `5` for the southern `6:11`.

Equal by construction to [`Z7_DELETED_DIGIT`](@ref) (z7.jl), which this module
derives independently from `vectors/pentagon_chains.csv`; the two are pinned
equal in `test/IGeo7/test_gbt_neighbors.jl` rather than one being defined as
the other, so the fit and the port stay separate evidence.
"""
const EXCLUDE_NEIGHBOURS = (2, 2, 2, 2, 2, 2, 5, 5, 5, 5, 5, 5)

"""
    ROTATIONS

`ROTATIONS[base+1]` — digit-frame rotation steps applied when a carry off
`base` crosses into a *polar* base (0 or 11), before the finer correction of
[`POLE_0_ROTATIONS`](@ref). Counted in multiplications by 5 mod 7; see
[`POW5_MOD7`](@ref).
"""
const ROTATIONS = (0, 5, 0, 1, 3, 4, 5, 4, 3, 1, 0, 0)

"""
    POLE_0_ROTATIONS

`POLE_0_ROTATIONS[row][col]` — extra digit-frame rotation steps for a carry
*out of* a polar base, indexed by the level-1 digit of the source (`row`) and
of the result (`col`), both `1:6`. The south pole reuses the north pole's
table under the digit reflection `d -> 7 - d`.
"""
const POLE_0_ROTATIONS = (
    (0, 1, 0, 1, 0, 2),
    (0, 0, 0, 0, 0, 0),
    (0, 0, 0, 3, 2, 2),
    (5, 5, 0, 0, 0, 0),
    (0, 1, 0, 0, 0, 0),
    (5, 0, 4, 0, 4, 0),
)

# GBT addition, clockwise (odd levels). `GBT_CW_0[a+1][b+1]` is the digit and
# `GBT_CW_1[a+1][b+1]` the carry of `a + b`; row/column 0 is the identity
# because digit 0 is the aggregate's center.
const GBT_CW_0 = (
    (0, 1, 2, 3, 4, 5, 6),
    (1, 4, 3, 6, 5, 2, 0),
    (2, 3, 1, 4, 6, 0, 5),
    (3, 6, 4, 5, 0, 1, 2),
    (4, 5, 6, 0, 2, 3, 1),
    (5, 2, 0, 1, 3, 6, 4),
    (6, 0, 5, 2, 1, 4, 3),
)

const GBT_CW_1 = (
    (0, 0, 0, 0, 0, 0, 0),
    (0, 1, 0, 1, 0, 5, 0),
    (0, 0, 2, 3, 0, 0, 2),
    (0, 1, 3, 3, 0, 0, 0),
    (0, 0, 0, 0, 4, 4, 6),
    (0, 5, 0, 0, 4, 5, 0),
    (0, 0, 2, 0, 6, 0, 6),
)

# GBT addition, counter-clockwise (even levels). The digit table is plain
# addition mod 7 — the counter-clockwise frame is the one in which the six
# directions are numbered consecutively — but the carries are not, so the pair
# still has to be looked up.
const GBT_CCW_0 = (
    (0, 1, 2, 3, 4, 5, 6),
    (1, 2, 3, 4, 5, 6, 0),
    (2, 3, 4, 5, 6, 0, 1),
    (3, 4, 5, 6, 0, 1, 2),
    (4, 5, 6, 0, 1, 2, 3),
    (5, 6, 0, 1, 2, 3, 4),
    (6, 0, 1, 2, 3, 4, 5),
)

const GBT_CCW_1 = (
    (0, 0, 0, 0, 0, 0, 0),
    (0, 1, 0, 3, 0, 1, 0),
    (0, 0, 2, 2, 0, 0, 6),
    (0, 3, 2, 3, 0, 0, 0),
    (0, 0, 0, 0, 4, 5, 4),
    (0, 1, 0, 0, 5, 5, 0),
    (0, 0, 6, 0, 4, 0, 6),
)

"""
    POW5_MOD7

`POW5_MOD7[n+1] = 5^n mod 7` for the `n ∈ 0:6` the rotation tables can ask
for. 5 generates `(Z/7)^*`, so multiplying every digit by `5^n` rotates the
six directions `n` steps while fixing the center digit 0 — which is what
re-seating a digit string into an adjacent face's frame is.
"""
const POW5_MOD7 = (1, 5, 4, 6, 2, 3, 1)

# --- the step --------------------------------------------------------------

# GBT addition at level `k`: chirality alternates on absolute level parity,
# matching `_encode_lattice`'s `isodd(k)`.
@inline function _gbt_add(k::Int, a::Int, b::Int)
    return isodd(k) ?
           (@inbounds(GBT_CW_1[a+1][b+1]), @inbounds(GBT_CW_0[a+1][b+1])) :
           (@inbounds(GBT_CCW_1[a+1][b+1]), @inbounds(GBT_CCW_0[a+1][b+1]))
end

# Multiply digits `from:to` by `mult` mod 7, i.e. rotate that stretch of the
# digit string into another face's frame. Digit 0 is fixed, so the pentagon
# chain of a prefix survives the rotation.
@inline function _rotate_digits(z::UInt64, from::Int, to::Int, mult::Int)
    @inbounds for k in from:to
        z = _z7_set_digit(z, k, (_z7_digit(z, k) * mult) % 7)
    end
    return z
end

# Level of the first non-zero digit, or `res + 1` when every active digit is
# zero (a pentagon). Used to find the digit that decides the exclusion-zone
# rotation below.
@inline function _first_nonzero_level(z::UInt64)
    masked = z & Z7_PAD_MASK
    masked == zero(UInt64) && return Z7_MAX_RESOLUTION + 1
    return (leading_zeros(masked) - 4) ÷ 3 + 1
end

# One neighbor: add direction `dir` at the finest level and ripple the carry
# up. Returns the index and the carry that survived past level 1 (0 when the
# neighbor stayed inside the base cell).
@inline function _gbt_step(z7::UInt64, dir::Int, res::Int)
    carry, digit = _gbt_add(res, _z7_digit(z7, res), dir)
    out = _z7_set_digit(z7, res, digit)
    if carry != 0
        @inbounds for k in (res-1):-1:1
            carry, digit = _gbt_add(k, _z7_digit(z7, k), carry)
            out = _z7_set_digit(out, k, digit)
            carry == 0 && break
        end
    end
    return out, carry
end

# Re-seat a neighbor that carried off level 1 into the base cell it landed in,
# rotating its digit string into that face's frame.
@inline function _gbt_cross_base(nb::UInt64, z7::UInt64, base::Int, carry::Int, res::Int)
    new_base = @inbounds BASE_CELL_NEIGHBOURS[base+1][carry]
    nb = (nb & Z7_PAD_MASK) | (UInt64(new_base) << Z7_BASE_SHIFT)

    # Into a pole: the source face's own twist, plus one more step when the
    # crossing leaves by a digit adjacent to the pole's seam.
    if new_base == 0 || new_base == Z7_NUM_BASES - 1
        rot = @inbounds ROTATIONS[base+1]
        d1 = _z7_digit(z7, 1)
        (d1 == 1 || d1 == 6) && (rot += 1)
        rot > 0 && (nb = _rotate_digits(nb, 1, res, @inbounds POW5_MOD7[rot+1]))
    end

    # Out of a pole: the correction depends on both level-1 digits, and the
    # south pole is the north pole reflected.
    if base == 0 || base == Z7_NUM_BASES - 1
        row = _z7_digit(z7, 1)
        col = _z7_digit(nb, 1)
        if base == Z7_NUM_BASES - 1
            row = 7 - row
            col = 7 - col
        end
        if 1 <= row <= 6 && 1 <= col <= 6
            rot = @inbounds POLE_0_ROTATIONS[row][col]
            rot > 0 && (nb = _rotate_digits(nb, 1, res, @inbounds POW5_MOD7[rot+1]))
        end
    end
    return nb
end

"""
    _cell_neighbors(z7) -> SmallList{6,UInt64}

Canonical ids of the cells sharing an edge with `z7`, ascending: 6 for a
hexagon, 5 for a pentagon. Pure integer GBT digit arithmetic — no geometry, no
floating point, no allocation (see the block comment above for the derivation
and the provenance of the tables).

Exactly agrees with [`_cell_neighbors_geometric`](@ref), the oracle-validated
lattice implementation, over the sweep recorded above; it is the
implementation the [`DiscreteGlobalGrids.cell_neighbors`](@ref) kernel wiring
calls.

Throws [`InvalidZ7Error`](@ref) for invalid ids and for resolution-20 ids
(valid for prefix arithmetic, no geometry — hence no neighbors).
"""
function _cell_neighbors(z7::UInt64)
    res = _geometry_checked(z7)
    base = z7_base_cell(z7)
    excluded = @inbounds EXCLUDE_NEIGHBOURS[base+1]
    out = Helpers.empty_small_list(Val(6), zero(UInt64))

    # Resolution 0 is the base cell itself: no digits to add into, so the
    # neighbors are just its five distinct entries in the adjacency table.
    if res == 0
        @inbounds for dir in 1:6
            dir == excluded && continue
            nb = BASE_CELL_NEIGHBOURS[base+1][dir]
            out = Helpers.small_push(out, (UInt64(nb) << Z7_BASE_SHIFT) | Z7_PAD_MASK)
        end
        return Helpers.small_sort(out)
    end

    # A pentagon has no edge in its base's deleted direction, so that step is
    # skipped rather than taken and deduplicated. Away from the chain the
    # deleted digit is a real direction again, and one of the six results can
    # land on it — the exclusion-zone rotation below moves such a result back
    # into the frame its own prefix implies.
    pentagon = z7_is_pentagon(z7)
    mult = 0
    if !pentagon
        zone = _z7_digit(z7, _first_nonzero_level(z7))
        mult = (zone * 5) % 7 == excluded ? 5 : (zone * 3) % 7 == excluded ? 3 : 0
    end

    for dir in 1:6
        pentagon && dir == excluded && continue
        nb, carry = _gbt_step(z7, dir, res)
        carry != 0 && (nb = _gbt_cross_base(nb, z7, base, carry, res))
        if mult != 0
            lvl = _first_nonzero_level(nb)
            lvl <= res && _z7_digit(nb, lvl) == excluded &&
                (nb = _rotate_digits(nb, lvl, res, mult))
        end
        out = Helpers.small_push(out, nb)
    end
    return Helpers.small_sort(out)
end
