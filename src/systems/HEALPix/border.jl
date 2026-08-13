# ---------------------------------------------------------------------------
# The subtree rim, without a single neighbour query
#
# Ported from the pre-redesign `src/HEALPix/HealpixKernel.jl`. The derivation
# below is that file's, re-checked against the new interface; the code is the
# same Morton-quadrant recursion.
#
# WHY THE SUBTREE IS A SQUARE BLOCK. A nested id is
# `face * nside^2 + morton(ix, iy)` (`xyf_to_nested`), and refining Δ levels
# scales the lattice by `s = 2^Δ` on the SAME face: pixel `(ix, iy, face)` at
# `level` covers exactly `[ix*s, (ix+1)*s) x [iy*s, (iy+1)*s)` of the leaf
# lattice. Bit-interleaving is positional, so that scaling is a shift of the
# Morton code:
#
#     morton(ix*s + dx, iy*s + dy) = morton(ix, iy) * s^2 + morton(dx, dy)
#
# — the high bits stay the parent's code, the low `2Δ` bits are free. Two
# consequences, and the whole method is built on them:
#
#   * the subtree is the contiguous id range `[id * 4^Δ, (id+1) * 4^Δ)`, which
#     is exactly what `descendant_range` returns; and
#   * a descendant's OFFSET within that range *is* `morton(dx, dy)`, so
#     ascending id order over the subtree is Morton order over the block, and
#     the operation reduces to "emit the perimeter of an `s x s` square in
#     Morton order".
#
# WHY THE RIM IS THAT PERIMETER. The neighbourhood here is the 3x3 lattice
# block, so a descendant at offset `(dx, dy)` has every neighbour inside the
# subtree iff `0 < dx < s-1` and `0 < dy < s-1`. The rim is therefore
# `dx ∈ {0, s-1} || dy ∈ {0, s-1}`, of size `4s - 4`.
#
# Three things could break that argument. None does:
#
#   * *Corner vs edge adjacency.* The 8-neighbourhood gives the IDENTICAL rim
#     to a 4-neighbourhood: every perimeter cell already has an EDGE neighbour
#     outside, and every interior cell has all 8 inside. So this answers the
#     same question under `Vertex()` and under `Edge()`.
#   * *Face seams.* A lattice step off a face edge wraps through `NB_FACEARRAY`,
#     which could in principle deposit the neighbour back inside the block. It
#     cannot: no entry of any non-centre row of that table maps a face to
#     itself, so a wrapped neighbour always lands on a different base face — and
#     a subtree block never leaves its own face.
#   * *The 24 degenerate pixels.* The pixels with only seven neighbours sit on
#     the degree-3 vertices of the base tiling, which in lattice terms are face
#     CORNERS. A face corner is a corner of whatever block contains it, and a
#     block corner has 5 of its 8 neighbours outside, so losing the one missing
#     diagonal still leaves 4. The margin, not luck, is why the missing
#     neighbour cannot un-rim it.
#
# WHY THE WALK RECURSES. Walking the perimeter geometrically (along the bottom,
# up the right, back along the top) is NOT monotone in Morton order and would
# need a sort. `_rim_walk!` instead recurses in Morton-quadrant order: the four
# quadrants of an `s x s` square occupy consecutive id blocks at offsets
# `q * (s/2)^2`, so visiting `q = 0, 1, 2, 3` and pruning the quadrants that
# inherit none of the parent square's exposed sides emits in strictly ascending
# id order by construction. The pruning is what keeps it O(rim).
#
# STATUS: T7 added the generic `subtree_border` hook to `src/interface/`, and
# this walk is now HEALPix's method on it — the signature grew a `connectivity`
# keyword it ignores (see the docstring) and nothing else changed.
# ---------------------------------------------------------------------------

# Extended, not shadowed: `HEALPix.subtree_border` and
# `DiscreteGlobalGrids.subtree_border` are the same function, so generic code
# gets the Morton walk without knowing HEALPix has one.
import ..DiscreteGlobalGrids: subtree_border

# Which sides of a square are exposed to the outside of the subtree, one bit per
# side of the face-local lattice. A quadrant inherits only the sides it actually
# shares with its parent square — the low-x half can be exposed at `dx = 0` but
# never at `dx = s-1` — which is the pruning rule in one line.
const _RIM_XMIN = 0x1
const _RIM_XMAX = 0x2
const _RIM_YMIN = 0x4
const _RIM_YMAX = 0x8

"""
    _rim_walk!(out, k, base, code, sz, mask) -> k

Append the perimeter of the `sz x sz` sub-square whose lower-left corner has
within-subtree Morton code `code`, restricted to the sides named by `mask`, to
`out` starting at `out[k + 1]`; return the new fill mark. Cells are emitted as
`base + code'`, `base` being the subtree's first nested id.

Ascending in the id, because Morton quadrants are visited in id order; O(cells
emitted), because a quadrant inheriting no exposed side is skipped whole rather
than descended and filtered. Never called with `mask == 0`, which is exactly
that skipped case.
"""
function _rim_walk!(out::Vector{Int64}, k::Int, base::Int64, code::Int64,
        sz::Int64, mask::UInt8)
    if sz == 1
        @inbounds out[k += 1] = base + code
        return k
    end
    half = sz >> 1
    quarter = half * half
    for q in Int64(0):Int64(3)
        # `q`'s low bit picks the x half, its high bit the y half (the Morton
        # convention `nested_to_xyf` inverts: even bits -> ix, odd -> iy), so a
        # quadrant keeps the parent's XMIN side only if it is the low-x half,
        # and so on. `m == 0` means every side of this quadrant is interior.
        m = ((q & 1) == 0 ? (mask & _RIM_XMIN) : (mask & _RIM_XMAX)) |
            ((q >> 1) == 0 ? (mask & _RIM_YMIN) : (mask & _RIM_YMAX))
        m == 0x0 && continue
        k = _rim_walk!(out, k, base, code + q * quarter, half, m)
    end
    return k
end

"""
    subtree_border(sys::HEALPixSystem, c::LevelIndex, leaf_level::Integer; connectivity = Vertex()) -> Vector{LevelIndex}

HEALPix's method on the package's [`subtree_border`](@ref) generic.

The **rim** of `c`'s subtree at `leaf_level`: the descendants that have at
least one neighbour outside the subtree, ascending by canonical id.

`4 * 2^Δ - 4` cells for `Δ = leaf_level - level(c) > 0`, and `[c]` at `Δ == 0`.
Θ(rim) time and one allocation — the rim is read straight off the leaf lattice,
with no neighbour query at any level.

`connectivity` is accepted for the generic's signature and does not change the
answer. A subtree is a square block of the face lattice, so a leaf is on the rim
exactly when it sits on one of the block's four sides — and a cell on a side has
a neighbour outside the block under edge adjacency already. `Vertex()` adds
diagonal neighbours, which can only be outside the block for cells that are
already on a side. See the file header.

Positions, if wanted, are `rawid + 1` — the same convention
[`descendant_range`](@ref) returns directly.
"""
function subtree_border(sys::HEALPixSystem, c::DGG.LevelIndex, leaf_level::Integer;
        connectivity::DGG.Connectivity=DGG.Vertex())
    # `descendant_range` runs FIRST and is the only guard needed: it raises the
    # `leaf_level < level(c)` and `> max_level` `ArgumentError`s, and its `lo`
    # is the subtree's first position. Nothing below re-derives either.
    r = DGG.descendant_range(sys, c, leaf_level)
    lo = Int64(first(r)) - 1                 # back to the 0-based nested id
    delta = Int(leaf_level) - DGG.level(c)
    _checked_index(HEALPixGrid(DGG.level(c)), c)
    delta == 0 && return [DGG.LevelIndex(leaf_level, lo)]

    s = Int64(1) << delta
    out = Vector{Int64}(undef, 4 * s - 4)
    k = _rim_walk!(out, 0, lo, Int64(0), s,
        _RIM_XMIN | _RIM_XMAX | _RIM_YMIN | _RIM_YMAX)
    # Deliberately NOT an `@assert`: assertions are elided under
    # `--check-bounds=no`-style flags, and this is the one invariant here whose
    # violation would be SILENT — a short walk leaves `undef` slots, i.e.
    # arbitrary Int64s handed back as cell ids. One comparison against Θ(2^Δ)
    # work, so it stays on unconditionally.
    k == length(out) || error(
        "HEALPix subtree_border filled $k of $(length(out)) rim slots for \
         $c at leaf level $(Int(leaf_level))")
    return [DGG.LevelIndex(leaf_level, id) for id in out]
end
