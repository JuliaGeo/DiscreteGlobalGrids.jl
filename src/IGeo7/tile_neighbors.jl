# ---------------------------------------------------------------------------
# Edge neighbors inside a subtree, by integer digit arithmetic
#
# `_cell_neighbors` (grid.jl) answers for *any* cell, so it decodes to the
# Eisenstein lattice, steps, and hands the result back through `dev_to_xyz` and
# a full `_xyz_to_z7` re-decode — floating-point angles, the cone cut, the
# pentagon collapse, base-cell crossing. All of that machinery is about leaving
# the neighborhood a subtree defines, and inside one tile none of it applies:
# every cell shares the root's prefix, hence its base cell.
#
# What is left is exact and small. The Horner encode splits at the root:
#
#     L(cell) = X_p · L(prefix)  +  L(suffix)
#
# because level `k`'s contribution is `sigma(d_k)` carried through the product
# of the chirality matrices above it, so the prefix digits factor out as a
# constant `X_p · L(prefix)` shared by the whole tile. An edge neighbor is one
# Eisenstein unit away, and that unit lands entirely on `L(suffix)`:
#
#     neighbor:  L(suffix) + u,   u in UNITS
#
# Peeling `r - p` digits back off the sum inverts the Horner (`res_*` picks the
# digit, `div_*` removes it). The residual is the punchline: it is zero exactly
# when the neighbor kept the root's prefix, i.e. when it is still in the tile.
# Membership is not a separate test — it falls out of computing the address.
#
# For a hexagon root the suffix digits are a plain base-7 numeral and the
# subtree holds `7^d` cells, so the digit string *is* the position. A pentagon
# root is missing its base's deleted digit at every level it stays a pentagon
# (`p(d) = (5·7^d + 1)/6` cells, not `7^d`), the numeral is no longer base-7,
# and this construction does not apply — `neighbor_stepper` keeps the generic
# path there. The failure is total rather than subtle (every position is
# wrong, since the radix is wrong), which is why the guard is a construction
# check and not a per-cell branch.
# ---------------------------------------------------------------------------

"""
    Z7TileNeighborStepper(tile)

Edge neighbors inside an IGEO7 subtree by exact integer arithmetic on the Z7
digit suffix — no geometry, no floating point, no memory. See the block comment
above for the derivation.

Requires a non-pentagon root; [`DGG.neighbor_stepper`](@ref) checks that and
falls back to the generic stepper otherwise.
"""
struct Z7TileNeighborStepper{T} <: DGG.AbstractNeighborStepper
    tile::T
    p::Int          # root level: the number of prefix digits held constant
    r::Int          # leaf level
end

# `POW7[d+1] == 7^d` is already this module's own (grid.jl), where it sizes a
# hexagon prefix's subtree — the same quantity, for the same reason.

# Suffix numeral -> lattice point. The chirality alternates on *absolute* level
# parity (`_encode_lattice`'s `isodd(k)`), so the loop is over true levels
# `p+1:r`, not over `1:d` — a tile rooted at an odd level and one rooted at an
# even level do not share this map.
@inline function _suffix_lattice(s::Int, p::Int, r::Int)
    a = 0
    b = 0
    @inbounds for k in (p + 1):r
        a, b = isodd(k) ? (2a + b, 3b - a) : (3a - b, a + 2b)   # mul_cbar / mul_c
        d = (s ÷ POW7[r - k + 1]) % 7
        if d != 0
            e = SIGMA_AB[d]
            a += e[1]
            b += e[2]
        end
    end
    return (a, b)
end

# Lattice point -> suffix numeral, or `-1` when the point leaves the subtree.
# Peels fine to coarse; the residual after `r - p` digits is zero exactly when
# the prefix survived.
@inline function _lattice_suffix(a::Int, b::Int, p::Int, r::Int)
    s = 0
    @inbounds for k in r:-1:(p + 1)
        d = isodd(k) ? RES_TO_DIGIT_CBAR[mod(a + 2b, 7) + 1] :
                       RES_TO_DIGIT_C[mod(a + 4b, 7) + 1]
        if d != 0
            e = SIGMA_AB[d]
            a -= e[1]
            b -= e[2]
        end
        # exact division: the digit above removed the residue
        a, b = isodd(k) ? ((3a - b) ÷ 7, (a + 2b) ÷ 7) : ((2a + b) ÷ 7, (3b - a) ÷ 7)
        s += d * POW7[r - k + 1]
    end
    return ((a == 0) & (b == 0)) ? s : -1
end

@inline function DGG.step_neighbors(s::Z7TileNeighborStepper, i::Int)
    p, r = s.p, s.r
    a, b = _suffix_lattice(i - 1, p, r)
    out = SmallVector{6,Int}()
    @inbounds for j in 1:6
        u = UNITS[j]
        t = _lattice_suffix(a + u[1], b + u[2], p, r)
        # The six unit directions are not in id order, so insert rather than
        # push: the stepper contract is ascending, as `cell_neighbors` is.
        t < 0 || (out = DGG._insert_sorted(out, t + 1))
    end
    return out
end

# The guard the block comment above argues for. A pentagon root makes the
# suffix numeral non-base-7 and every position wrong, so it is decided once,
# here, rather than branched on per cell. A non-pentagon root cannot have a
# pentagon descendant — a pentagon's ancestors are all pentagons — so this one
# check covers the whole tile.
function DGG.neighbor_stepper(t::DGG.DGGSSubtreeIds{DGG.IGEO7DGGS})
    z7_is_pentagon(UInt64(t.root_id)) && return DGG.GenericNeighborStepper(t)
    return Z7TileNeighborStepper(t, t.root_level, t.level)
end
