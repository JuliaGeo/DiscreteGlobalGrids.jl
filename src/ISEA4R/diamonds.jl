# ---------------------------------------------------------------------------
# The ten ISEA4R diamonds — the layout table
#
# An icosahedron's twenty faces pair up into ten rhombi ("diamonds"), each two
# faces glued along a shared edge. That pairing, the numbering of the ten, and
# the orientation of the `(x, y)` square inside each of them is the *layout*:
# the one thing an ISEA4R grid needs that the Snyder machinery in `ISEA` does
# not already provide. This file derives it, snaps the resulting affine
# coefficients onto exact values, and asserts the result against pinned
# literals at package load.
#
# ## Provenance — read this before claiming interoperability
#
# This numbering, pairing, and in-diamond orientation is this package's
# canonical choice, derived at build time from `ISEA.FACE_TRIPLES` /
# `ISEA.NBRS_CCW` / `ISEA.NEIGHBORS` in the submodule's standard ISEA
# placement (identity `Orientation`), anchored on the vertex pair `(0, 11)`.
# There is no external oracle pinning it: identifier compatibility with any
# external ISEA4R/ISEA9R product, DGGAL included, is deliberately not claimed
# and must not be inferred. The layout coincides in shape with
# SphericalSpatialTrees.jl's `ISEACircleTree` (10 × 2^r × 2^r), but the diamond
# numbering and per-diamond axis orientation have not been cross-pinned against
# SST either. Anyone needing DGGAL/SST or any other external identifier interop
# must first pin a permutation against fixtures.
#
# That is also why this module lives in `src/ISEA4R/` rather than inside
# `src/ISEA/`: see the `ISEA4R` module docstring.
#
# ## The derivation rule
#
# Three hardcoded convention choices, and everything else follows:
#
#   1. The anchor vertex pair is base `0` and its antipode base `11`.
#   2. *Northern* diamond `d ∈ 0:4` is anchored on dev-slot `d` of base 0: its
#      upper face is `(0, n_d, n_{d+1})` where `n_d = NBRS_CCW[1][d+1]` (the
#      counterclockwise ring of base 0, indices mod 5), and its lower face is
#      the unique *other* face carrying the edge `(n_d, n_{d+1})`.
#   3. *Southern* diamond `5 + d` is the one whose top apex is `n_d`: its upper
#      face is the unique face holding `n_d` and both of `n_d`'s neighbours on
#      the southern ring (`NEIGHBORS[n_d+1] ∩ NEIGHBORS[12]`), and its lower
#      face joins that pair to base 11.
#
# The five northern uppers/lowers are faces `{0,1,2,3,4,5,6,8,11,12}` and the
# southern ones `{7,9,10,13,14,15,16,17,18,19}`; that the twenty faces are used
# exactly once is asserted, not assumed.
#
# ## The chart corners
#
# Per diamond the four corners are `verts = (v00, v10, v11, v01)`, sitting at
# chart `(0,0)`, `(1,0)`, `(1,1)`, `(0,1)`. The glued icosahedron edge — the
# *seam* — is the MAIN DIAGONAL `y == x`, with endpoints `v00` and `v11`; `v01`
# is the upper face's remaining corner and `v10` the lower face's. Northern
# diamonds all have `v01 == 0`, southern ones all have `v10 == 11`.
#
# Which of the two seam endpoints is `v00` is not free: exactly one of the two
# assignments makes the chart positively oriented (counterclockwise seen from
# outside), and the build picks it by testing, then asserts the result against
# the pinned literals below.
# ---------------------------------------------------------------------------

"""
    Diamond

One of the ten ISEA4R diamonds: two icosahedron faces glued along a seam, with
the affine maps that carry the `[0, 1]²` chart square onto their two Snyder
planes.

- `upper`, `lower`  the two 0-based face indices (`ISEA.FACES`); `upper` carries
                    the half `y >= x` and `lower` the half `y < x`
- `verts`           the four corner base ids `(v00, v10, v11, v01)` at chart
                    `(0,0)`, `(1,0)`, `(1,1)`, `(0,1)`; the seam is `(v00, v11)`
- `cP0, aP, bP`     upper-half affine: `u_P(x, y) = cP0 + x·aP + y·bP`, a map
                    into face `upper`'s planar frame
- `cQ0, aQ, bQ`     lower-half affine: `u_Q(x, y) = cQ0 + x·aQ + y·bQ`, into
                    face `lower`'s frame

Both halves agree on the seam (`u_P(t, t)` and `u_Q(t, t)` are the same
spherical point to `< 1e-13` rad — the two planar segments are the two
developments of the same icosahedron edge), so the pair defines one continuous
chart of the square. Each affine has `imag(conj(a) * b) == 2π/5` *exactly*
after the snapping below: the chart is equal-area, and a chart rectangle of
area `A` covers solid angle `A · 4π/10`.

See [`DIAMONDS`](@ref) for the table and [`xyd_to_point`](@ref) for the chart.
"""
struct Diamond
    upper::Int
    lower::Int
    verts::NTuple{4,Int}
    cP0::ComplexF64
    aP::ComplexF64
    bP::ComplexF64
    cQ0::ComplexF64
    aQ::ComplexF64
    bQ::ComplexF64
end

"""
    _PINNED_DIAMOND_VERTS

The corner tuple `(v00, v10, v11, v01)` of each of the ten diamonds, written
out. [`_make_diamonds`](@ref) derives the same table from `ISEA`'s icosahedron
tables and `@assert`s equality against this literal, so drift in an upstream
`ISEA` table fails at package load rather than silently renumbering a shipped
grid. `test/ISEA4R/test_diamonds.jl` restates the same numbers independently.
"""
const _PINNED_DIAMOND_VERTS = (
    (1, 6, 2, 0), (2, 7, 3, 0), (3, 8, 4, 0), (4, 9, 5, 0), (5, 10, 1, 0),
    (10, 11, 6, 1), (6, 11, 7, 2), (7, 11, 8, 3), (8, 11, 9, 4), (9, 11, 10, 5),
)

"""
    _PINNED_DIAMOND_FACES

The `(upper, lower)` 0-based face-index pair of each of the ten diamonds; the
face half of the pin described on [`_PINNED_DIAMOND_VERTS`](@ref). The twenty
entries are the twenty faces, each exactly once.
"""
const _PINNED_DIAMOND_FACES = (
    (1, 4), (5, 11), (6, 12), (2, 8), (0, 3),
    (7, 13), (10, 16), (17, 19), (14, 18), (9, 15),
)

"""
    _CIS30

`cis(k · 30°)` for `k ∈ 0:11`, written with exact components drawn from
`{0, ±1/2, ±√3/2, ±1}` rather than evaluated by `cis`. Used to snap the affine
edge vectors (see [`_snap_edge`](@ref)) — the same snap-and-assert discipline
`ISEA.snyder.jl` applies to its `DevSlot` rotations, and for the same reason:
the exact components remove FP angle noise from a quantity that is known
analytically.
"""
const _CIS30 = (
    complex(1.0, 0.0), complex(SQRT3 / 2, 0.5), complex(0.5, SQRT3 / 2),
    complex(0.0, 1.0), complex(-0.5, SQRT3 / 2), complex(-SQRT3 / 2, 0.5),
    complex(-1.0, 0.0), complex(-SQRT3 / 2, -0.5), complex(-0.5, -SQRT3 / 2),
    complex(0.0, -1.0), complex(0.5, -SQRT3 / 2), complex(SQRT3 / 2, -0.5),
)

"""
    _diamond_corner(f, v) -> ComplexF64

Planar corner position of base vertex `v` on face `f` (build-time helper).

A three-line local clone of `ISEA.snyder.jl`'s `_corner_pos` written against
the exported [`ISEA.FACES`](@ref) table: `_corner_pos` is an `ISEA` build
helper and deliberately private, and reaching into another module's private
build surface is exactly the coupling this package avoids.
"""
function _diamond_corner(f::Int, v::Int)
    fc = FACES[f+1]
    i = findfirst(==(v), fc.verts)
    @assert i !== nothing "vertex $v is not a corner of face $f"
    return fc.corner[i]
end

"""
    _snap_edge(c) -> ComplexF64

Snap a raw face-edge vector onto `L_PLANE · cis(30k°)` with exact components.

Face corners sit at `R_EA · cis(0°/±120°)`, so the difference of any two is an
icosahedron edge in the plane: magnitude exactly `L_PLANE = √3·R_EA`, direction
an exact multiple of 30°. The raw difference carries FP noise on both; this
replaces it with the analytic value and `@assert`s that the raw angle really was
within `1e-9` degrees of a multiple of 30° (measured worst: `2.8e-14`). Same
mechanism, same justification, as the `DevSlot` rotation snap in
`ISEA.snyder.jl`.

The payoff is not cosmetic: after snapping, `imag(conj(aP) * bP)` is
`2π/5` to the *last bit* on all ten diamonds, so the equal-area statement the
chart makes is exact rather than approximate.
"""
function _snap_edge(c::ComplexF64)
    a = rad2deg(angle(c))
    k = round(Int, a / 30)
    @assert abs(a - 30k) < 1e-9 "edge direction $a° is not a multiple of 30°"
    @assert abs(abs(c) - L_PLANE) < 1e-9 "edge length $(abs(c)) is not L_PLANE"
    return L_PLANE * @inbounds _CIS30[mod(k, 12)+1]
end

"""
    _make_diamonds() -> NTuple{10, Diamond}

Derive the ten-diamond table by the rule in this file's header, snap the affine
coefficients, and assert the result — against the pinned literals, against the
face partition, and against the geometry itself. Runs once, at package load.

The assertions, in the order they fire:

1. **Structure.** Every derived face triple exists, the ten `(upper, lower)`
   pairs use each of the twenty faces exactly once, and the seam pair is a real
   icosahedron edge.
2. **Orientation.** Exactly one of the two seam-endpoint assignments gives
   `imag(conj(aP) * bP) > 0` (counterclockwise seen from outside); the other is
   its mirror. Both halves then carry `|det| == 2π/5`, i.e. `4π/10` — the ten
   diamonds account for the whole sphere with nothing left over.
3. **Corners.** The six affine corner predictions (three per half) reproduce the
   planar corner positions to `< 1e-12` (measured `2.5e-16`), and the chart
   corners land on their `ISEA.VERTICES` entries to `< 1e-13` rad (measured
   `1.5e-15`).
4. **Seam.** Upper and lower halves agree along `y == x` to `< 1e-13` rad
   (measured `5.1e-15` over a dense sweep; the coarse build-time check samples
   five parameters).
5. **The pin.** `verts` and `(upper, lower)` equal
   [`_PINNED_DIAMOND_VERTS`](@ref) / [`_PINNED_DIAMOND_FACES`](@ref).
"""
function _make_diamonds()
    ring = NBRS_CCW[1]                       # counterclockwise ring of base 0
    south = NEIGHBORS[12]                    # ring of the antipodal base 11

    face_of(tri) = let k = findfirst(t -> issetequal(t, tri), FACE_TRIPLES)
        @assert k !== nothing "no icosahedron face with corners $tri"
        k - 1
    end

    out = Vector{Diamond}(undef, 10)
    for d in 0:9
        # --- 1. structure: the two faces, the seam pair, and the two apexes ---
        if d < 5                                             # northern
            nd = ring[d+1]
            nd1 = ring[mod(d + 1, 5)+1]
            upper = face_of((0, nd, nd1))
            others = [k - 1 for k in 1:20
                      if nd in FACE_TRIPLES[k] && nd1 in FACE_TRIPLES[k] && k - 1 != upper]
            @assert length(others) == 1 "edge ($nd, $nd1) is not on exactly two faces"
            lower = others[1]
            s1, s2 = nd, nd1
            v01 = 0
            v10 = only(setdiff(FACE_TRIPLES[lower+1], (nd, nd1)))
        else                                                 # southern
            nd = ring[d-5+1]
            pair = sort!(collect(intersect(NEIGHBORS[nd+1], south)))
            @assert length(pair) == 2 "base $nd has $(length(pair)) southern-ring neighbours"
            upper = face_of((nd, pair[1], pair[2]))
            lower = face_of((pair[1], pair[2], 11))
            s1, s2 = pair[1], pair[2]
            v01 = nd
            v10 = 11
        end
        @assert abs(vdot(VERTICES[s1+1], VERTICES[s2+1]) - ADJ_DOT) < 1e-9 "seam ($s1, $s2) is not an icosahedron edge"

        # --- 2. orientation picks which seam endpoint is v00 ---
        chosen = nothing
        for (v00, v11) in ((s1, s2), (s2, s1))
            aP = _snap_edge(_diamond_corner(upper, v11) - _diamond_corner(upper, v01))
            bP = _snap_edge(_diamond_corner(upper, v01) - _diamond_corner(upper, v00))
            imag(conj(aP) * bP) > 0 || continue
            @assert chosen === nothing "both seam orders are positively oriented (d=$d)"
            aQ = _snap_edge(_diamond_corner(lower, v10) - _diamond_corner(lower, v00))
            bQ = _snap_edge(_diamond_corner(lower, v11) - _diamond_corner(lower, v10))
            chosen = Diamond(upper, lower, (v00, v10, v11, v01),
                _diamond_corner(upper, v00), aP, bP,
                _diamond_corner(lower, v00), aQ, bQ)
        end
        @assert chosen !== nothing "neither seam order is positively oriented (d=$d)"
        dm = chosen::Diamond
        # Planar rhombus area, both halves: 4π/10 of the sphere, ten of them 4π.
        @assert imag(conj(dm.aP) * dm.bP) == 2pi / 5 "upper half of diamond $d does not carry 4π/10"
        @assert imag(conj(dm.aQ) * dm.bQ) == 2pi / 5 "lower half of diamond $d does not carry 4π/10"

        # --- 3. corners: the six affine predictions, then the sphere ---
        v00, _v10, v11, _v01 = dm.verts
        @assert abs(dm.cP0 - _diamond_corner(upper, v00)) < 1e-12
        @assert abs(dm.cP0 + dm.bP - _diamond_corner(upper, v01)) < 1e-12
        @assert abs(dm.cP0 + dm.aP + dm.bP - _diamond_corner(upper, v11)) < 1e-12
        @assert abs(dm.cQ0 - _diamond_corner(lower, v00)) < 1e-12
        @assert abs(dm.cQ0 + dm.aQ - _diamond_corner(lower, v10)) < 1e-12
        @assert abs(dm.cQ0 + dm.aQ + dm.bQ - _diamond_corner(lower, v11)) < 1e-12
        for ((x, y), v) in (((0.0, 0.0), dm.verts[1]), ((1.0, 0.0), dm.verts[2]),
            ((1.0, 1.0), dm.verts[3]), ((0.0, 1.0), dm.verts[4]))
            p = y >= x ? snyder_inv_xyz(dm.upper, dm.cP0 + x * dm.aP + y * dm.bP) :
                snyder_inv_xyz(dm.lower, dm.cQ0 + x * dm.aQ + y * dm.bQ)
            @assert deg2rad(angdist(p, VERTICES[v+1])) < 1e-13 "chart corner ($x, $y) of diamond $d misses vertex $v"
        end

        # --- 4. seam: the two halves are the two developments of one edge ---
        for t in (0.0, 0.25, 0.5, 0.75, 1.0)
            p = snyder_inv_xyz(dm.upper, dm.cP0 + t * dm.aP + t * dm.bP)
            q = snyder_inv_xyz(dm.lower, dm.cQ0 + t * dm.aQ + t * dm.bQ)
            @assert deg2rad(angdist(p, q)) < 1e-13 "diamond $d halves disagree on the seam at t=$t"
        end

        out[d+1] = dm
    end

    # --- 5. the pin, and the face partition ---
    faces = sort!(vcat([dm.upper for dm in out], [dm.lower for dm in out]))
    @assert faces == collect(0:19) "the ten diamonds do not partition the twenty faces"
    for d in 0:9
        @assert out[d+1].verts == _PINNED_DIAMOND_VERTS[d+1] "diamond $d corners drifted from the pinned table"
        @assert (out[d+1].upper, out[d+1].lower) == _PINNED_DIAMOND_FACES[d+1] "diamond $d faces drifted from the pinned table"
    end
    return ntuple(i -> out[i], 10)
end

"""
    DIAMONDS

The ten [`Diamond`](@ref)s, indexed by `diamond + 1` (diamond ids are 0-based
throughout this module, as face ids are in `HEALPix` and `S2`).

Derived at load time by [`_make_diamonds`](@ref) and asserted equal to
[`_PINNED_DIAMOND_VERTS`](@ref) / [`_PINNED_DIAMOND_FACES`](@ref):

| d | upper face   | lower face   | `(v00, v10, v11, v01)` | seam    |
|:--|:-------------|:-------------|:-----------------------|:--------|
| 0 | 1  `(0,1,2)` | 4  `(1,2,6)` | `(1, 6, 2, 0)`         | `(1,2)` |
| 1 | 5  `(0,2,3)` | 11 `(2,3,7)` | `(2, 7, 3, 0)`         | `(2,3)` |
| 2 | 6  `(0,3,4)` | 12 `(3,4,8)` | `(3, 8, 4, 0)`         | `(3,4)` |
| 3 | 2  `(0,4,5)` | 8  `(4,5,9)` | `(4, 9, 5, 0)`         | `(4,5)` |
| 4 | 0  `(0,1,5)` | 3 `(1,5,10)` | `(5, 10, 1, 0)`        | `(5,1)` |
| 5 | 7 `(1,6,10)` | 13 `(6,10,11)` | `(10, 11, 6, 1)`     | `(10,6)`|
| 6 | 10 `(2,6,7)` | 16 `(6,7,11)`  | `(6, 11, 7, 2)`      | `(6,7)` |
| 7 | 17 `(3,7,8)` | 19 `(7,8,11)`  | `(7, 11, 8, 3)`      | `(7,8)` |
| 8 | 14 `(4,8,9)` | 18 `(8,9,11)`  | `(8, 11, 9, 4)`      | `(8,9)` |
| 9 | 9 `(5,9,10)` | 15 `(9,10,11)` | `(9, 11, 10, 5)`     | `(9,10)`|

Diamonds `0:4` are northern (all have `v01 == 0`), `5:9` southern (all have
`v10 == 11`). **This numbering is a package convention with no external oracle
behind it** — see the provenance block at the top of this file and
`docs/design/isea4r_diamond_layout.md` before assuming DGGAL / SST / other
external identifier compatibility.
"""
const DIAMONDS = _make_diamonds()
