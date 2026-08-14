# Exact Eisenstein-lattice arithmetic for Z7 encoding and decoding.
# Level chirality alternates `cbar`, `c`; digit directions follow the fitted
# GBT cycle. Floating point is confined to the lattice/plane bridge.

# ---------------------------------------------------------------------------
# 1. Eisenstein integers
# ---------------------------------------------------------------------------

"""
    MAX_DIGITS

Number of Z7 digit levels the bit format can hold. Matches
`Z7_MAX_RESOLUTION` while keeping this integer layer self-contained.
"""
const MAX_DIGITS = 20

"`omega = e^{2 pi i/3}`, the primitive cube root of unity generating `Z[omega]`."
const OMEGA = complex(-0.5, sqrt(3.0) / 2)

"""
    ecpx(a, b) -> ComplexF64

Complex embedding of the Eisenstein integer `a + b*omega` (lattice -> plane).
"""
ecpx(a::Real, b::Real) = complex(Float64(a), 0.0) + Float64(b) * OMEGA
ecpx(z::Tuple{<:Real,<:Real}) = ecpx(z[1], z[2])

"""
    norm_eis(a, b) -> Integer

Multiplicative norm `|a + b*omega|² = a² - ab + b²`. Units have norm 1;
`norm_eis(3, 1) == norm_eis(2, -1) == 7`.
"""
norm_eis(a::Integer, b::Integer) = a * a - a * b + b * b

"""
    UNITS

The six units of `Z[omega]` in axial coordinates, in counterclockwise order
starting at `+1`; `UNITS[j+1]` has argument `60j` degrees. These are the six
hexagon neighbor directions.
"""
const UNITS = ((1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1))

"""
    ALPHA_DEG

The aperture-7 rotation angle `atand(√3/5) = 19.106605350869096` degrees:
`arg(c) = +ALPHA_DEG`, `arg(cbar) = -ALPHA_DEG`.
"""
const ALPHA_DEG = atand(sqrt(3.0) / 5)

# Multiplication/division/residue for the two chiralities. Matrix form
# (spec section 1.2): mult by a+b*omega is [a -b; b a-b] on axial columns, so
# c = 3+omega is [3 -1; 1 2] and cbar = 2-omega is [2 1; -1 3]; the inverses
# are the conjugate matrices divided by 7.

"Multiply the axial pair `(a, b)` by `c = 3 + omega`."
mul_c(a::Integer, b::Integer) = (3a - b, a + 2b)

"Multiply the axial pair `(a, b)` by `cbar = 2 - omega`."
mul_cbar(a::Integer, b::Integer) = (2a + b, 3b - a)

"""
    div_c(a, b)

Divide by `c = 3 + omega`, exactly: `cbar * (a + b*omega) / 7`. The caller
must have removed the digit first (`res_c(a, b) == 0`), otherwise the
division silently truncates.
"""
div_c(a::Integer, b::Integer) = ((2a + b) ÷ 7, (3b - a) ÷ 7)

"Divide by `cbar = 2 - omega`, exactly (see [`div_c`](@ref))."
div_cbar(a::Integer, b::Integer) = ((3a - b) ÷ 7, (a + 2b) ÷ 7)

"Residue of `a + b*omega` modulo `c = 3 + omega` in `0:6` (`omega ↦ -3 ≡ 4`)."
res_c(a::Integer, b::Integer) = mod(a + 4b, 7)

"Residue of `a + b*omega` modulo `cbar = 2 - omega` in `0:6` (`omega ↦ 2`)."
res_cbar(a::Integer, b::Integer) = mod(a + 2b, 7)

# residue -> unit index j (-1 marks the zero residue, i.e. digit 0). `{0}` plus
# the six units is a complete residue system mod chi (N(chi) = 7 is prime), so
# the table fills without gaps and each level's digit is unique.
function _make_res_to_j(reschi)
    t = fill(-1, 7)
    for j in 0:5
        t[reschi(UNITS[j+1]...)+1] = j
    end
    return ntuple(i -> t[i], 7)
end

"""
    RES_TO_J_C

Residue (`0:6`) to unit index `j` for chirality `c`; `-1` at residue 0 marks
the center digit. Convention-free: it depends only on `c`, not on the digit
labelling.
"""
const RES_TO_J_C = _make_res_to_j(res_c)

"Residue (`0:6`) to unit index `j` for chirality `cbar` (see [`RES_TO_J_C`](@ref))."
const RES_TO_J_CBAR = _make_res_to_j(res_cbar)

"""
    unitmul(a, b, g) -> (a, b)

Multiply the axial pair by the unit of argument `60g` degrees (`g` may be
negative). One step is multiplication by `1 + omega`: `(a, b) -> (a - b, a)`.
Used by the pentagon-collapse and cone-wrap rotations of grid.jl.
"""
function unitmul(a::Integer, b::Integer, g::Integer)
    for _ in 1:mod(g, 6)
        (a, b) = (a - b, a)
    end
    return (a, b)
end

# ---------------------------------------------------------------------------
# 2. Fitted IGEO7 digit conventions
# ---------------------------------------------------------------------------

"""
    chi_is_c(k) -> Bool

Chirality of level `k`: `true` for `chi_k = c = 3 + omega`, `false` for
`cbar = 2 - omega`. IGEO7 alternates, starting with `cbar` at level 1
**[fitted; see `spec/igeo7-geometry-diagnosis.md` §3 — the measured per-res
lattice angles alternate 49.1066°/30° in the face frame, matching the
`cbar`-first sequence to 3.8e-10 deg while every uniform-chirality candidate
misses by 21.8°; the res-5 dump pins the continuation]** — so even-res
lattices are edge-aligned with the res-0 lattice and odd-res lattices are
rotated by `+ALPHA_DEG` (counterclockwise seen from outside).
"""
chi_is_c(k::Integer) = iseven(k)

"""
    SIGMA_J

Digit `d = 1:6` to unit index `j` (dev-frame angle `60j` degrees before the
per-level rotation) **[fitted, A-gauge; see
`spec/igeo7-geometry-diagnosis.md` §4]**: `5->0°, 4->60°,
6->120°, 2->180°, 3->240°, 1->300°`. The GBT complement pairs
`(1,6), (2,5), (3,4)` are antipodal; up to the gauge rotation this is the
published GBT digit cycle `4,6,2,3,1,5` running counterclockwise seen from
outside the sphere (spec/z7-paper-spec.md §3.4).
"""
const SIGMA_J = (5, 3, 4, 1, 0, 2)

"Unit index `j = 0:5` to digit, the inverse of [`SIGMA_J`](@ref)."
const DIGIT_OF_J = ntuple(j -> findfirst(==(j - 1), SIGMA_J), 6)

"""
    sigma(d) -> NTuple{2,Int}

Axial offset contributed by digit `d ∈ 0:6` at its own level (digit 0 is the
center child, offset `(0, 0)`) **[fitted; see [`SIGMA_J`](@ref)]**.
"""
function sigma(d::Integer)
    d == 0 && return (0, 0)
    1 <= d <= 6 || throw(ArgumentError("Z7 digit must be in 0:6, got $d"))
    return @inbounds UNITS[SIGMA_J[d]+1]
end

# residue -> digit under the fitted sigma
_make_res_to_digit(rtj) = ntuple(i -> (j = rtj[i]; j < 0 ? 0 : DIGIT_OF_J[j+1]), 7)

"Residue (`0:6`) to Z7 digit for a level of chirality `c` **[fitted]**."
const RES_TO_DIGIT_C = _make_res_to_digit(RES_TO_J_C)

"Residue (`0:6`) to Z7 digit for a level of chirality `cbar` **[fitted]**."
const RES_TO_DIGIT_CBAR = _make_res_to_digit(RES_TO_J_CBAR)

"Multiply by the chirality of level `k` (see [`chi_is_c`](@ref))."
mulchi(a::Integer, b::Integer, k::Integer) = chi_is_c(k) ? mul_c(a, b) : mul_cbar(a, b)

"Exactly divide by the chirality of level `k` (digit must already be removed)."
divchi(a::Integer, b::Integer, k::Integer) = chi_is_c(k) ? div_c(a, b) : div_cbar(a, b)

"Residue modulo the chirality of level `k`, in `0:6`."
reschi(a::Integer, b::Integer, k::Integer) = chi_is_c(k) ? res_c(a, b) : res_cbar(a, b)

"Digit of a level-`k` residue, under the fitted digit map."
digit_of_res(r::Integer, k::Integer) =
    chi_is_c(k) ? (@inbounds RES_TO_DIGIT_C[r+1]) : (@inbounds RES_TO_DIGIT_CBAR[r+1])

# ---------------------------------------------------------------------------
# 3. Powers of the chirality sequence
# ---------------------------------------------------------------------------

function _make_p_r()
    ps = Vector{Tuple{Int64,Int64}}(undef, MAX_DIGITS + 1)
    ps[1] = (Int64(1), Int64(0))
    for k in 1:MAX_DIGITS
        ps[k+1] = mulchi(ps[k][1], ps[k][2], k)
    end
    return ntuple(i -> ps[i], MAX_DIGITS + 1)
end

"""
    P_R

`P_R[r+1] = prod_{k=1..r} chi_k`, the exact Eisenstein integer relating the
res-`r` lattice to the res-0 lattice: a cell with res-`r` axial coordinates
`X` sits at plane position `L * X / P_r`. With the fitted alternating
chirality `c*cbar = 7`, so `P_{2m} = 7^m` and `P_{2m+1} = 7^m * cbar` — the
accumulated rotation never drifts (design section 6: never accumulate the
angle in floating point).
"""
const P_R = _make_p_r()

"""
    ARGP_DEG

`arg(P_R[r+1])` in degrees, in `[0, 360)`: `0` for even `r`,
`360 − ALPHA_DEG = 340.8933946…` for odd `r`. Used by the pentagon-collapse
rule to convert Horner directions into dev-frame angles.
"""
const ARGP_DEG = ntuple(i -> mod(rad2deg(angle(ecpx(P_R[i]))), 360.0), MAX_DIGITS + 1)

# ---------------------------------------------------------------------------
# 4. Encode / decode
# ---------------------------------------------------------------------------

"""
    horner(digits) -> (a, b)

Horner evaluation of a Z7 digit string (coarse to fine, digits in `0:6`) into
exact axial coordinates on the resolution-`length(digits)` lattice of the
base's frame:

    X = sum_k sigma(d_k) * chi_{k+1} * ... * chi_r

Allocation-free for `NTuple` and `Vector` inputs; two integer multiplies and
an add per level. This is the *lattice* position — the pentagon-collapse
rotation and cone wrap that turn it into a physical dev-frame position live
in grid.jl.
"""
function horner(digits)
    a = Int64(0)
    b = Int64(0)
    k = 0
    for d in digits
        k += 1
        (a, b) = mulchi(a, b, k)
        u = sigma(d)
        a += u[1]
        b += u[2]
    end
    return (a, b)
end

"""
    decode_step(a, b, k) -> (digit, a′, b′)

One exact fine-to-coarse decode step at level `k`: read the digit from the
residue modulo `chi_k`, subtract it, and divide exactly by `chi_k`. Streaming
these from `k = r` down to `1` recovers the digit string; the final `(a, b)`
is the residual res-0 lattice point, which is `(0, 0)` exactly when the cell
belongs to this base's subtree.
"""
function decode_step(a::Integer, b::Integer, k::Integer)
    d = digit_of_res(reschi(a, b, k), k)
    u = sigma(d)
    (qa, qb) = divchi(a - u[1], b - u[2], k)
    return (d, qa, qb)
end

"""
    digit_decode!(out, a, b, r) -> (leftover_a, leftover_b)

Decode `r` digits of the lattice point `(a, b)` into `out[1:r]`
(coarse-to-fine) and return the residual res-0 lattice point. A nonzero
leftover means `(a, b)` is not in this base's planar subtree.
"""
function digit_decode!(out::AbstractVector{<:Integer}, a::Integer, b::Integer, r::Integer)
    for k in r:-1:1
        (d, a, b) = decode_step(a, b, k)
        @inbounds out[k] = d
    end
    return (a, b)
end

"""
    digit_decode(a, b, r) -> (digits, leftover_a, leftover_b)

Allocating convenience wrapper around [`digit_decode!`](@ref).
"""
function digit_decode(a::Integer, b::Integer, r::Integer)
    out = zeros(Int, r)
    (la, lb) = digit_decode!(out, a, b, r)
    return (out, la, lb)
end

# ---------------------------------------------------------------------------
# 5. Continuous position -> lattice point
# ---------------------------------------------------------------------------

"""
    hex_round(a, b) -> (Int64, Int64)

Round continuous axial coordinates to the containing hexagon's center.

Works in cube coordinates `(x, y, z) = (a, b-a, -b)` on the plane
`x + y + z = 0`, where the lattice is the integer triples and the six
neighbors are the permutations of `(1, -1, 0)`: round each coordinate, and if
the sum drifts to ±1 reset the component with the largest rounding error.
The result is the nearest lattice point (proof in
spec/aperture7-indexing-spec.md section 3.2); ties are broken
deterministically by that reset rule.
"""
function hex_round(at::Float64, bt::Float64)
    x, y, z = at, bt - at, -bt
    rx, ry, rz = round(x), round(y), round(z)
    if rx + ry + rz != 0
        ex, ey, ez = abs(rx - x), abs(ry - y), abs(rz - z)
        if ex >= ey && ex >= ez
            rx = -ry - rz
        elseif ey >= ez
            ry = -rx - rz
        else
            rz = -rx - ry
        end
    end
    return (Int64(rx), Int64(-rz))
end
