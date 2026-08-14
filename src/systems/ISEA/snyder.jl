# Snyder ISEA forward and inverse charts plus per-base development-frame slot
# maps. This layer depends only on the shared icosahedron geometry.

# ---------------------------------------------------------------------------
# Planar constants of the equal-area face
# ---------------------------------------------------------------------------

"""
    R_EA

Equal-area scale of an icosahedron face: `sqrt(4π / (15·√3))`, the planar
circumradius of a face triangle whose area is `4π/20` (unit sphere).
Vertices map to planar radius `R_EA` in their faces' frames.
"""
const R_EA = sqrt(4pi / (15 * sqrt(3.0)))

"""
    L_PLANE

Planar icosahedron edge length `L = √3 · R_EA` — the lattice scale: res-`r`
cell centers sit at `u = L · X / P_r` for Eisenstein integers `X` in the
dev frame **[design Section 11]**.
"""
const L_PLANE = sqrt(3.0) * R_EA

"`cos g = cot 60°·cot 36°`; `g = 37.37736814064969°` is the vertex→face-center arc."
const COS_LG = cotd(60.0) * cotd(36.0)

"`tan g = 3 − √5`."
const TAN_LG = tan(acos(COS_LG))

# ---------------------------------------------------------------------------
# Snyder face constants (spec/isea-projection-spec.md §3; G = 36°, θ = 30°
# for the icosahedron face; g the vertex→face-center arc via COS_LG)
# ---------------------------------------------------------------------------

"Snyder's `G = 36°` (face-corner angle of the plane triangle), radians."
const SNY_G = deg2rad(36.0)

"`sin G`, `cos G` with `G = 36°`."
const SNY_SING = sind(36.0)
const SNY_COSG = cosd(36.0)

"`cot θ` with Snyder's `θ = 30°` (per-face sector half-angle): `√3`."
const SNY_COTT = sqrt(3.0)

"One face sector, `2π/3` — the plane angle covered per icosa-face corner."
const SNY_SECTOR = 2pi / 3

"`R_EA²` (the equal-area face circumradius squared)."
const R_EA2 = R_EA * R_EA

const SQRT3 = sqrt(3.0)

"`tan² g` — used by the half-angle form of `sin(q/2)` (see [`_sin_half_q`](@ref))."
const TAN_LG_SQ = TAN_LG * TAN_LG

"Spec §7.1's area coefficient `R_EA²·sin 30° / 2 = R_EA²/4`."
const SNY_AG_COEF = R_EA2 / 4

"`cos 30° = √3/2` (exact-component sector rotations, `sin(· + 30°)` expansion)."
const COS_30 = SQRT3 / 2

"Unit rotations by `k·120°` (one face sector), with exact `±1/2, ±√3/2` parts."
const ROT_SECTOR = (complex(1.0, 0.0), complex(-0.5, COS_30), complex(-0.5, -COS_30))

"rotate the unit complex `e = cis(θ)` by `k` face sectors: `cis(θ + k·2π/3`)"
@inline rot_sector(e::ComplexF64, k::Int) = e * (@inbounds ROT_SECTOR[mod(k, 3)+1])

"""
    _sin_half_q(den) -> Float64

`sin(q/2)` for Snyder's `q = atan(tan g, den)` with `den = cos Az + sin Az·cot θ`
(> 0 on the face, so `q ∈ (0, π/2)`): `cos q = den / √(tan²g + den²)` and the
half-angle identity `sin(q/2) = √((1 − cos q)/2)` — no trig calls
**[first-principles identities; spec §6/§7 define q]**.
"""
@inline function _sin_half_q(den::Float64)
    cq = den / sqrt(TAN_LG_SQ + den * den)
    return sqrt((1 - cq) / 2)
end

# ---------------------------------------------------------------------------
# Dev-frame angular bookkeeping
# ---------------------------------------------------------------------------

"""
    DEV_CONE_DEG

Total planar angle of the development around a base vertex, `300°`: the five
faces at a vertex carry `5·72° = 360°` of spherical angle but only
`5·60° = 300°` of planar angle — the `60°` deficit is why the res-0 cells
are pentagons **[spec/isea-projection-spec.md 8.5]**. Canonical dev-frame
positions live on `[0°, 300°)`.
"""
const DEV_CONE_DEG = 300.0

"""
    DEV_CUT_GUARD_DEG

Planar angles above this are read as the FP shadow of dev angle 0 (a planar
`-ε`) rather than as a continuation past the cone cut. `345°` is the guard
used by the decoder's collapse-order heuristic; [`dev_slot_index`](@ref)
uses `330°`, the midpoint of the missing 60° wedge **[fitted guard bands,
design Section 7]**.
"""
const DEV_CUT_GUARD_DEG = 345.0

"""
    dev_angle_deg(u, cut = DEV_CUT_GUARD_DEG) -> Float64

Planar angle of `u` in degrees, wrapped to `[0, 360)` and then folded to a
small negative value above `cut`, so that the FP shadow of dev angle 0
(a planar `-ε`, which `mod` turns into `360 − ε`) is not mistaken for a
continuation past the cone cut **[design Section 7 guard bands]**.
"""
@inline function dev_angle_deg(u::ComplexF64, cut::Float64=DEV_CUT_GUARD_DEG)
    p = mod(rad2deg(angle(u)), 360.0)
    return p > cut ? p - 360.0 : p
end

# ---------------------------------------------------------------------------
# Faces and face frames (spec §5; triples hard-coded from table T2)
# ---------------------------------------------------------------------------

"The 20 icosahedron faces as base-id triples **[spec/isea-projection-spec.md
§5.2 T2]**; face index is 0-based position in this tuple."
const FACE_TRIPLES = (
    (0, 1, 5), (0, 1, 2), (0, 4, 5), (1, 5, 10), (1, 2, 6),
    (0, 2, 3), (0, 3, 4), (1, 6, 10), (4, 5, 9), (5, 9, 10),
    (2, 6, 7), (2, 3, 7), (3, 4, 8), (6, 10, 11), (4, 8, 9),
    (9, 10, 11), (6, 7, 11), (3, 7, 8), (8, 9, 11), (7, 8, 11),
)

"""
    Face

One icosahedron face of the Snyder chart:

- `verts`   the three base ids of its corners (ascending)
- `c`       unit face center
- `u`, `w`  orthonormal tangent frame at `c`; planar angle 0 points at the
            lowest-numbered corner, `w = c × u` (CCW seen from outside)
- `corner`  planar corner positions `R_EA·cis(0/±120°)`, aligned with `verts`

**[spec/isea-projection-spec.md §5.3–§5.5]**
"""
struct Face
    verts::NTuple{3,Int}
    c::NTuple{3,Float64}
    u::NTuple{3,Float64}
    w::NTuple{3,Float64}
    corner::NTuple{3,ComplexF64}
end

function _make_faces()
    faces = Vector{Face}(undef, 20)
    for (fi, tri) in enumerate(FACE_TRIPLES)
        c = vnormalize(vadd(vadd(VERTICES[tri[1]+1], VERTICES[tri[2]+1]),
            VERTICES[tri[3]+1]))
        v0 = minimum(tri)
        t = vsub(VERTICES[v0+1], vscale(c, vdot(VERTICES[v0+1], c)))
        u = vnormalize(t)
        w = vcross(c, u)
        corner = ntuple(3) do i
            v = VERTICES[tri[i]+1]
            a = atand(vdot(v, w), vdot(v, u))
            ia = round(Int, a / 120) * 120
            @assert abs(a - ia) < 1e-9 "face $(fi-1) corner angle $a not a multiple of 120"
            R_EA * cis(deg2rad(Float64(ia)))
        end
        faces[fi] = Face(tri, c, u, w, corner)
    end
    return ntuple(i -> faces[i], 20)
end

"The 20 [`Face`](@ref)s, indexed by `face + 1`."
const FACES = _make_faces()

"planar corner position of vertex `v` on face `f` (build-time helper)"
function _corner_pos(f::Int, v::Int)
    fc = FACES[f+1]
    i = findfirst(==(v), fc.verts)
    @assert i !== nothing "vertex $v not on face $f"
    return fc.corner[i]
end

# ---------------------------------------------------------------------------
# Snyder ISEA forward / inverse (spec §6, §7)
# ---------------------------------------------------------------------------

"""
    snyder_fwd(p) -> (face, w)

Map grid-frame unit vector `p` to its containing icosahedron face and planar
Snyder coordinate. Ties choose the lowest face index. The map is equal-area
within each face; face centers map to zero and vertices to the face triangle's
corners.
"""
function snyder_fwd(p::NTuple{3,Float64})
    f = 0
    best = -Inf
    @inbounds for k in 0:19
        d = vdot(FACES[k+1].c, p)
        if d > best
            best = d
            f = k
        end
    end
    return (f, snyder_fwd_face(f, p))
end

"""Map `p` into a specified Snyder face chart, including boundary ties."""
function snyder_fwd_face(f::Int, p::NTuple{3,Float64})
    0 <= f < 20 || throw(ArgumentError("Snyder face must lie in 0:19"))
    fc = @inbounds FACES[f+1]
    shz = vnorm(vsub(p, fc.c)) / 2           # sin(z/2) = half the chord to c
    shz <= 0.5e-15 && return complex(0.0, 0.0)   # z < 1e-15: the center
    Az = atan(vdot(p, fc.w), vdot(p, fc.u))
    k = floor(Int, Az / SNY_SECTOR)
    Azs = Az - k * SNY_SECTOR
    ss, cs = sincos(Azs)
    H = acos(clamp(ss * SNY_SING * COS_LG - cs * SNY_COSG, -1.0, 1.0))
    AG = max(Azs + SNY_G + H - pi, 0.0)
    sa = 2 * AG                              # ∝ sin Az′s
    ca = R_EA2 - sa * SNY_COTT               # ∝ cos Az′s
    h = sqrt(sa * sa + ca * ca)
    dp = R_EA * h / (ca + sa * SNY_COTT)     # R_EA / (cos Az′s + sin Az′s·cot θ)
    rho = dp * shz / _sin_half_q(cs + ss * SNY_COTT)
    return rho * rot_sector(complex(ca, sa) / h, k)
end

"""
    snyder_inv_xyz(face, w) -> NTuple{3,Float64}

Map planar Snyder coordinate `w` on `face` back to a grid-frame unit vector.
The azimuth solve uses Newton iteration capped at 10 steps, which spec §7.1
bounds at 4. Coordinates beyond the face triangle are valid for
development-frame fringe cells.
"""
function snyder_inv_xyz(f::Int, w::ComplexF64)
    fc = @inbounds FACES[f+1]
    rho = abs(w)
    rho == 0.0 && return fc.c
    Azp = angle(w)
    k = floor(Int, Azp / SNY_SECTOR)
    Azps = Azp - k * SNY_SECTOR
    sp, cp = sincos(Azps)
    AG = SNY_AG_COEF * sp / (sp * COS_30 + cp * 0.5)   # /sin(Azps + 30°), expanded
    Az = Azps
    sA, cA = sp, cp
    for _ in 1:10
        x = clamp(sA * SNY_SING * COS_LG - cA * SNY_COSG, -1.0, 1.0)
        H = acos(x)
        F = Az + SNY_G + H - pi - AG
        Fp = 1 - (cA * SNY_SING * COS_LG + sA * SNY_COSG) / sqrt(1 - x * x)
        dlt = -F / Fp
        Az += dlt
        abs(dlt) <= 1e-7 && break                # applied step: O(δ²) ≈ 1e-14 left
        sA, cA = sincos(Az)
    end
    sA, cA = sincos(Az)
    dp = R_EA / (cp + sp * SNY_COTT)
    t = clamp((rho / dp) * _sin_half_q(cA + sA * SNY_COTT), -1.0, 1.0)
    cz = 1 - 2 * t * t                           # cos z,  z = 2·asin t
    sz = 2 * t * sqrt(1 - t * t)                 # sin z
    e = rot_sector(complex(cA, sA), k)           # cis(Az + k·2π/3)
    dir = vadd(vscale(fc.u, real(e)), vscale(fc.w, imag(e)))
    return vadd(vscale(fc.c, cz), vscale(dir, sz))
end

# ---------------------------------------------------------------------------
# Dev-frame slot maps: base vertex 300° wedge <-> per-face planar triangles
# ---------------------------------------------------------------------------

"""
    DevSlot

One dev-frame slot of a base vertex: dev angles `[60j, 60j+60)` map rigidly
onto face `f`'s planar triangle by `w = cb + rot·u`, where `cb` is the
base's planar corner position on `f` and `rot = cis(φ)` a multiple of 30°.
`irot = cis(−φ)` is the inverse (decode) rotation. Corner-anchored and a
multiple of 60° between adjacent slots, so the maps are exact
lattice-to-lattice **[fitted; see `spec/igeo7-geometry-diagnosis.md` §4]**.
"""
struct DevSlot
    f::Int
    cb::ComplexF64
    rot::ComplexF64
    irot::ComplexF64
end

function _make_dev_slots()
    slotmaps = Vector{NTuple{5,DevSlot}}(undef, NBASE)
    soffull = Vector{NTuple{20,Int}}(undef, NBASE)
    for b in 0:11
        slots = Vector{DevSlot}(undef, 5)
        sof = fill(-1, 20)
        for j in 0:4
            nj = NBRS_CCW[b+1][j+1]
            nn = NBRS_CCW[b+1][mod(j + 1, 5)+1]
            f = findfirst(t -> b in t && nj in t && nn in t, FACE_TRIPLES) - 1
            cb = _corner_pos(f, b)
            cnj = _corner_pos(f, nj)
            phi = mod(rad2deg(angle(cnj - cb)) - 60.0 * j, 360.0)
            # corner angles are multiples of 120°, so edge directions are
            # multiples of 30° — snap and assert (removes FP angle noise)
            iphi = round(Int, phi / 30) * 30
            @assert abs(phi - iphi) < 1e-9 "slot rotation $phi not a multiple of 30 (b=$b j=$j)"
            rot = cis(deg2rad(Float64(mod(iphi, 360))))
            # consistency: the CCW-next neighbor must land on the third corner
            cnn = _corner_pos(f, nn)
            pred = cb + rot * (L_PLANE * cis(deg2rad(60.0 * (j + 1))))
            @assert abs(pred - cnn) < 1e-12 "dev->face orientation broken (b=$b j=$j)"
            slots[j+1] = DevSlot(f, cb, rot, conj(rot))
            sof[f+1] = j
        end
        slotmaps[b+1] = ntuple(i -> slots[i], 5)
        soffull[b+1] = ntuple(i -> sof[i], 20)
    end
    return ntuple(i -> slotmaps[i], NBASE), ntuple(i -> soffull[i], NBASE)
end

const _DEV_SLOT_TABLES = _make_dev_slots()

"`DEV_SLOTS[b+1][j+1]`: dev-frame slot `j` of base `b` (see [`DevSlot`](@ref))."
const DEV_SLOTS = _DEV_SLOT_TABLES[1]

"""
    DEV_SLOT_OF_FACE

`DEV_SLOT_OF_FACE[b+1][f+1]` = slot index `j` of face `f` around base `b`
(`-1` when `b` is not a corner of `f`) — the decode-side inverse of
[`DEV_SLOTS`](@ref).
"""
const DEV_SLOT_OF_FACE = _DEV_SLOT_TABLES[2]

"""
    dev_slot_index(u) -> Int

Dev-frame slot `j ∈ 0:4` of the dev position `u`: `floor(angle/60°)`, with a
planar angle above `330°` read as the FP shadow of dev angle 0 (canonical
dev positions live on `[0°, 300°)`; see [`DEV_CUT_GUARD_DEG`](@ref))
**[fitted guard; see `spec/igeo7-geometry-diagnosis.md` §4]**.

Computed by half-plane tests instead of `atan`: each 60° boundary ray is the
line `y = ±√3·x` (or the x-axis), so the sector falls out of sign
comparisons — same slot for every input, no trig.
"""
@inline function dev_slot_index(u::ComplexF64)
    x, y = reim(u)
    if y > 0.0 || (y == 0.0 && x >= 0.0)         # angle ∈ [0°, 180°)
        y < SQRT3 * x && return 0                #   [0°, 60°)
        return y > -SQRT3 * x ? 1 : 2            #   [60°, 120°) / [120°, 180°)
    end
    SQRT3 * x < y && return 3                    # [180°, 240°)
    return x <= -SQRT3 * y ? 4 : 0               # [240°, 330°] / (330°, 360°): shadow of 0
end

"""
    dev_to_xyz(base, u) -> NTuple{3,Float64}

Dev-frame position `u` at `base`'s vertex to the unit sphere: slot lookup,
rigid map into the slot's face plane, Snyder inverse. Fringe positions
(beyond the base's own pentagon) are ordinary dev positions — the dev frame
covers the full 5-face neighborhood, so no transport step exists.
"""
@inline function dev_to_xyz(base::Int, u::ComplexF64)
    s = @inbounds DEV_SLOTS[base+1][dev_slot_index(u)+1]
    return snyder_inv_xyz(s.f, s.cb + s.rot * u)
end

"""
    face_to_dev(base, f, w) -> ComplexF64

Planar face position `w` on face `f` into `base`'s dev frame (the decode
direction of [`dev_to_xyz`](@ref)'s rigid map). `base` must be a corner
of `f`.
"""
@inline function face_to_dev(base::Int, f::Int, w::ComplexF64)
    j = @inbounds DEV_SLOT_OF_FACE[base+1][f+1]
    s = @inbounds DEV_SLOTS[base+1][j+1]
    return s.irot * (w - s.cb)
end
