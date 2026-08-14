# grid.jl — the composed grid: encode (cell id -> center/boundary/area) and
# decode (lon/lat -> cell id), spec/igeo7-geometry-diagnosis.md §6–§7 /
# spec/design.md sections 11–12.
#
# This layer glues the four modules below it together (design Section 2):
#   icosahedron.jl  vertex/neighbor tables, `Orientation`
#   z7.jl           the UInt64 cell id, digit slots, validation
#   engine.jl       exact Eisenstein lattice arithmetic + fitted digit tables
#   snyder.jl        the per-face Snyder ISEA plane + dev-frame slot maps
#
# Ported from the verified clean-provenance prototype fit
# (spec/igeo7-geometry-diagnosis.md §4), which reproduces every oracle cell
# center at res 1..5 (all 196,080) within the dumps' own print noise and
# decodes each to the exact z7 string. Provenance of each rule is cited
# inline.
#
# The geometry pipeline:
#
#   Z7 digits --Horner--> res-r lattice point X in the base's dev frame
#            --pentagon collapse + cone wrap--> physical lattice point
#            --u = L*X/P_r--> dev-frame position
#            --slot map + Snyder inverse--> unit sphere --Orientation--> world
#
# and the decoder runs it backwards — Snyder forward into the containing
# face, the face's three corner bases nearest-first through the slot maps,
# hex rounding in place of the lattice evaluation, and a strict re-encode
# acceptance test in place of any nearest-center arbitration. The dev frame
# covers each base's full 5-face neighborhood, so fringe cells are ordinary
# dev positions: no rim transport exists anywhere.

# ---------------------------------------------------------------------------
# Constant tables (design Section 6: hoist every per-call lookup)
# ---------------------------------------------------------------------------

"Unit rotations by ±60°, the cone-wrap step (wrapping ±300° of dev angle)."
const CIS_P60 = cis(deg2rad(60.0))
const CIS_M60 = cis(-deg2rad(60.0))

"""
    CELL_SCALE

`CELL_SCALE[r+1] = L / P_r :: ComplexF64` — the res-`r` lattice-to-dev
scale: a cell with physical lattice coordinates `X` sits at dev position
`u = ecpx(X) * CELL_SCALE[r+1]` **[design Section 11]**.
"""
const CELL_SCALE = ntuple(i -> L_PLANE / ecpx(P_R[i]), MAX_DIGITS + 1)

"`INV_CELL_SCALE[r+1] = P_r / L`, the dev-to-lattice scale (decode side)."
const INV_CELL_SCALE = ntuple(i -> ecpx(P_R[i]) / L_PLANE, MAX_DIGITS + 1)

"""
    THETA_DIR

`THETA_DIR[m+1][d]` = dev direction (degrees in `[0, 360)`) of digit `d`'s
Horner contribution at level `m`: `mod(60·σ(d) − arg P_m, 360)` — the digit's
unit direction rotated into the res-`m` lattice **[a7; fitted]**.
Level 0 (`m = 0`) is the res-0 lattice itself (used by pentagon corners).
"""
const THETA_DIR = ntuple(mi -> ntuple(d -> mod(60.0 * SIGMA_J[d] - ARGP_DEG[mi], 360.0), 6),
    MAX_DIGITS + 1)

"""
    CORNER_OFFS

Hexagon corner offsets in lattice units: corner `j` of the cell at lattice
point `X` is `(X + (e_j + e_{j+1})/3) * L/P_r` — the circumradius of a unit
hexagon is `1/√3` and its corners bisect adjacent neighbor directions
**[design Section 4.1; first-principles hex geometry]**.
"""
const CORNER_OFFS = ntuple(i -> (ecpx(UNITS[i]) + ecpx(UNITS[mod1(i + 1, 6)])) / 3, 6)

"Closed-form cell areas: `HEX_AREA[r+1] = 4πR²/(10·7^r)` m² **[contract]**."
const HEX_AREA = ntuple(i -> 4pi * R_AUTHALIC^2 / (10 * 7.0^(i - 1)), MAX_DIGITS)

"`SIGMA_AB[d]` = axial offset of digit `d ∈ 1:6` (fused `UNITS[SIGMA_J[d]+1]`)."
const SIGMA_AB = ntuple(d -> UNITS[SIGMA_J[d]+1], 6)

"FP guard band (degrees) for dev angles on the cone cut **[design Section 7]**."
const CUT_EPS_DEG = 1e-9

"""
    TIE_EPS

Rounding-tie band in lattice units (design Sections 7 and 10 item 4): a
point within this margin of a cell-boundary bisector counts as a tie. Dev
position noise is ~1e-15 plane units = `1e-15·7^(r/2)/L ≈ 1e-7` lattice
units at res 19, so `1e-6` covers every resolution with an order of margin
while staying far below any geometrically meaningful distance. Used only by
the decoder's last-resort tie fallback (see [`lonlat_to_z7`](@ref)).
"""
const TIE_EPS = 1e-6

# ---------------------------------------------------------------------------
# Encode: digits -> physical lattice point (Horner + pentagon collapse +
# cone wrap)
# ---------------------------------------------------------------------------

"""
    _encode_lattice(z, base, res) -> (a, b)

Physical res-`res` lattice point of cell `z` in its base's dev frame,
streaming the digits straight out of the `UInt64` (no buffer):

1. Horner: `X = Σ σ(d_k)·χ_{k+1}…χ_r` (2 integer multiplies + add per level)
   **[a7]**;
2. pentagon collapse **[fitted; see `spec/igeo7-geometry-diagnosis.md`
   §4]**: at the level `m` of the first nonzero digit, the deleted
   digit's 60° sector is excised and every kept subtree whose level-`m`
   direction is CCW-past the deleted direction rotates by −60° (an exact
   unit rotation);
3. cone wrap: the whole-subtree representative is wrapped into dev angle
   `[0°, 300°)` — wrapping by ±300° of angle is a unit rotation by ∓60°.

Returns `(0, 0)` exactly for pentagons (all-zero digits).
"""
function _encode_lattice(z::UInt64, base::Int, res::Int)
    a = Int64(0)
    b = Int64(0)
    m = 0
    dm = 0
    for k in 1:res
        # mulchi(a, b, k) with the chirality branch and digit tables inlined
        # (identical arithmetic to engine.jl's mul_cbar/mul_c/sigma:
        # cbar at odd levels, c at even)
        (a, b) = isodd(k) ? (2a + b, 3b - a) : (3a - b, a + 2b)
        d = _z7_digit(z, k)
        if d != 0
            e = @inbounds SIGMA_AB[d]
            a += e[1]
            b += e[2]
            if m == 0
                m = k
                dm = d
            end
        end
    end
    m == 0 && return (Int64(0), Int64(0))
    thc = @inbounds THETA_DIR[m+1][dm]
    thdel = @inbounds THETA_DIR[m+1][z7_deleted_digit(base)]
    shift = thc > thdel ? -1 : 0
    (a, b) = unitmul(a, b, shift)
    thslot = mod(thc + 60.0 * shift, 360.0)
    # wrap the whole-subtree representative into [0, 300) of dev angle:
    psi = mod(rad2deg(angle(ecpx(a, b))) - (@inbounds ARGP_DEG[res+1]), 360.0)
    dd = mod(psi - thslot + 180.0, 360.0) - 180.0    # cell angle rel. its slot
    psit = thslot + dd
    if psit < 0.0
        (a, b) = unitmul(a, b, -1)
    elseif psit >= DEV_CONE_DEG
        (a, b) = unitmul(a, b, 1)
    end
    return (a, b)
end

# ---------------------------------------------------------------------------
# Shared validation
# ---------------------------------------------------------------------------

"""
    _geometry_checked(z) -> res

Validate a cell id for a geometry operation: structurally valid Z7 index
(`is_valid_z7`) at a resolution with geometry (`0:19`); returns the
resolution. Res-20 ids are valid for prefix arithmetic but have no geometry
**[contract]**.

Both failure edges throw an [`InvalidZ7Error`](@ref) (`:invalid_index` /
`:res20_geometry`) directly. The message is built only if the error is
printed, so this validator (inlined into every geometry and hierarchy entry
point) carries no string machinery.
"""
@inline function _geometry_checked(z::UInt64)
    is_valid_z7(z) || throw(InvalidZ7Error(:invalid_index, z, 0, 0))
    res = _z7_leading_resolution(z)
    res <= MAX_RESOLUTION ||
        throw(InvalidZ7Error(:res20_geometry, z, Z7_MAX_RESOLUTION, MAX_RESOLUTION))
    return res
end

# ---------------------------------------------------------------------------
# Public encode API (design Section 4.1)
# ---------------------------------------------------------------------------

"grid-frame unit vector of the cell center (validated id, res <= 19)"
function _cell_center_xyz(z::UInt64, res::Int)
    base = z7_base_cell(z)
    (a, b) = _encode_lattice(z, base, res)
    (a == 0 && b == 0) && return vertex(base)        # pentagon: the vertex
    u = ecpx(a, b) * (@inbounds CELL_SCALE[res+1])
    return dev_to_xyz(base, u)
end

"""
    cell_center(z7; orientation=ORIENT_IDENTITY) -> (lon, lat)

Center of cell `z7` in degrees, `(lon, lat)` order **[contract]**. Pentagon
centers are exactly the icosahedron vertices; hexagon centers are the cell's
physical lattice point in its base's dev frame, pulled back through the
slot map and Snyder inverse. Under the default identity orientation this is
the standard ISEA placement — validated at 100% exact z7 decode on all
196,080 oracle cell centers, res 1–5. The `orientation` rotation is applied
to the output (`from_grid`, design Section 5).

Throws [`InvalidZ7Error`](@ref) for invalid ids and for resolution-20 ids
(valid for prefix arithmetic, no geometry).
"""
function cell_center(z7::UInt64; orientation::Orientation=ORIENT_IDENTITY)
    res = _geometry_checked(z7)
    return xyz_to_lonlat(from_grid(orientation, _cell_center_xyz(z7, res)))
end
cell_center(z7::Unsigned; kwargs...) = cell_center(UInt64(z7); kwargs...)

"""
    cell_boundary_cartesian(z7; closed_ring=true, orientation=ORIENT_IDENTITY)
        -> Vector{NTuple{3,Float64}}

Boundary ring of cell `z7` as unit-sphere `xyz` tuples: 6 corners for a
hexagon, 5 for a pentagon; `closed_ring` repeats the first corner at the end
**[contract]**. Rings wind counterclockwise (seen from outside) under the
identity orientation.

Hexagon corners are the lattice-space midpoint constructions
`(X + (e_j + e_{j+1})/3)·L/P_r`, each canonicalized across the cone cut by
the cell center's branch (fringe corners are ordinary dev positions — no
transport). Pentagon corners sit at radius `s_r/√3` at dev angles
`θ_del + 60k − 30` — the bisectors between the pentagon's five ring slots in
the res-`r` lattice **[fitted; design Section 4.1]**.
"""
function cell_boundary_cartesian(z7::UInt64; closed_ring::Bool=true,
    orientation::Orientation=ORIENT_IDENTITY)
    res = _geometry_checked(z7)
    base = z7_base_cell(z7)
    pent = z7_is_pentagon(z7)
    n = pent ? 5 : 6
    out = Vector{NTuple{3,Float64}}(undef, n + (closed_ring ? 1 : 0))
    if pent
        # corner directions live in the res-r lattice (rotation arg P_r)
        thdel = THETA_DIR[res+1][z7_deleted_digit(base)]
        rho = abs(CELL_SCALE[res+1]) / SQRT3
        for k in 1:5
            ang = thdel + 60.0 * k - 30.0            # bisector between slots
            u = rho * cis(deg2rad(mod(ang, DEV_CONE_DEG)))
            out[k] = from_grid(orientation, dev_to_xyz(base, u))
        end
    else
        (a, b) = _encode_lattice(z7, base, res)
        c0 = ecpx(a, b)
        thc = mod(rad2deg(angle(c0)) - ARGP_DEG[res+1], 360.0)  # center dev angle
        scale = CELL_SCALE[res+1]
        for j in 0:5
            u = (c0 + CORNER_OFFS[j+1]) * scale
            psi = mod(rad2deg(angle(u)), 360.0)
            if psi >= DEV_CONE_DEG - CUT_EPS_DEG
                # corner rep is across the cut: canonicalize by the branch
                # the cell center sits on (cone angle ψ−300 CCW / ψ−60 CW)
                u *= (thc > 150.0 ? CIS_P60 : CIS_M60)
            end
            out[j+1] = from_grid(orientation, dev_to_xyz(base, u))
        end
    end
    closed_ring && (out[end] = out[1])
    return out
end
cell_boundary_cartesian(z7::Unsigned; kwargs...) =
    cell_boundary_cartesian(UInt64(z7); kwargs...)

"""
    cell_boundary(z7; closed_ring=true, orientation=ORIENT_IDENTITY)
        -> Vector{NTuple{2,Float64}}

[`cell_boundary_cartesian`](@ref) in `(lon, lat)` degrees **[contract]**.
"""
function cell_boundary(z7::UInt64; closed_ring::Bool=true,
    orientation::Orientation=ORIENT_IDENTITY)
    cart = cell_boundary_cartesian(z7; closed_ring, orientation)
    return NTuple{2,Float64}[xyz_to_lonlat(p) for p in cart]
end
cell_boundary(z7::Unsigned; kwargs...) = cell_boundary(UInt64(z7); kwargs...)

"""
    cell_area(z7) -> Float64

Cell area in square metres, closed form **[contract]**: every hexagon at
resolution `r` has exactly `4πR²/(10·7^r)` and every pentagon exactly `5/6`
of that (12 pentagons + `10·7^r − 10` hexagons tile the authalic sphere).
The Snyder chart is exactly equal-area, so the closed form is the geometry
(design Section 4.1).
"""
function cell_area(z7::UInt64)
    res = _geometry_checked(z7)
    hex = @inbounds HEX_AREA[res+1]
    return z7_is_pentagon(z7) ? 5 * hex / 6 : hex
end
cell_area(z7::Unsigned) = cell_area(UInt64(z7))

# ---------------------------------------------------------------------------
# Decode (spec/igeo7-geometry-diagnosis.md §6; ported from the verified
# prototype fit of that document's §4)
# ---------------------------------------------------------------------------

"""
    _try_decode(base, u, res[, tie, gs]) -> (ok, z)

Decode the dev-frame position `u` (at `base`'s vertex) at resolution `res`:

1. hex-round `u·P_r/L` to the nearest lattice point **[a7]**;
2. canonicalize a rounded point whose dev angle lands in `[300°, 360°)`
   across the cone cut — which continuation branch (`ψ−300` vs `ψ−60`) is
   decided by the side the *measured* angle came from (> 150° or not)
   **[fitted, geometric]**;
3. un-collapse: try the unit rotations `g ∈ (0, 1, −1, 2, −2)` — `g` absorbs
   the pentagon-collapse rotation and the cut wrap. For each, a pure-Horner
   digit decode streams digits directly into the Z7 word (residue mod 7 →
   digit, subtract, exact divide by χ) **[a7]**; a candidate is accepted iff
   the leftover is exactly `(0, 0)`, the digit string is valid (first nonzero
   digit ≠ the base's deleted digit), and a strict re-encode
   ([`_encode_lattice`](@ref)) reproduces the canonical rounded point.
   Global-lattice consistency makes at most one `g` accept, so the search is
   equivalent to computing `g` directly and doubles as a self-check.

Returns `ok = false` when no candidate accepts — the owner then roots in
another corner base of the containing face (see [`lonlat_to_z7`](@ref)).

With `tie = true` (the last-resort pass), a failing primary candidate falls
back to the runner-up lattice points whose Voronoi margin is within
[`TIE_EPS`](@ref) of the primary, in ascending-margin order: on an exact
cell-boundary tie (e.g. a point *on* an icosahedron edge midpoint) FP noise
can flip the rounding to a lattice point outside every reachable subtree in
every base; the runner-up is then the equally-near owner. Behavior on ties
is deterministic but not bit-promised **[design Section 10 item 4]**.
"""
function _try_decode(base::Int, u::ComplexF64, res::Int, tie::Bool=false,
    gs::Tuple{Vararg{Int}}=(0, 1, -1, 2, -2))
    zc = u * (@inbounds INV_CELL_SCALE[res+1])
    x = real(zc)
    y = imag(zc)
    (a, b) = hex_round(x + y / SQRT3, 2 * y / SQRT3)
    del = z7_deleted_digit(base)
    zbase = (UInt64(base) << Z7_BASE_SHIFT) | _z7_tail_mask(res)
    ok, z = _decode_candidate(base, u, res, a, b, del, zbase, gs)
    (ok || !tie) && return (ok, z)
    # tie fallback: neighbors of the rounded point, nearest bisector first
    dc = zc - ecpx(a, b)
    tried = 0
    for _ in 1:6
        bj = 0
        bm = TIE_EPS                     # only genuine ties qualify
        for j in 1:6
            (tried >> j) & 1 == 1 && continue
            mj = (abs2(dc - ecpx(UNITS[j])) - abs2(dc)) / 2
            if mj < bm
                bm = mj
                bj = j
            end
        end
        bj == 0 && break
        tried |= 1 << bj
        e = @inbounds UNITS[bj]
        ok, z = _decode_candidate(base, u, res, a + e[1], b + e[2], del, zbase, gs)
        ok && return (true, z)
    end
    return (false, zbase)
end

"cut-canonicalize one rounded lattice point and run the `g` search on it"
function _decode_candidate(base::Int, u::ComplexF64, res::Int,
    a::Int64, b::Int64, del::Int, zbase::UInt64,
    gs::Tuple{Vararg{Int}}=(0, 1, -1, 2, -2))
    if !(a == 0 && b == 0)
        px = mod(rad2deg(angle(ecpx(a, b))) - (@inbounds ARGP_DEG[res+1]), 360.0)
        if px >= DEV_CONE_DEG - CUT_EPS_DEG
            pu = mod(rad2deg(angle(u)), 360.0)
            (a, b) = unitmul(a, b, pu > 150.0 ? 1 : -1)
        end
    end
    for g in gs
        (ca, cb) = unitmul(a, b, g)
        z = zbase
        m = 0
        dm = 0
        # engine.jl's decode_step chain with the chirality branch and digit
        # tables inlined (identical arithmetic: res_cbar/res_c residue,
        # sigma subtraction, exact div_cbar/div_c — cbar at odd levels,
        # c at even), streaming digits directly into the Z7 word
        for k in res:-1:1
            if isodd(k)
                d = @inbounds RES_TO_DIGIT_CBAR[mod(ca + 2cb, 7)+1]
                if d != 0
                    e = @inbounds SIGMA_AB[d]
                    ca -= e[1]
                    cb -= e[2]
                    z |= UInt64(d) << _z7_shift(k)
                    m = k
                    dm = d
                end
                (ca, cb) = ((3ca - cb) ÷ 7, (ca + 2cb) ÷ 7)
            else
                d = @inbounds RES_TO_DIGIT_C[mod(ca + 4cb, 7)+1]
                if d != 0
                    e = @inbounds SIGMA_AB[d]
                    ca -= e[1]
                    cb -= e[2]
                    z |= UInt64(d) << _z7_shift(k)
                    m = k
                    dm = d
                end
                (ca, cb) = ((2ca + cb) ÷ 7, (3cb - ca) ÷ 7)
            end
        end
        (ca == 0 && cb == 0) || continue             # not in this base's subtree
        (m == 0 || dm != del) || continue            # deleted pentagon chain
        _encode_lattice(z, base, res) == (a, b) || continue  # strict re-encode
        return (true, z)
    end
    return (false, zbase)
end

"""
    _g_pair(base, u) -> (g1, g2)

Predicted un-collapse rotation order for the decoder's fast pass: a point's
dev angle relative to the base's deleted-digit direction (`60·σ(del)`, the
collapse boundary) predicts the accepting `g`. Pure ordering heuristic — a
misprediction costs one wasted decode + re-encode and never changes the
accepted candidate (at most one `(base, g)` accepts, pinned by the suites).
"""
@inline function _g_pair(base::Int, u::ComplexF64)
    thdel = 60.0 * (@inbounds SIGMA_J[z7_deleted_digit(base)])
    g = dev_angle_deg(u) > thdel ? 1 : 0
    return (g, 1 - g)
end

"""
    _corner_bases(p, f) -> NTuple{3,Int}

The three corner bases of Snyder face `f`, nearest-first to the unit vector
`p` (descending vertex dot product) — the decode candidate order
**[diagnosis §6: containing face + its 3 corners]**. Allocation-free 3-sort.
"""
@inline function _corner_bases(p::NTuple{3,Float64}, f::Int)
    v1, v2, v3 = (@inbounds FACES[f+1]).verts
    d1 = vdot(p, @inbounds VERTICES[v1+1])
    d2 = vdot(p, @inbounds VERTICES[v2+1])
    d3 = vdot(p, @inbounds VERTICES[v3+1])
    if d1 < d2
        (v1, v2), (d1, d2) = (v2, v1), (d2, d1)
    end
    if d1 < d3
        (v1, v3), (d1, d3) = (v3, v1), (d3, d1)
    end
    d2 < d3 && ((v2, v3) = (v3, v2))
    return (v1, v2, v3)
end

"""
    lonlat_to_z7(lon, lat, res; orientation=ORIENT_IDENTITY) -> UInt64

Z7 cell id of the res-`res` cell whose published polygon contains the point
(degrees, `(lon, lat)` order). The Snyder forward map picks the containing
face; its three corner bases are tried nearest-first through the dev-frame
slot maps and [`_try_decode`](@ref) — first with the predicted `g` pair,
then the full `g` search, then the rounding-tie fallback. Exactly one owner
accepts (global-lattice consistency); there is no nearest-center arbitration
anywhere. Validated at 100% exact decode on all 196,080 oracle cell centers,
res 1–5.

The search passes are result-neutral ordering (the at-most-one-accepts
invariant): pass 1 is the predicted `(g, 1−g)` pair over the three corner
bases, pass 2 the full `g` search (insurance + self-check), pass 3 the
rounding-tie fallback — kept last so ordinary points are never affected.

Throws [`InvalidZ7Error`](@ref) (`:resolution_range`) for `res ∉ 0:19`.
"""
function lonlat_to_z7(lon::Real, lat::Real, res::Integer;
    orientation::Orientation=ORIENT_IDENTITY)
    0 <= res <= MAX_RESOLUTION || throw(InvalidZ7Error(
        :resolution_range, zero(UInt64), _z7_int(res), MAX_RESOLUTION))
    return _xyz_to_z7(to_grid(orientation, lonlat_to_xyz(Float64(lon), Float64(lat))), Int(res))
end

"""
    _xyz_to_z7(p, res) -> UInt64

[`lonlat_to_z7`](@ref)'s decode body on a *grid-frame* unit vector: the
containing face, its corner bases nearest-first, the three result-neutral
search passes. Factored out so grid-frame producers — the neighbor step in
[`_cell_neighbors`](@ref) hands over `dev_to_xyz` output directly — skip the
degrees round trip and the orientation rotation.
"""
function _xyz_to_z7(p::NTuple{3,Float64}, r::Int)
    f, w = snyder_fwd(p)
    bs = _corner_bases(p, f)
    for i in 1:3
        b = @inbounds bs[i]
        u = face_to_dev(b, f, w)
        ok, z = _try_decode(b, u, r, false, _g_pair(b, u))
        ok && return z
    end
    for tie in (false, true)
        for i in 1:3
            b = @inbounds bs[i]
            u = face_to_dev(b, f, w)
            ok, z = _try_decode(b, u, r, tie)
            ok && return z
        end
    end
    throw(ErrorException(
        "IGeo7 internal error: no cell accepted grid-frame point $p at res $r"))
end

"""
    lonlat_to_cell(lon, lat, res; orientation=ORIENT_IDENTITY) -> UInt64

Alias of [`lonlat_to_z7`](@ref): the Z7 `UInt64` *is* the cell id
(design Section 3).
"""
lonlat_to_cell(lon::Real, lat::Real, res::Integer; kwargs...) =
    lonlat_to_z7(lon, lat, res; kwargs...)

# ---------------------------------------------------------------------------
# Edge neighbors
#
# A res-`r` cell's region is the Voronoi hexagon of its lattice point, so its
# edge neighbors sit at exactly the six Eisenstein units away — the same
# first-principles hex geometry the corner construction rests on. The step is
# taken on the *physical* (post-collapse, cone-wrapped) lattice point that
# `_encode_lattice` returns, a representative crossing the cone cut is
# canonicalized by the center's branch exactly as `cell_boundary_cartesian`
# canonicalizes corner representatives, and the resulting position — the
# neighbor's own center, exact lattice arithmetic through the exact slot
# maps — is handed to the standard decoder. No new convention enters: every
# step is a fitted-and-validated piece of the existing pipeline, and the
# strict re-encode acceptance inside `_try_decode` rejects any candidate that
# is not exactly the cell standing at that position.
#
# Pentagons need no special geometry: at the cone apex the six unit
# directions cover 360° of raw angle folded onto the 300° cone, so exactly
# two of them land on the same physical ring slot (for either wrap branch the
# folded direction coincides with another unit's slot) and the pentagon's
# five neighbors fall out of id-deduplication.
# ---------------------------------------------------------------------------

"""
    _cell_neighbors_directional(z7) -> SmallList{6,UInt64}

Canonical ids of the cells sharing an edge with `z7`, in the counterclockwise
order of `UNITS`: 6 for a hexagon, 5 for a pentagon. Exact lattice
arithmetic for the neighbor positions (see the block comment above); the
position-to-id step is the decoder validated at 100% exact decode on all
196,080 oracle cell centers.

Throws [`InvalidZ7Error`](@ref) for invalid ids and for resolution-20 ids
(valid for prefix arithmetic, no geometry — hence no neighbors).
"""
function _cell_neighbors_directional(z7::UInt64)
    res = _geometry_checked(z7)
    base = z7_base_cell(z7)
    (a, b) = _encode_lattice(z7, base, res)
    c0 = ecpx(a, b)
    scale = @inbounds CELL_SCALE[res+1]
    # Which branch a cut-crossing representative wraps to is decided by the
    # side the center sits on (`cell_boundary_cartesian`'s rule; a neighbor
    # step never swings more than the 180° that rule discriminates). At the
    # apex there is no side and no need for one: both branches fold the
    # crossing direction onto another unit's slot, and the duplicate id is
    # dropped below.
    thc = mod(rad2deg(angle(c0)) - (@inbounds ARGP_DEG[res+1]), 360.0)
    out = Helpers.empty_small_list(Val(6), zero(UInt64))
    for j in 1:6
        e = @inbounds UNITS[j]
        u = (c0 + ecpx(e)) * scale
        psi = mod(rad2deg(angle(u)), 360.0)
        if psi >= DEV_CONE_DEG - CUT_EPS_DEG
            u *= (thc > 150.0 ? CIS_P60 : CIS_M60)
        end
        z = _xyz_to_z7(dev_to_xyz(base, u), res)
        seen = false
        for k in 1:length(out)
            (@inbounds out[k]) == z && (seen = true; break)
        end
        seen || (out = Helpers.small_push(out, z))
    end
    return out
end

"""
    _cell_neighbors(z7) -> SmallList{6,UInt64}

Canonical ids of the cells sharing an edge with `z7`, ascending. This is the
system-kernel order; [`_cell_neighbors_directional`](@ref) exposes the native
counterclockwise construction order.
"""
_cell_neighbors(z7::UInt64) = Helpers.small_sort(_cell_neighbors_directional(z7))

# ---------------------------------------------------------------------------
# Dense full-world indexing, hierarchy wrappers and introspection
# (design Sections 3, 4.3; contract "Native API")
#
# The canonical full-world order is **ascending encoded cell id**. Because the
# base cell occupies bits [63:60] and the padding sentinel `7` sorts after
# every active digit, ascending `UInt64` order *is* (base, digit-string)
# lexicographic order, so the dense index of a cell is a positional rank in a
# mixed-radix digit string:
#
#   * a hexagon prefix has all seven child digits, so its subtree at depth `d`
#     holds `7^d` cells;
#   * a pentagon prefix (all digits so far zero) is missing its base's deleted
#     digit at every level for which it is still a pentagon, so its subtree
#     holds `p(d) = (5·7^d + 1) / 6` cells (1, 6, 41, 288, … — the 6 and 41
#     are pinned by the contract);
#   * hence `num_cells(r) = 12·p(r) = 10·7^r + 2` **[contract]**.
#
# `cell_to_index` walks the digit string once, adding the subtree size of
# every earlier sibling; `index_to_cell` runs the same walk backwards. Both
# are O(res) integer ops with no allocation and no division by a runtime
# value other than the precomputed table entries.
# ---------------------------------------------------------------------------

"`POW7[d+1] = 7^d` — subtree size of a *hexagon* prefix at depth `d`."
const POW7 = ntuple(i -> Int64(7)^(i - 1), MAX_RESOLUTION + 1)

"""
    PENT_COUNT

`PENT_COUNT[d+1] = p(d) = (5·7^d + 1) / 6` — subtree size of a *pentagon*
prefix at depth `d` **[contract: 6 at d = 1, 41 at d = 2]**. Recurrence
`p(d) = 6·7^(d-1) + p(d-1)`: a pentagon has six children, five of them
hexagons (full `7^(d-1)` subtrees) and one — digit 0 — a pentagon again.
"""
const PENT_COUNT = ntuple(i -> (5 * Int64(7)^(i - 1) + 1) ÷ 6, MAX_RESOLUTION + 1)

"`NUM_CELLS[r+1] = 12·p(r) = 10·7^r + 2` **[contract, oracle: num_cells.csv]**."
const NUM_CELLS = ntuple(i -> 10 * Int64(7)^(i - 1) + 2, MAX_RESOLUTION + 1)

"""
    num_cells(res) -> Int64

Number of cells covering the whole world at resolution `res`:
`10·7^res + 2` **[contract; test/IGeo7/vectors/num_cells.csv]** — twelve
pentagons and `10·7^res − 10` hexagons. Fits `Int64` through `res = 19`
(1.14e17).

Throws [`InvalidZ7Error`](@ref) (`:resolution_range`) for
`res ∉ 0:$MAX_RESOLUTION`.
"""
@inline function num_cells(res::Integer)
    0 <= res <= MAX_RESOLUTION || throw(InvalidZ7Error(
        :resolution_range, zero(UInt64), _z7_int(res), MAX_RESOLUTION))
    return @inbounds NUM_CELLS[Int(res)+1]
end

"""
    res0_cells() -> SmallList{12,UInt64}

The twelve resolution-0 cells (one pentagon per icosahedron vertex) in
ascending id order **[contract]**. The res-0 id of base `b` is
`(b << 60) | 0x0fff…f`, so ascending id order is ascending base order and the
list is exactly dense indexes `1:12`.
"""
function res0_cells()
    out = Helpers.empty_small_list(Val(Z7_NUM_BASES), zero(UInt64))
    for base in 0:(Z7_NUM_BASES-1)
        out = Helpers.small_push(out, (UInt64(base) << Z7_BASE_SHIFT) | Z7_PAD_MASK)
    end
    return out
end

"""
    cell_to_index(id) -> Int

One-based rank of `id` among all cells of its own resolution, in the
canonical full-world order (ascending encoded id) **[contract, design
Section 4.3]**.

The walk is `base·p(res)` for the whole-base blocks that precede it, plus,
for each level `k`, the total size of the subtrees of the siblings that sort
before digit `d_k`. While the prefix is still all-zero (still a pentagon)
digit 0 leads to a pentagon subtree of size `p(res−k)` and the base's deleted
digit is absent; afterwards every sibling subtree has `7^(res−k)` cells.
O(res), allocation-free.

Throws [`InvalidZ7Error`](@ref) for a structurally invalid id or one at
resolution 20 (no geometry, hence no dense index).
"""
function cell_to_index(z7::UInt64)
    res = _geometry_checked(z7)
    base = z7_base_cell(z7)
    deleted = @inbounds Z7_DELETED_DIGIT[base+1]
    rank = base * (@inbounds PENT_COUNT[res+1])
    pentagon = true
    for k in 1:res
        d = _z7_digit(z7, k)
        depth = res - k
        if pentagon
            if d != 0
                # digit 0 (the pentagon child) precedes every other digit;
                # the hexagon siblings below `d` are 1:(d-1) minus the
                # deleted digit, which `is_valid_z7` has ruled out for `d`.
                nhex = d - 1 - (deleted < d ? 1 : 0)
                rank += (@inbounds PENT_COUNT[depth+1]) + nhex * (@inbounds POW7[depth+1])
                pentagon = false
            end
        else
            rank += d * (@inbounds POW7[depth+1])
        end
    end
    return Int(rank) + 1
end
cell_to_index(z7::Unsigned) = cell_to_index(UInt64(z7))

"""
    index_to_cell(index, res) -> UInt64

Inverse of [`cell_to_index`](@ref): the `index`-th cell of resolution `res`
in canonical order. Greedy inverse of the rank walk — peel the base block,
then at each level take the pentagon child while the remainder fits its
subtree, otherwise divide by `7^(res−k)` to pick the hexagon sibling
(re-inserting the deleted digit's gap). O(res), allocation-free.

Throws `BoundsError` for `index ∉ 1:num_cells(res)` **[contract — the
`BoundsError` is part of the API shape and is kept]** and
[`InvalidZ7Error`](@ref) (`:resolution_range`, via `num_cells`) for
`res ∉ 0:$MAX_RESOLUTION`.
"""
function index_to_cell(index::Integer, res::Integer)
    n = num_cells(res)                       # also validates `res`
    1 <= index <= n || throw(BoundsError(Base.OneTo(n), index))
    r = Int(res)
    remainder = Int64(index) - 1
    block = @inbounds PENT_COUNT[r+1]
    base = Int(remainder ÷ block)
    remainder -= base * block
    deleted = @inbounds Z7_DELETED_DIGIT[base+1]
    z = (UInt64(base) << Z7_BASE_SHIFT) | Z7_PAD_MASK
    pentagon = true
    for k in 1:r
        depth = r - k
        width = @inbounds POW7[depth+1]
        digit = 0
        if pentagon
            centre = @inbounds PENT_COUNT[depth+1]
            if remainder >= centre
                remainder -= centre
                q = remainder ÷ width
                remainder -= q * width
                # ascending hexagon siblings are 1:6 with `deleted` removed
                digit = Int(q) + 1
                digit >= deleted && (digit += 1)
                pentagon = false
            end
        else
            q = remainder ÷ width
            remainder -= q * width
            digit = Int(q)
        end
        z = _z7_set_digit(z, k, digit)
    end
    return z
end

"""
    lonlat_to_index(lon, lat, res; orientation=ORIENT_IDENTITY) -> Int

Dense full-world index of the res-`res` cell containing `(lon, lat)`
(degrees) — [`lonlat_to_z7`](@ref) composed with [`cell_to_index`](@ref),
allocation-free **[contract]**.
"""
lonlat_to_index(lon::Real, lat::Real, res::Integer; kwargs...) =
    cell_to_index(lonlat_to_z7(lon, lat, res; kwargs...))

# --- hierarchy wrappers (geometry-side names) ------------------------------

"""
    cell_to_parent(id) -> UInt64
    cell_to_parent(id, res) -> UInt64

Ancestor of `id` at resolution `res` (default: one coarser) — [`z7_parent`](@ref)
with the geometry-side validation of the cell id **[contract]**.

Throws [`InvalidZ7Error`](@ref) for an invalid id (`:invalid_index`), for a
resolution-20 id (`:res20_geometry`, no geometry), for
`res ∉ 0:get_resolution(id)` (`:parent_res`) and for a resolution-0 cell
(`:no_parent`).
"""
@inline function cell_to_parent(z7::UInt64, res::Integer)
    current = _geometry_checked(z7)
    0 <= res <= current ||
        throw(InvalidZ7Error(:parent_res, z7, _z7_int(res), current))
    return z7 | _z7_tail_mask(res)
end

@inline function cell_to_parent(z7::UInt64)
    current = _geometry_checked(z7)
    current > 0 || throw(InvalidZ7Error(:no_parent, z7, 0, 0))
    return z7 | _z7_tail_mask(current - 1)
end

@inline cell_to_parent(z7::Unsigned) = cell_to_parent(UInt64(z7))
@inline cell_to_parent(z7::Unsigned, res::Integer) = cell_to_parent(UInt64(z7), res)

"""
    cell_to_children(id) -> SmallList{7,UInt64}

Immediate children of `id` in ascending id order — seven for a hexagon, six
for a pentagon **[contract]**. [`z7_children`](@ref) with the geometry-side
validation; resolution-$MAX_RESOLUTION cells throw [`InvalidZ7Error`](@ref)
(`:no_child_geometry`) because their children would have no geometry.
"""
function cell_to_children(z7::UInt64)
    res = _geometry_checked(z7)
    res < MAX_RESOLUTION ||
        throw(InvalidZ7Error(:no_child_geometry, z7, 0, MAX_RESOLUTION))
    return z7_children(z7)
end

"""
    cell_to_children(id, res) -> SmallList{1,UInt64} / SmallList{7,UInt64} / Vector{UInt64}

*All* descendants of `id` at resolution `res`, ascending **[contract]**. The
return type follows the depth so the shallow cases stay allocation-free, as
the spatial-tree wiring expects:

| depth `res − get_resolution(id)` | result |
|:--|:--|
| 0 | `SmallList{1,UInt64}` holding `id` itself |
| 1 | `SmallList{7,UInt64}` (6 entries for a pentagon) |
| ≥ 2 | `Vector{UInt64}` of `7^depth` (hexagon) or `p(depth)` (pentagon) ids |

Enumeration is digit-lexicographic depth-first, which *is* ascending id order
(a deeper level only touches lower bits), so the result needs no sort.

Throws [`InvalidZ7Error`](@ref) (`:descendant_res`) when
`res < get_resolution(id)` or `res > $MAX_RESOLUTION`.
"""
function cell_to_children(z7::UInt64, res::Integer)
    own = _geometry_checked(z7)
    own <= res <= MAX_RESOLUTION ||
        throw(InvalidZ7Error(:descendant_res, z7, _z7_int(res), own))
    target = Int(res)
    target == own && return Helpers.small_push(
        Helpers.empty_small_list(Val(1), zero(UInt64)), z7)
    target == own + 1 && return z7_children(z7)
    out = Vector{UInt64}(undef, _subtree_count(z7, own, target))
    _fill_descendants!(out, 1, z7, own, target)
    return out
end

cell_to_children(z7::Unsigned) = cell_to_children(UInt64(z7))
cell_to_children(z7::Unsigned, res::Integer) = cell_to_children(UInt64(z7), res)

"""
    _subtree_count(z7, own, target) -> Int

Number of descendants of `z7` (at resolution `own`) at resolution `target`:
`7^d` for a hexagon prefix, `p(d)` for a pentagon prefix (design Section 4.3).
"""
@inline function _subtree_count(z7::UInt64, own::Int, target::Int)
    depth = target - own
    active = Z7_PAD_MASK ⊻ _z7_tail_mask(own)
    return (z7 & active) == zero(UInt64) ? (@inbounds PENT_COUNT[depth+1]) :
           (@inbounds POW7[depth+1])
end

"""
    _fill_descendants!(out, pos, z, res, target) -> next_pos

Digit-lexicographic DFS writing every res-`target` descendant of `z` into
`out` starting at `pos`. Recursion depth ≤ 20; the pentagon test is redone at
each node because a prefix stops being a pentagon at its first nonzero digit.
"""
function _fill_descendants!(out::Vector{UInt64}, pos::Int, z::UInt64, res::Int, target::Int)
    if res == target
        @inbounds out[pos] = z
        return pos + 1
    end
    deleted = @inbounds Z7_DELETED_DIGIT[z7_base_cell(z)+1]
    active = Z7_PAD_MASK ⊻ _z7_tail_mask(res)
    pentagon = (z & active) == zero(UInt64)
    shift = _z7_shift(res + 1)
    cleared = z & ~(UInt64(7) << shift)
    for digit in 0:6
        pentagon && digit == deleted && continue
        pos = _fill_descendants!(out, pos, cleared | (UInt64(digit) << shift),
            res + 1, target)
    end
    return pos
end

# --- subtree border --------------------------------------------------------
#
# Which descendants of a cell touch a cell outside it, from the digit string
# alone — no neighbor query, no geometry.
#
# The mechanism is the refinement itself. A child sits at `chi_k * x + sigma(d)`
# (engine.jl's Horner step), and for any unit `u` the neighbor `child + u` has
# parent `x + delta` with `N(delta) <= 1`, because `|sigma(d) + u - sigma(d')|`
# is at most 3 and `sqrt(3) < sqrt(7)`. So a cell's children only ever touch
# children of that cell or of its six edge neighbors. Two things follow:
# `delta(0, u) = 0` for every unit, so the center child is enclosed by its own
# siblings and no border suffix contains a zero digit anywhere; and a border
# cell descends only from border cells, so the digit search prunes.
#
# What a border cell carries is *which* directions it is exposed in, and those
# sets are always contiguous arcs of the six units. A state is therefore the arc
# `{s, s+1, ..., s+L-1}` (mod 6), where `L` is exactly the number of border
# children the cell has — 2, 3 or 4, plus `L = 6` for the subtree root, which is
# exposed all round. `_border_step` is that arc's transition table.
#
# The table is *level-parity*-dependent, not fixed: refinement alternates
# chirality (`chi_is_c`, engine.jl), multiplying by `c = 3 + omega` at even
# levels and by `cbar` at odd ones, so the two tables are mirror images and a
# cell's border suffixes depend on the resolutions its digits sit at, not only
# on the depth below it. The digit alphabet itself never rotates — `SIGMA_J` is
# one table for every level, and the ALPHA_DEG turn lives in the frame (`P_R`,
# `THETA_DIR`), never in digit arithmetic.
#
# The census follows from the same table: an arc of length `L` has `L` border
# children, of which one is a 2-arc, one a 3-arc and the rest 4-arcs, so with
# `n_L` the count of each at a depth, `n_2' = n_3' = n_2 + n_3 + n_4` and
# `n_4' = n_3 + 2 * n_4`. Hence `n_4 - n_2` is invariant — 6 for a hexagon root,
# 5 for a pentagon, the discrete turning number of the closed rim — and
# `B(d) = 3 * B(d-1) + (n_4 - n_2)`, which is `_border_count`.

"""
    _border_step(state, digit, level) -> NTuple{2,Int}

State a `digit` child at absolute resolution `level` inherits from a cell in
`state`, or `(0, 0)` when that child is interior to the subtree. A state is
`(L, s)`: the arc of exposed directions `s, s+1, ..., s+L-1` (mod 6), in
`SIGMA_J`'s unit indices. `L == 6` is the subtree root, the one state with no
arc ends.
"""
@inline function _border_step(state::NTuple{2,Int}, digit::Int, level::Int)
    L, s = state
    (L == 0 || digit == 0) && return (0, 0)
    t = @inbounds SIGMA_J[digit]
    o = mod(t - s, 6)
    o < L || return (0, 0)
    if iseven(level)                         # chi_is_c(level): this level turns +ALPHA
        L < 6 && o == 0 && return (3, mod(t - 1, 6))
        L < 6 && o == L - 1 && return (2, mod(t - 2, 6))
        return (4, mod(t - 2, 6))
    else                                     # the mirror table
        L < 6 && o == L - 1 && return (3, mod(t - 1, 6))
        L < 6 && o == 0 && return (2, mod(t + 1, 6))
        return (4, mod(t - 1, 6))
    end
end

"""
    _border_count(z7, own, target) -> Int

Number of border descendants of `z7` (at resolution `own`) at resolution
`target`, from the census recurrence above: `3^(d+1) - 3` for a hexagon,
`5 * (3^d - 1) / 2` for a pentagon, `d = target - own`. Used to size the
output; it is a count, so a wrong one would cost a reallocation, not a result.
"""
@inline function _border_count(z7::UInt64, own::Int, target::Int)
    depth = target - own
    depth == 0 && return 1
    active = Z7_PAD_MASK ⊻ _z7_tail_mask(own)
    p3 = 3^depth
    return (z7 & active) == zero(UInt64) ? (5 * (p3 - 1)) ÷ 2 : 3 * p3 - 3
end

"""
    border_descendants(id, res) -> Vector{UInt64}

Descendants of `id` at resolution `res` that share an edge with a cell outside
`id`'s subtree — the subtree's rim — in ascending id order **[contract]**.
`res == get_resolution(id)` returns `[id]`, whose whole neighborhood is outside
its own subtree.

Decided from the Z7 digits alone (the block comment above), so the cost is
`O(result)`: the rim of a hexagon subtree holds `3^(d+1) - 3` cells at depth
`d`, against the `7^d` the subtree holds. Enumeration is digit-lexicographic
depth-first, which *is* ascending id order, so the result needs no sort.

Throws [`InvalidZ7Error`](@ref) (`:descendant_res`) when
`res < get_resolution(id)` or `res > $MAX_RESOLUTION`, as
[`cell_to_children`](@ref) does.
"""
function border_descendants(z7::UInt64, res::Integer)
    own = _geometry_checked(z7)
    own <= res <= MAX_RESOLUTION ||
        throw(InvalidZ7Error(:descendant_res, z7, _z7_int(res), own))
    target = Int(res)
    out = Vector{UInt64}()
    if target == own
        push!(out, z7)
        return out
    end
    sizehint!(out, _border_count(z7, own, target))
    return _fill_border!(out, z7, own, target, (6, 0))
end

border_descendants(z7::Unsigned, res::Integer) = border_descendants(UInt64(z7), res)

"""
    _fill_border!(out, z, res, target, state) -> out

Digit-lexicographic DFS over the border automaton, appending every res-`target`
border descendant of `z` to `out`. The pentagon test is redone at each node for
the same reason [`_fill_descendants!`](@ref) redoes it — a prefix stops being a
pentagon at its first nonzero digit — even though a border suffix's first digit
is already nonzero, so only the root node can ever be one.
"""
function _fill_border!(out::Vector{UInt64}, z::UInt64, res::Int, target::Int,
    state::NTuple{2,Int})
    deleted = @inbounds Z7_DELETED_DIGIT[z7_base_cell(z)+1]
    active = Z7_PAD_MASK ⊻ _z7_tail_mask(res)
    pentagon = (z & active) == zero(UInt64)
    shift = _z7_shift(res + 1)
    cleared = z & ~(UInt64(7) << shift)
    for digit in 0:6
        pentagon && digit == deleted && continue
        child_state = _border_step(state, digit, res + 1)
        child_state[1] == 0 && continue
        child = cleared | (UInt64(digit) << shift)
        if res + 1 == target
            push!(out, child)
        else
            _fill_border!(out, child, res + 1, target, child_state)
        end
    end
    return out
end

# --- introspection ---------------------------------------------------------

"""
    get_resolution(id) -> Int

Resolution `0:$MAX_RESOLUTION` of a cell id — [`z7_resolution`](@ref) with the
geometry-side validation **[contract]**. Throws [`InvalidZ7Error`](@ref) for a
structurally invalid id or one at resolution 20.
"""
@inline get_resolution(z7::UInt64) = _geometry_checked(z7)
@inline get_resolution(z7::Unsigned) = get_resolution(UInt64(z7))

"""
    is_valid_cell(id) -> Bool

`true` when `id` is a structurally valid Z7 index ([`is_valid_z7`](@ref)) at a
resolution that has geometry (`0:$MAX_RESOLUTION`) **[contract]**. Never
throws — resolution-20 indexes are valid Z7 but are not cells.
"""
@inline function is_valid_cell(z7::UInt64)
    is_valid_z7(z7) || return false
    return _z7_leading_resolution(z7) <= MAX_RESOLUTION
end
@inline is_valid_cell(z7::Unsigned) = is_valid_cell(UInt64(z7))

"""
    is_pentagon(id) -> Bool

`true` for the twelve pentagons of the id's resolution (every active digit
zero) — [`z7_is_pentagon`](@ref) with the geometry-side validation
**[contract]**.
"""
@inline function is_pentagon(z7::UInt64)
    res = _geometry_checked(z7)
    active = Z7_PAD_MASK ⊻ _z7_tail_mask(res)
    return (z7 & active) == zero(UInt64)
end
@inline is_pentagon(z7::Unsigned) = is_pentagon(UInt64(z7))

"""
    cell_to_z7(id) -> UInt64
    z7_to_cell(z7) -> UInt64

Validation-only identities: the Z7 `UInt64` **is** the cell id in this
implementation (design Section 3), so both directions just check the id and
return it. They exist for API parity with wrappers whose native ids differ
from Z7, and they are the documented place where a resolution-20 index is
rejected **[contract]**.
"""
@inline function cell_to_z7(cell::UInt64)
    _geometry_checked(cell)
    return cell
end
@inline cell_to_z7(cell::Unsigned) = cell_to_z7(UInt64(cell))

@doc (@doc cell_to_z7)
@inline function z7_to_cell(z7::UInt64)
    _geometry_checked(z7)
    return z7
end
@inline z7_to_cell(z7::Unsigned) = z7_to_cell(UInt64(z7))

# ---------------------------------------------------------------------------
# Deliberately absent: `serialize` / `deserialize` of cell ids beyond the Z7
# string and hex codecs in z7.jl. The sealed oracle exposes no observable
# serialization semantics that could be fitted from the recorded oracle
# vectors, and
# guessing one would violate the provenance rule in CLEANROOM.md. Callers that
# need a stable on-disk form use `z7_to_hex` / `z7_from_hex` (fixed 16-char,
# order-preserving) or the raw `UInt64`.
# ---------------------------------------------------------------------------
