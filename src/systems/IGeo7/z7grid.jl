# Z7 geometry and dense indexing. Encoding maps digits through the Eisenstein
# lattice and Snyder chart; decoding reverses that pipeline and verifies each
# candidate by strict re-encoding.

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

Return the physical level-`res` lattice point in the base's development frame.
Digits are accumulated by Horner evaluation, the deleted pentagon sector is
collapsed, and the result is wrapped into the `[0°, 300°)` cone. An all-zero
pentagon prefix returns `(0, 0)`.
"""
@inline function _encode_lattice(z::UInt64, base::Int, res::Int)
    (a, b, _) = _encode_lattice_rot(z, base, res)
    return (a, b)
end

"""
    _encode_lattice_rot(z, base, res) -> (a, b, g)

[`_encode_lattice`](@ref) with its own frame rotation exposed: `g` is the net
number of 60° [`unitmul`](@ref) steps the collapse and the cone wrap applied to
the Horner accumulator, so the returned representative sits `g` sixths of a turn
counterclockwise of where the raw digit sum put it.

That number is what relates the two descriptions of a cell's six neighbour
directions. A digit's own direction is [`SIGMA_J`](@ref)'s, fixed; this
function's `(a, b)` is in the rotated frame, so a neighbour reached by adding
digit `d` sits at unit index `mod(SIGMA_J[d] + g, 6)`. That is the whole bridge
between the GBT digit step of `gbt.jl` and this file's counterclockwise unit
order — see `_cell_neighbors_ccw`. An all-zero prefix (a pentagon) never rotates:
`g == 0`.
"""
function _encode_lattice_rot(z::UInt64, base::Int, res::Int)
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
    m == 0 && return (Int64(0), Int64(0), 0)
    # collapse at level `m`, the first nonzero digit: the deleted digit's 60°
    # sector is excised, so a subtree CCW-past it rotates back by 60° [fitted].
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
        shift -= 1
    elseif psit >= DEV_CONE_DEG
        (a, b) = unitmul(a, b, 1)
        shift += 1
    end
    return (a, b, shift)
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

Center of `z7` as `(longitude, latitude)` in degrees. Pentagon centers are
icosahedron vertices; hexagon centers are obtained from the lattice point via
the slot map and Snyder inverse. `orientation` rotates the result.

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

Boundary of `z7` as unit-sphere `xyz` tuples: six corners for a hexagon and
five for a pentagon. Rings wind counterclockwise (seen from outside the sphere)
under the identity orientation; `closed_ring=true` repeats the first corner.
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

Decode development-frame position `u` at resolution `res`. The method rounds
to the nearest lattice point, canonicalizes the cone cut, tries the possible
pentagon-collapse rotations, and accepts only a valid digit string that
strictly re-encodes to the rounded point. Return `ok=false` if this base has no
owner. With `tie=true`, equally near lattice points within [`TIE_EPS`](@ref)
are tried in ascending-margin order.
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

Return the Z7 cell containing `(lon, lat)` in degrees at resolution `res`.
The containing face's corner bases are tried nearest-first, with the rounding
tie fallback applied only after ordinary decoding fails.

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

# Edge neighbors are the six Eisenstein-unit steps from the physical lattice
# point. Cone wrapping uses the center branch; pentagon duplicates are removed.

"""
    _cell_neighbors(z7) -> SmallList{6,UInt64}

Canonical ids of the cells sharing an edge with `z7`, ascending: 6 for a
hexagon and 5 for a pentagon. `gbt.jl`'s digit step answers this; the sort is
over [`_cell_neighbors_ccw`](@ref)'s rotational order.

Throws [`InvalidZ7Error`](@ref) for invalid ids and for resolution-20 ids
(valid for prefix arithmetic, no geometry — hence no neighbors).
"""
function _cell_neighbors(z7::UInt64)
    return Helpers.small_sort(_cell_neighbors_ccw(z7))
end

"""
    _cell_neighbors_ccw_geometric(z7) -> SmallList{6,UInt64}

The neighbours by geometry: step one Eisenstein unit from the cell's physical
lattice point in each of the six unit directions, project the result through
`dev_to_xyz`, and re-decode it with [`_xyz_to_z7`](@ref). Counterclockwise from
the development frame's `+1` direction, by construction — the loop *is* the
unit order.

Nothing on the hot path calls this — a floating-point round trip with a
three-base search in it, per neighbour, is roughly forty times the cost of
`gbt.jl`'s digit arithmetic. It stays because it is the **differential oracle**:
it answers adjacency from this package's own oracle-validated lattice and
decoder, sharing no reasoning with the ported digit kernel, so agreement between
the two is real evidence rather than a restatement. Testset `9b` of
`test/systems/IGeo7/runtests.jl` is that comparison.
"""
function _cell_neighbors_ccw_geometric(z7::UInt64)
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

# ---------------------------------------------------------------------------
# Dense order is ascending encoded identifier, equivalent to lexicographic
# `(base, digits)` order. Hexagon subtrees contain `7^d` cells and pentagon
# subtrees contain `(5*7^d + 1)/6`; rank and inverse-rank are `O(res)`.
# ---------------------------------------------------------------------------

"`POW7[d+1] = 7^d` — subtree size of a *hexagon* prefix at depth `d`."
const POW7 = ntuple(i -> Int64(7)^(i - 1), MAX_RESOLUTION + 1)

"""
    PENT_COUNT

`PENT_COUNT[d+1] = p(d) = (5·7^d + 1) / 6` — subtree size of a *pentagon*
prefix at depth `d`. Recurrence
`p(d) = 6·7^(d-1) + p(d-1)`: a pentagon has six children, five of them
hexagons (full `7^(d-1)` subtrees) and one — digit 0 — a pentagon again.
"""
const PENT_COUNT = ntuple(i -> (5 * Int64(7)^(i - 1) + 1) ÷ 6, MAX_RESOLUTION + 1)

"`NUM_CELLS[r+1] = 12·p(r) = 10·7^r + 2`."
const NUM_CELLS = ntuple(i -> 10 * Int64(7)^(i - 1) + 2, MAX_RESOLUTION + 1)

"""
    num_cells(res) -> Int64

Number of cells covering the whole world at resolution `res`:
`10·7^res + 2` — twelve
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

One-based rank of `id` among cells at its resolution in ascending encoded-id
order. The mixed-radix prefix walk is `O(res)` and allocation-free.

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

Return the `index`-th cell at resolution `res`, inverse to
[`cell_to_index`](@ref). The greedy mixed-radix walk is `O(res)` and
allocation-free.

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

All descendants of `id` at resolution `res`, in ascending identifier order.
The return type depends on depth:

| depth `res − get_resolution(id)` | result |
|:--|:--|
| 0 | `SmallList{1,UInt64}` holding `id` itself |
| 1 | `SmallList{7,UInt64}` (6 entries for a pentagon) |
| ≥ 2 | `Vector{UInt64}` of `7^depth` (hexagon) or `p(depth)` (pentagon) ids |

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

# The subtree-border automaton tracks each cell's contiguous arc of exposed
# lattice directions. Its mirrored transition table follows level chirality,
# prunes interior children, and enumerates only the `O(rim)` result. Pruning is
# sound because a child's neighbors descend only from its parent or the parent's
# six edge neighbors (`|sigma(d) + u - sigma(d')| <= 3 < sqrt 7`). Census: an arc
# of length `L` has exactly `L` border children — one 2-arc, one 3-arc, the rest
# 4-arcs — so `n_4 - n_2` is invariant (6 for a hexagon root, 5 for a pentagon)
# and `B(d) = 3*B(d-1) + (n_4 - n_2)`, which is `_border_count`.

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
the subtree, in ascending identifier order. `res == get_resolution(id)` returns
`[id]`. The digit automaton runs in `O(result)`.

Throws [`InvalidZ7Error`](@ref) (`:descendant_res`) when
`res < get_resolution(id)` or `res > $MAX_RESOLUTION`, as
[`cell_to_children`](@ref) does.
"""
function border_descendants(z7::UInt64, res::Integer)
    own = _geometry_checked(z7)
    own <= res <= MAX_RESOLUTION ||
        throw(InvalidZ7Error(:descendant_res, z7, _z7_int(res), own))
    # Through the guard, not a comprehension: a comprehension over a `HasLength`
    # walk preallocates from the count and would leave an `undef` tail if the two
    # ever disagreed. `collect_subtree` is the one place that comparison lives.
    return [c.id for c in DGG.collect_subtree(Z7RimEngine(z7, own, Int(res)))]
end

border_descendants(z7::Unsigned, res::Integer) = border_descendants(UInt64(z7), res)

# The automaton as a resumable walk. `next` is the digit to try when the frame is
# next on top, so a frame retires at `next > 6`; `full` marks a branch the
# automaton pruned, under which every cell is interior. Resolution is not stored
# — frame `k` sits at `root + k - 1` — and neither is the id: one scalar `z`
# carries the current path, because descending overwrites exactly one digit slot
# and slots below the target keep the 7 padding that never gets written.
struct Z7Frame
    L::Int8
    s::Int8
    next::Int8
    full::Bool
end

const _Z7_STACK_CAP = MAX_RESOLUTION + 1
const Z7Stack = Helpers.SmallList{_Z7_STACK_CAP,Z7Frame}

@inline _z7_empty_stack() = Helpers.empty_small_list(Val(_Z7_STACK_CAP),
    Z7Frame(Int8(0), Int8(0), Int8(0), false))

"""
    Z7Walk(z, stack)

The walk state: the current path's id and the frame stack, both inline.
"""
struct Z7Walk
    z::UInt64
    stack::Z7Stack
end

@inline _z7_root_walk(e) = Z7Walk(e.z,
    Helpers.small_push(_z7_empty_stack(), Z7Frame(Int8(6), Int8(0), Int8(0), false)))

@inline _z7_bump(f::Z7Frame) = Z7Frame(f.L, f.s, f.next + Int8(1), f.full)

# The pentagon test is redone at each node for the same reason
# [`_fill_descendants!`](@ref) redoes it — a prefix stops being a pentagon at its
# first nonzero digit. `res` masks off the stale digits the scalar `z` still
# carries from a deeper branch, so it reads the same value the recursion did.
@inline function _z7_node(z::UInt64, res::Int)
    active = Z7_PAD_MASK ⊻ _z7_tail_mask(res)
    return (@inbounds Z7_DELETED_DIGIT[z7_base_cell(z)+1]),
        (z & active) == zero(UInt64)
end

"""
    Z7RimEngine(z, res, target)

Digit-lexicographic walk over the border automaton, yielding every res-`target`
rim descendant of `z` as a `Z7Cell`, ascending, in `O(depth)` memory.
"""
struct Z7RimEngine
    z::UInt64
    res::Int
    target::Int
end

Base.eltype(::Type{Z7RimEngine}) = Z7Cell
Base.IteratorSize(::Type{Z7RimEngine}) = Base.HasLength()
Base.length(e::Z7RimEngine) = Int(_border_count(e.z, e.res, e.target))

function Base.iterate(e::Z7RimEngine)
    e.target == e.res && return (Z7Cell(e.z), Z7Walk(e.z, _z7_empty_stack()))
    return iterate(e, _z7_root_walk(e))
end

Base.iterate(e::Z7RimEngine, w::Z7Walk) = _z7_rim_advance(e.res, e.target, w)

# The walk itself, taking the two numbers it reads rather than an engine: the
# seeded arc engine below runs the same automaton from a different root frame,
# and this is the whole of what the two share.
function _z7_rim_advance(res0::Int, target::Int, w::Z7Walk)
    z = w.z
    st = w.stack
    while !isempty(st)
        k = length(st)
        f = @inbounds st[k]
        if f.next > Int8(6)
            st = Helpers.small_pop(st)
            continue
        end
        digit = Int(f.next)
        st = Helpers.small_setlast(st, _z7_bump(f))
        res = res0 + k - 1
        deleted, pentagon = _z7_node(z, res)
        pentagon && digit == deleted && continue
        child = _border_step((Int(f.L), Int(f.s)), digit, res + 1)
        child[1] == 0 && continue
        shift = _z7_shift(res + 1)
        z = (z & ~(UInt64(7) << shift)) | (UInt64(digit) << shift)
        res + 1 == target && return (Z7Cell(z), Z7Walk(z, st))
        st = Helpers.small_push(st,
            Z7Frame(Int8(child[1]), Int8(child[2]), Int8(0), false))
    end
    return nothing
end

"""
    Z7ArcEngine(z, res, target, L, s)

The same automaton entered at an arbitrary arc: the res-`target` descendants of
`z` reachable along the exposed directions `s, s+1, …, s+L-1 (mod 6)`, ascending,
in `O(depth)` memory. [`Z7RimEngine`](@ref) is this with `(L, s) == (6, 0)`,
which is the one state the census recurrence describes — so this engine declares
`SizeUnknown()` and no `length`.

The seed cell's own pentagon deletion needs no field here, unlike H3's: `_z7_node`
re-derives it from `z` at every node, and at the root frame that is exactly "is
`z` a pentagon at resolution `res`". A calibrated arc IS seeded at cells that are
pentagons, so this is load-bearing rather than incidental.
"""
struct Z7ArcEngine
    z::UInt64
    res::Int
    target::Int
    L::Int8
    s::Int8
end

Base.eltype(::Type{Z7ArcEngine}) = Z7Cell
Base.IteratorSize(::Type{Z7ArcEngine}) = Base.SizeUnknown()

Base.iterate(e::Z7ArcEngine) = iterate(e, Z7Walk(e.z,
    Helpers.small_push(_z7_empty_stack(), Z7Frame(e.L, e.s, Int8(0), false))))

Base.iterate(e::Z7ArcEngine, w::Z7Walk) = _z7_rim_advance(e.res, e.target, w)

"""
    Z7InteriorEngine(z, res, target)

The complement, off the same automaton: where the rim walk prunes, the whole
branch below is interior, so this one descends it in full instead of dropping
it. No rim membership is ever tested or stored.
"""
struct Z7InteriorEngine
    z::UInt64
    res::Int
    target::Int
end

Base.eltype(::Type{Z7InteriorEngine}) = Z7Cell
Base.IteratorSize(::Type{Z7InteriorEngine}) = Base.HasLength()
function Base.length(e::Z7InteriorEngine)
    e.target == e.res && return 0
    return Int(_descendant_count(e.z, e.res, e.target)) -
           Int(_border_count(e.z, e.res, e.target))
end

function Base.iterate(e::Z7InteriorEngine)
    e.target == e.res && return nothing
    return iterate(e, _z7_root_walk(e))
end

function Base.iterate(e::Z7InteriorEngine, w::Z7Walk)
    z = w.z
    st = w.stack
    while !isempty(st)
        k = length(st)
        f = @inbounds st[k]
        if f.next > Int8(6)
            st = Helpers.small_pop(st)
            continue
        end
        digit = Int(f.next)
        st = Helpers.small_setlast(st, _z7_bump(f))
        res = e.res + k - 1
        deleted, pentagon = _z7_node(z, res)
        pentagon && digit == deleted && continue
        child = f.full ? (0, 0) : _border_step((Int(f.L), Int(f.s)), digit, res + 1)
        shift = _z7_shift(res + 1)
        z = (z & ~(UInt64(7) << shift)) | (UInt64(digit) << shift)
        if child[1] != 0
            res + 1 == e.target && continue     # a rim cell: not ours
            st = Helpers.small_push(st,
                Z7Frame(Int8(child[1]), Int8(child[2]), Int8(0), false))
            continue
        end
        res + 1 == e.target && return (Z7Cell(z), Z7Walk(z, st))
        st = Helpers.small_push(st, Z7Frame(Int8(0), Int8(0), Int8(0), true))
    end
    return nothing
end

"""
    _descendant_count(z7, own, target) -> Int

Number of resolution-`target` descendants of `z7`: `7^d` for a hexagon, and the
all-zero chain plus five branches, `1 + 5*(7^d - 1)/6`, for a pentagon.
"""
@inline function _descendant_count(z7::UInt64, own::Int, target::Int)
    depth = target - own
    active = Z7_PAD_MASK ⊻ _z7_tail_mask(own)
    p7 = 7^depth
    return (z7 & active) == zero(UInt64) ? 1 + 5 * (p7 - 1) ÷ 6 : p7
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

# Stable storage uses the raw `UInt64` or the fixed-width hexadecimal codec.
