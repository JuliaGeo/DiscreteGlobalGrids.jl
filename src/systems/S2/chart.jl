# ---------------------------------------------------------------------------
# S2's six charts map `[0, 1]²` to the sphere through quadratic `ST → UV`,
# cube-face `UV → XYZ`, and radial normalization. Chart lines map to great-circle
# arcs, so four corners exactly describe each cell. Face frames, projection, and
# Hilbert tables follow the published s2geometry conventions; no s2geometry code
# or library is used. Those conventions are `s2coords.h` and `s2cell_id.h`, as
# documented at https://s2geometry.io/devguide/s2cell_hierarchy.html.
#
# Index conventions:
#
# - `ix`, `iy` — 0-based face-local lattice coordinates in `0:nside-1`. `ix`
#   runs along the `s`/`u` axis and `iy` along the `t`/`v` axis; continuous
#   chart coordinates are `s = ix/nside`, `t = iy/nside`.
# - `face` — 0-based, `0:5`, in s2geometry's numbering: face `f` has its normal
#   on axis `f % 3`, positive for `f < 3` and negative for `f >= 3`.
# - Row-major id — 0-BASED: `face * nside² + iy * nside + ix`.
# - Hilbert id — 0-BASED: `face * nside² + hilbert_position`, which at
#   `nside = 2^level` is exactly the SCAFFOLD ORDINAL this system's canonical
#   `LevelIndex` carries (`face * 4^level + hilbert_position`).
#
# Both ids are 0-based; complete-grid position is `id + 1`. Argument order is
# `(ix, iy, face, nside)` and inverse codecs take the id first.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Face frames
#
# Right-handed face frames `(û, v̂, ŵ)`, with `û × v̂ == ŵ`. Face `f` has its
# outward normal on axis `f % 3`, positive for `f < 3`.
# ---------------------------------------------------------------------------

"""
    FACE_U_AXIS[face + 1] -> NTuple{3, Int}

Unit vector `û` of face `face` (0-based, `0:5`): the direction the chart's `u`
coordinate advances in. Together with [`FACE_V_AXIS`](@ref) and
[`FACE_NORMAL`](@ref) it forms a right-handed orthonormal frame, `û × v̂ == ŵ`.
"""
const FACE_U_AXIS = ((0, 1, 0), (-1, 0, 0), (-1, 0, 0), (0, 0, -1), (0, 0, -1), (0, 1, 0))

"""
    FACE_V_AXIS[face + 1] -> NTuple{3, Int}

Unit vector `v̂` of face `face` (0-based, `0:5`): the direction the chart's `v`
coordinate advances in. See [`FACE_U_AXIS`](@ref).
"""
const FACE_V_AXIS = ((0, 0, 1), (0, 0, 1), (0, -1, 0), (0, -1, 0), (1, 0, 0), (1, 0, 0))

"""
    FACE_NORMAL[face + 1] -> NTuple{3, Int}

Outward normal `ŵ` of face `face` (0-based, `0:5`) — the centre of the face,
and a signed coordinate axis. The six entries are the six signed axes, each
exactly once: `+x, +y, +z, -x, -y, -z`, which is s2geometry's face numbering.
"""
const FACE_NORMAL = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (-1, 0, 0), (0, -1, 0), (0, 0, -1))

"""
    face_uv_to_xyz(face, u, v) -> NTuple{3, Float64}

Lift face-local chart coordinates `(u, v) ∈ [-1, 1]²` onto the plane of cube
face `face` (0-based, `0:5`), returning an *unnormalised* point — the frame sum
`ŵ + u û + v v̂` of [`FACE_U_AXIS`](@ref) / [`FACE_V_AXIS`](@ref) /
[`FACE_NORMAL`](@ref).

The explicit switch preserves shared seam coordinates bit-identically by using
literal `±1.0` and exact negation rather than arithmetic accumulation.

!!! note "Signed zero is the one exception"
    Seam coordinates may differ only by `0.0` versus `-0.0`. They compare with
    `==` but not `isequal`; normalize signed zero before using coordinates as
    `Dict` or `Set` keys.
"""
function face_uv_to_xyz(face::Integer, u::Float64, v::Float64)
    face == 0 && return (1.0, u, v)
    face == 1 && return (-u, 1.0, v)
    face == 2 && return (-u, -v, 1.0)
    face == 3 && return (-1.0, -v, -u)
    face == 4 && return (v, -1.0, -u)
    return (v, u, -1.0)
end

"""
    xyz_to_face(p) -> Int

The 0-based cube face selected by the signed largest-magnitude component of `p`.

**Ties are broken toward the lower axis**, `x` before `y` before `z`, and a
component of zero counts as positive. This deterministically assigns cube-edge
and cube-corner points and is the face tie rule used by [`cellat`](@ref).
"""
function xyz_to_face(p)
    ax, ay, az = abs(p[1]), abs(p[2]), abs(p[3])
    axis = ax >= ay ? (ax >= az ? 1 : 3) : (ay >= az ? 2 : 3)
    return p[axis] >= 0 ? axis - 1 : axis + 2
end

# ---------------------------------------------------------------------------
# ST <-> UV: the quadratic transform
#
# s2geometry's quadratic projection places the dyadic `(s, t)` lattice on
# canonical S2 boundaries and gives an approximately 2.08× within-level area
# spread.
# ---------------------------------------------------------------------------

"""
    st_to_uv(s) -> Float64

s2geometry's `S2_QUADRATIC_PROJECTION`, mapping a chart coordinate
`s ∈ [0, 1]` to the cube-face coordinate `u ∈ [-1, 1]`:

```
u = (4s² - 1) / 3              for s >= 1/2
u = (1 - 4(1 - s)²) / 3        otherwise
```

It is strictly increasing, C¹ at `s = 1/2`, and exactly odd there in floating
point: `st_to_uv(1 - s) == -st_to_uv(s)`. Values at `0`, `1/2`, and `1` are
exactly `-1`, `0`, and `1`.

Exact oddness lets seam topology use signed-permutation integer lattice
arithmetic without reprojection. The transform is not equal-area.

See [`uv_to_st`](@ref) for the inverse.
"""
st_to_uv(s::Real) = (sf = Float64(s); sf >= 0.5 ? (4 * sf^2 - 1) / 3 : (1 - 4 * (1 - sf)^2) / 3)

"""
    uv_to_st(u) -> Float64

Inverse of [`st_to_uv`](@ref): cube-face coordinate `u ∈ [-1, 1]` back to chart
coordinate `s ∈ [0, 1]`.

```
s = sqrt(1 + 3u) / 2           for u >= 0
s = 1 - sqrt(1 - 3u) / 2       otherwise
```

Used by [`point_to_xyf`](@ref) in `xyz → face/uv → st → lattice` conversion.
"""
uv_to_st(u::Real) = (uf = Float64(u); uf >= 0 ? 0.5 * sqrt(1 + 3 * uf) : 1 - 0.5 * sqrt(1 - 3 * uf))

"""
    stf_to_point(s, t, face) -> GO.UnitSphericalPoint

Evaluate the S2 chart of `face` (0-based, `0:5`) at continuous face coordinates
`(s, t) ∈ [0, 1]²`, returning the corresponding point on the unit sphere:
quadratic [`st_to_uv`](@ref) per axis, then [`face_uv_to_xyz`](@ref), then
radial normalisation.

Defined for all real `s, t` in the unit square, independent of lattice depth.

The gnomonic projection maps chart lines to great-circle arcs, so a four-corner
ring exactly represents a cell. Normalization is well-conditioned because one
coordinate is always `±1` and the norm lies in `[1, √3]`.
"""
function stf_to_point(s::Real, t::Real, face::Integer)
    x, y, z = face_uv_to_xyz(Int(face), st_to_uv(s), st_to_uv(t))
    n = sqrt(x * x + y * y + z * z)
    return GO.UnitSphericalPoint(x / n, y / n, z / n)
end

"""
    point_to_xyf(p, nside) -> (ix, iy, face)

The face-local lattice cell of resolution `nside` containing the unit-sphere
point `p`: the chart's analytic inverse.

It applies `xyz → face/uv → st → lattice`, using [`xyz_to_face`](@ref),
[`uv_to_st`](@ref), and `floor`.

**Ties.** `floor` assigns chart cut lines to the higher-side cell; cube edges
and corners use [`xyz_to_face`](@ref)'s face rule.

Clamping corrects sub-ulp excursions outside `[0, 1]` on a face border.
"""
function point_to_xyf(p, nside::Integer)
    face = xyz_to_face(p)
    u_axis = FACE_U_AXIS[face + 1]
    v_axis = FACE_V_AXIS[face + 1]
    w_axis = FACE_NORMAL[face + 1]
    w = p[1] * w_axis[1] + p[2] * w_axis[2] + p[3] * w_axis[3]
    u = (p[1] * u_axis[1] + p[2] * u_axis[2] + p[3] * u_axis[3]) / w
    v = (p[1] * v_axis[1] + p[2] * v_axis[2] + p[3] * v_axis[3]) / w
    n = Int64(nside)
    ix = clamp(floor(Int64, uv_to_st(u) * n), Int64(0), n - 1)
    iy = clamp(floor(Int64, uv_to_st(v) * n), Int64(0), n - 1)
    return (Int(ix), Int(iy), face)
end

# ---------------------------------------------------------------------------
# Cell geometry
#
# The corner order is counter-clockwise in `(s, t)`. Monotone `st_to_uv` and the
# right-handed face frames preserve that orientation when viewed from outside.
# ---------------------------------------------------------------------------

"""
    cell_corners(ix, iy, face, nside) -> NTuple{4, GO.UnitSphericalPoint}

The four corners of cell `(ix, iy)` on `face` at resolution `nside`, as lattice
points evaluated through [`stf_to_point`](@ref).

Corners are emitted counter-clockwise as seen from outside the sphere,
in the order `((ix+1)/n, (iy+1)/n)`, `(ix/n, (iy+1)/n)`, `(ix/n, iy/n)`,
`((ix+1)/n, iy/n)`. This winding is required by spherical clipping: a clockwise
ring clips to empty rather than erroring.

Shared lattice corners are bit-identical except for signed zero; see
[`face_uv_to_xyz`](@ref). Since cell edges are great-circle arcs, no
densification is required.
"""
function cell_corners(ix::Integer, iy::Integer, face::Integer, nside::Integer)
    n = nside
    return (
        stf_to_point((ix + 1) / n, (iy + 1) / n, face),
        stf_to_point( ix      / n, (iy + 1) / n, face),
        stf_to_point( ix      / n,  iy      / n, face),
        stf_to_point((ix + 1) / n,  iy      / n, face),
    )
end

"""
    cell_center(ix, iy, face, nside) -> GO.UnitSphericalPoint

The center of cell `(ix, iy)` on `face` at resolution `nside`: the chart
evaluated at the lattice cell's midpoint `((ix + 0.5)/nside, (iy + 0.5)/nside)`.

This matches `S2CellId::ToPoint`. It is not the spherical quadrilateral's area
centroid.
"""
cell_center(ix::Integer, iy::Integer, face::Integer, nside::Integer) =
    stf_to_point((ix + 0.5) / nside, (iy + 0.5) / nside, face)

# ---------------------------------------------------------------------------
# Row-major order
#
# Lexicographic `nside × nside` layout with `ix` fastest, then `iy`, then face.
# It supports any `nside >= 1` but is not the canonical nested system id.
# ---------------------------------------------------------------------------

"""
    xyf_to_rowmajor(ix, iy, face, nside) -> Int64

0-based row-major id of cell `(ix, iy)` on `face` at resolution `nside`:
`face * nside² + iy * nside + ix`. Valid for any `nside >= 1`.

Inputs are assumed lattice-valid (`0 <= ix, iy < nside`, `0 <= face <= 5`) and
are not checked — garbage in, garbage out. [`rowmajor_to_xyf`](@ref) *is*
checked, being the direction a user-supplied id enters through.
"""
xyf_to_rowmajor(ix::Integer, iy::Integer, face::Integer, nside::Integer) =
    Int64(face) * Int64(nside)^2 + Int64(iy) * Int64(nside) + Int64(ix)

"""
    rowmajor_to_xyf(q, nside) -> (ix, iy, face)

Inverse of [`xyf_to_rowmajor`](@ref): face-local lattice coordinates and 0-based
face of the 0-based row-major id `q`. Valid for any `nside >= 1`; throws an
`ArgumentError` for an id outside `0:6nside²-1`.
"""
function rowmajor_to_xyf(q::Integer, nside::Integer)
    npface = Int64(nside)^2
    qi = Int64(q)
    0 <= qi < 6npface || throw(ArgumentError(
        "row-major id $qi out of range for nside=$nside (expected 0:$(6npface - 1))"))
    face, r = divrem(qi, npface)
    iy, ix = divrem(r, Int64(nside))
    return (Int(ix), Int(iy), Int(face))
end

# ---------------------------------------------------------------------------
# Hilbert order
#
# Canonical within-face Hilbert order. The tables are s2geometry's `kIJtoPos` / `kPosToIJ` /
# `kPosToOrientation`, indexed by a two-bit orientation state built from
# `SWAP_MASK` and `INVERT_MASK`. At each level the quadrant `ij = 2i + j`
# (`i` = the `ix`/s-axis bit, high; `j` = the `iy`/t-axis bit, low) is looked up
# to a position along the curve, the position contributes two bits, and the
# orientation is advanced by XOR. Odd faces start swapped, which is what chains
# the six per-face curves into a single closed loop over the whole sphere.
# ---------------------------------------------------------------------------

"""
    SWAP_MASK

Orientation bit meaning "the `i` and `j` axes are exchanged" in the Hilbert
state machine (s2geometry's `kSwapMask`). See [`IJ_TO_POS`](@ref).
"""
const SWAP_MASK = 1

"""
    INVERT_MASK

Orientation bit meaning "both axes are reflected" in the Hilbert state machine
(s2geometry's `kInvertMask`). See [`IJ_TO_POS`](@ref).
"""
const INVERT_MASK = 2

"""
    IJ_TO_POS[orientation + 1][ij + 1] -> Int

Quadrant `ij = 2i + j` (`i` the `ix`/s-axis bit, `j` the `iy`/t-axis bit) to its
0-based position along the Hilbert curve, under the current `orientation` in
`0:3` (`SWAP_MASK | INVERT_MASK`). s2geometry's `kIJtoPos`.
"""
const IJ_TO_POS = ((0, 1, 3, 2), (0, 3, 1, 2), (2, 3, 1, 0), (2, 1, 3, 0))

"""
    POS_TO_IJ[orientation + 1][pos + 1] -> Int

Inverse of [`IJ_TO_POS`](@ref): position along the curve back to the quadrant
`ij = 2i + j`. s2geometry's `kPosToIJ`.
"""
const POS_TO_IJ = ((0, 1, 3, 2), (0, 2, 3, 1), (3, 2, 0, 1), (3, 1, 0, 2))

"""
    POS_TO_ORIENTATION[pos + 1] -> Int

Orientation delta XORed into the state after descending into the sub-cell at
0-based curve position `pos`. s2geometry's `kPosToOrientation`.
"""
const POS_TO_ORIENTATION = (SWAP_MASK, 0, 0, SWAP_MASK | INVERT_MASK)

"""
    xyf_to_hilbert(ix, iy, face, nside) -> Int64

0-based Hilbert id of cell `(ix, iy)` on `face` at resolution `nside`, which
must be a power of two — throws `ArgumentError` otherwise.

The id is `face * nside² + hilbert_position`, i.e. the **scaffold ordinal**
`face * 4^level + hilbert_position` that this system's canonical
[`LevelIndex`](@ref) carries. Within a face the position is built two bits per
level, most significant first, by walking [`IJ_TO_POS`](@ref) and advancing the
orientation state through [`POS_TO_ORIENTATION`](@ref); odd faces start at
`SWAP_MASK`.

# Nesting

Dropping the low two position bits steps one level up:
`xyf_to_hilbert(ix >> 1, iy >> 1, face, nside >> 1)`'s position is
`xyf_to_hilbert(ix, iy, face, nside)`'s position `>> 2`. Because the face term
`face * nside²` divides through by 4 as well, the same statement holds of the
whole ordinal, and the system's parent/children arithmetic is exactly `÷ 4` and
`4p + k`, making subtrees contiguous.

!!! note "Alignment with native `s2_cellid` is intended, not yet verified"
    The tables and odd-face initial swap follow s2geometry conventions, but
    native S2 cell-id compatibility has not been fixture-verified. No native
    `s2_cellid` reindexing scheme is provided.
"""
function xyf_to_hilbert(ix::Integer, iy::Integer, face::Integer, nside::Integer)
    ispow2(nside) || throw(ArgumentError(
        "the Hilbert index is only defined for nside = 2^k, got nside=$nside; \
         use `xyf_to_rowmajor` for arbitrary nside"))
    (0 <= ix < nside && 0 <= iy < nside) || throw(ArgumentError(
        "lattice coordinates ($ix, $iy) out of range for nside=$nside (expected 0:$(nside - 1))"))
    0 <= face <= 5 || throw(ArgumentError("face $face out of range (expected 0:5)"))
    k = trailing_zeros(Int64(nside))                # nside == 2^k
    # Odd faces start swapped; that is what joins the six per-face curves into
    # one continuous loop over the sphere.
    orientation = isodd(face) ? SWAP_MASK : 0
    pos = Int64(0)
    for b in (k - 1):-1:0
        ij = 2 * ((Int64(ix) >> b) & 1) + ((Int64(iy) >> b) & 1)
        p = IJ_TO_POS[orientation + 1][ij + 1]
        pos = (pos << 2) | p
        orientation ⊻= POS_TO_ORIENTATION[p + 1]
    end
    return Int64(face) * Int64(nside)^2 + pos
end

"""
    hilbert_to_xyf(h, nside) -> (ix, iy, face)

Inverse of [`xyf_to_hilbert`](@ref): face-local lattice coordinates and 0-based
face of the 0-based Hilbert id (scaffold ordinal) `h`. `nside` must be a power
of two — throws `ArgumentError` otherwise, and for an id outside `0:6nside²-1`.

Splits off the face and consumes two position bits per level through
[`POS_TO_IJ`](@ref), rebuilding `ix` and `iy`.

!!! note "Alignment with native `s2_cellid` is intended, not yet verified"
    See [`xyf_to_hilbert`](@ref) — the same caveat applies verbatim.
"""
function hilbert_to_xyf(h::Integer, nside::Integer)
    ispow2(nside) || throw(ArgumentError(
        "the Hilbert index is only defined for nside = 2^k, got nside=$nside; \
         use `rowmajor_to_xyf` for arbitrary nside"))
    npface = Int64(nside)^2
    hid = Int64(h)
    0 <= hid < 6npface || throw(ArgumentError(
        "Hilbert id $hid out of range for nside=$nside (expected 0:$(6npface - 1))"))
    face, pos = divrem(hid, npface)
    k = trailing_zeros(Int64(nside))
    orientation = isodd(face) ? SWAP_MASK : 0
    ix = Int64(0)
    iy = Int64(0)
    for b in (k - 1):-1:0
        p = Int((pos >> (2b)) & 3)
        ij = POS_TO_IJ[orientation + 1][p + 1]
        ix = (ix << 1) | ((ij >> 1) & 1)
        iy = (iy << 1) | (ij & 1)
        orientation ⊻= POS_TO_ORIENTATION[p + 1]
    end
    return (Int(ix), Int(iy), Int(face))
end
