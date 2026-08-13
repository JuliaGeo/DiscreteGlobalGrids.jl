# icosahedron.jl — the IGEO7 icosahedron: sphere helpers, the standard ISEA
# vertex table, neighbor rings, per-base counterclockwise rings and the
# `Orientation` mechanism (spec/design.md sections 2, 5 and 12).
#
# Provenance
#   * vertex layout (standard ISEA placement, closed form): golden-rectangle
#     icosahedron rotated +11.25° about z, base-id assignment per
#     spec/isea-projection-spec.md §4.2 / table T1 and spec/z7-paper-spec.md
#     §1.2 [published]; locked against the recorded res-0 oracle centers
#     (whose own print noise is ≤ 5.92e-9 deg) and cell-level against the
#     res-1..5 oracle center dumps (spec/igeo7-geometry-diagnosis.md §2).
#   * authalic radius -> spec/interface-contract.md  [contract]
#   * adjacency / antipodal structure
#                     -> spec/aperture7-indexing-spec.md section 5.1
#   * per-base reference edge (dev-frame angle 0)
#                     -> fitted, spec/igeo7-geometry-diagnosis.md §4 (the
#                        A-gauge); validated on all 196,080 oracle cell
#                        centers res 1-5
#
# Everything here lives in the *grid frame* — the standard ISEA placement
# itself: the `Orientation` rotation is applied to world points on the way in
# and inverted on the way out, so the tables below never move (design
# section 5). The identity orientation therefore returns coordinates in the
# standard ISEA frame directly — the frame the recorded oracle dumps are in.
#
# The planar constants of the chart (`R_EA`, `L_PLANE`, ...) belong to
# snyder.jl; this file is pure spherical geometry.

# ---------------------------------------------------------------------------
# Sphere helpers (degrees on the boundary, radians internally)
# ---------------------------------------------------------------------------

"""
    lonlat_to_xyz(lon, lat) -> NTuple{3,Float64}

Unit vector of a spherical (lon, lat) pair in degrees, right-handed with
`+z` at the north pole and `+x` at (0, 0). Lon/lat are plain spherical
coordinates — IGEO7 applies no geodetic/authalic latitude conversion
**[contract; confirmed against every oracle dump]**.
"""
function lonlat_to_xyz(lon::Real, lat::Real)
    cl = cosd(lat)                 # hoisted: LLVM does not CSE the two cosd(lat)
    return (cl * cosd(lon), cl * sind(lon), sind(lat))
end

"""
    xyz_to_lonlat(v) -> (lon, lat)

Inverse of [`lonlat_to_xyz`](@ref); `lon ∈ (-180, 180]`, `lat ∈ [-90, 90]`.
"""
xyz_to_lonlat(v::NTuple{3,Float64}) =
    (atand(v[2], v[1]), atand(v[3], hypot(v[1], v[2])))

vdot(a::NTuple{3,Float64}, b::NTuple{3,Float64}) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
vcross(a::NTuple{3,Float64}, b::NTuple{3,Float64}) =
    (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])
vnorm(a::NTuple{3,Float64}) = sqrt(vdot(a, a))
vscale(a::NTuple{3,Float64}, s::Real) = (a[1] * s, a[2] * s, a[3] * s)
vadd(a::NTuple{3,Float64}, b::NTuple{3,Float64}) = (a[1] + b[1], a[2] + b[2], a[3] + b[3])
vsub(a::NTuple{3,Float64}, b::NTuple{3,Float64}) = (a[1] - b[1], a[2] - b[2], a[3] - b[3])
vnormalize(a::NTuple{3,Float64}) = vscale(a, 1 / vnorm(a))

"""
    angdist(a, b) -> Float64

Angular distance in degrees between two unit vectors, via `atan2` of the
cross/dot pair (accurate at both ends of the range).
"""
angdist(a::NTuple{3,Float64}, b::NTuple{3,Float64}) = atand(vnorm(vcross(a, b)), vdot(a, b))

# ---------------------------------------------------------------------------
# Exact constants
# ---------------------------------------------------------------------------

"""
    R_AUTHALIC

WGS84 authalic Earth radius in metres **[contract]**. It enters the grid only
through areas: `hex_area(r) = 4πR²/(10·7^r)`, `pentagon = 5/6 · hex`; the
res-0 cell area `4πR²/12 = 4.250546847700739e13 m²` matches the recorded
oracle areas exactly.
"""
const R_AUTHALIC = 6371007.180918475

"`cos` of the icosahedron edge arc, `1/√5` (arc `acosd(1/√5) = 63.4349...°`)."
const ADJ_DOT = 1 / sqrt(5.0)

"""
    NBASE

Number of base cells / icosahedron vertices. Same fact as z7.jl's
`Z7_NUM_BASES`, deliberately restated there because `z7.jl` must stay
standalone-includable (pure integer layer, see its file header and
`test/test_z7.jl`); the two are pinned equal by `test/test_icosahedron.jl`
("cross-file constant consistency").
"""
const NBASE = 12

# The per-base deleted (pentagon-chain) digit table has ONE definition,
# `Z7_DELETED_DIGIT` in z7.jl, reached through `z7_deleted_digit(base)`; it is
# a digit-alphabet fact, not spherical geometry, and z7.jl owns it.

"""
    ISEA_LON0

Longitude of the standard ISEA IGEO7 vertex 0, `11.25` degrees
(`= π/16` rad, the value of PROJ's `ISEA_STD_LONG`)
**[spec/isea-projection-spec.md §4.1]**.
"""
const ISEA_LON0 = 11.25

"""
    ISEA_LAT_HI

Latitude of the standard ISEA IGEO7 vertex 0, `atand(φ) =
58.282525588538995` degrees (PROJ's `ISEA_STD_LAT`)
**[spec/isea-projection-spec.md §4.1]**. The 12 vertices sit at latitudes
`±ISEA_LAT_HI`, `±(90 - ISEA_LAT_HI)` and `0`.
"""
const ISEA_LAT_HI = atand((1 + sqrt(5.0)) / 2)

# ---------------------------------------------------------------------------
# Vertex table and adjacency (standard ISEA placement, closed form)
# ---------------------------------------------------------------------------

const _PHI = (1 + sqrt(5.0)) / 2

# Prime-frame golden-rectangle icosahedron; base-id order per
# spec/isea-projection-spec.md §4.2 / table T1.
const _PRIME = (
    (1.0, 0.0, _PHI), (-1.0, 0.0, _PHI), (0.0, -_PHI, 1.0),
    (_PHI, -1.0, 0.0), (_PHI, 1.0, 0.0), (0.0, _PHI, 1.0),
    (-_PHI, -1.0, 0.0), (0.0, -_PHI, -1.0), (1.0, 0.0, -_PHI),
    (0.0, _PHI, -1.0), (-_PHI, 1.0, 0.0), (-1.0, 0.0, -_PHI),
)

function _make_vertices()
    s = sqrt(_PHI + 2)
    c, sn = cosd(ISEA_LON0), sind(ISEA_LON0)
    ntuple(NBASE) do i
        x, y, z = _PRIME[i]
        x, y, z = x / s, y / s, z / s
        (c * x - sn * y, sn * x + c * y, z)
    end
end

"""
    VERTICES

The 12 icosahedron vertices (= res-0 pentagon centers) as unit vectors in the
grid frame, indexed by `base + 1` — the standard ISEA placement: vertex 0
at `(11.25°E, 58.28252559°N)`, azimuth 0 **[spec/isea-projection-spec.md
§4.2 T1; locked against the recorded res-0 oracle centers]**.
"""
const VERTICES = _make_vertices()

"""
    vertex(base) -> NTuple{3,Float64}

Unit vector of base cell `base`'s icosahedron vertex (grid frame).
"""
function vertex(base::Integer)
    0 <= base <= 11 || throw(ArgumentError("base cell must be in 0:11, got $base"))
    return @inbounds VERTICES[Int(base)+1]
end

# Adjacency: two vertices are neighbors iff their dot product is 1/√5.
function _make_neighbors()
    ntuple(NBASE) do i
        b = i - 1
        ns = [j for j in 0:11 if j != b && abs(vdot(VERTICES[i], VERTICES[j+1]) - ADJ_DOT) < 1e-9]
        @assert length(ns) == 5
        ntuple(k -> ns[k], 5)          # ascending base index
    end
end

"""
    NEIGHBORS

Neighbor ring of each base cell in *ascending base order*, indexed by
`base + 1` (matches the adjacency table of
`spec/aperture7-indexing-spec.md` section 5.1). For the geometric ring —
counterclockwise from the reference edge — use [`NBRS_CCW`](@ref).
"""
const NEIGHBORS = _make_neighbors()

"""
    nearest_vertex(p) -> Int

Base cell whose icosahedron vertex is closest to the unit vector `p`
(grid frame). Allocation-free: 12 dot products and an argmax.
"""
function nearest_vertex(p::NTuple{3,Float64})
    best = 0
    bestd = -Inf
    @inbounds for b in 0:11
        d = vdot(p, VERTICES[b+1])
        if d > bestd
            bestd = d
            best = b
        end
    end
    return best
end

nearest_vertex(lon::Real, lat::Real) = nearest_vertex(lonlat_to_xyz(lon, lat))

# ---------------------------------------------------------------------------
# Per-base counterclockwise rings (the dev-frame gauge)
# ---------------------------------------------------------------------------

"""
    REFERENCE_EDGE

Per-base reference edge: the neighbor base that sits at dev-frame angle 0 in
each base's vertex chart **[fitted; see `spec/igeo7-geometry-diagnosis.md` §4
(the A-gauge, whose pentagon-collapse rule is grid.jl's verbatim); validated
on all 196,080 oracle cell centers res 1–5]**. One global
chirality/digit map fits all 12 bases with these references — there are no
mirrored bases.
"""
const REFERENCE_EDGE = (1, 10, 6, 7, 8, 9, 11, 11, 11, 11, 11, 8)

function _make_ccw_ring(b::Int, ref::Int)
    Vb = VERTICES[b+1]
    tangent(n) = vnormalize(vsub(VERTICES[n+1], vscale(Vb, vdot(VERTICES[n+1], Vb))))
    ut = tangent(ref)
    wt = vcross(Vb, ut)
    angs = map(NEIGHBORS[b+1]) do n
        n == ref && return 0.0
        t = tangent(n)
        mod(atand(vdot(t, wt), vdot(t, ut)), 360.0)
    end
    order = sortperm(collect(angs))
    # the five neighbors form a perfect 72-degree comb starting at `ref`
    for i in 1:5
        @assert abs(angs[order[i]] - 72.0 * (i - 1)) < 1e-8
    end
    return ntuple(i -> NEIGHBORS[b+1][order[i]], 5)
end

"""
    NBRS_CCW

Per base, the five neighbors counterclockwise (seen from outside) starting
at [`REFERENCE_EDGE`](@ref): `NBRS_CCW[b+1][j+1]` sits at azimuth `72j`° —
dev-frame slot `j` (snyder.jl's slot tables key on this).
"""
const NBRS_CCW = ntuple(i -> _make_ccw_ring(i - 1, REFERENCE_EDGE[i]), NBASE)

# ---------------------------------------------------------------------------
# Orientation (design section 5)
# ---------------------------------------------------------------------------

const _IDENTITY_R = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)

"""
    Orientation(R)

Rigid rotation of the sphere taking **world** coordinates to the **grid
frame** in which the tables of this file live, stored row-major as an
`NTuple{9,Float64}` (no StaticArrays dependency: the package core is
stdlib-only).

Decode applies `R` to the input point ([`to_grid`](@ref)); encode applies
`R'` to the output ([`from_grid`](@ref)). Everything else — vertex table,
charts, lattices — is orientation-independent.

`R` is the *only* argument: the `identity` field beside it is a cached
`R == I` predicate, not an independent fact, so it is computed here and there
is no constructor that takes it (see the note on the struct).

[`ORIENT_IDENTITY`](@ref) is the identity: the standard ISEA IGEO7
placement, the default for every geometric entry point. Pass a custom
`Orientation` to run the same grid in a rotated placement.
"""
struct Orientation
    R::NTuple{9,Float64}
    identity::Bool
    # `identity` is a fast path, not a degree of freedom: [`to_grid`](@ref) and
    # [`from_grid`](@ref) *skip the rotation entirely* when it is set, so an
    # `Orientation(R, true)` built beside a non-identity `R` would return every
    # point unrotated — a silently wrong grid placement, not an error. The flag
    # is therefore derived here and nowhere else; this inner constructor is the
    # only way to make an `Orientation`, including for `ORIENT_IDENTITY` below.
    Orientation(R::NTuple{9,Float64}) = new(R, R == _IDENTITY_R)
end

"""
    ORIENT_IDENTITY

Identity orientation — the standard ISEA IGEO7 placement of
[`VERTICES`](@ref) (vertex 0 at 11.25°E, 58.28252559°N). Default for every
geometric entry point; boundary rings wind counterclockwise seen from
outside under it.
"""
const ORIENT_IDENTITY = Orientation(_IDENTITY_R)

"""
    to_grid(orientation, p) -> NTuple{3,Float64}

Rotate a world-frame unit vector into the grid frame (`R * p`). Returns `p`
itself for the identity orientation.
"""
function to_grid(o::Orientation, p::NTuple{3,Float64})
    o.identity && return p
    R = o.R
    return (R[1] * p[1] + R[2] * p[2] + R[3] * p[3],
        R[4] * p[1] + R[5] * p[2] + R[6] * p[3],
        R[7] * p[1] + R[8] * p[2] + R[9] * p[3])
end

"""
    from_grid(orientation, p) -> NTuple{3,Float64}

Rotate a grid-frame unit vector back into world coordinates (`R' * p`).
"""
function from_grid(o::Orientation, p::NTuple{3,Float64})
    o.identity && return p
    R = o.R
    return (R[1] * p[1] + R[4] * p[2] + R[7] * p[3],
        R[2] * p[1] + R[5] * p[2] + R[8] * p[3],
        R[3] * p[1] + R[6] * p[2] + R[9] * p[3])
end
