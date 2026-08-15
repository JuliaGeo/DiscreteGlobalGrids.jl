# Standard ISEA icosahedron, adjacency tables, spherical helpers, and frame
# orientation. Tables remain in the grid frame; `Orientation` rotates inputs
# and outputs without changing them.

# ---------------------------------------------------------------------------
# Sphere helpers (degrees on the boundary, radians internally)
# ---------------------------------------------------------------------------

"""
    lonlat_to_xyz(lon, lat) -> NTuple{3,Float64}

Unit vector of a spherical (lon, lat) pair in degrees, right-handed with
`+z` at the north pole and `+x` at (0, 0). Lon/lat are plain spherical
coordinates; no geodetic/authalic latitude conversion is applied.
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

WGS84 authalic Earth radius in metres. It enters the grid only through areas:
`hex_area(r) = 4πR²/(10·7^r)` and `pentagon = 5/6 * hex`.
"""
const R_AUTHALIC = 6371007.180918475

"`cos` of the icosahedron edge arc, `1/√5` (arc `acosd(1/√5) = 63.4349...°`)."
const ADJ_DOT = 1 / sqrt(5.0)

"""
    NBASE

Number of base cells and icosahedron vertices.
"""
const NBASE = 12

# Pentagon-chain digit deletion is defined by `Z7_DELETED_DIGIT` in `z7.jl`.

"""
    ISEA_LON0

Longitude of standard ISEA vertex 0: `11.25` degrees (`π/16` radians).
"""
const ISEA_LON0 = 11.25

"""
    ISEA_LAT_HI

Latitude of standard ISEA vertex 0: `atand(φ) = 58.282525588538995` degrees.
The 12 vertices sit at latitudes
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
at `(11.25°E, 58.28252559°N)`, azimuth 0.
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

Neighbor base at development-frame angle zero for each base. The same chirality
and digit map applies to all bases.
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

Row-major rotation from world coordinates to the ISEA grid frame.
[`to_grid`](@ref) applies `R`; [`from_grid`](@ref) applies `R'`. The cached
`identity` flag is derived from `R`.
"""
struct Orientation
    R::NTuple{9,Float64}
    identity::Bool
    # Cache whether `R` is exactly the identity rotation.
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
