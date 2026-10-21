# ---------------------------------------------------------------------------
# Geometry and location
#
# Geometry methods preserve libh3's cell definitions while converting results
# to unit-sphere points.
# ---------------------------------------------------------------------------

"""
    cell_boundary(::H3System, c::H3Cell) -> Helpers.SmallList{6,UnitSphericalPoint}

The cell's six true vertices (five at a pentagon) on the unit sphere,
implicitly closed and counter-clockwise seen from outside, joined by great
circles: libh3's `cellToVertexes` points, read through `vertexToLatLng` so
that two cells sharing a vertex get bit-identical coordinates.

libh3's own `cellToBoundary` adds *distortion vertices* where a class III
(odd-level) cell crosses an icosahedron edge, bending that edge in the
gnomonic projection. This ring leaves them out: every cell is one convex
spherical polygon (GridInterface's `ConvexCells`), at the cost of an edge
that differs from libh3's by the distortion: a straddling cell's area differs
from libh3's `cellAreaRads2` by exactly that sliver, up to about 12 % at every
class III level (the level still sums to `4π`).

The ring is rotated so that corner `k` and corner `k + 1` are the two
vertices shared with Edge neighbour `k` of [`one_ring`](@ref) — the alignment
GridInterface asks of a counter-clockwise grid.
"""
function cell_boundary(::H3System, c::H3Cell)
    ids = _aligned_vertex_ids(c)
    out = DGG.Helpers.empty_small_list(Val(6), USPoint(1.0, 0.0, 0.0))
    for v in ids
        v == 0 && continue
        x, y, z = H3Native.vertex_to_cartesian(v)
        out = DGG.Helpers.small_push(out, USPoint(x, y, z))
    end
    return out
end

# The cell's vertex ids, rotated so the first two are the ones it shares with
# the first neighbour of `one_ring` (both walks are counter-clockwise from
# outside, so the rest then line up edge for edge). Absent slots are `0`.
@inline function _aligned_vertex_ids(c::H3Cell)
    ids = H3Native.cell_to_vertexes(c.id)
    first_nb = _first_neighbor(c)
    first_nb == 0 && return ids
    n = _nvertex_ids(ids)
    start = _shared_edge_start(ids, H3Native.cell_to_vertexes(first_nb), n)
    start <= 1 && return ids
    return ntuple(k -> k > n ? H3Native.H3Index(0) : ids[mod1(start + k - 1, n)], Val(6))
end

@inline _nvertex_ids(ids::NTuple{6,UInt64}) = ids[6] == 0 ? 5 : 6

# The slot `k` for which vertices `k` and `k + 1` (cyclically) both belong to
# the neighbour, or 0 when no edge is shared.
@inline function _shared_edge_start(ids::NTuple{6,UInt64}, nids::NTuple{6,UInt64}, n::Int)
    for k in 1:n
        a = ids[k]
        b = ids[k == n ? 1 : k + 1]
        _has_vertex(nids, a) && _has_vertex(nids, b) && return k
    end
    return 0
end

@inline _has_vertex(ids::NTuple{6,UInt64}, v::UInt64) =
    ids[1] == v || ids[2] == v || ids[3] == v || ids[4] == v || ids[5] == v || ids[6] == v

# The first cell of `one_ring`, by the same rule it uses.
@inline function _first_neighbor(c::H3Cell)
    ring = one_ring(H3LevelGrid(H3System(), level(c)), c, Edge())
    return isempty(ring) ? H3Native.H3Index(0) : first(ring).id
end

"""
    cell_centroid(::H3System, c::H3Cell) -> UnitSphericalPoint

The cell centre from libh3's `cellToLatLng`, converted to a unit-sphere point.
It is not the mean of the boundary vertices.
"""
function cell_centroid(::H3System, c::H3Cell)
    x, y, z = H3Native.cell_center_cartesian(c.id)
    return USPoint(x, y, z)
end

# ===========================================================================
# Location
# ===========================================================================

"""
    cellat(grid::H3LevelGrid, p::UnitSphericalPoint) -> H3Cell
    cellat(grid::H3LevelGrid, lon::Real, lat::Real) -> H3Cell

The cell containing a point, computed by libh3's O(1) `latLngToCell` inverse.

Never `nothing`: H3 tessellates the whole sphere, so every point is in a cell.

Shared-edge ties follow libh3's deterministic rule.

The `(lon, lat)` overload takes degrees.
"""
function cellat(grid::H3LevelGrid, p::GO.UnitSphericalPoint)
    lon = atand(p[2], p[1])
    lat = asind(clamp(p[3], -1.0, 1.0))
    return H3Cell(H3Native.lonlat_to_cell(lon, lat, grid.level))
end

cellat(grid::H3LevelGrid, lon::Real, lat::Real) =
    H3Cell(H3Native.lonlat_to_cell(lon, lat, grid.level))

# GridInterface's position form: libh3's location, repaired to the straight
# ring `cell_boundary` reports. At class III (odd) levels libh3 locates by its
# bent boundary, which differs from the straight one only in the sliver between
# a cell edge and its distortion vertex on an icosahedron edge. A point farther
# than one cell edge from every icosahedron edge cannot be in a sliver and pays
# 20 dot products; a closer one pays one containment test, and in a sliver the
# neighbours' too.

# The icosahedron's 20 face centres (normalised sums of three adjacent pentagon
# centres) and the chord between adjacent ones: the great circle of their
# shared edge is the bisector plane of the two centres.
const _ICOSA_FACE_CENTERS = let v = [H3Native.cell_center_cartesian(p) for p in H3Native.get_pentagons(0)]
    d(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
    fs = NTuple{3,Float64}[]
    for a in 1:12, b in a+1:12, c in b+1:12
        d(v[a], v[b]) > 0.4 && d(v[a], v[c]) > 0.4 && d(v[b], v[c]) > 0.4 || continue
        s = v[a] .+ v[b] .+ v[c]
        push!(fs, s ./ sqrt(d(s, s)))
    end
    length(fs) == 20 || error("expected 20 icosahedron faces, found $(length(fs))")
    Tuple(fs)
end
const _ADJ_FACE_CHORD = 0.7136441795461800
# An upper bound on an H3 edge's length in radians times sqrt(7)^res
# (0.2206 measured at res 1, 3 and 5).
const _EDGE_LENGTH_BOUND = 0.25

# Whether `p` may lie within one cell edge of an icosahedron edge at `res`: the
# gap between the nearest and second-nearest face centre bounds the distance
# to the nearest face edge from below.
@inline function _near_icosa_edge(p, res::Int)
    d1 = d2 = -2.0
    for f in _ICOSA_FACE_CENTERS
        d = f[1] * p[1] + f[2] * p[2] + f[3] * p[3]
        if d > d1
            d2 = d1; d1 = d
        elseif d > d2
            d2 = d
        end
    end
    return (d1 - d2) / _ADJ_FACE_CHORD < _EDGE_LENGTH_BOUND / sqrt(7.0)^res
end

# Closed containment in the straight ring of `id` (the half-plane test
# GridInterface applies to `ConvexCells`, with its tolerance).
function _in_straight_ring(id, p)
    ids = H3Native.cell_to_vertexes(id)
    n = _nvertex_ids(ids)
    a1 = H3Native.vertex_to_cartesian(ids[1])
    a = a1
    for k in 1:n
        b = k == n ? a1 : H3Native.vertex_to_cartesian(ids[k + 1])
        nab = (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])
        nab[1] * p[1] + nab[2] * p[2] + nab[3] * p[3] <
            -1e-12 * sqrt(nab[1]^2 + nab[2]^2 + nab[3]^2) && return false
        a = b
    end
    return true
end

function GridInterface.cellat(grid::H3LevelGrid, p)
    q = GridInterface.check_point(GridInterface.best_manifold(grid), p)
    c = cellat(grid, q)
    if isodd(grid.level) && _near_icosa_edge(q, grid.level) && !_in_straight_ring(c.id, q)
        for nb in one_ring(grid, c, Edge())
            _in_straight_ring(nb.id, q) && return DGG.localindex(grid, nb)
        end
    end
    return DGG.localindex(grid, c)
end
