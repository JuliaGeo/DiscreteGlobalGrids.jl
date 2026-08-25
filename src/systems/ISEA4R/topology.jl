# Integer topology for axis-aligned diamond cells. Axis offsets are edge
# neighbors and diagonal offsets are vertex-only neighbors. Border crossings pair
# reversed oriented edges; corner fans handle valence-3 vertices and the
# valence-5 vertices 0 and 11. Results wind counterclockwise from `(+1,0)`; both
# affine halves and `snyder_inv_xyz` preserve orientation, so chart-CCW is CCW
# seen from outside the sphere, which is the order the interface asks for.

# Border edge slots, in boundary-walk order. `S, E, N, W` are the chart directions
# `y = 0`, `x = 1`, `y = 1`, `x = 0` — chart directions, not compass ones: a
# diamond is not aligned with anything on the globe.
const EDGE_S = 1
const EDGE_E = 2
const EDGE_N = 3
const EDGE_W = 4

# Corner slots, in the same boundary-walk order: chart `(0,0)`, `(1,0)`,
# `(1,1)`, `(0,1)`, i.e. the `Diamond.verts` slots.
const CORNER_00 = 1
const CORNER_10 = 2
const CORNER_11 = 3
const CORNER_01 = 4

"""
    _edge_pair(d, e) -> (a, b)

The ordered base-vertex pair of border edge slot `e` of diamond `d`, in the
counterclockwise boundary walk of the chart square.
"""
function _edge_pair(d::Int, e::Int)
    v00, v10, v11, v01 = @inbounds DIAMONDS[d+1].verts
    e == EDGE_S && return (v00, v10)
    e == EDGE_E && return (v10, v11)
    e == EDGE_N && return (v11, v01)
    return (v01, v00)
end

"""
    EDGE_NEIGHBORS

`EDGE_NEIGHBORS[d + 1][e] == (d', e')`: the diamond and border edge slot on the
other side of slot `e` of diamond `d`.

Derived by matching reversed ordered vertex pairs and asserted at load time to
be an involution over all forty slots.
"""
const EDGE_NEIGHBORS = let
    out = Vector{NTuple{4,Tuple{Int,Int}}}(undef, 10)
    for d in 0:9
        out[d+1] = ntuple(4) do e
            a, b = _edge_pair(d, e)
            hits = [(d2, e2) for d2 in 0:9 for e2 in 1:4 if _edge_pair(d2, e2) == (b, a)]
            @assert length(hits) == 1 "border edge ($a, $b) of diamond $d has $(length(hits)) partners, not 1"
            hits[1]
        end
    end
    for d in 0:9, e in 1:4
        d2, e2 = out[d+1][e]
        @assert out[d2+1][e2] == (d, e) "the border-edge pairing is not involutive at ($d, $e)"
        @assert (d2, e2) != (d, e) "border edge slot $e of diamond $d is paired with itself"
    end
    ntuple(i -> out[i], 10)
end

"""
    _corner_vertex(d, c) -> Int

The base vertex id at corner slot `c` of diamond `d` — `DIAMONDS[d+1].verts[c]`,
named for what it is.
"""
_corner_vertex(d::Int, c::Int) = @inbounds DIAMONDS[d+1].verts[c]

# The two border edges meeting at corner slot `c`, in counterclockwise sweep order:
# the one whose direction PRECEDES the corner's diagonal and the one that
# FOLLOWS it. Corner `c` sits between edge slots `c-1` and `c` in the boundary
# walk, and the walk is counterclockwise, so those are exactly the two.
@inline _corner_edge_before(c::Int) = mod1(c - 1, 4)
@inline _corner_edge_after(c::Int) = c

"""
    CORNER_FANS

Other diamond corners at the same icosahedron vertex as `(d,c)`, in
counterclockwise fan order. The endpoints are reached by axis offsets; interior
entries are the additional diagonal neighbors at vertices 0 and 11.
"""
const CORNER_FANS = let
    out = Vector{NTuple{4,Vector{Tuple{Int,Int}}}}(undef, 10)
    for d in 0:9
        out[d+1] = ntuple(4) do c
            v = _corner_vertex(d, c)
            fan = Tuple{Int,Int}[]
            cur_d, via = d, _corner_edge_before(c)
            while true
                nd, ne = EDGE_NEIGHBORS[cur_d+1][via]
                # The corner of `nd` on the same icosahedron vertex. The border
                # edge we arrived through has `v` as an endpoint, so `nd` really
                # does have a corner there, and a diamond's four corners are
                # four distinct vertices, so it is unique.
                nc = findfirst(==(v), DIAMONDS[nd+1].verts)
                @assert nc !== nothing "vertex $v is not a corner of diamond $nd"
                (nd, nc) == (d, c) && break
                push!(fan, (nd, nc))
                # Leave through the OTHER border edge at this corner.
                via = ne == _corner_edge_before(nc) ? _corner_edge_after(nc) :
                      _corner_edge_before(nc)
                cur_d = nd
                @assert length(fan) <= 5 "the fan around vertex $v does not close"
            end
            @assert length(fan) in (2, 4) "vertex $v carries $(length(fan) + 1) diamond-corners"
            @assert allunique(fan) "the fan around vertex $v repeats a diamond-corner"
            # The walk STARTS across the preceding border edge by construction; that
            # it also ENDS across the following one is the substantive claim, and
            # it is what makes the interior entries exactly the cells no axis
            # offset reaches.
            @assert fan[end] == let nd = first(EDGE_NEIGHBORS[d+1][_corner_edge_after(c)])
                (nd, findfirst(==(v), DIAMONDS[nd+1].verts))
            end "the fan around vertex $v does not end across the following border edge"
            fan
        end
    end
    # Every vertex's fans agree on how many diamond-corners it carries, and the
    # forty corner-slots account for the twelve vertices exactly.
    counts = zeros(Int, 12)
    for d in 0:9, c in 1:4
        counts[_corner_vertex(d, c)+1] += 1
    end
    @assert counts == [5, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 5] "diamond-corner valences drifted: $counts"
    ntuple(i -> out[i], 10)
end

"""
    _edge_cell(e, j, nside) -> (ix, iy)

The border cell of diamond-local slot `e` at boundary-walk index `j`, using
`S: j=ix`, `E: j=iy`, `N: j=nside-1-ix`, and `W: j=nside-1-iy`.
"""
@inline function _edge_cell(e::Int, j::Int, nside::Int)
    e == EDGE_S && return (j, 0)
    e == EDGE_E && return (nside - 1, j)
    e == EDGE_N && return (nside - 1 - j, nside - 1)
    return (0, nside - 1 - j)
end

"""
    _corner_cell(c, nside) -> (ix, iy)

The corner cell of a diamond at corner slot `c`.
"""
@inline function _corner_cell(c::Int, nside::Int)
    c == CORNER_00 && return (0, 0)
    c == CORNER_10 && return (nside - 1, 0)
    c == CORNER_11 && return (nside - 1, nside - 1)
    return (0, nside - 1)
end

# Which corner slot the extended lattice position `(xx, yy)` names, when both
# coordinates are off the lattice. Only the four diagonal offsets at the four
# corner cells reach here, so the sign of each coordinate decides it outright.
@inline function _corner_slot(xx::Int, yy::Int)
    xx < 0 && return yy < 0 ? CORNER_00 : CORNER_01
    return yy < 0 ? CORNER_10 : CORNER_11
end

"""
    NEIGHBOR_OFFSETS

The eight lattice offsets `(dx, dy)`, **counterclockwise in the chart plane
starting at `(+1, 0)`**: E, NE, N, NW, W, SW, S, SE in chart directions.

Slots `1, 3, 5, 7` — the four axis offsets — are the [`Edge`](@ref) subset:
those share a whole cell edge, while the diagonals share only a lattice point.
"""
const NEIGHBOR_OFFSETS = ((1, 0), (1, 1), (0, 1), (-1, 1),
    (-1, 0), (-1, -1), (0, -1), (1, -1))

"""
    _offset_slots(connectivity) -> Tuple

The indices of [`NEIGHBOR_OFFSETS`](@ref) that `connectivity` admits, in
counterclockwise order from the same start.

`Vertex()` is all eight; `Edge()` is the four axis offsets `1, 3, 5, 7`.
Restricting a cycle cannot change its winding, so the `Edge()` order is
counterclockwise for free.
"""
_offset_slots(::DGG.Vertex) = (1, 2, 3, 4, 5, 6, 7, 8)
_offset_slots(::DGG.Edge) = (1, 3, 5, 7)

"""
    lattice_neighbors(ix, iy, diamond, nside, connectivity) -> SmallVector{9,NTuple{3,Int}}

Return `(ix,iy,diamond)` neighbors counterclockwise from chart direction
`(+1,0)`. Vertex connectivity yields 8 cells in the interior, 7 at valence-3
corners, and 9 at vertices 0 and 11; edge connectivity yields at most 4.
Integer table lookup is allocation-free, and duplicates at `nside == 1` are
removed without changing order.
"""
function lattice_neighbors(ix::Integer, iy::Integer, diamond::Integer,
        nside::Integer, connectivity::DGG.Connectivity = DGG.Vertex())
    n = Int(nside)
    x = Int(ix)
    y = Int(iy)
    d = Int(diamond)
    out = SmallVector{9,NTuple{3,Int}}()
    for m in _offset_slots(connectivity)
        dx, dy = @inbounds NEIGHBOR_OFFSETS[m]
        xx = x + dx
        yy = y + dy
        inx = 0 <= xx < n
        iny = 0 <= yy < n
        if inx && iny
            out = _pushunique(out, (xx, yy, d))
        elseif inx || iny
            # One border crossing: read the extended square's edge index, mirror it
            # into the paired slot.
            e, j = xx >= n ? (EDGE_E, yy) :
                   xx < 0 ? (EDGE_W, n - 1 - yy) :
                   yy >= n ? (EDGE_N, n - 1 - xx) : (EDGE_S, xx)
            d2, e2 = @inbounds EDGE_NEIGHBORS[d+1][e]
            jx, jy = _edge_cell(e2, n - 1 - j, n)
            out = _pushunique(out, (jx, jy, d2))
        else
            # An icosahedron vertex: everything in the fan except the two ends,
            # which the axis offsets on either side of this diagonal already
            # emitted.
            c = _corner_slot(xx, yy)
            fan = @inbounds CORNER_FANS[d+1][c]
            for i in 2:(length(fan) - 1)
                d2, c2 = fan[i]
                cx, cy = _corner_cell(c2, n)
                out = _pushunique(out, (cx, cy, d2))
            end
        end
    end
    return out
end

@inline _pushunique(out::SmallVector, v) = v in out ? out : SmallCollections.push(out, v)
