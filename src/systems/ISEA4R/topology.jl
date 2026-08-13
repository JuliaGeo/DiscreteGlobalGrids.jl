# ---------------------------------------------------------------------------
# The rhombus lattice one-ring, and the ten-diamond seam topology
#
# Inside a diamond the neighbourhood is the plainest thing in this package: the
# chart cell of `(ix, iy)` is the axis-aligned square `[ix, ix+1] × [iy, iy+1]`,
# so the four axis offsets share a whole cell EDGE and the four diagonal offsets
# share a single lattice POINT. (Note how this differs from HEALPix, where the
# pixel is a diamond rotated 45° against the lattice and the pairing is the
# other way round. Same lattice, opposite reading — do not transplant one file's
# `Edge()` slot tuple into the other.)
#
# What needs deriving is what happens at a diamond's rim. Everything below comes
# out of one observation about the layout table:
#
#   > A diamond's four chart edges are four icosahedron edges, named by their
#   > endpoint vertex pairs; the SEAM (the main diagonal `y == x`) is a fifth,
#   > and it is interior. Ten diamonds contribute forty rim edge-slots and ten
#   > seams — exactly the icosahedron's thirty edges, the ten seams used once
#   > and the twenty others twice.
#
# So the diamond-to-diamond adjacency is read straight off the vertex pairs,
# with no geometry and no fitting.
#
# ## Edge slots, and why the pairing needs no orientation flag
#
# Number a diamond's rim edges by the counterclockwise boundary walk of its
# chart square, `(0,0) → (1,0) → (1,1) → (0,1) → (0,0)`:
#
#     slot 1  S  y = 0   from v00 to v10      slot 3  N  y = 1   from v11 to v01
#     slot 2  E  x = 1   from v10 to v11      slot 4  W  x = 0   from v01 to v00
#
# Read as ORDERED vertex pairs those are `(v00,v10)`, `(v10,v11)`, `(v11,v01)`,
# `(v01,v00)`. Two diamonds sharing an edge traverse it in opposite directions —
# the ten charts are consistently oriented (each is counterclockwise seen from
# outside, which `diamonds.jl` asserts), and that is exactly what an orientable
# surface's face boundaries do. So the neighbour of `(d, slot)` is the unique
# `(d', slot')` whose ordered pair is the REVERSE, and the parameter along the
# edge always runs backwards across the join. There is no per-edge "reversed"
# bit to fit or to get wrong; the [`EDGE_NEIGHBORS`](@ref) build asserts the
# reversal is unique and involutive rather than assuming it.
#
# ## Crossing an edge
#
# Give each rim cell of slot `e` an EDGE INDEX `j ∈ 0:nside-1`, its position
# along the boundary walk:
#
#     S: j = ix     E: j = iy     N: j = nside-1-ix     W: j = nside-1-iy
#
# and let `_edge_cell` invert that. Then crossing is `j' = nside - 1 - j` into
# the paired slot — the reversal above, on the lattice.
#
# The rule that makes this cover the diagonal offsets too is to read an
# out-of-range offset as the extended chart square `[xx, xx+1] × [yy, yy+1]`
# beyond the rim and ask which edge index IT spans:
#
#     xx == nside  ->  slot 2 (E), j = yy          yy == nside  ->  slot 3 (N), j = nside-1-xx
#     xx == -1     ->  slot 4 (W), j = nside-1-yy  yy == -1     ->  slot 1 (S), j = xx
#
# At a diagonal offset that lands one step further along the rim than the axis
# offset does, which is precisely the corner-only neighbour. One rule, both
# cases, and it is the same rule read forwards.
#
# ## The twelve icosahedron vertices
#
# An offset with BOTH coordinates out of range can only happen at the four
# corner cells of a diamond, and it points at an icosahedron vertex. That is
# where the lattice stops being a lattice.
#
# Count the diamond-corners at each vertex: vertex 0 is `v01` of all five
# northern diamonds and vertex 11 is `v10` of all five southern ones, so those
# two carry FIVE diamond-corners each; the remaining ten vertices carry THREE.
# (Forty corner-slots, five plus five plus ten times three. The five faces
# meeting at a vertex are shared out as 60° or 120° per diamond-corner depending
# on whether the vertex is a seam endpoint of that diamond, which is why the
# counts differ while the 300° cone angle does not.)
#
# The cells touching a vertex are exactly the corner cells of those
# diamond-corners — no other cell has the vertex as one of its four corners. Of
# the `k-1` others, TWO are already reached by the axis offsets across the two
# rim edges meeting at that corner, so the diagonal offset yields `k-3`: NONE at
# a valence-3 vertex, and TWO at vertices 0 and 11.
#
# Hence the neighbour counts: 8 in the interior, 7 at a valence-3 corner cell,
# and **9** at the twenty corner cells sitting on vertices 0 and 11 — which is
# why `max_neighbors(sys, Vertex())` is 9 and not 8.
#
# ## Rotational order
#
# The offsets are emitted counterclockwise in the `(x, y)` chart plane starting
# at `(+1, 0)`, and the chart is orientation-preserving onto the sphere seen
# from outside (`diamonds.jl` asserts `imag(conj(a)·b) > 0` on all twenty
# half-maps, and `snyder_inv_xyz` preserves orientation from a face's `(u, w)`
# frame), so chart-counterclockwise IS counterclockwise seen from outside. That
# is the argument; `test/systems/ISEA4R/runtests.jl` measures it, because the
# argument has been wrong before in this package.
#
# The two cells a vertex offset yields must be wound too. They are emitted in
# the order of the FAN WALK around the vertex, which [`CORNER_FANS`](@ref)
# builds by repeatedly crossing rim edges: starting from the rim edge that precedes the
# diagonal in the counterclockwise sweep and ending on the one that follows it,
# so the extras fall between their two axis neighbours in the same sweep.
# ---------------------------------------------------------------------------

# Rim edge slots, in boundary-walk order. `S, E, N, W` are the chart directions
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

The ordered base-vertex pair of rim edge slot `e` of diamond `d`, in the
counterclockwise boundary walk of the chart square. See the file header.
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

`EDGE_NEIGHBORS[d + 1][e] == (d', e')`: the diamond and rim edge slot on the
other side of slot `e` of diamond `d`.

Derived by matching REVERSED ordered vertex pairs (see the file header), and
asserted at load time to be a well-defined involution over all forty slots.
"""
const EDGE_NEIGHBORS = let
    out = Vector{NTuple{4,Tuple{Int,Int}}}(undef, 10)
    for d in 0:9
        out[d+1] = ntuple(4) do e
            a, b = _edge_pair(d, e)
            hits = [(d2, e2) for d2 in 0:9 for e2 in 1:4 if _edge_pair(d2, e2) == (b, a)]
            @assert length(hits) == 1 "rim edge ($a, $b) of diamond $d has $(length(hits)) partners, not 1"
            hits[1]
        end
    end
    for d in 0:9, e in 1:4
        d2, e2 = out[d+1][e]
        @assert out[d2+1][e2] == (d, e) "the rim-edge pairing is not involutive at ($d, $e)"
        @assert (d2, e2) != (d, e) "rim edge slot $e of diamond $d is paired with itself"
    end
    ntuple(i -> out[i], 10)
end

"""
    _corner_vertex(d, c) -> Int

The base vertex id at corner slot `c` of diamond `d` — `DIAMONDS[d+1].verts[c]`,
named for what it is.
"""
_corner_vertex(d::Int, c::Int) = @inbounds DIAMONDS[d+1].verts[c]

# The two rim edges meeting at corner slot `c`, in counterclockwise sweep order:
# the one whose direction PRECEDES the corner's diagonal and the one that
# FOLLOWS it. Corner `c` sits between edge slots `c-1` and `c` in the boundary
# walk, and the walk is counterclockwise, so those are exactly the two.
@inline _corner_edge_before(c::Int) = mod1(c - 1, 4)
@inline _corner_edge_after(c::Int) = c

"""
    CORNER_FANS

`CORNER_FANS[d + 1][c]`: the OTHER diamond-corners sitting on the same
icosahedron vertex as corner slot `c` of diamond `d`, as `(d', c')` pairs, in
counterclockwise fan order about that vertex.

The walk starts by crossing the rim edge that precedes the corner's diagonal
direction and ends on the one that follows it, so the first and last entries are
always the two corners the AXIS offsets already reach and the interior entries —
two of them at vertices 0 and 11, none anywhere else — are what the diagonal
offset yields. The build asserts exactly that, plus that the walk
closes on its starting corner and visits every diamond-corner of the vertex
once.
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
                # The corner of `nd` on the same icosahedron vertex. The rim
                # edge we arrived through has `v` as an endpoint, so `nd` really
                # does have a corner there, and a diamond's four corners are
                # four distinct vertices, so it is unique.
                nc = findfirst(==(v), DIAMONDS[nd+1].verts)
                @assert nc !== nothing "vertex $v is not a corner of diamond $nd"
                (nd, nc) == (d, c) && break
                push!(fan, (nd, nc))
                # Leave through the OTHER rim edge at this corner.
                via = ne == _corner_edge_before(nc) ? _corner_edge_after(nc) :
                      _corner_edge_before(nc)
                cur_d = nd
                @assert length(fan) <= 5 "the fan around vertex $v does not close"
            end
            @assert length(fan) in (2, 4) "vertex $v carries $(length(fan) + 1) diamond-corners"
            @assert allunique(fan) "the fan around vertex $v repeats a diamond-corner"
            # The walk STARTS across the preceding rim edge by construction; that
            # it also ENDS across the following one is the substantive claim, and
            # it is what makes the interior entries exactly the cells no axis
            # offset reaches.
            @assert fan[end] == let nd = first(EDGE_NEIGHBORS[d+1][_corner_edge_after(c)])
                (nd, findfirst(==(v), DIAMONDS[nd+1].verts))
            end "the fan around vertex $v does not end across the following rim edge"
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

The rim cell of diamond-local slot `e` at boundary-walk position `j` — the
inverse of the edge-index table in the file header (`S: j = ix`, `E: j = iy`,
`N: j = nside-1-ix`, `W: j = nside-1-iy`).
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
those share a whole cell edge, the diagonals share a single lattice point. See
the file header on why this is the opposite reading from HEALPix's.
"""
const NEIGHBOR_OFFSETS = ((1, 0), (1, 1), (0, 1), (-1, 1),
    (-1, 0), (-1, -1), (0, -1), (1, -1))

"""
    _offset_slots(connectivity) -> Tuple

The positions of [`NEIGHBOR_OFFSETS`](@ref) that `connectivity` admits, in
counterclockwise order from the same start.

`Vertex()` is all eight; `Edge()` is the four axis offsets `1, 3, 5, 7`.
Restricting a cycle cannot change its winding, so the `Edge()` order is
counterclockwise for free.
"""
_offset_slots(::DGG.Vertex) = (1, 2, 3, 4, 5, 6, 7, 8)
_offset_slots(::DGG.Edge) = (1, 3, 5, 7)

"""
    lattice_neighbors(ix, iy, diamond, nside, connectivity) -> SmallVector{9,NTuple{3,Int}}

The `(ix, iy, diamond)` neighbours of a cell, in **counterclockwise rotational
order seen from outside the sphere, starting at the `(+1, 0)` chart direction**.

Eight cells in a diamond's interior, seven at a corner cell on one of the ten
valence-3 icosahedron vertices, nine at a corner cell on vertex 0 or vertex 11
(where the diagonal offset yields two cells rather than none) — see the file
header for the count.

Pure integer arithmetic on the layout tables: no geometry, no floating point,
and no allocation — the result is a `SmallCollections.SmallVector` sized by the
nine of `max_neighbors`. Duplicates are removed keeping the FIRST occurrence,
which preserves the cycle; they can only arise at `nside == 1`, where a whole
diamond is one cell.
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
            # One rim crossing: read the extended square's edge index, mirror it
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
