# ===========================================================================
# Analytic vertices, for GridInterface's dual verbs.
#
# A HEALPix pixel's corners are N, W, S, E (the order of `cell_corners`), and
# its compass neighbours come from `nested_neighbors` in the slots SW, W, NW,
# N, NE, E, SE, S (1..8), `-1` where a degree-3 vertex has no diagonal
# neighbour. Every vertex is a lattice point `(face, x, y)` with `x, y ∈ 0:N`,
# `N = nside`, and the faces' shared edges give each one a single owner, so
# the vertices have dense numbers without a table.
# ===========================================================================

using SmallCollections: push, pushfirst

# Around corner k (N, W, S, E) of pixel p, counter-clockwise from p: the
# compass slots of the three other cells (one is missing at a degree-3 vertex).
const _HP_FAN = ((5, 4, 3), (3, 2, 1), (1, 8, 7), (7, 6, 5))
# Lattice offset of corner k (N, W, S, E) from the pixel's (ix, iy).
const _HP_COFF = ((1, 1), (0, 1), (0, 0), (1, 0))

GridInterface.vertices_are_analytic(::HEALPixLevelGrid) = true
GridInterface.vertex_numbering_is_analytic(::HEALPixLevelGrid) = true

@inline function _hp_nb(g::HEALPixLevelGrid, i::Integer)
    1 <= i <= DGG.ncells(g) || throw(BoundsError(g, i))
    return nested_neighbors(Int(i) - 1, g.level)
end

# The positions in the given compass slots, skipping absent ones and `i`.
@inline function _hp_pick(::Val{C}, nb, i::Int, slots::NTuple{K,Int}) where {C,K}
    out = SmallVector{C,Int}()
    for m in slots
        q = nb[m]
        q < 0 && continue
        p = Int(q) + 1
        (p == i || p in out) && continue
        out = push(out, p)
    end
    return out
end

@inline function _hp_corner(g::HEALPixLevelGrid, i::Int, k::Int)
    nside = 1 << g.level
    ix, iy, f = nested_to_xyf(i - 1, nside)
    return pixel_corners(ix, iy, f, nside)[k]
end

@inline function _hp_fan(g::HEALPixLevelGrid, i::Int, k::Int)
    1 <= k <= 4 || throw(BoundsError(1:4, k))
    nb = _hp_nb(g, i)
    return pushfirst(_hp_pick(Val(4), nb, i, _HP_FAN[k]), i)
end

@inline _d2(a, b) = (a[1] - b[1])^2 + (a[2] - b[2])^2 + (a[3] - b[3])^2

# The key of corner k of cell i: owned by the lowest position of its fan, as
# that cell's corner — found on the lattice on the same face, by the nearest
# corner across a face edge.
function GridInterface.vertex_key(g::HEALPixLevelGrid, i::Integer, k::Integer)
    fan = _hp_fan(g, Int(i), Int(k))
    o = minimum(fan)
    o == i && return GridInterface.VertexKey(Int(i), Int(k))
    nside = 1 << g.level
    ix, iy, f = nested_to_xyf(Int(i) - 1, nside)
    ox, oy, fo = nested_to_xyf(o - 1, nside)
    if fo == f          # same face: the corner index follows from the lattice
        vx, vy = ix + _HP_COFF[k][1], iy + _HP_COFF[k][2]
        for j in 1:4
            (ox + _HP_COFF[j][1], oy + _HP_COFF[j][2]) == (vx, vy) &&
                return GridInterface.VertexKey(o, j)
        end
    end
    p = pixel_corners(ix, iy, f, nside)[k]       # across a face edge: nearest corner
    cs = pixel_corners(ox, oy, fo, nside)
    j = 1; d = _d2(cs[1], p)
    for jj in 2:4
        dj = _d2(cs[jj], p)
        dj < d && ((j, d) = (jj, dj))
    end
    return GridInterface.VertexKey(o, j)
end
GridInterface.vertex_cells(g::HEALPixLevelGrid, key::GridInterface.VertexKey) =
    _hp_fan(g, key.cell, key.corner)
GridInterface.vertex_point(g::HEALPixLevelGrid, key::GridInterface.VertexKey) =
    _hp_corner(g, key.cell, key.corner)

# Dense vertex numbers, O(1). Lattice point (f, x, y), x, y ∈ 0:N, N = nside.
# Every face owns x ∈ 0:N-1, y ∈ 1:N (its SW and NW edges and W corner, not
# its NE or SE edge): each shared edge is NE/SE on one side and SW/NW on the
# other (face tables in `neighbors.jl`), so each non-pole vertex has one owner.
# The poles (N corner of faces 0:3, S corner of 8:11) are left over and
# numbered last.
#   v = f N² + (y - 1) N + x + 1,   north pole 12N² + 1,   south pole 12N² + 2

# Owner of lattice point (f, x, y), or face -1 / -2 for the north / south pole.
# Crossing NE (x = N) or SE (y = 0): north f ↔ (f+1)%4 via (N,y) ↔ (y,N) and
# ↔ 4+(f+1)%4 via (x,0) ↔ (x,N); equatorial 4+j ↔ north j via (N,y) ↔ (0,y) and
# ↔ 8+j via (x,0) ↔ (x,N); south 8+j ↔ 4+(j+1)%4 via (N,y) ↔ (0,y) and
# ↔ 8+(j+1)%4 via (x,0) ↔ (0,x). An E corner (N, 0) crosses twice.
@inline function _hp_owner(f::Int, x::Int, y::Int, N::Int)
    r, j1 = f >> 2, (f + 1) & 3
    if x == N
        if r == 0
            y == N && return (-1, 0, 0)
            return (j1, y, N)
        end
        y > 0 && return (r == 1 ? f & 3 : 4 + j1, 0, y)
        return (r == 1 ? 4 + j1 : 8 + j1, 0, N)
    elseif y == 0
        r == 0 && return (4 + j1, x, N)
        r == 1 && return (f + 4, x, N)
        x == 0 && return (-2, 0, 0)
        return (8 + j1, 0, x)
    end
    return (f, x, y)
end

# Euler: V = F + 2 for a quad tessellation of the sphere.
function GridInterface.nvertices(g::HEALPixLevelGrid)
    N = 1 << g.level
    return 12 * N * N + 2
end

function GridInterface.vertex_key(g::HEALPixLevelGrid, v::Integer)
    lv = g.level; N = 1 << lv; nf = N * N
    1 <= v <= 12nf + 2 || throw(BoundsError(1:12nf+2, v))
    v > 12nf && return v == 12nf + 1 ? GridInterface.vertex_key(g, nf, 1) :
        GridInterface.vertex_key(g, 8nf + 1, 3)
    w = Int(v) - 1
    f, r = w >> 2lv, w & (nf - 1)
    x, y = r & (N - 1), r >> lv + 1
    # The owner point is the W corner of (x, y-1).
    return GridInterface.vertex_key(g, xyf_to_nested(x, y - 1, f, N) + 1, 2)
end

function GridInterface.vertex_index(g::HEALPixLevelGrid, key::GridInterface.VertexKey)
    lv = g.level; N = 1 << lv
    k = key.corner
    1 <= k <= 4 || throw(BoundsError(1:4, k))
    ix, iy, f = nested_to_xyf(key.cell - 1, N)
    fo, x, y = _hp_owner(f, ix + _HP_COFF[k][1], iy + _HP_COFF[k][2], N)
    fo < 0 && return 12 * N * N - fo
    return (fo << 2lv) + ((y - 1) << lv) + x + 1
end
