# HEALPix 3×3 lattice neighbourhood in compass order
# `SW, W, NW, N, NE, E, SE, S`. Because pixels are diamond-shaped in this
# lattice, one-axis offsets share edges and two-axis offsets share corners;
# `Edge()` therefore selects positions 1, 3, 5, and 7. The reference tuple is
# clockwise on the sphere, so `_neighbor_cycle` reverses it. Degree-3 vertices
# have one missing neighbour, and level-0 results require deduplication.
#
# The offset and face-transition tables below are transcribed from the HEALPix
# C++ `neighbors()` (healpix_base.cc). At level 0 every offset wraps through
# them and all twelve base pixels come out with exactly six distinct
# neighbours; the degree-3 "seven" case starts at level 1.

const NB_XOFFSET = (-1, -1, 0, 1, 1, 1, 0, -1)
const NB_YOFFSET = (0, 1, 1, 1, 0, -1, -1, -1)
const NB_FACEARRAY = (
    (8, 9, 10, 11, -1, -1, -1, -1, 10, 11, 8, 9),   # S
    (5, 6, 7, 4, 8, 9, 10, 11, 9, 10, 11, 8),       # SE
    (-1, -1, -1, -1, 5, 6, 7, 4, -1, -1, -1, -1),   # E
    (4, 5, 6, 7, 11, 8, 9, 10, 11, 8, 9, 10),       # SW
    (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11),         # center
    (1, 2, 3, 0, 0, 1, 2, 3, 5, 6, 7, 4),           # NE
    (-1, -1, -1, -1, 7, 4, 5, 6, -1, -1, -1, -1),   # W
    (3, 0, 1, 2, 3, 0, 1, 2, 4, 5, 6, 7),           # NW
    (2, 3, 0, 1, -1, -1, -1, -1, 0, 1, 2, 3),       # N
)
const NB_SWAPARRAY = (
    (0, 0, 3), (0, 0, 6), (0, 0, 0), (0, 0, 5), (0, 0, 0),
    (5, 0, 0), (0, 0, 0), (6, 0, 0), (3, 0, 0),
)

"""
    nested_neighbors(pix, level) -> NTuple{8,Int64}

Return the eight lattice neighbours of 0-based nested pixel `pix` in compass
order `SW, W, NW, N, NE, E, SE, S`. Missing degree-3 neighbours are `-1`;
level-0 results may repeat and must be deduplicated by callers.
"""
function nested_neighbors(pix::Integer, level::Integer)
    nside = Int64(1) << Int(level)
    x, y, f = nested_to_xyf(pix, nside)
    return ntuple(8) do m
        xx = Int64(x) + NB_XOFFSET[m]
        yy = Int64(y) + NB_YOFFSET[m]
        nbnum = 5                                   # 1-based "center" row
        if xx < 0
            xx += nside; nbnum -= 1
        elseif xx >= nside
            xx -= nside; nbnum += 1
        end
        if yy < 0
            yy += nside; nbnum -= 3
        elseif yy >= nside
            yy -= nside; nbnum += 3
        end
        nf = NB_FACEARRAY[nbnum][f + 1]
        nf < 0 && return Int64(-1)
        # The face transition may mirror or transpose the lattice; these are the
        # reference implementation's swap bits, applied before re-encoding.
        bits = NB_SWAPARRAY[nbnum][(f >> 2) + 1]
        (bits & 1) != 0 && (xx = nside - xx - 1)
        (bits & 2) != 0 && (yy = nside - yy - 1)
        (bits & 4) != 0 && ((xx, yy) = (yy, xx))
        return Int64(xyf_to_nested(xx, yy, nf, nside))
    end
end

"""
    _neighbor_cycle(connectivity) -> Tuple

Positions in [`nested_neighbors`](@ref) ordered counter-clockwise from `SW` as
seen from outside the sphere. `Vertex()` yields `SW, S, SE, E, NE, N, NW, W`;
`Edge()` yields `SW, SE, NE, NW`. Callers omit absent neighbours.
"""
_neighbor_cycle(::DGG.Vertex) = (1, 8, 7, 6, 5, 4, 3, 2)
_neighbor_cycle(::DGG.Edge) = (1, 7, 5, 3)
