# ---------------------------------------------------------------------------
# The 3x3 lattice neighbourhood, on the nested id
#
# Transcription of the HEALPix C++ `neighbors()` (healpix_base.cc), lifted here
# from the pre-redesign `src/HEALPix/HealpixLookups.jl`. Two things changed in
# the lift, and both are corrections rather than tidying:
#
#   * it no longer rides Healpix.jl. The old version called
#     `Healpix.pix2xyfNest` / `Healpix.xyf2pixNest` to get in and out of the
#     lattice, which put a C-library-shaped dependency in the middle of an
#     otherwise pure-Julia system and — more to the point — put the neighbour
#     tables in the LOOKUP layer, which the old kernel then reached back into.
#     The chart's own `nested_to_xyf` / `xyf_to_nested` are the same codec, so
#     the tables belong next to them.
#   * the results are DEDUPLICATED. At `nside == 1` (level 0) every offset
#     wraps through the face tables, and distinct offsets can land on the same
#     base pixel; the interface contract says a neighbour list has no
#     duplicates, so `nested_neighbors` is filtered by the callers below rather
#     than handed out raw.
#
# ## Geometry of the neighbourhood, and what `Edge()` means here
#
# `nested_neighbors` answers in the reference implementation's compass order,
#
#     SW, W, NW, N, NE, E, SE, S
#
# which pairs with the lattice offsets `(dx, dy)`
#
#     (-1, 0) (-1,+1) (0,+1) (+1,+1) (+1, 0) (+1,-1) (0,-1) (-1,-1)
#
# A HEALPix pixel is a diamond on the lattice: `pixel_corners` emits its
# corners as `(north, west, south, east)` at `(ix+1,iy+1)`, `(ix,iy+1)`,
# `(ix,iy)`, `(ix+1,iy)`. So an offset with exactly ONE non-zero component
# shares a whole pixel EDGE, and an offset with two shares only the single
# lattice point that is one of the four named corners — which is why the
# compass label of a diagonal offset is exactly the corner it touches
# (`(+1,+1)` is labelled `N` and touches the north corner, and so on).
#
# Therefore, on this grid:
#
#   * `Vertex()` (Moore, the interface default) is all eight, and
#   * `Edge()` (von Neumann) is the four one-component offsets — the compass
#     directions SW, NW, NE, SE, at odd positions `1, 3, 5, 7` of the tuple.
#
# The four DROPPED entries are the diagonals W, N, E, S. Read that pairing
# twice before touching it: the entries whose labels look like edge directions
# are the corner-only ones, because the pixel is rotated 45 degrees against the
# lattice. Getting it backwards produces a plausible-looking neighbour set that
# fails the conformance suite's `Edge() ⊆ Vertex()` and symmetry clauses.
#
# The 24 pixels sitting on a degree-3 vertex of the base tiling have only seven
# neighbours; the table yields `-1` for the missing one and it is skipped.
#
# Level 0 is its own case: with `nside == 1` a base pixel IS a face, every
# offset wraps through the face tables, and ALL TWELVE base pixels come out
# with six neighbours — the polar faces (0-3, 8-11) lose the `E` and `W`
# diagonals, the equatorial faces (4-7) lose `N` and `S`. So the base tiling is
# 6-regular, and the "24 pixels with seven" rule starts at level 1.
#
# ## Rotational order
#
# The interface's canonical `neighbors` order is ROTATIONAL: counter-clockwise
# seen from OUTSIDE the sphere, from a documented starting direction. The
# compass tuple above is a natural cycle for that, but it runs the WRONG WAY —
# its offsets sweep 180°, 135°, 90°, ... in the `(dx, dy)` lattice plane, i.e.
# clockwise, and the chart is orientation-PRESERVING (the `(x, y)` lattice maps
# to `(φ, z)` with positive Jacobian, and `(φ, z)` is right-handed seen from
# outside — which is also why `pixel_corners`' N, W, S, E order is CCW). So the
# tuple order is clockwise on the sphere and the cycle below REVERSES it.
#
# This is measured, not just argued: summing the signed azimuth steps of the
# neighbour centres around the cell centre gives exactly -2π for the tuple
# order at every level and on polar and equatorial faces alike, and +2π for the
# cycle below. `test/systems/HEALPix/runtests.jl` keeps that assertion running.
# ---------------------------------------------------------------------------

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

The eight lattice neighbours of 0-based nested pixel `pix` at refinement
`level`, in compass order `SW, W, NW, N, NE, E, SE, S`. Entries that do not
exist — the one missing diagonal at each of the 24 pixels on a degree-3
base-tiling vertex — are `-1`.

Pure integer arithmetic on the chart codec: no geometry, no tables beyond the
face-transition ones above, no Healpix.jl. The result may contain repeats at
`level == 0`, where every offset wraps; callers deduplicate.
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

The positions of [`nested_neighbors`](@ref)' compass tuple, in **counter-
clockwise order seen from outside the sphere, starting at `SW`**.

`Vertex()` gives the full cycle `SW, S, SE, E, NE, N, NW, W` — the compass
tuple reversed (see the file header: the tuple itself runs clockwise), rotated
so that `SW`, its own first entry, stays first.

`Edge()` is that same cycle restricted to its four members, order preserved:
`SW, SE, NE, NW`. Restricting a cycle cannot change its winding, so the
`Edge()` order is CCW for free.

Absent neighbours — the missing diagonal at the 24 degree-3 pixels, and the two
missing entries at every level-0 base pixel — simply drop out of the cycle.
"""
_neighbor_cycle(::DGG.Vertex) = (1, 8, 7, 6, 5, 4, 3, 2)
_neighbor_cycle(::DGG.Edge) = (1, 7, 5, 3)
