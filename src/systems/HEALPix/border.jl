# A nested id is
# `face * nside^2 + morton(ix, iy)` (`xyf_to_nested`), and refining Δ levels
# scales the lattice by `s = 2^Δ` on the SAME face: pixel `(ix, iy, face)` at
# `level` covers exactly `[ix*s, (ix+1)*s) x [iy*s, (iy+1)*s)` of the leaf
# lattice. Bit-interleaving is positional, so that scaling is a shift of the
# Morton code:
#
#     morton(ix*s + dx, iy*s + dy) = morton(ix, iy) * s^2 + morton(dx, dy)
#
# Thus the high bits remain the parent code and the low `2Δ` bits are free:
#
#   * the subtree is the contiguous id range `[id * 4^Δ, (id+1) * 4^Δ)`, which
#     is exactly what `descendant_range` returns; and
#   * a descendant's OFFSET within that range *is* `morton(dx, dy)`, so
#     ascending id order over the subtree is Morton order over the block, and
#     the operation reduces to "emit the perimeter of an `s x s` square in
#     Morton order".
#
# In the 3x3 lattice neighbourhood, a descendant has all neighbours inside iff
# `0 < dx < s-1` and `0 < dy < s-1`. The rim is
# `dx ∈ {0, s-1} || dy ∈ {0, s-1}`, of size `4s - 4`.
#
# Edge and vertex connectivity yield the same rim; seam crossings leave the
# subtree face (no non-centre row of `NB_FACEARRAY` maps a face to itself), and
# missing diagonals at degree-3 vertices do not affect rim membership.
#
# That leaves nothing HEALPix-specific in the walk itself: it is the package's
# `SquareRimEngine` under `MortonCurve`, shared with ISEA4R (same curve) and S2
# (Hilbert), which descends the quadtree in curve order and prunes a quadrant
# inheriting none of the parent square's exposed sides.

# Extended, not shadowed: `HEALPix.subtree_border` and
# `DiscreteGlobalGrids.subtree_border` are the same function, so generic code
# gets the Morton walk without knowing HEALPix has one.
import ..DiscreteGlobalGrids: subtree_border

# `descendant_range` runs FIRST and is the only level guard needed: it raises the
# `target < level(c)` and `> max_level` `ArgumentError`s, and its `lo` is the
# subtree's first nested id. Nothing below re-derives either.
function _healpix_square(sys::HEALPixSystem, c::DGG.LevelIndex, target::Int)
    r = DGG.descendant_range(sys, c, target)
    _checked_index(c)
    lo = Int64(first(r)) - 1                 # back to the 0-based nested id
    return lo, Int64(1) << (target - DGG.level(c))
end

function DGG.rim_engine(sys::HEALPixSystem, c::DGG.LevelIndex, target::Int,
        connectivity::DGG.Connectivity)
    lo, side = _healpix_square(sys, c, target)
    return DGG.SquareRimEngine(DGG.MortonCurve(), lo, target, side, 0x0)
end

function DGG.interior_engine(sys::HEALPixSystem, c::DGG.LevelIndex, target::Int,
        connectivity::DGG.Connectivity)
    lo, side = _healpix_square(sys, c, target)
    return DGG.SquareInteriorEngine(DGG.MortonCurve(), lo, target, side, 0x0)
end

# The halo — the outside face of the same boundary — is the width-1 band
# around the block, walked lazily by the package's face-quadtree descent. Away
# from the face edge that band is entirely in-face, where adjacency is the plain
# 3×3 lattice (every irregular neighbourhood in `NB_FACEARRAY` sits on a face
# edge), so the band IS the halo, minus its four corners under `Edge()`. Flush
# with the edge it crosses the seam — `NB_SWAPARRAY` territory, ids on other
# faces — and `square_halo_engine` derives the candidate rectangles by asking
# `neighbors` about a few rim cells and filtering every candidate with the
# native one-ring. Neither case needs a seam table here.
#
# `nested_to_xyf` of the block's first id names its lattice origin, because
# min-Morton is min-corner. (That step is Morton-specific and does NOT
# generalise: see S2, where the Hilbert curve's first id is any corner of the
# block and the origin has to come from the parent's own lattice coordinates
# instead.)
#
# `side == 1` is depth zero, which the generic engine answers with the cell's
# own one-ring — exact at the degree-3 vertices, where a band of one is not.
function DGG.halo_engine(sys::HEALPixSystem, c::DGG.LevelIndex, target::Int,
        connectivity::DGG.Connectivity)
    # Before `_healpix_square`, whose own guard is `descendant_range`'s and says
    # something else about the same user error. The halo verb's wording is one
    # wording on every system; see `check_halo_level`.
    DGG.check_halo_level(sys, c, target)
    lo, side = _healpix_square(sys, c, target)
    side == 1 && return DGG.generic_halo_engine(sys, c, target, connectivity)
    n = _nside(target)
    ix, iy, face = nested_to_xyf(lo, n)
    return DGG.square_halo_engine(sys, DGG.MortonCurve(), c, target, connectivity,
        Int64(ix), Int64(iy), side, Int64(face), n)
end

# The three lines the shared square walk needs from a system. `_nside` is the
# face's lattice side at a level, and the Morton curve's state never changes, so
# every face root is read under `0x0`.
DGG.lattice_decode(sys::HEALPixSystem, c::DGG.LevelIndex) =
    nested_to_xyf(c.index, _nside(DGG.level(c)))

DGG.lattice_cell(sys::HEALPixSystem, l::Int, ix::Integer, iy::Integer,
    face::Integer) = DGG.LevelIndex(l, xyf_to_nested(ix, iy, face, _nside(l)))

DGG.face_orientation(sys::HEALPixSystem, face::Integer) = 0x0

"""
    subtree_border(sys::HEALPixSystem, c::LevelIndex, leaf_level::Integer; connectivity = Vertex()) -> Vector{LevelIndex}

The **rim** of `c`'s subtree at `leaf_level`: the descendants that have at
least one neighbour outside the subtree, ascending by canonical id.

`4 * 2^Δ - 4` cells for `Δ = leaf_level - level(c) > 0`, and `[c]` at `Δ == 0`.
Θ(rim) time and one allocation — the rim is read straight off the leaf lattice,
with no neighbour query at any level. `collect` of
[`EdgeCellIterator`](@ref), which is the same walk without the allocation.

`connectivity` does not change the result: edge adjacency already identifies
every boundary cell of the square lattice block.

Positions are `rawid + 1`.
"""
subtree_border(sys::HEALPixSystem, c::DGG.LevelIndex, leaf_level::Integer;
    connectivity::DGG.Connectivity=DGG.Vertex()) =
    DGG.collect_subtree(DGG.EdgeCellIterator(sys, c, leaf_level; connectivity))
