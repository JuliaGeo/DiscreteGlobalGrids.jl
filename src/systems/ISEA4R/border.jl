# A subtree is a square lattice block, and a Morton id refines by
# `id * 4^Δ + morton(dx, dy)` (see `descendant_range`), so a descendant's offset
# within the subtree's contiguous id range IS its Morton code in the block.
# Ascending id order over the subtree is therefore Morton order over the block,
# and the rim is the block's perimeter: `4*2^Δ - 4` cells.
#
# A boundary cell always has an AXIS neighbour outside the block — across a
# diamond rim as readily as within one — so neither connectivity nor a block
# sitting against a seam needs a special case.
#
# The walk is the package's `SquareRimEngine` under `MortonCurve`, shared with
# HEALPix (same curve) and S2 (Hilbert): descend the quadtree in curve order and
# prune a quadrant that inherits none of the parent square's exposed sides. That
# emits in ascending id by construction, so the perimeter is no longer built and
# sorted — the cost drops from `O(rim log rim)` to `O(rim)`, in `O(Δ)` memory.

function DGG.rim_engine(sys::ISEA4RSystem, c::DGG.LevelIndex, target::Int,
        connectivity::DGG.Connectivity)
    lo, side = _isea4r_square(sys, c, target)
    return DGG.SquareRimEngine(DGG.MortonCurve(), lo, target, side, 0x0)
end

function DGG.interior_engine(sys::ISEA4RSystem, c::DGG.LevelIndex, target::Int,
        connectivity::DGG.Connectivity)
    lo, side = _isea4r_square(sys, c, target)
    return DGG.SquareInteriorEngine(DGG.MortonCurve(), lo, target, side, 0x0)
end

# `descendant_range` is the only level guard needed — it raises both
# `ArgumentError`s — and its `lo` is the subtree's first Morton id.
function _isea4r_square(sys::ISEA4RSystem, c::DGG.LevelIndex, target::Int)
    r = DGG.descendant_range(sys, c, target)
    _checked_index(c)
    lo = Int64(first(r)) - 1
    return lo, Int64(1) << (target - DGG.level(c))
end

# The halo — the outside face of the same boundary — is the width-1 band
# around the block, walked lazily by the package's diamond-quadtree descent.
# Away from the diamond edge that band is entirely in-diamond, where adjacency
# is the plain 3×3 lattice (every irregular neighbourhood — the corner fans, the
# valence-3 corners — sits on a diamond edge), so the band IS the halo, minus
# its four corners under `Edge()`. Flush with the edge it crosses
# `EDGE_NEIGHBORS` onto other diamonds, and at an icosahedral vertex
# `CORNER_FANS` onto a third: `square_halo_engine` derives those candidates by
# asking `neighbors` about a few rim cells, and filters every one of them with
# the native one-ring. No seam table is read here.
#
# `morton_to_xyd` of the block's first id names its lattice origin: min-Morton
# is min-corner. (Morton-specific, and not true of S2's Hilbert curve — see its
# `border.jl`.) `side == 1` is depth zero, which the generic engine answers with
# the cell's own one-ring — exact at vertices 0 and 11, where a nine-cell
# neighbourhood is not a band of one.
function DGG.halo_engine(sys::ISEA4RSystem, c::DGG.LevelIndex, target::Int,
        connectivity::DGG.Connectivity)
    lo, side = _isea4r_square(sys, c, target)
    side == 1 && return DGG.generic_halo_engine(sys, c, target, connectivity)
    n = _nside(target)
    ix, iy, diamond = morton_to_xyd(lo, n)
    return DGG.square_halo_engine(sys, DGG.MortonCurve(), c, target, connectivity,
        Int64(ix), Int64(iy), side, Int64(diamond), n)
end

# The three lines the shared square walk needs from a system. The Morton curve's
# state never changes, so every diamond root is read under `0x0`.
DGG.lattice_decode(sys::ISEA4RSystem, c::DGG.LevelIndex) =
    morton_to_xyd(c.index, _nside(DGG.level(c)))

DGG.lattice_cell(sys::ISEA4RSystem, l::Int, ix::Integer, iy::Integer,
    diamond::Integer) = DGG.LevelIndex(l, xyd_to_morton(ix, iy, diamond, _nside(l)))

DGG.face_orientation(sys::ISEA4RSystem, diamond::Integer) = 0x0

"""
    subtree_border(ISEA4RSystem(), c, l; connectivity = Vertex()) -> Vector{LevelIndex}

Return the boundary of the `2^Δ × 2^Δ` descendant block in ascending order:
`4*2^Δ - 4` cells. Connectivity does not change this rim. Cost is `O(rim)`.
`collect` of [`EdgeCellIterator`](@ref), which is the same walk lazily.

`subtree_border(sys, c, level(c))` is `[c]`; `l < level(c)` throws an
`ArgumentError`, uniformly with every other system.
"""
DGG.subtree_border(sys::ISEA4RSystem, c::DGG.LevelIndex, l::Integer;
    connectivity::DGG.Connectivity = DGG.Vertex()) =
    DGG.collect_subtree(DGG.EdgeCellIterator(sys, c, l; connectivity))
