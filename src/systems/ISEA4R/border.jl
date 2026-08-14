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
# around the block, walked lazily by the package's diamond-quadtree descent
# wherever the block is nowhere flush with its diamond's edge: a non-flush
# block's halo is entirely in-diamond, where adjacency is the plain 3×3 lattice
# (every irregular neighbourhood — the corner fans, the valence-3 corners —
# sits on a diamond edge), so the band IS the halo, minus its four corners
# under `Edge()`. A flush block's halo crosses `EDGE_NEIGHBORS` onto other
# diamonds and takes the generic materialising engine instead. `morton_to_xyd`
# of the block's first id names its lattice origin: min-Morton is min-corner.
function DGG.neighbor_engine(sys::ISEA4RSystem, c::DGG.LevelIndex, target::Int,
        connectivity::DGG.Connectivity)
    lo, side = _isea4r_square(sys, c, target)
    if side > 1
        n = _nside(target)
        ix, iy, diamond = morton_to_xyd(lo, n)
        x0, y0 = Int64(ix), Int64(iy)
        if 1 <= x0 && x0 + side <= n - 1 && 1 <= y0 && y0 + side <= n - 1
            return DGG.SquareHaloEngine(DGG.MortonCurve(), Int64(diamond) * n * n,
                target, n, 0x0, x0, y0, side, connectivity isa DGG.Vertex)
        end
    end
    return DGG.generic_neighbor_engine(sys, c, target, connectivity)
end

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
